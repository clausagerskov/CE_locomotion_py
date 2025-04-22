import math
import random
import time
import sys
import os
import pickle
import threading
from enum import Enum
from typing import List, Callable, Optional, Any, Tuple, Union
from loguru import logger
# Import TVector from the other file
from vector_matrix import TVector
from random_state_legacy import RandomState

# Constants from C++
MIN_SEARCH_VALUE = -1.0
MAX_SEARCH_VALUE = 1.0

# Threading configuration
THREADED_SEARCH = False  # Set to False to disable threading
THREAD_COUNT = 16      # Number of threads for evaluation

# --- Helper Functions ---

def clip(x: float, min_val: float, max_val: float) -> float:
    """Clips a value to be within [min_val, max_val]."""
    return max(min_val, min(x, max_val))

def MapSearchParameter(x: float, min_val: float, max_val: float,
                        clipmin: float = -1.0e99, clipmax: float = 1.0e99) -> float:
    """Maps a search parameter from [-1, 1] to [min_val, max_val]."""
    m = (max_val - min_val) / (MAX_SEARCH_VALUE - MIN_SEARCH_VALUE)
    b = min_val - m * MIN_SEARCH_VALUE
    return clip(m * x + b, clipmin, clipmax)

def InverseMapSearchParameter(x: float, min_val: float, max_val: float) -> float:
    """Maps a model parameter from [min_val, max_val] back to [-1, 1]."""
    m = (MAX_SEARCH_VALUE - MIN_SEARCH_VALUE) / (max_val - min_val)
    b = MIN_SEARCH_VALUE - m * min_val
    return m * x + b

# --- Enums ---

class TSelectionMode(Enum):
    FITNESS_PROPORTIONATE = 1
    RANK_BASED = 2

class TReproductionMode(Enum):
    HILL_CLIMBING = 1
    GENETIC_ALGORITHM = 2

class TCrossoverMode(Enum):
    UNIFORM = 1
    TWO_POINT = 2

# --- Threading Helper Class ---
# Not strictly necessary in Python as args can be passed directly,
# but kept for structural similarity.
class PopRangeSpec:
    def __init__(self, search_instance, start, end):
        self.search = search_instance
        self.start = start
        self.end = end

# --- TSearch Class ---

# Define type hints for the function pointers
EvaluationFuncType = Callable[[TVector[float], RandomState], float]
BestActionFuncType = Callable[[int, TVector[float]], None]
PopStatsDisplayFuncType = Callable[[int, float, float, float], None]
SearchTerminateFuncType = Callable[[int, float, float, float], bool]
SearchResultsDisplayFuncType = Callable[["TSearch"], None] # Forward ref

class TSearch:
    def __init__(self, vectorSize: int = 0, EvalFn: Optional[EvaluationFuncType] = None):
        # Internal State
        self.rs = RandomState() # Master random state
        self.RandomStates: TVector[RandomState] = TVector(1, 0) # Per-individual random states
        self.Gen: int = 0
        self.SearchInitialized: bool = False
        self.Population: TVector[TVector[float]] = TVector(1, 0)
        self.Perf: TVector[float] = TVector(1, 0) # Performance scores
        self.fitness: TVector[float] = TVector(1, 0) # Scaled/ranked fitness
        self.UpdateBestFlag: bool = False
        self.bestVector: TVector[float] = TVector(1, 0)
        self.BestPerf: float = -1.0e99 # Initialize very low
        self.MinPerf: float = 0.0
        self.MaxPerf: float = 0.0
        self.AvgPerf: float = 0.0
        self.PerfVar: float = 0.0

        # Search Modes (Defaults)
        self.SelectMode: TSelectionMode = TSelectionMode.RANK_BASED
        self.RepMode: TReproductionMode = TReproductionMode.GENETIC_ALGORITHM
        self.CrossMode: TCrossoverMode = TCrossoverMode.TWO_POINT

        # Search Parameters (Defaults)
        self._vectorSize: int = 0 # Use property for access
        self.MaxGens: int = 0
        self.EFraction: float = 0.0 # Elitist Fraction
        self.MaxExpOffspring: float = 1.1
        self.MutationVar: float = 1.0
        self.CrossProb: float = 0.0
        self.crossTemplate: TVector[int] = TVector(1, 0)
        self.crossPoints: TVector[int] = TVector(1, 0)
        self.ConstraintVector: TVector[int] = TVector(1, 0) # 1=constrained, 0=unconstrained
        self.ReEvalFlag: bool = False # Re-evaluate elites/parents?
        self.CheckpointInt: int = 0 # Checkpoint interval (0 = disabled)

        # Function Pointers (Callables)
        self.EvaluationFunction: Optional[EvaluationFuncType] = EvalFn
        self.BestActionFunction: Optional[BestActionFuncType] = None
        self.PopulationStatisticsDisplayFunction: Optional[PopStatsDisplayFuncType] = None
        self.SearchTerminationFunction: Optional[SearchTerminateFuncType] = None
        self.SearchResultsDisplayFunction: Optional[SearchResultsDisplayFuncType] = None

        # Initialize vector size if provided
        if vectorSize > 0:
            self.SetVectorSize(vectorSize)
        # Set default population size (example)
        self.SetPopulationSize(1)


    # --- Basic Accessors ---
    def VectorSize(self) -> int:
        return self._vectorSize

    def SetVectorSize(self, NewSize: int):
        if NewSize <= 0:
            raise ValueError(f"Invalid vector size: {NewSize}")
        self._vectorSize = NewSize
        # Resize population vectors
        pop_size = self.PopulationSize()
        new_pop = TVector[TVector[float]](1, pop_size)
        for i in range(1, pop_size + 1):
             # Create new vectors; preserving old data is complex, often not needed here
            new_pop[i] = TVector[float](1, NewSize)
            if i <= self.Population.Size() and self.Population[i].Size() == NewSize:
                 # If sizes match, maybe copy? C++ version implicitly reallocates.
                 # Let's re-init, consistent with C++ destructor/constructor effect.
                 pass # Keep newly created vector
            elif i <= self.Population.Size():
                 # Sizes differ, cannot directly copy
                 pass # Keep newly created vector
        self.Population = new_pop

        # Adjust bestVector
        self.bestVector.SetBounds(1, NewSize) # Preserves data if size increases/decreases
        # Reset crossover template and points
        v = TVector[int](1, NewSize)
        for i in range(1, NewSize + 1):
            v[i] = i
        self.SetCrossoverTemplate(v) # This also sets crossPoints
        # Reset constraint vector
        self.ConstraintVector = TVector[int](1, NewSize)
        self.ConstraintVector.FillContents(1) # Default: constrained


    def SetRandomSeed(self, seed: int):
        self.rs.seed(seed)
        # Also re-seed individual random states if they exist
        for i in range(1, self.RandomStates.Size() + 1):
            # Generate seeds deterministically from the main seed for reproducibility
             individual_seed = self.rs.UniformRandomInteger(0, 2**32 - 1)
             self.RandomStates[i].seed(individual_seed)
        # Reset main generator after seeding individuals
        self.rs.seed(seed)


    # --- Search Mode Accessors ---
    def SelectionMode(self) -> TSelectionMode: return self.SelectMode
    def SetSelectionMode(self, newmode: TSelectionMode): self.SelectMode = newmode
    def ReproductionMode(self) -> TReproductionMode: return self.RepMode
    def SetReproductionMode(self, newmode: TReproductionMode): self.RepMode = newmode
    def CrossoverMode(self) -> TCrossoverMode: return self.CrossMode
    def SetCrossoverMode(self, newmode: TCrossoverMode): self.CrossMode = newmode

    # --- Search Parameter Accessors ---
    def PopulationSize(self) -> int: return self.Population.Size()

    def SetPopulationSize(self, NewSize: int):
        if NewSize <= 0:
             raise ValueError(f"Invalid population size: {NewSize}")
        current_size = self.Population.Size()

        # Resize main population storage
        # Create new individual vectors only if needed
        new_pop = TVector[TVector[float]](1, NewSize)
        for i in range(1, NewSize + 1):
             if i <= current_size:
                 new_pop[i] = self.Population[i]
                 # Ensure individual vector size is correct if vectorSize changed prior
                 if new_pop[i].Size() != self._vectorSize:
                      new_pop[i].SetBounds(1, self._vectorSize)
             else:
                 new_pop[i] = TVector[float](1, self._vectorSize)
        self.Population = new_pop

        # Resize Perf and fitness vectors
        self.Perf.SetBounds(1, NewSize)
        self.fitness.SetBounds(1, NewSize)

        # Resize and initialize RandomStates
        old_rs_size = self.RandomStates.Size()
        new_rs = TVector[RandomState](1, NewSize)
        # Preserve old random states if possible
        for i in range(1, NewSize + 1):
            if i <= old_rs_size:
                new_rs[i] = self.RandomStates[i]
            else:
                # Seed new random states deterministically
                individual_seed = self.rs.UniformRandomInteger(0, 2**32 - 1)
                new_rs[i] = RandomState(individual_seed)
        self.RandomStates = new_rs


    def MaxGenerations(self) -> int: return self.MaxGens
    def SetMaxGenerations(self, NewMax: int):
        if NewMax < 0: raise ValueError(f"Invalid MaxGenerations: {NewMax}")
        self.MaxGens = NewMax

    def ElitistFraction(self) -> float: return self.EFraction
    def SetElitistFraction(self, NewFraction: float):
        if not (0.0 <= NewFraction <= 1.0): raise ValueError(f"Invalid ElitismFraction: {NewFraction}")
        self.EFraction = NewFraction

    def MaxExpectedOffspring(self) -> float: return self.MaxExpOffspring
    def SetMaxExpectedOffspring(self, NewVal: float):
        if not (1.0 <= NewVal <= 2.0): raise ValueError(f"Invalid MaxExpectedOffspring: {NewVal}")
        self.MaxExpOffspring = NewVal

    def MutationVariance(self) -> float: return self.MutationVar
    def SetMutationVariance(self, NewVariance: float):
        if NewVariance <= 0.0: raise ValueError(f"Invalid MutationVariance: {NewVariance}")
        self.MutationVar = NewVariance

    def CrossoverProbability(self) -> float: return self.CrossProb
    def SetCrossoverProbability(self, NewProb: float):
        if not (0.0 <= NewProb <= 1.0): raise ValueError(f"Invalid CrossoverProbability: {NewProb}")
        self.CrossProb = NewProb

    def CrossoverTemplate(self) -> TVector[int]: return self.crossTemplate
    def SetCrossoverTemplate(self, NewTemplate: TVector[int]):
        if NewTemplate.Size() != self._vectorSize:
            raise ValueError(f"Invalid vector size for CrossoverTemplate: {NewTemplate.Size()}")
        # Validate format (monotonically increasing integers starting from 1)
        x = 1
        valid = True
        if NewTemplate.Size() > 0 and NewTemplate[NewTemplate.LowerBound()] != 1:
            valid = False
        else:
            for i in range(NewTemplate.LowerBound(), NewTemplate.UpperBound() + 1):
                if NewTemplate[i] == x:
                    continue
                elif NewTemplate[i] == x + 1:
                    x += 1
                else:
                    valid = False
                    break
        if not valid:
            raise ValueError(f"Invalid format for CrossoverTemplate: {NewTemplate}")

        # Use copy assignment for TVector
        self.crossTemplate = TVector(NewTemplate.LowerBound(), NewTemplate.UpperBound(), NewTemplate.to_list()) # Deep copy

        # Modify CrossoverPoints appropriately
        num_points = x
        self.crossPoints = TVector[int](1, num_points)
        self.crossPoints[1] = NewTemplate.LowerBound() # Assumes template is 1-based
        current_point = 1
        for i in range(NewTemplate.LowerBound() + 1, NewTemplate.UpperBound() + 1):
            if NewTemplate[i] != NewTemplate[i-1]:
                current_point += 1
                self.crossPoints[current_point] = i

    def CrossoverPoints(self) -> TVector[int]: return self.crossPoints
    def SetCrossoverPoints(self, NewPoints: TVector[int]):
         # Validation similar to C++
        if NewPoints.Size() < 1 or NewPoints[NewPoints.LowerBound()] != self._vectorSize > 0: # Check first point is start
            raise ValueError(f"Invalid format for Crossover Points: {NewPoints}")
        last_point = 0
        valid = True
        for i in range(NewPoints.LowerBound(), NewPoints.UpperBound() + 1):
            if NewPoints[i] > last_point and NewPoints[i] <= self._vectorSize + NewPoints.LowerBound() - 1: # Adjust for bounds
                 last_point = NewPoints[i]
            else:
                 valid = False
                 break
        if not valid:
             raise ValueError(f"Invalid format for Crossover Points: {NewPoints}")

        self.crossPoints = TVector(NewPoints.LowerBound(), NewPoints.UpperBound(), NewPoints.to_list()) # Deep copy

        # Modify CrossoverTemplate appropriately
        self.crossTemplate = TVector[int](1, self._vectorSize) # Assuming 1-based template
        template_val = 1
        point_idx = NewPoints.LowerBound()
        next_point_start = NewPoints[point_idx] if point_idx <= NewPoints.UpperBound() else self._vectorSize + 1

        for i in range(1, self._vectorSize + 1):
            if point_idx < NewPoints.UpperBound() and i >= NewPoints[point_idx + 1]:
                point_idx += 1
                template_val += 1
                next_point_start = NewPoints[point_idx]
            self.crossTemplate[i] = template_val


    def SearchConstraint(self) -> TVector[int]: return self.ConstraintVector
    def SetSearchConstraint(self, Constraint: Union[TVector[int], int]):
        if isinstance(Constraint, int):
            self.ConstraintVector.FillContents(Constraint)
        elif isinstance(Constraint, TVector):
            if Constraint.Size() != self._vectorSize:
                raise ValueError("Invalid vector size for SearchConstraint")
            # Perform deep copy if necessary, assuming TVector assignment is shallow
            self.ConstraintVector = TVector(Constraint.LowerBound(), Constraint.UpperBound(), Constraint.to_list())
        else:
            raise TypeError("Constraint must be an int or TVector[int]")

    def ReEvaluationFlag(self) -> bool: return self.ReEvalFlag
    def SetReEvaluationFlag(self, flag: bool): self.ReEvalFlag = bool(flag)

    def CheckpointInterval(self) -> int: return self.CheckpointInt
    def SetCheckpointInterval(self, NewInterval: int):
        if NewInterval < 0: raise ValueError(f"Invalid CheckpointInterval: {NewInterval}")
        self.CheckpointInt = NewInterval

    # --- Function Pointer Accessors ---
    def SetEvaluationFunction(self, EvalFn: EvaluationFuncType): self.EvaluationFunction = EvalFn
    def SetBestActionFunction(self, BestFn: BestActionFuncType): self.BestActionFunction = BestFn
    def SetPopulationStatisticsDisplayFunction(self, DisplayFn: PopStatsDisplayFuncType): self.PopulationStatisticsDisplayFunction = DisplayFn
    def SetSearchTerminationFunction(self, TerminationFn: SearchTerminateFuncType): self.SearchTerminationFunction = TerminationFn
    def SetSearchResultsDisplayFunction(self, DisplayFn: SearchResultsDisplayFuncType): self.SearchResultsDisplayFunction = DisplayFn

    # --- Status Accessors ---
    def Generation(self) -> int: return self.Gen
    def Individual(self, i: int) -> TVector[float]: return self.Population[i]
    def Fitness(self, i: int) -> float: return self.fitness[i] # Scaled/Ranked
    def Performance(self, i: int) -> float: return self.Perf[i] # Raw score
    def BestPerformance(self) -> float: return self.BestPerf
    def BestIndividual(self) -> TVector[float]: return self.bestVector

    # --- Control ---
    def InitializeSearch(self):
        """Initializes a new search."""
        self.Gen = 0
        logger.info("Randomize Population")
        self.RandomizePopulation()
        self.BestPerf = -1.0e99 # Reset best performance
        self.SearchInitialized = True

    def ExecuteSearch(self):
        """Executes a search from scratch."""
        self._DoSearch(ResumeFlag=False)

    def ResumeSearch(self):
        """Resumes a search from a checkpoint file."""
        try:
            self.ReadCheckpointFile()
            logger.info(f"Resuming search from generation {self.Gen}")
            self._DoSearch(ResumeFlag=True)
        except FileNotFoundError:
            logger.info("Checkpoint file 'search.cpt' not found. Starting new search.")
            self.ExecuteSearch()
        except Exception as e:
            logger.info(f"Error reading checkpoint file: {e}. Starting new search.")
            self.ExecuteSearch()

    # --- Input and Output ---
    def WriteCheckpointFile(self, filename: str = "search.cpt"):
        """Saves the current state of the search to a file using pickle."""
        checkpoint_data = {
            'vectorSize': self._vectorSize,
            'populationSize': self.PopulationSize(),
            'Gen': self.Gen,
            'MaxGens': self.MaxGens,
            'rs_state': self.rs.getstate(),
            'SelectMode': self.SelectMode,
            'RepMode': self.RepMode,
            'CrossMode': self.CrossMode,
            'SearchInitialized': self.SearchInitialized,
            'ReEvalFlag': self.ReEvalFlag,
            'CheckpointInt': self.CheckpointInt,
            'ConstraintVector': self.ConstraintVector, # Save TVector directly
            'MutationVar': self.MutationVar,
            'CrossProb': self.CrossProb,
            'crossTemplate': self.crossTemplate, # Save TVector directly
            'EFraction': self.EFraction,
            'MaxExpOffspring': self.MaxExpOffspring,
            'BestPerf': self.BestPerf,
            'bestVector': self.bestVector, # Save TVector directly
            'Population': self.Population, # Save TVector[TVector]
            'Perf': self.Perf, # Save TVector directly
            'RandomStates': self.RandomStates # Save TVector[RandomState]
            # Note: Function pointers (callables) are NOT saved.
            # Note: Fitness is not saved as it's recalculated.
        }
        try:
            with open(filename, 'wb') as f:
                pickle.dump(checkpoint_data, f)
            #print(f"Checkpoint saved to {filename} at generation {self.Gen}")
        except Exception as e:
            logger.info(f"Error writing checkpoint file {filename}: {e}")

    def ReadCheckpointFile(self, filename: str = "search.cpt"):
        """Loads the state of the search from a file using pickle."""
        with open(filename, 'rb') as f:
            checkpoint_data = pickle.load(f)

        # Restore state carefully, checking types and sizes
        self.SetVectorSize(checkpoint_data['vectorSize'])
        self.SetPopulationSize(checkpoint_data['populationSize']) # Must be after SetVectorSize

        self.Gen = checkpoint_data['Gen']
        self.SetMaxGenerations(checkpoint_data['MaxGens'])
        self.rs.setstate(checkpoint_data['rs_state'])
        self.SetSelectionMode(checkpoint_data['SelectMode'])
        self.SetReproductionMode(checkpoint_data['RepMode'])
        self.SetCrossoverMode(checkpoint_data['CrossMode'])
        self.SearchInitialized = checkpoint_data['SearchInitialized']
        self.SetReEvaluationFlag(checkpoint_data['ReEvalFlag'])
        self.SetCheckpointInterval(checkpoint_data['CheckpointInt'])
        self.SetSearchConstraint(checkpoint_data['ConstraintVector']) # Assumes TVector loaded correctly
        self.SetMutationVariance(checkpoint_data['MutationVar'])
        self.SetCrossoverProbability(checkpoint_data['CrossProb'])
        self.SetCrossoverTemplate(checkpoint_data['crossTemplate']) # Assumes TVector loaded correctly
        self.SetElitistFraction(checkpoint_data['EFraction'])
        self.SetMaxExpectedOffspring(checkpoint_data['MaxExpOffspring'])
        self.BestPerf = checkpoint_data['BestPerf']
        self.bestVector = checkpoint_data['bestVector'] # Assumes TVector loaded correctly
        self.Population = checkpoint_data['Population'] # Assumes TVector[TVector] loaded correctly
        self.Perf = checkpoint_data['Perf'] # Assumes TVector loaded correctly
        self.RandomStates = checkpoint_data['RandomStates'] # Assumes TVector[RandomState] loaded correctly

        # Validate loaded data if necessary (e.g., check sizes)
        if self.Population.Size() != self.PopulationSize() or \
           self.Perf.Size() != self.PopulationSize() or \
           self.RandomStates.Size() != self.PopulationSize() or \
           (self.Population.Size() > 0 and self.Population[1].Size() != self.VectorSize()):
            raise ValueError("Checkpoint data size mismatch.")

        # Fitness needs to be recalculated if needed before reproduction
        self.fitness.SetBounds(1, self.PopulationSize())


    # --- Private Helper Methods ---

    def _DoSearch(self, ResumeFlag: bool):
        """The main internal search loop."""
        if not self.SearchInitialized and not ResumeFlag:
            logger.info("Initialize Search")
            self.InitializeSearch()

        if self.EvaluationFunction is None:
            raise ValueError("Error: NULL evaluation function")

        if not ResumeFlag:
            logger.info("Evaluating initial population...")
            self.EvaluatePopulation()
            self.BestPerf = -1.0e99 # Ensure the first calculated best is stored
            self.UpdateBestFlag = False # Reset flag
        else:
            # Resuming: BestPerf and bestVector are already loaded.
            # Fitness needs recalculation if reproduction follows immediately.
            # C++ version does this implicitly by calling UpdatePopStats -> Sort -> UpdateFitness.
             pass # No immediate action needed here.

        logger.info("Updating population statistics...")
        self.UpdatePopulationStatistics()
        self.DisplayPopulationStatistics()

        if self.UpdateBestFlag and self.BestActionFunction is not None:
            self.BestActionFunction(self.Gen, self.bestVector)

        while not self.SearchTerminated():
            self.Gen += 1
            self.UpdateBestFlag = False # Reset for this generation

            logger.info(f"--- Generation {self.Gen} ---")
            logger.info("Reproducing population...")
            self._ReproducePopulation() # This calls EvaluatePopulation internally

            logger.info("Updating population statistics...")
            self.UpdatePopulationStatistics()
            self.DisplayPopulationStatistics()

            if self.UpdateBestFlag and self.BestActionFunction is not None:
                self.BestActionFunction(self.Gen, self.bestVector)

            if self.CheckpointInt > 0 and (self.Gen % self.CheckpointInt) == 0:
                 logger.info(f"Writing checkpoint at generation {self.Gen}...")
                 self.WriteCheckpointFile()

        logger.info("Search finished.")
        self.DisplaySearchResults()

    def _EqualVector(self, v1: TVector[float], v2: TVector[float]) -> bool:
        """Compares two TVectors for equality (using TVector's __eq__)."""
        return v1 == v2

    def RandomizeVector(self, v: TVector[float]):
        """Fills a vector with random values in the search range."""
        for i in range(v.LowerBound(), v.UpperBound() + 1):
            v[i] = self.rs.UniformRandom(MIN_SEARCH_VALUE, MAX_SEARCH_VALUE)

    def RandomizePopulation(self):
        """Randomizes all individuals in the population."""
        for i in range(1, self.Population.Size() + 1):
            logger.info(f"Randomize Vector {i}/{self.Population.Size()}")
            self.RandomizeVector(self.Population[i])

    def EvaluateVector(self, Vector: TVector[float], ind_rs: RandomState) -> float:
        """Evaluates a single individual's vector."""
        if self.EvaluationFunction is None:
            raise ValueError("Evaluation function not set.")
        perf = self.EvaluationFunction(Vector, ind_rs)
        return max(0.0, perf) # Ensure non-negative performance


    # --- Threaded Evaluation ---
    def _EvaluatePopulationRange(self, start_idx: int, end_idx: int, results: list, lock: threading.Lock):
        """Target function for evaluation threads."""
        local_results = []
        for i in range(start_idx, end_idx + 1):
            logger.info(f"  Evaluating individual {i}/{self.PopulationSize()}...")
            try:
                # Fetch the correct random state for this individual
                individual_rs = self.RandomStates[i]
                perf = self.EvaluateVector(self.Population[i], individual_rs)
                local_results.append((i, perf))
            except Exception as e:
                logger.info(f"Error evaluating individual {i}: {e}")
                local_results.append((i, 0.0)) # Assign zero fitness on error

        # Safely append results to the shared list
        with lock:
            results.extend(local_results)

    def EvaluatePopulation(self, start: int = 1):
        """Evaluates the population from index start to end."""
        if start > self.PopulationSize():
            return # Nothing to evaluate

        num_to_evaluate = self.PopulationSize() - start + 1
        if num_to_evaluate <= 0:
            return

        if THREADED_SEARCH and num_to_evaluate > 1 and THREAD_COUNT > 1:
            threads = []
            results = []
            lock = threading.Lock()
            individuals_per_thread = math.ceil(num_to_evaluate / THREAD_COUNT)

            for i in range(THREAD_COUNT):
                thread_start_idx = start + i * individuals_per_thread
                thread_end_idx = min(start + (i + 1) * individuals_per_thread - 1, self.PopulationSize())

                if thread_start_idx > thread_end_idx: # No individuals left for this thread
                    continue

                # Create PopRangeSpec equivalent info as args
                thread_args = (thread_start_idx, thread_end_idx, results, lock)
                thread = threading.Thread(target=self._EvaluatePopulationRange, args=thread_args)
                threads.append(thread)
                thread.start()

            # Wait for all threads to complete
            for thread in threads:
                thread.join()

            # Update Perf vector from results (results may be out of order)
            for idx, perf_val in results:
                self.Perf[idx] = perf_val

        else: # Serial evaluation
            for i in range(start, self.PopulationSize() + 1):
                try:
                     individual_rs = self.RandomStates[i]
                     logger.info(f"EvaluateVector {i}")
                     self.Perf[i] = self.EvaluateVector(self.Population[i], individual_rs)
                except Exception as e:
                     print(f"Error evaluating individual {i}: {e}")
                     self.Perf[i] = 0.0


    def SortPopulation(self):
        """Sorts the population based on performance (descending)."""
        pop_size = self.PopulationSize()
        if pop_size < 2:
            return

        # Create pairs of (performance, index) for sorting
        # Use negative performance for descending sort
        perf_index_pairs = [(-self.Perf[i], i) for i in range(1, pop_size + 1)]
        perf_index_pairs.sort()

        # Create new sorted population and performance vectors
        sorted_pop = TVector[TVector[float]](1, pop_size)
        sorted_perf = TVector[float](1, pop_size)
        # Keep track of original random states to move them accordingly
        sorted_rs = TVector[RandomState](1, pop_size)

        for new_idx, (_, old_idx) in enumerate(perf_index_pairs, start=1):
            sorted_pop[new_idx] = self.Population[old_idx]
            sorted_perf[new_idx] = self.Perf[old_idx]
            sorted_rs[new_idx] = self.RandomStates[old_idx]

        # Replace old vectors with sorted ones
        self.Population = sorted_pop
        self.Perf = sorted_perf
        self.RandomStates = sorted_rs

    # --- Fitness Scaling ---
    def _LinearScaleFactor(self, min_p, max_p, avg_p, FMultiple):
        """Calculates linear scaling factor (from C++)."""
        if max_p == min_p: # Avoid division by zero if all performances are equal
            return 0.0
        # Check if scaled min > 0
        scaled_min_threshold = (FMultiple * avg_p - max_p) / (FMultiple - 1.0) if FMultiple != 1.0 else -float('inf')
        if min_p > scaled_min_threshold:
            delta = max_p - avg_p
            return (FMultiple - 1.0) * avg_p / delta if delta > 0.0 else 0.0
        else:
            delta = avg_p - min_p
            return avg_p / delta if delta > 0.0 else 0.0


    def UpdatePopulationFitness(self):
        """Calculates normalized fitness based on performance and selection mode."""
        pop_size = self.PopulationSize()
        if pop_size == 0: return

        self.SortPopulation() # Population is now sorted by performance (desc)

        if self.SelectMode == TSelectionMode.FITNESS_PROPORTIONATE:
            # Recalculate stats on potentially *unsorted* original performance if needed
            # But C++ version uses stats calculated *before* sort for scaling.
            # Let's stick to C++ logic: use pre-calculated MinPerf, MaxPerf, AvgPerf
            m = self._LinearScaleFactor(self.MinPerf, self.MaxPerf, self.AvgPerf, self.MaxExpOffspring)
            total_fitness = 0.0
            raw_fitness = [0.0] * pop_size # Temp list for raw scaled fitness

            for i in range(1, pop_size + 1):
                # Perf[i] is now the performance of the i-th *ranked* individual
                # We need the original performance associated with this individual *before* scaling
                # This is complex. Let's assume scaling happens *after* sorting on the sorted Perf values.
                # This matches Goldberg's description less but fits the C++ code structure better.
                f_raw = m * (self.Perf[i] - self.AvgPerf) + self.AvgPerf
                f_raw = max(0.0, f_raw) # Ensure non-negative fitness
                raw_fitness[i-1] = f_raw
                total_fitness += f_raw

            if total_fitness <= 0.0: # Handle case where all fitness is zero
                 # Assign equal probability if total is zero
                 for i in range(1, pop_size + 1):
                      self.fitness[i] = 1.0 / pop_size if pop_size > 0 else 0.0
            else:
                 for i in range(1, pop_size + 1):
                      self.fitness[i] = raw_fitness[i-1] / total_fitness

        elif self.SelectMode == TSelectionMode.RANK_BASED:
             # Fitness based on rank (index 'i' after sorting)
             if pop_size == 1:
                  self.fitness[1] = 1.0
             else:
                 for i in range(1, pop_size + 1):
                      # Formula from Goldberg/Mitchell, adapted from C++
                      rank_fitness = (self.MaxExpOffspring + (2.0 - 2.0 * self.MaxExpOffspring) * ((i - 1.0) / (pop_size - 1.0)))
                      self.fitness[i] = rank_fitness / pop_size # Normalize

        else:
            raise ValueError("Invalid selection mode")

        # Optional: Verify sum(fitness) is close to 1.0
        # print(f"Fitness sum: {sum(self.fitness[i] for i in range(1, pop_size + 1))}")

    # --- Reproduction ---
    def _ReproducePopulation(self):
        """Creates the next generation."""
        if self.RepMode == TReproductionMode.HILL_CLIMBING:
            logger.info("Reproduce Population using Hill Climbing")
            self._ReproducePopulationHillClimbing()
        elif self.RepMode == TReproductionMode.GENETIC_ALGORITHM:
            logger.info("Reproduce Population using Genetic Algorithm")
            self._ReproducePopulationGeneticAlgorithm()
        else:
            raise ValueError("Invalid reproduction mode")

    def _StochasticUniversalSampling(self, num_to_select: int) -> List[int]:
        """Performs Stochastic Universal Sampling based on self.fitness. (Corrected)"""
        pop_size = self.PopulationSize()
        selected_indices = []
        if pop_size == 0 or num_to_select <= 0: # Added check for num_to_select <= 0
            return selected_indices

        # Ensure fitness vector bounds are valid for population size
        if self.fitness.Size() != pop_size or self.fitness.LowerBound() != 1:
             raise RuntimeError(f"Fitness vector bounds ({self.fitness.LowerBound()}..{self.fitness.UpperBound()}) mismatch population size {pop_size}")

        # Assuming self.fitness is already calculated and normalized (sums to ~1.0)
        # You could add a check here:
        # fitness_sum = sum(self.fitness[i] for i in range(1, pop_size + 1))
        # if not math.isclose(fitness_sum, 1.0, abs_tol=1e-5):
        #     print(f"Warning: Fitness sum is {fitness_sum}, should be close to 1.0", file=sys.stderr)

        if num_to_select == 0 : return [] # Handle edge case

        step_size = 1.0 / num_to_select
        start_pointer = self.rs.UniformRandom(0.0, step_size)
        pointers = [start_pointer + i * step_size for i in range(num_to_select)] # i = 0, 1, ... num_to_select-1

        cumulative_fitness = 0.0
        current_member_idx = 1 # Start checking from the first member (index 1)

        for pointer in pointers:
            # Move forward through members until cumulative fitness exceeds the pointer
            while cumulative_fitness < pointer:
                 # Check if we have exceeded the number of members
                 if current_member_idx > pop_size:
                      # This means pointers sum > total fitness (shouldn't happen if fitness sums to 1)
                      print(f"Warning: SUS pointer {pointer} exceeded total cumulative fitness. Index {current_member_idx} > Pop size {pop_size}.", file=sys.stderr)
                      # Fallback: select the last member
                      current_member_idx = pop_size
                      break # Exit inner while loop

                 # Add fitness of the current member
                 cumulative_fitness += self.fitness[current_member_idx]

                 # If adding this member's fitness *now* makes cumulative >= pointer,
                 # then this is the member we select. But the standard way is to increment
                 # *after* adding, so the index selects the member whose range *contains* the pointer.
                 # Let's stick to the standard: increment index only if cumulative is still less.
                 if cumulative_fitness >= pointer:
                     break # Found the member, its index is current_member_idx

                 current_member_idx += 1 # Move to check the next member in the next inner loop iteration

            # Append the index of the member whose fitness range contained the pointer
            # If the loop exited because current_member_idx > pop_size, this will append the fallback value (pop_size)
            selected_indices.append(current_member_idx)

            # --- Crucial: Do NOT reset cumulative_fitness or current_member_idx ---
            # They carry over for the next pointer comparison.

        # Final check: Ensure all selected indices are within valid bounds [1, pop_size]
        for idx in selected_indices:
             if not (1 <= idx <= pop_size):
                  raise RuntimeError(f"SUS generated invalid index {idx}. Bounds [1, {pop_size}]")

        return selected_indices


    def _ReproducePopulationHillClimbing(self):
        pop_size = self.PopulationSize()
        if pop_size == 0: return

        # Calculate fitness (requires sorted population)
        self.UpdatePopulationFitness() # Sorts population and sets fitness

        # Parent selection (essentially selects everyone, fitness doesn't matter here)
        # We just need a copy of the current population and performance
        ParentPopulation = TVector[TVector[float]](1, pop_size)
        ParentPerf = TVector[float](1, pop_size)
        ParentRS = TVector[RandomState](1, pop_size) # Keep track of random states
        for i in range(1, pop_size + 1):
             # Use deep copies for vectors
             ParentPopulation[i] = TVector(self.Population[i].LowerBound(), self.Population[i].UpperBound(), self.Population[i].to_list())
             ParentPerf[i] = self.Perf[i]
             ParentRS[i] = self.RandomStates[i] # Keep same random state object associated


        # If ReEvalFlag is set, re-evaluate parents (though they were just evaluated)
        # C++ version does this, maybe for consistency or if eval is stochastic?
        if self.ReEvalFlag:
            logger.info("Re-evaluating parents (ReEvalFlag=True)...")
            self.EvaluatePopulation() # Re-evaluates self.Population
            # Update ParentPerf with potentially new scores
            for i in range(1, pop_size + 1):
                ParentPerf[i] = self.Perf[i]
            self.BestPerf = -1.0e99 # Reset best since parents might change perf

        # Create children by mutation
        # Mutate the *current* population vectors in-place
        logger.info("Mutating parents to create children...")
        for i in range(1, pop_size + 1):
             # Ensure we use the individual's random state for mutation consistency
            # This is different from C++, which used the global rs.
            # Let's stick to the global rs for mutation as in C++.
            self.MutateVector(self.Population[i]) # Mutates in place

        # Evaluate children
        logger.info("Evaluating children...")
        self.EvaluatePopulation() # Evaluates the modified self.Population

        # Compare children with parents
        logger.info("Performing selection (Parent vs Child)...")
        for i in range(1, pop_size + 1):
            if ParentPerf[i] > self.Perf[i]:
                # Restore parent
                self.Population[i] = ParentPopulation[i] # Restore deep copy
                self.Perf[i] = ParentPerf[i]
                # Restore random state association (although it wasn't used for eval here)
                self.RandomStates[i] = ParentRS[i]

        # Population now contains the winners of parent vs child comparisons


    def _ReproducePopulationGeneticAlgorithm(self):
        pop_size = self.PopulationSize()
        if pop_size == 0: return

        # Calculate fitness (requires sorted population)
        self.UpdatePopulationFitness() # Sorts population and sets fitness

        # Elitism: Copy top individuals directly
        ElitePop = int(round(self.EFraction * pop_size))
        NewPopulation = TVector[TVector[float]](1, pop_size)
        NewPerf = TVector[float](1, pop_size)
        NewRS = TVector[RandomState](1, pop_size) # Track random states too

        logger.info(f"Selecting {ElitePop} elite individuals...")
        for i in range(1, ElitePop + 1):
             # Deep copy elite individuals
             NewPopulation[i] = TVector(self.Population[i].LowerBound(), self.Population[i].UpperBound(), self.Population[i].to_list())
             NewPerf[i] = self.Perf[i] # Store elite performance
             NewRS[i] = self.RandomStates[i] # Keep elite random state

        # Select parents for the rest of the population using SUS
        num_parents_needed = pop_size - ElitePop
        if num_parents_needed > 0:
            logger.info(f"Selecting {num_parents_needed} parents using SUS...")
            # SUS uses self.fitness which is based on the sorted Population
            parent_indices = self._StochasticUniversalSampling(num_parents_needed)

            # Create the pool of selected parents (non-elite part)
            parent_pool = []
            for idx in parent_indices:
                # We need the individual from the *original sorted* population at index 'idx'
                parent_pool.append({
                     'vector': TVector(self.Population[idx].LowerBound(), self.Population[idx].UpperBound(), self.Population[idx].to_list()), # Deep copy
                     'rs': self.RandomStates[idx] # Associate random state
                     })

            # Shuffle the parent pool for random pairing in crossover
            # --- Manual Shuffle using self.rs (Fisher-Yates) ---
            # Replace self.rs.shuffle(parent_pool) with this:
            if self.CrossProb > 0 and len(parent_pool) > 1: # Only shuffle if crossover possible and more than 1 parent
                 print("Shuffling non-elite parents...")
                 n = len(parent_pool)
                 for i in range(n - 1, 0, -1): # Iterate from n-1 down to 1
                      # Choose random index j from [0, i] inclusive
                      j = self.rs.UniformRandomInteger(0, i) # Use our custom generator
                      # Swap elements at i and j
                      parent_pool[i], parent_pool[j] = parent_pool[j], parent_pool[i]
            # --- End of Manual Shuffle ---

            # Perform Crossover and Mutation
            child_idx = ElitePop + 1
            parent_idx = 0
            logger.info("Performing crossover and mutation...")
            while child_idx <= pop_size:
                if parent_idx >= len(parent_pool):
                     logger.info("Warning: Not enough parents generated by SUS?", file=sys.stderr)
                     # Fill remaining spots with mutations of random parents? Or elites?
                     # Let's mutate a random parent from the pool
                     if not parent_pool: # Should not happen if num_parents_needed > 0
                         rand_parent_vec = self.Population[self.rs.UniformRandomInteger(1, ElitePop)] if ElitePop > 0 else TVector(1, self.VectorSize())
                         rand_parent_rs = self.RandomStates[self.rs.UniformRandomInteger(1, ElitePop)] if ElitePop > 0 else RandomState()
                     else:
                         rand_parent_data = self.rs.choice(parent_pool)
                         rand_parent_vec = rand_parent_data['vector']
                         rand_parent_rs = rand_parent_data['rs']

                     NewPopulation[child_idx] = TVector(rand_parent_vec.LowerBound(), rand_parent_vec.UpperBound(), rand_parent_vec.to_list()) # Copy
                     self.MutateVector(NewPopulation[child_idx])
                     NewRS[child_idx] = rand_parent_rs # Assign RS (may be duplicated)
                     child_idx += 1
                     continue


                # Check for crossover
                if self.rs.ran1() < self.CrossProb and parent_idx + 1 < len(parent_pool):
                    # Perform crossover
                    p1_data = parent_pool[parent_idx]
                    p2_data = parent_pool[parent_idx + 1]
                    # Create children as copies initially
                    child1_vec = TVector(p1_data['vector'].LowerBound(), p1_data['vector'].UpperBound(), p1_data['vector'].to_list())
                    child2_vec = TVector(p2_data['vector'].LowerBound(), p2_data['vector'].UpperBound(), p2_data['vector'].to_list())

                    if self.CrossMode == TCrossoverMode.UNIFORM:
                        self.UniformCrossover(child1_vec, child2_vec)
                    elif self.CrossMode == TCrossoverMode.TWO_POINT:
                        self.TwoPointCrossover(child1_vec, child2_vec)
                    else:
                        raise ValueError("Invalid crossover mode")

                    # Check if children are identical to parents after crossover
                    if self._EqualVector(child1_vec, p1_data['vector']):
                        self.MutateVector(child1_vec)
                    if self._EqualVector(child2_vec, p2_data['vector']):
                        self.MutateVector(child2_vec)

                    # Add children to new population
                    NewPopulation[child_idx] = child1_vec
                    NewRS[child_idx] = p1_data['rs'] # Assign parent1's RS
                    child_idx += 1
                    if child_idx <= pop_size:
                        NewPopulation[child_idx] = child2_vec
                        NewRS[child_idx] = p2_data['rs'] # Assign parent2's RS
                        child_idx += 1

                    parent_idx += 2
                else:
                    # Perform mutation only
                    p_data = parent_pool[parent_idx]
                    child_vec = TVector(p_data['vector'].LowerBound(), p_data['vector'].UpperBound(), p_data['vector'].to_list()) # Copy
                    self.MutateVector(child_vec) # Mutate the copy
                    NewPopulation[child_idx] = child_vec
                    NewRS[child_idx] = p_data['rs'] # Assign parent's RS
                    child_idx += 1
                    parent_idx += 1

        # Replace old population with the new one
        self.Population = NewPopulation
        self.RandomStates = NewRS

        # Evaluate the new population
        # If ReEvalFlag is true, evaluate everyone (including elites)
        # Otherwise, only evaluate the newly created children
        eval_start_index = 1 if self.ReEvalFlag else ElitePop + 1
        if eval_start_index <= pop_size:
             logger.info(f"Evaluating new individuals from index {eval_start_index}...")
             self.EvaluatePopulation(start=eval_start_index)
        else:
             logger.info("No new individuals to evaluate.")

        # If elites were not re-evaluated, their performance is already known
        if not self.ReEvalFlag:
             for i in range(1, ElitePop + 1):
                  self.Perf[i] = NewPerf[i] # Restore known elite performance


    # --- Genetic Operators ---
    def MutateVector(self, v: TVector[float]):
        """Applies Gaussian mutation to a vector."""
        if self._vectorSize == 0: return

        # GaussianRandom(mean, stddev) -> MutationVar is variance in C++
        std_dev = math.sqrt(self.MutationVar)
        magnitude = self.rs.GaussianRandom(0.0, std_dev)

        # Random unit vector (simplified: generate components and normalize)
        temp_vector = [self.rs.GaussianRandom(0.0, 1.0) for _ in range(self._vectorSize)]
        norm = math.sqrt(sum(x*x for x in temp_vector))
        if norm == 0: norm = 1.0 # Avoid division by zero

        unit_vector = [x / norm for x in temp_vector]

        # Apply mutation
        for i in range(1, self._vectorSize + 1):
            v[i] += magnitude * unit_vector[i-1] # Adjust index for unit_vector
            # Apply constraints
            if self.ConstraintVector[i]: # If constraint is active (1)
                 v[i] = clip(v[i], MIN_SEARCH_VALUE, MAX_SEARCH_VALUE)
            # else: value can go out of bounds (matches C++ comment)

    def UniformCrossover(self, v1: TVector[float], v2: TVector[float]):
        """Performs modular uniform crossover."""
        num_modules = self.crossPoints.Size()
        if num_modules < 1: return # Need at least one module defined by points

        # Iterate through modules defined by crossPoints
        for module_idx in range(1, num_modules + 1):
            start_gene_idx = self.crossPoints[module_idx]
            # End index is start of next point - 1, or vector end if last module
            end_gene_idx = self.crossPoints[module_idx + 1] - 1 if module_idx < num_modules else self._vectorSize

            if self.rs.ran1() < 0.5: # Swap this module
                for j in range(start_gene_idx, end_gene_idx + 1):
                    # Swap elements using Python's tuple unpacking
                    v1[j], v2[j] = v2[j], v1[j]


    def TwoPointCrossover(self, v1: TVector[float], v2: TVector[float]):
        """Performs modular two-point crossover."""
        num_modules = self.crossPoints.Size()
        if num_modules < 2: return # Need at least two points to define a segment

        # Choose two distinct module indices (1-based)
        i1 = self.rs.UniformRandomInteger(1, num_modules)
        i2 = i1
        while i2 == i1:
            i2 = self.rs.UniformRandomInteger(1, num_modules)

        # Ensure i1 < i2
        if i1 > i2:
            i1, i2 = i2, i1

        # Get the actual gene indices corresponding to the start of these modules
        start_gene_idx = self.crossPoints[i1]
        # End gene index is the start of module i2, minus 1
        end_gene_idx = self.crossPoints[i2] - 1

        # Swap the genes between these points
        for i in range(start_gene_idx, end_gene_idx + 1):
             v1[i], v2[i] = v2[i], v1[i]


    # --- Statistics and Display ---
    def UpdatePopulationStatistics(self):
        """Updates performance stats: Min, Max, Avg, Var, Best."""
        pop_size = self.PopulationSize()
        if pop_size == 0:
            self.MinPerf = self.MaxPerf = self.AvgPerf = self.PerfVar = 0.0
            self.UpdateBestFlag = False # Should not update if no population
            return

        best_index_in_current_pop = 1
        current_max_perf = -float('inf')
        current_min_perf = float('inf')
        total_perf = 0.0

        for i in range(1, pop_size + 1):
            perf = self.Perf[i]
            if perf > current_max_perf:
                current_max_perf = perf
                best_index_in_current_pop = i
            if perf < current_min_perf:
                current_min_perf = perf
            total_perf += perf

        self.MaxPerf = current_max_perf
        self.MinPerf = current_min_perf
        self.AvgPerf = total_perf / pop_size if pop_size > 0 else 0.0

        # Clip AvgPerf (as in C++)
        self.AvgPerf = max(self.MinPerf, min(self.MaxPerf, self.AvgPerf))

        # Variance
        if pop_size > 1:
            sum_sq_diff = sum((self.Perf[i] - self.AvgPerf)**2 for i in range(1, pop_size + 1))
            self.PerfVar = sum_sq_diff / (pop_size - 1)
        else:
            self.PerfVar = 0.0

        # Update overall best if improved or ReEvalFlag is set
        # Note: C++ updates if MaxPerf > BestPerf. Use >= for robustness?
        # Stick to > to match C++.
        if self.MaxPerf > self.BestPerf or self.ReEvalFlag:
             # Check if the actual best vector changed, even if perf is same (due to ReEvalFlag)
             potential_new_best = self.Population[best_index_in_current_pop]
             if self.MaxPerf > self.BestPerf or not self._EqualVector(self.bestVector, potential_new_best):
                  self.UpdateBestFlag = True
                  self.BestPerf = self.MaxPerf
                  # Deep copy the best vector
                  self.bestVector = TVector(potential_new_best.LowerBound(), potential_new_best.UpperBound(), potential_new_best.to_list())


    def DisplayPopulationStatistics(self):
        """Calls the registered display function or prints defaults."""
        if self.PopulationStatisticsDisplayFunction is not None:
            self.PopulationStatisticsDisplayFunction(self.Gen, self.BestPerf, self.AvgPerf, self.PerfVar)
        else:
            # Default print format
            logger.info(f"Gen {self.Gen}: Best={self.BestPerf:.6f}, Avg={self.AvgPerf:.6f}, Var={self.PerfVar:.6f}")

    def SearchTerminated(self) -> bool:
        """Checks if the search termination conditions are met."""
        term_gen = self.Gen >= self.MaxGens
        term_func = False
        if self.SearchTerminationFunction is not None:
            term_func = self.SearchTerminationFunction(self.Gen, self.BestPerf, self.AvgPerf, self.PerfVar)
        return term_gen or term_func

    def DisplaySearchResults(self):
        """Calls the registered results display function."""
        if self.SearchResultsDisplayFunction is not None:
            self.SearchResultsDisplayFunction(self)
        else:
            # Default: print best performance and vector
            logger.info("\nSearch Results:")
            logger.info(f"Best Performance: {self.BestPerf}")
            # print(f"Best Vector: {self.bestVector}") # Can be long