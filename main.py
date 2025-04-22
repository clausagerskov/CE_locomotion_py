import math
import time
import sys
import os
import contextlib # For redirecting stdout

# Import necessary components from our other files
from vector_matrix import TVector
from t_search import (
    TSearch, TSelectionMode, TReproductionMode, TCrossoverMode,
    MapSearchParameter, clip
)
from loguru import logger
from worm import Worm
from random_state_legacy import RandomState
from worm_body import WormBody, InitializeBodyConstants

# --- Constants ---
PRINTTOFILE = True # Set to False to print to console
SKIP_STEPS = 10

# Integration parameters (match C++)
DURATION = 1.0 # Use float for time
TRANSIENT = 0.5
STEP_SIZE = 0.005
N_CURVS = 23 # Seems unused in provided C++ main, but kept for context

# Fitness trajectory (match C++)
AVG_SPEED = 0.0001 # 0.00022 in C++ comment, but 0.0001 used in code
BBC_FIT = AVG_SPEED * DURATION

# Genotype -> Phenotype Mapping Ranges (match C++)
BIAS_RANGE = 16.0
SC_RANGE = 16.0 # Self-Connection Range
CS_RANGE = 16.0 # Chemical Synapse Range
ES_RANGE = 2.0  # Electrical Synapse Range (Gap Junction)
SR_MAX = 200.0  # Stretch Receptor Max Gain
NMJ_MAX = 0.8
NMJ_MIN = 0.0

# Indices for Phenotype vector (1-based like C++)
# These MUST match the order in GenPhenMapping
SR_A = 1
SR_B = 2
# Bias: 3, 4, 5
# Self connections: 6, 7, 8
# DA, DB (Excitatory): 9, 10
# VD (Inhibitory): 11, 12
# Interunit Gap Junctions: 13, 14
# Excitatory NMJ: 15, 16
# Inhibitory NMJ: 17

# Size of genotype (phenotype has same size after mapping)
VECT_SIZE = 17

# --- Genotype-Phenotype Mapping ---
def GenPhenMapping(gen: TVector[float], phen: TVector[float]):
    """Maps genotype from [-1, 1] to phenotype parameters."""
    if gen.Size() != VECT_SIZE or phen.Size() != VECT_SIZE:
        raise ValueError("Genotype/Phenotype vector size mismatch")

    # Parameters for the Stretch Receptors
    phen[SR_A] = MapSearchParameter(gen[SR_A], 0.0, SR_MAX)
    phen[SR_B] = MapSearchParameter(gen[SR_B], 0.0, SR_MAX)

    k = 3 # Start index for biases (1-based)
    # Bias (3 parameters)
    for _ in range(3):
        phen[k] = MapSearchParameter(gen[k], -BIAS_RANGE, BIAS_RANGE)
        k += 1
    # Self connections (3 parameters)
    for _ in range(3):
        phen[k] = MapSearchParameter(gen[k], -SC_RANGE, SC_RANGE)
        k += 1
    # DA, DB Chemical synapses (excitatory) (2 parameters)
    for _ in range(2):
        phen[k] = MapSearchParameter(gen[k], 0.0, CS_RANGE)
        k += 1
    # VD Chemical synapses (Inhibitory) (2 parameters)
    for _ in range(2):
        phen[k] = MapSearchParameter(gen[k], -CS_RANGE, 0.0)
        k += 1
    # Interunits Gap junctions (2 parameters)
    for _ in range(2):
        phen[k] = MapSearchParameter(gen[k], 0.0, ES_RANGE)
        k += 1
    # Excitatory NMJ Weight (2 parameters)
    for _ in range(2):
        phen[k] = MapSearchParameter(gen[k], NMJ_MIN, NMJ_MAX)
        k += 1
    # Inhibitory NMJ Weight (1 parameter)
    for _ in range(1):
        phen[k] = MapSearchParameter(gen[k], -NMJ_MAX, -NMJ_MIN)
        k += 1

    if k - 1 != VECT_SIZE:
         logger.info(f"Warning: GenPhenMapping processed {k-1} parameters, expected {VECT_SIZE}", file=sys.stderr)


# --- Fitness Function ---
def Evaluation(v_genotype: TVector[float], rs: RandomState, direction: int) -> float:
    """Evaluates a genotype for forward (1) or backward (-1) locomotion."""
    distancetravelled = 0.0
    xt = 0.0
    yt = 0.0

    # Genotype-Phenotype Mapping
    phenotype = TVector[float](1, VECT_SIZE)
    GenPhenMapping(v_genotype, phenotype)

    # Create Worm instance
    w = Worm(phenotype, 1) # Second arg seems unused in C++ Worm constructor call
    w.InitializeState(rs)

    # Set Command Interneuron Activation based on direction
    if direction == 1: # Forward
        w.AVA_output = w.AVA_inact # Typically 0 or low
        w.AVB_output = w.AVB_act   # Typically 1 or high
    else: # Backward (-1)
        w.AVA_output = w.AVA_act   # Typically 1 or high
        w.AVB_output = w.AVB_inact # Typically 0 or low

    # --- Simulation ---
    # Transient phase
    current_time = 0.0
    while current_time <= TRANSIENT:
        w.Step(STEP_SIZE, 1) # Second arg seems unused
        current_time += STEP_SIZE

    xt = w.CoMx()
    yt = w.CoMy()
    oxt = xt # Original position X
    oyt = yt # Original position Y

    # Run phase
    current_time = 0.0
    step_counter_eval = 0
    while current_time <= DURATION:
        w.Step(STEP_SIZE, 1)
        step_counter_eval += 1
        if step_counter_eval % 500 == 0:
            logger.info(f"    Eval step {step_counter_eval}/{int(DURATION/STEP_SIZE)} (t={current_time:.2f})")
        # Current and past centroid position
        xtp = xt
        ytp = yt
        xt = w.CoMx()
        yt = w.CoMy()

        # Integration error check (like C++)
        if math.isnan(xt) or math.isnan(yt):
            logger.info("Warning: NaN detected in worm position during evaluation.", file=sys.stderr)
            return 0.0 # Penalize NaN

        # Calculate displacement in this step
        delta_x = xt - xtp
        delta_y = yt - ytp
        step_distance = math.sqrt(delta_x**2 + delta_y**2)

        # Check for excessive movement (instability)
        if step_distance > 10 * AVG_SPEED * STEP_SIZE:
             logger.info("Warning: Excessive movement detected during evaluation.", file=sys.stderr)
             return 0.0 # Penalize instability

        # Velocity Fitness Calculation (from C++)
        body_orientation = w.Orientation()
        if step_distance > 1e-9: # Avoid atan2(0,0)
            movement_orientation = math.atan2(delta_y, delta_x)
        else:
            movement_orientation = body_orientation # Assume aligned if no movement

        angle_diff = movement_orientation - body_orientation
        # Normalize angle difference to [-pi, pi]
        angle_diff = (angle_diff + math.pi) % (2 * math.pi) - math.pi

        # temp determines if movement counts positively or negatively based on direction
        cos_anglediff = math.cos(angle_diff)
        if direction == 1: # Forward
            temp = 1.0 if cos_anglediff > 0.0 else -1.0
        else: # Backward (-1)
            temp = -1.0 if cos_anglediff > 0.0 else 1.0

        distancetravelled += temp * step_distance
        current_time += STEP_SIZE

    # --- Fitness Calculation ---
    fxt = w.CoMx() # Final position X
    fyt = w.CoMy() # Final position Y
    # Final displacement distance
    # distance = math.sqrt((oxt - fxt)**2 + (oyt - fyt)**2)

    # Fitness A (displacement-based, commented out in C++)
    # fitA = 1.0 - (abs(BBC_FIT - distance) / BBC_FIT) if BBC_FIT > 0 else 0.0
    # fitA = max(0.0, fitA)

    # Fitness B (path-based, used in C++)
    fitB = 1.0 - (abs(BBC_FIT - distancetravelled) / BBC_FIT) if BBC_FIT > 0 else 0.0
    fitB = max(0.0, fitB)

    return fitB


def EvaluationFunction(v: TVector[float], rs: RandomState) -> float:
    """
    The main evaluation function called by TSearch.
    Handles modifications to the genotype for directed evaluation.
    """
    # Store original SR gene values
    sra_gene = v[SR_A]
    srb_gene = v[SR_B]

    # --- Forward Evaluation ---
    # Modify genotype temporarily: Turn OFF SR_A for forward eval
    v[SR_A] = -1.0 # MapSearchParameter(-1.0, 0.0, SR_MAX) -> 0.0
    v[SR_B] = srb_gene # Keep original SR_B gene value
    fitnessForward = Evaluation(v, rs, 1)

    # --- Backward Evaluation (Commented out as in C++) ---
    # v[SR_A] = sra_gene # Restore SR_A gene
    # v[SR_B] = -1.0     # Turn OFF SR_B for backward eval
    # fitnessBackward = Evaluation(v, rs, -1)

    # Restore original genotype
    v[SR_A] = sra_gene
    v[SR_B] = srb_gene

    # Return combined fitness or just forward/backward
    # return (fitnessForward + fitnessBackward) / 2.0
    return fitnessForward
    # return fitnessBackward


# --- Plotting / Saving Traces ---
def save_traces(v_genotype: TVector[float], rs: RandomState):
    """Simulates the best genotype and saves traces to files."""
    logger.info("\nGenerating traces for the best individual...")

    # Ensure output directories exist if needed (not strictly necessary for current files)
    # os.makedirs("output_data", exist_ok=True)

    #try:
    with open("curv.dat", "w") as curvfile, \
            open("body.dat", "w") as bodyfile, \
            open("act.dat", "w") as actfile, \
            open("phenotype.dat", "w") as phenfile:

        # Header for phenotype file
        phenfile.write(f"# Phenotype Parameters for Best Genotype\n")
        phenfile.write(f"# Genotype: {v_genotype}\n")

        # Genotype-Phenotype Mapping
        phenotype = TVector[float](1, VECT_SIZE)
        GenPhenMapping(v_genotype, phenotype)

        # Get mapped SR gains
        sra_pheno = phenotype[SR_A]
        srb_pheno = phenotype[SR_B]

        # Create Worm instance
        w = Worm(phenotype, 1)
        w.DumpParams(phenfile) # Save phenotype details

        w.InitializeState(rs)

        total_sim_time = 0.0

        # --- Phase 1: Forward Movement ---
        logger.info("Simulating Phase 1: Forward...")
        w.sr.SR_A_gain = 0.0 # SR_A off
        w.sr.SR_B_gain = srb_pheno # SR_B on (using mapped value)
        w.AVA_output = w.AVA_inact # Forward command
        w.AVB_output = w.AVB_act

        current_time = 0.0
        sim_end_time = TRANSIENT + DURATION
        while current_time <= sim_end_time:
            logger.info("Step")
            w.Step(STEP_SIZE, 1)
            logger.info("DumpyBodyState")
            w.DumpBodyState(bodyfile, SKIP_STEPS)
            logger.info("DumpCurvature")
            w.DumpCurvature(curvfile, SKIP_STEPS)
            logger.info("DumpActState")
            w.DumpActState(actfile, SKIP_STEPS)
            current_time += STEP_SIZE
        total_sim_time += current_time

        # --- Phase 2: SR Off ---
        logger.info("Simulating Phase 2: SR Off...")
        w.sr.SR_A_gain = 0.0
        w.sr.SR_B_gain = 0.0
        # Keep command neurons as they were? C++ code doesn't change them here.

        current_time = 0.0
        sim_end_time = 12.0
        while current_time <= sim_end_time:
                w.Step(STEP_SIZE, 1)
                w.DumpBodyState(bodyfile, SKIP_STEPS)
                w.DumpCurvature(curvfile, SKIP_STEPS)
                w.DumpActState(actfile, SKIP_STEPS)
                current_time += STEP_SIZE
        total_sim_time += current_time

        # --- Phase 3: Backward Movement ---
        logger.info("Simulating Phase 3: Backward...")
        w.sr.SR_A_gain = sra_pheno # SR_A on
        w.sr.SR_B_gain = 0.0      # SR_B off
        w.AVA_output = w.AVA_act  # Backward command
        w.AVB_output = w.AVB_inact

        current_time = 0.0
        sim_end_time = 20.0
        while current_time <= sim_end_time:
                w.Step(STEP_SIZE, 1)
                w.DumpBodyState(bodyfile, SKIP_STEPS)
                w.DumpCurvature(curvfile, SKIP_STEPS)
                w.DumpActState(actfile, SKIP_STEPS)
                current_time += STEP_SIZE
        total_sim_time += current_time

        # --- Phase 4: SR Off ---
        logger.info("Simulating Phase 4: SR Off...")
        w.sr.SR_A_gain = 0.0
        w.sr.SR_B_gain = 0.0
        # Keep command neurons as they were.

        current_time = 0.0
        sim_end_time = 12.0
        while current_time <= sim_end_time:
                w.Step(STEP_SIZE, 1)
                w.DumpBodyState(bodyfile, SKIP_STEPS)
                w.DumpCurvature(curvfile, SKIP_STEPS)
                w.DumpActState(actfile, SKIP_STEPS)
                current_time += STEP_SIZE
        total_sim_time += current_time

        logger.info(f"Trace simulation complete. Total time: {total_sim_time:.2f}s")

    #except IOError as e:
        #logger.info(f"Error opening or writing trace files: {e}", file=sys.stderr)
    #except Exception as e:
        #logger.info(f"An error occurred during save_traces: {e}", file=sys.stderr)

    return 0 # Match C++ return type

# --- Display Functions ---
def EvolutionaryRunDisplay(Generation: int, BestPerf: float, AvgPerf: float, PerfVar: float):
    """Callback function to display progress during evolution."""
    # Output format matches C++ default for fitness.dat
    logger.info(f"{BestPerf:.10f} {AvgPerf:.10f} {PerfVar:.10f}")
    sys.stdout.flush() # Ensure output is written immediately, especially if redirected

def ResultsDisplay(s: TSearch):
    """Callback function to display results at the end of the search."""
    bestVector = s.BestIndividual()
    logger.info("\nEvolution Finished.")
    logger.info(f"Best Performance: {s.BestPerformance()}")
    # print(f"Best Genotype Vector: {bestVector}") # Can be very long

    # Save the best genotype vector to a file
    try:
        with open("best.gen.dat", "w") as f:
            # Use TVector's __str__ method which prints space-separated values
            f.write(f"{bestVector}\n")
        logger.info("Best genotype saved to best.gen.dat")
    except IOError as e:
        logger.info(f"Error writing best genotype file: {e}", file=sys.stderr)


# --- Main Program ---
if __name__ == "__main__":
    logger.info(f"Evolving Locomotion - Python Port")
    logger.info(f"Vector Size: {VECT_SIZE}")

    # --- Random Seed Setup ---
    # Use time and optional command-line argument for seed
    random_seed = int(time.time())
    if len(sys.argv) == 2:
        try:
            random_seed += int(sys.argv[1])
            logger.info(f"Using provided offset {sys.argv[1]} for random seed.")
        except ValueError:
            logger.info(f"Warning: Invalid command-line argument '{sys.argv[1]}'. Using time-based seed only.")

    logger.info(f"Setting master random seed to: {random_seed}")
    # Save the seed to a file
    try:
        with open("seed.dat", "w") as seedfile:
            seedfile.write(str(random_seed) + "\n")
    except IOError as e:
        logger.info(f"Warning: Could not write seed file: {e}", file=sys.stderr)

    # --- Initialize Search ---
    s = TSearch(VECT_SIZE)
    s.SetRandomSeed(random_seed) # Seed master and individual generators

    # --- Configure Search Parameters (match C++) ---
    s.SetPopulationStatisticsDisplayFunction(EvolutionaryRunDisplay)
    s.SetSearchResultsDisplayFunction(ResultsDisplay)
    s.SetSelectionMode(TSelectionMode.RANK_BASED)
    s.SetReproductionMode(TReproductionMode.HILL_CLIMBING)
    s.SetPopulationSize(1) # Default was 1, C++ uses 96
    s.SetMaxGenerations(2) # C++ uses 1000, use 10 for quick test
    s.SetMutationVariance(0.05) # C++ uses 0.1, changed to 0.05 in comments? Use 0.05
    s.SetCrossoverProbability(0.5)
    s.SetCrossoverMode(TCrossoverMode.UNIFORM) # C++ uses TWO_POINT, changed to UNIFORM? Use UNIFORM
    s.SetMaxExpectedOffspring(1.1)
    s.SetElitistFraction(0.02)
    s.SetSearchConstraint(1) # 1 means constrain to [-1, 1]
    s.SetReEvaluationFlag(False) # C++ uses 0

    # --- File Redirection ---
    output_target = None
    if PRINTTOFILE:
        try:
            output_target = open("fitness.dat", "w")
            logger.info("Redirecting evolution statistics output to fitness.dat")
        except IOError as e:
            logger.info(f"Warning: Could not open fitness.dat for writing: {e}. Outputting to console.", file=sys.stderr)
            output_target = sys.stdout # Fallback to console
    else:
        output_target = sys.stdout # Print to console

    # --- Run Evolution ---
    # Initialize Body Constants
    logger.info("Initializing Body Constants...")
    InitializeBodyConstants()

    s.SetEvaluationFunction(EvaluationFunction)

    logger.info("\nStarting evolutionary search...")
    start_time = time.time()

    # Use context manager for potential file redirection
    with contextlib.redirect_stdout(output_target):
        s.ExecuteSearch() # Run the evolution

    end_time = time.time()
    logger.info(f"\nEvolution took {end_time - start_time:.2f} seconds.")

    # Close the fitness file if it was opened
    if PRINTTOFILE and output_target is not None and output_target is not sys.stdout:
        output_target.close()
        logger.info("Closed fitness.dat")


    # --- Post-Evolution: Save Traces ---
    # Re-seed a random generator for trace generation consistency if needed,
    # or just use a new one. C++ used a new time-based seed.
    trace_rs = RandomState()
    trace_seed = int(time.time())
    trace_rs.seed(trace_seed)
    logger.info(f"\nUsing new random seed {trace_seed} for trace generation.")

    #try:
    # Load the best genotype from the file saved by ResultsDisplay
    with open("best.gen.dat", "r") as bestfile:
        line = bestfile.readline()
        best_gen_data = [float(x) for x in line.strip().split()]
        best_genotype = TVector[float](1, VECT_SIZE, best_gen_data)

    # Generate and save traces
    save_traces(best_genotype, trace_rs)
    logger.info("Trace generation finished.")

    #except FileNotFoundError:
    logger.info("Error: best.gen.dat not found. Cannot generate traces.", file=sys.stderr)
    #except Exception as e:
    logger.info(f"An error occurred during trace generation setup: {e}", file=sys.stderr)

    logger.info("\nProgram finished.")
