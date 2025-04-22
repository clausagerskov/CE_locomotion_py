import math
import sys
import pickle
from typing import List, Tuple, Any, TextIO

# Assuming vector_matrix.py is in the same directory or accessible
from vector_matrix import TVector

# --- Global defines from random.h ---
IA = 16807
IM = 2147483647
AM = (1.0 / IM)
IQ = 127773
IR = 2836
NTAB = 32
NDIV = (1 + (IM - 1) // NTAB) # Use // for integer division
EPS = 1.2e-7
RNMX = (1.0 - EPS)


class RandomState:
    """
    Python port of the C++ RandomState class using the ran1 algorithm
    from Numerical Recipes for pseudo-random number generation.
    Intended for reproducibility with the original C++ simulation.
    """

    def __init__(self, seed_val: int = 0):
        """Initializes the random state."""
        # State variables for ran1
        self.idum: int = 0
        self.iy: int = 0
        self.iv: List[int] = [0] * NTAB # Internal table for ran1

        # State variable for Gaussian generation (Box-Muller)
        self.gaussian_flag: bool = False
        self.gX1: float = 0.0
        self.gX2: float = 0.0

        # Store the initial seed
        self.seed: int = 0 # Will be set by SetRandomSeed

        self.SetRandomSeed(seed_val)
        self.gaussian_flag = False # Ensure flag is reset on init

    def SetRandomSeed(self, s: int):
        """
        Seeds the ran1 random number generator.
        Adapted from Numerical Recipes.
        """
        self.seed = s # Store the seed used
        self.idum = int(s) # Ensure it's an integer

        # Original C++ code check: if (idum == 0) idum = 1;
        # Let's replicate, although 0 seed is often discouraged.
        if self.idum == 0:
             self.idum = 1
             print("Warning: Random seed 0 provided, using 1 instead.", file=sys.stderr)
        elif self.idum < 0:
             # ran1 expects positive seed, NR uses negative to re-init.
             # Let's just take absolute value for simplicity here.
             self.idum = abs(self.idum)
             print(f"Warning: Negative seed {s} provided, using absolute value {self.idum}.", file=sys.stderr)


        # Initialization loop from ran1
        for j in range(NTAB + 7, -1, -1): # Loop from NTAB+7 down to 0
            k = self.idum // IQ # Integer division
            self.idum = IA * (self.idum - k * IQ) - IR * k
            if self.idum < 0:
                self.idum += IM
            if j < NTAB:
                self.iv[j] = self.idum

        self.iy = self.iv[0]
        self.gaussian_flag = False # Reset Gaussian flag when seeding

    def GetRandomSeed(self) -> int:
        """Returns the initial seed used to initialize the generator."""
        return self.seed

    def ran1(self) -> float:
        """
        Generates a uniform random deviate between 0.0 and 1.0 (exclusive of 1.0).
        Based on ran1 from Numerical Recipes.
        """
        # Ensure generator is initialized if called before SetRandomSeed
        # The check `if not self.iy:` from C++ isn't directly portable as iy can be 0.
        # Rely on __init__ having called SetRandomSeed.

        if self.idum <= 0 or self.iy < 0:
             raise ValueError(f"Internal PRNG state invalid: idum={self.idum}, iy={self.iy}. Must call SetRandomSeed with positive seed.")

        k = self.idum // IQ
        self.idum = IA * (self.idum - k * IQ) - IR * k
        if self.idum < 0:
            self.idum += IM

        j = self.iy // NDIV # Index into the shuffle table
        if not (0 <= j < NTAB):
            # This indicates a problem with the state or constants
             raise IndexError(f"Shuffle table index j={j} out of range [0, {NTAB-1}] (iy={self.iy}, NDIV={NDIV})")

        self.iy = self.iv[j]
        self.iv[j] = self.idum

        temp = AM * self.iy # Scale to [0, 1)
        # Ensure result is capped slightly below 1.0 as per RNMX
        return min(RNMX, temp)


    def UniformRandom(self, min_val: float, max_val: float) -> float:
        """
        Returns a uniformly-distributed random double between min_val and max_val.
        [min_val, max_val) - Interval likely includes min, excludes max.
        """
        return (max_val - min_val) * self.ran1() + min_val

    def UniformRandomInteger(self, min_val: int, max_val: int) -> int:
        """
        Returns a uniformly-distributed random integer between min_val and max_val (inclusive).
        """
        if min_val > max_val:
             raise ValueError(f"Minimum value {min_val} cannot be greater than maximum value {max_val}")
        # Formula from C++ uses floor(0.5 + UniformRandom(min-0.5, max+0.5))
        # This aims for equal probability bins centered on integers.
        # Python's random.randint(a, b) includes both endpoints.
        # Let's use the direct approach:
        return min_val + int(self.ran1() * (max_val - min_val + 1))
        # Alternative C++ logic port:
        # return math.floor(0.5 + self.UniformRandom(float(min_val) - 0.5, float(max_val) + 0.5))


    def GenerateNormals(self):
        """
        Generates two normally-distributed random variables gX1 and gX2
        using the Box-Muller transform. Should only be called internally.
        """
        while True:
            v1 = self.UniformRandom(-1.0, 1.0)
            v2 = self.UniformRandom(-1.0, 1.0)
            s = v1 * v1 + v2 * v2
            if 0.0 < s < 1.0: # Check s is within (0, 1)
                break

        d = math.sqrt(-2.0 * math.log(s) / s)
        self.gX1 = v1 * d
        self.gX2 = v2 * d
        self.gaussian_flag = True # Set flag indicating normals are ready


    def GaussianRandom(self, mean: float, variance: float) -> float:
        """
        Generates a Gaussian random variable with specified mean and variance.
        Uses the Box-Muller method, generating pairs of deviates.
        Args:
            mean: The mean of the Gaussian distribution.
            variance: The variance (sigma^2) of the distribution.
        """
        if variance < 0:
            raise ValueError("Variance cannot be negative")

        if not self.gaussian_flag:
            # If flag is not set, generate a new pair of normals
            self.GenerateNormals()
            # Return the first normal from the pair, scaled and shifted
            std_dev = math.sqrt(variance)
            # Keep flag set, as gX2 is still available
            return std_dev * self.gX1 + mean
        else:
            # If flag is set, use the second normal from the previous pair
            std_dev = math.sqrt(variance)
            self.gaussian_flag = False # Reset the flag
            return std_dev * self.gX2 + mean


    def RandomUnitVector(self, v: TVector[float]):
        """
        Generates a random unit vector by filling v with Gaussian deviates
        and normalizing. Modifies v in-place.
        """
        r_sq = 0.0
        if v.Size() == 0: return # Cannot normalize zero-size vector

        for i in range(v.LowerBound(), v.UpperBound() + 1):
            val = self.GaussianRandom(0.0, 1.0) # Mean 0, Variance 1
            v[i] = val
            r_sq += val * val

        if r_sq == 0.0:
             # Handle rare case where all elements are zero
             # Set to an arbitrary unit vector, e.g., [1, 0, 0, ...]
             v[v.LowerBound()] = 1.0
             for i in range(v.LowerBound() + 1, v.UpperBound() + 1):
                  v[i] = 0.0
             return

        r = math.sqrt(r_sq)
        inv_r = 1.0 / r
        for i in range(v.LowerBound(), v.UpperBound() + 1):
            v[i] *= inv_r


    def ProbabilisticChoice(self, prob: float) -> int:
        """Returns 1 with probability prob, otherwise returns 0."""
        if not (0.0 <= prob <= 1.0):
             raise ValueError("Probability must be between 0.0 and 1.0")
        return 1 if self.UniformRandom(0.0, 1.0) < prob else 0


    # --- Input/Output and State Management ---

    def getstate(self) -> Tuple:
        """Returns the internal state of the generator as a tuple."""
        # Make sure to return copies, especially for the list 'iv'
        return (
            self.seed,
            self.idum,
            self.iy,
            list(self.iv), # Return a copy of the list
            self.gaussian_flag,
            self.gX1,
            self.gX2
        )

    def setstate(self, state: Tuple):
        """Sets the internal state of the generator from a tuple."""
        if len(state) != 7:
            raise ValueError("Invalid state tuple length")
        try:
            s, idm, iy_val, iv_list, g_flag, gx1_val, gx2_val = state
            self.seed = int(s)
            self.idum = int(idm)
            self.iy = int(iy_val)
            if len(iv_list) != NTAB:
                raise ValueError(f"Invalid length for iv list in state (expected {NTAB})")
            self.iv = list(iv_list) # Set as a copy
            self.gaussian_flag = bool(g_flag)
            self.gX1 = float(gx1_val)
            self.gX2 = float(gx2_val)
        except (TypeError, ValueError) as e:
             raise ValueError(f"Invalid data type or value in state tuple: {e}")


    def WriteRandomState(self, ofs: TextIO):
        """Writes the generator state to a text stream (space-separated)."""
        state_tuple = self.getstate()
        # Manually format like C++ version
        ofs.write(f"{state_tuple[0]} {state_tuple[1]} {state_tuple[2]} ") # seed idum iy
        ofs.write(f"{int(state_tuple[4])} {state_tuple[5]} {state_tuple[6]} ") # gaussian_flag gX1 gX2 (write bool as 0/1)
        ofs.write(" ".join(map(str, state_tuple[3]))) # The iv list
        ofs.write("\n") # Match C++ endl

    def ReadRandomState(self, ifs: TextIO):
        """Reads the generator state from a text stream."""
        line = ifs.readline()
        if not line:
            raise EOFError("Could not read random state line from stream")
        parts = line.split()
        # Expected number of parts: 3 + 3 + NTAB
        expected_parts = 6 + NTAB
        if len(parts) < expected_parts:
            raise ValueError(f"Incorrect number of items on random state line (got {len(parts)}, expected {expected_parts})")

        try:
            s = int(parts[0])
            idm = int(parts[1])
            iy_val = int(parts[2])
            g_flag = bool(int(parts[3]))
            gx1_val = float(parts[4])
            gx2_val = float(parts[5])
            iv_list = [int(p) for p in parts[6:6+NTAB]]
            self.setstate((s, idm, iy_val, iv_list, g_flag, gx1_val, gx2_val))
        except (ValueError, IndexError) as e:
            raise ValueError(f"Error parsing random state from text stream: {e}")

    # --- Binary I/O (Optional, for exact C++ file compatibility) ---
    # Note: Requires careful handling of integer sizes (e.g., long vs int)
    # Using pickle via getstate/setstate is generally easier in Python.

    def BinaryWriteRandomState(self, bofs): # bofs is a binary file stream opened for writing
        """Writes state in a binary format potentially compatible with C++."""
        # Using pickle is easier for pure Python save/load
        pickle.dump(self.getstate(), bofs)

    def BinaryReadRandomState(self, bifs): # bifs is a binary file stream opened for reading
        """Reads state from a binary format."""
        state_tuple = pickle.load(bifs)
        self.setstate(state_tuple)