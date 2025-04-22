import math
import random
from typing import TextIO

# Assuming vector_matrix.py is in the same directory or accessible
from vector_matrix import TVector
from stretch_receptor import StretchReceptor

from worm_body import WormBody, N_SEGMENTS, N_RODS
from muscles import Muscles
from nervous_system import NervousSystem
from loguru import logger
# --- Constants (mirrored from C++) ---
# PI defined using math.pi later where needed

# Parameters
N_MUSCLES = 24           # Number of muscles alongside the body
N_UNITS = 10             # Number of neural units in VNC
N_NEURONSPERUNIT = 6     # Number of neurons in a VNC neural unit (6 neurons)
# H_NEURONSPERUNIT = 3   # Half for DV symmetry (unused in provided Worm.cpp/h)

N_STRETCHREC = N_UNITS   # Number of stretch receptors
T_MUSCLE = 0.1           # Muscle time constant

# N_MUSCLEPERNU = 4 # C++ comment, not used directly in Worm.cpp code logic shown

# Motoneuron name conventions (1-based indices within a unit)
DA = 1
DB = 2
DD = 3
VD = 4
VA = 5
VB = 6

# Body segment name conventions (1-based indices)
HEAD = 1
TAIL = N_SEGMENTS # Depends on the actual N_SEGMENTS

class PlaceholderNervousSystem:
    def __init__(self):
        self.num_neurons = N_UNITS * N_NEURONSPERUNIT
        self.outputs = TVector(1, self.num_neurons, [0.0] * self.num_neurons)
        self.states = TVector(1, self.num_neurons, [0.0] * self.num_neurons)
        self.biases = TVector(1, self.num_neurons, [0.0] * self.num_neurons)
        self.taus = TVector(1, self.num_neurons, [1.0] * self.num_neurons)
        self.chem_weights = {} # Store as dict: (pre, post): weight
        self.elec_weights = {} # Store as dict: (n1, n2): weight


    def SetCircuitSize(self, num_neurons, state_size, output_size): # state/output size unused here
        print(f"PlaceholderNervousSystem: SetCircuitSize({num_neurons}, ...) called")
        # Reinitialize if size changes (though it shouldn't in this context)
        self.num_neurons = num_neurons
        self.outputs = TVector(1, self.num_neurons, [0.0] * self.num_neurons)
        self.states = TVector(1, self.num_neurons, [0.0] * self.num_neurons)
        self.biases = TVector(1, self.num_neurons, [0.0] * self.num_neurons)
        self.taus = TVector(1, self.num_neurons, [1.0] * self.num_neurons)
        self.chem_weights = {}
        self.elec_weights = {}

    def SetNeuronBias(self, i, bias):
        self.biases[i] = bias

    def SetNeuronTimeConstant(self, i, tau):
         self.taus[i] = tau

    def SetChemicalSynapseWeight(self, pre, post, weight):
        self.chem_weights[(pre, post)] = weight

    def SetElectricalSynapseWeight(self, n1, n2, weight):
         # Store electrical weights uniquely (e.g., always store with lower index first)
         if n1 > n2: n1, n2 = n2, n1
         self.elec_weights[(n1, n2)] = weight

    def SetNeuronExternalInput(self, i, input_val):
        # Placeholder: Could store this to influence EulerStep
        pass

    def SetNeuronOutput(self, i, output_val):
         if 1 <= i <= self.num_neurons:
              self.outputs[i] = output_val
              self.states[i] = math.atanh(output_val) # Approx inverse sigmoid if needed

    def EulerStep(self, StepSize):
        # Simulate minimal random activity
        for i in range(1, self.num_neurons + 1):
            # Very simple Euler step placeholder - applies bias and decays state
            # A real CTRNN step is much more complex (sum inputs, etc.)
            dydt = (-self.states[i] + self.biases[i]) / self.taus[i] # Simplified update
            self.states[i] += dydt * StepSize
            self.outputs[i] = math.tanh(self.states[i]) # Sigmoid output
            # Add noise
            self.outputs[i] = max(0, min(1, self.outputs[i] + random.gauss(0, 0.05)*StepSize)) # Clip to [0,1] like sigmoid
        pass

    def NeuronOutput(self, i):
        if 1 <= i <= self.num_neurons:
            return self.outputs[i]
        return 0.0 # Default if out of bounds

    def NeuronState(self, i):
         if 1 <= i <= self.num_neurons:
              return self.states[i]
         return 0.0

    def NeuronBias(self, i): return self.biases[i]
    def NeuronTimeConstant(self, i): return self.taus[i]

    def ChemicalSynapseWeight(self, pre, post):
        return self.chem_weights.get((pre, post), 0.0)

    def ElectricalSynapseWeight(self, n1, n2):
        if n1 > n2: n1, n2 = n2, n1
        return self.elec_weights.get((n1, n2), 0.0)

# --- Helper Function ---
def nn(neuronNumber: int, unitNumber: int) -> int:
  """Calculates the flat index for a neuron given its type and unit number."""
  # Assumes neuronNumber is 1-based index within unit
  # Assumes unitNumber is 1-based index
  return neuronNumber + ((unitNumber - 1) * N_NEURONSPERUNIT)

# --- Worm Class ---
class Worm:
    def __init__(self, v_phenotype: TVector[float], output_unused: float):
        """
        Worm constructor.
        Initializes Body, Muscles, NervousSystem, StretchReceptors,
        and sets up connections based on the phenotype vector 'v_phenotype'.
        """
        # Instantiate components (using placeholders for now)
        self.b = WormBody()
        self.m = Muscles()
        self.n = NervousSystem()
        self.sr = StretchReceptor()

        self.t = 0.0 # Simulation time
        self.step_count = 0 # For dump skipping logic

        # --- Initialize components based on phenotype vector ---
        # Muscles
        self.m.SetMuscleParams(N_MUSCLES, T_MUSCLE)

        # Nervous system size
        self.n.SetCircuitSize(N_UNITS * N_NEURONSPERUNIT, 3, 2) # State/output sizes unused in placeholder

        # Stretch receptor gains
        # Indices SR_A=1, SR_B=2 assumed from main.py
        self.sr.SetStretchReceptorParams(N_SEGMENTS, N_STRETCHREC, v_phenotype[1], v_phenotype[2])

        # --- Set up Nervous System Connections ---
        for u in range(1, N_UNITS + 1):
            # Neuron indices within this unit
            da = nn(DA, u)
            db = nn(DB, u)
            dd = nn(DD, u)
            vd = nn(VD, u)
            va = nn(VA, u)
            vb = nn(VB, u)

            # Indices for neurons in the *next* unit (for inter-unit connections)
            # Handle boundary condition for the last unit
            ddNext = nn(DD, u + 1) if u < N_UNITS else -1 # Use -1 or some invalid index
            vdNext = nn(VD, u + 1) if u < N_UNITS else -1
            vbNext = nn(VB, u + 1) if u < N_UNITS else -1
            dbNext = nn(DB, u + 1) if u < N_UNITS else -1

            # Biases (indices 3, 4, 5 from main.py)
            self.n.SetNeuronBias(da, v_phenotype[3])
            self.n.SetNeuronBias(va, v_phenotype[3])
            self.n.SetNeuronBias(db, v_phenotype[4])
            self.n.SetNeuronBias(vb, v_phenotype[4])
            self.n.SetNeuronBias(dd, v_phenotype[5])
            self.n.SetNeuronBias(vd, v_phenotype[5])

            # Time-constants fixed to 1.0
            for i in range(1, N_NEURONSPERUNIT + 1):
                self.n.SetNeuronTimeConstant(nn(i, u), 1.0)

            # Self-connections (indices 6, 7, 8)
            self.n.SetChemicalSynapseWeight(da, da, v_phenotype[6])
            self.n.SetChemicalSynapseWeight(va, va, v_phenotype[6])
            self.n.SetChemicalSynapseWeight(db, db, v_phenotype[7])
            self.n.SetChemicalSynapseWeight(vb, vb, v_phenotype[7])
            self.n.SetChemicalSynapseWeight(dd, dd, v_phenotype[8])
            self.n.SetChemicalSynapseWeight(vd, vd, v_phenotype[8])

            # Cross-connections Intraunit
            # Excitatory Chemical Synapses (indices 9, 10)
            self.n.SetChemicalSynapseWeight(da, vd, v_phenotype[9])
            self.n.SetChemicalSynapseWeight(va, dd, v_phenotype[9])
            self.n.SetChemicalSynapseWeight(vb, dd, v_phenotype[10])
            self.n.SetChemicalSynapseWeight(db, vd, v_phenotype[10])

            # Inhibitory Chemical Synapses (indices 11, 12)
            self.n.SetChemicalSynapseWeight(vd, va, v_phenotype[11])
            self.n.SetChemicalSynapseWeight(dd, da, v_phenotype[11])
            self.n.SetChemicalSynapseWeight(vd, vb, v_phenotype[12])
            self.n.SetChemicalSynapseWeight(dd, db, v_phenotype[12])

            # Electrical Synapse Intersegment connections (indices 13, 14)
            if u < N_UNITS:
                # Check if next indices are valid before setting
                if ddNext != -1: self.n.SetElectricalSynapseWeight(dd, ddNext, v_phenotype[13])
                if vdNext != -1: self.n.SetElectricalSynapseWeight(vd, vdNext, v_phenotype[13])
                if vbNext != -1: self.n.SetElectricalSynapseWeight(vb, vbNext, v_phenotype[14])
                if dbNext != -1: self.n.SetElectricalSynapseWeight(db, dbNext, v_phenotype[14])

        # --- Neuromuscular Junction (NMJ) Weights ---
        # Excitatory VNC NMJ Weight (indices 15, 16)
        self.NMJ_DA = v_phenotype[15]
        self.NMJ_VA = v_phenotype[15]
        self.NMJ_DB = v_phenotype[16]
        self.NMJ_VB = v_phenotype[16]

        # Inhibitory VNC NMJ Weight (index 17)
        self.NMJ_DD = v_phenotype[17]
        self.NMJ_VD = v_phenotype[17]

        # --- Command Interneuron Outputs ---
        self.AVA_output = 0.0 # Set during simulation step based on direction
        self.AVB_output = 0.0

        # Define states for active/inactive command neurons (used in main.py/Evaluation)
        # These could potentially be parameters themselves.
        self.AVA_act = 1.0
        self.AVA_inact = 0.0
        self.AVB_act = 1.0
        self.AVB_inact = 0.0

    def InitializeState(self, rs: random.Random):
        """Initializes the state of the worm components."""
        self.t = 0.0
        self.step_count = 0
        # Initialize neuron outputs (specific pattern from C++)
        for u in range(1, N_UNITS + 1):
            # Dorsal neurons
            self.n.SetNeuronOutput(nn(DA, u), 0.1)
            self.n.SetNeuronOutput(nn(DB, u), 0.1)
            self.n.SetNeuronOutput(nn(DD, u), 0.9)
            # Ventral neurons
            self.n.SetNeuronOutput(nn(VA, u), 0.9)
            self.n.SetNeuronOutput(nn(VB, u), 0.9)
            self.n.SetNeuronOutput(nn(VD, u), 0.1)
            # Note: C++ code used RandomizeCircuitState commented out,
            # using specific init values instead. We follow the specific values.

        self.b.InitializeBodyState()
        self.m.InitializeMuscleState()

    def Step(self, StepSize: float, output_unused: float):
        """Performs one simulation step."""
        # --- Update Body ---
        self.b.StepBody(StepSize)

        # --- Set input to Stretch Receptors from Body ---
        for i in range(1, N_SEGMENTS + 1):
            dorsal_len = self.b.DorsalSegmentLength(i)
            ventral_len = self.b.VentralSegmentLength(i)
            rest_len = self.b.RestingLength(i)

            # Calculate relative stretch/compression
            ds = (dorsal_len - rest_len) / rest_len if rest_len != 0 else 0.0
            vs = (ventral_len - rest_len) / rest_len if rest_len != 0 else 0.0

            # --- Apply SR Transformation (Ported from C++ #defines) ---
            # NOTE: These depend on configuration. Defaulting to LINEAR (no transformation)
            #       as no defines were active in the provided C++ header snippet.
            #       To enable these, set corresponding boolean flags (e.g., self.SR_TRANS_STRETCH)
            #       in __init__ based on some configuration parameter.

            # if self.SR_TRANS_STRETCH: # Example flag
            #     ds = max(0.0, ds)
            #     vs = max(0.0, vs)
            # if self.SR_TRANS_CONTRACT: # Example flag
            #     ds = min(0.0, ds)
            #     vs = min(0.0, vs)
            # if self.SR_TRANS_ABS: # Example flag
            #     ds = abs(ds)
            #     vs = abs(vs)
            # if self.SR_TRANS_NEG: # Example flag
            #     ds = -ds
            #     vs = -vs

            self.sr.SetDorsalInput(i, ds)
            self.sr.SetVentralInput(i, vs)

        # --- Update Stretch Receptors ---
        self.sr.Update()

        # --- Set input to Nervous System (VNC) from SRs and Command Neurons ---
        # To A_class motorneurons
        for i in range(1, N_UNITS + 1):
            self.n.SetNeuronExternalInput(nn(DA, i), self.sr.A_D_sr(i) + self.AVA_output)
            self.n.SetNeuronExternalInput(nn(VA, i), self.sr.A_V_sr(i) + self.AVA_output)
        # To B_class motorneurons
        for i in range(1, N_UNITS + 1):
            self.n.SetNeuronExternalInput(nn(DB, i), self.sr.B_D_sr(i) + self.AVB_output)
            self.n.SetNeuronExternalInput(nn(VB, i), self.sr.B_V_sr(i) + self.AVB_output)

        # --- Update Nervous System ---
        self.n.EulerStep(StepSize)

        # --- Set input to Muscles from VNC motorneurons ---
        # Calculate combined input signal per unit first
        dorsalInput = TVector(1, N_UNITS)
        ventralInput = TVector(1, N_UNITS)
        for i in range(1, N_UNITS + 1):
            d_input = (self.NMJ_DA * self.n.NeuronOutput(nn(DA, i)) +
                       self.NMJ_DB * self.n.NeuronOutput(nn(DB, i)) +
                       self.NMJ_DD * self.n.NeuronOutput(nn(DD, i)))
            v_input = (self.NMJ_VD * self.n.NeuronOutput(nn(VD, i)) +
                       self.NMJ_VA * self.n.NeuronOutput(nn(VA, i)) +
                       self.NMJ_VB * self.n.NeuronOutput(nn(VB, i)))
            dorsalInput[i] = d_input
            ventralInput[i] = v_input

        # --- Apply inputs to muscles based on complex innervation pattern ---
        # Muscles 1-3 (Innervated by unit 1)
        for mi in range(1, 3 + 1):
            self.m.SetVentralMuscleInput(mi, ventralInput[1])
            self.m.SetDorsalMuscleInput(mi, dorsalInput[1])

        # Muscle 4 (Innervated by units 1 + 2)
        mi = 4
        self.m.SetVentralMuscleInput(mi, (ventralInput[1] + ventralInput[2]))
        self.m.SetDorsalMuscleInput(mi, (dorsalInput[1] + dorsalInput[2]))

        # Muscle 5 (Innervated by unit 2)
        mi = 5
        self.m.SetVentralMuscleInput(mi, ventralInput[2])
        self.m.SetDorsalMuscleInput(mi, dorsalInput[2])

        # Muscles 6-19 (Alternating innervation pattern)
        mt = 2 # Index of the *first* innervating unit for the current pair
        for mi in range(6, 19 + 1):
            self.m.SetVentralMuscleInput(mi, (ventralInput[mt] + ventralInput[mt + 1]))
            self.m.SetDorsalMuscleInput(mi, (dorsalInput[mt] + dorsalInput[mt + 1]))
            # Increment unit index 'mt' every two muscles (i.e., when mi is odd)
            if mi % 2 != 0: # Check if mi is odd (7, 9, ..., 19)
                mt += 1

        # Muscle 20 (Innervated by unit 9)
        mi = 20
        self.m.SetVentralMuscleInput(mi, ventralInput[9])
        self.m.SetDorsalMuscleInput(mi, dorsalInput[9])

        # Muscle 21 (Innervated by units 9 + 10)
        mi = 21
        self.m.SetVentralMuscleInput(mi, (ventralInput[9] + ventralInput[10]))
        self.m.SetDorsalMuscleInput(mi, (dorsalInput[9] + dorsalInput[10]))

        # Muscles 22-24 (Innervated by unit 10)
        for mi in range(22, 24 + 1):
            self.m.SetVentralMuscleInput(mi, ventralInput[10])
            self.m.SetDorsalMuscleInput(mi, dorsalInput[10])

        # --- Update Muscle activation ---
        self.m.EulerStep(StepSize)

        # --- Set input to Mechanical Body from Muscle forces ---
        # Segments 1, 2 (Innervated by muscle 1)
        activation1_d = self.m.DorsalMuscleOutput(1) / 2.0
        activation1_v = self.m.VentralMuscleOutput(1) / 2.0
        self.b.SetDorsalSegmentActivation(1, activation1_d)
        self.b.SetVentralSegmentActivation(1, activation1_v)
        self.b.SetDorsalSegmentActivation(2, activation1_d)
        self.b.SetVentralSegmentActivation(2, activation1_v)

        # Segments 3 to N_segments-2 (Innervated by two muscles)
        for i in range(3, N_SEGMENTS - 2 + 1):
            # Muscle index calculation (matches C++)
            mi = (i - 1) // 2 # Integer division gives the index of the first muscle
            act_d = (self.m.DorsalMuscleOutput(mi) + self.m.DorsalMuscleOutput(mi + 1)) / 2.0
            act_v = (self.m.VentralMuscleOutput(mi) + self.m.VentralMuscleOutput(mi + 1)) / 2.0
            self.b.SetDorsalSegmentActivation(i, act_d)
            self.b.SetVentralSegmentActivation(i, act_v)

        # Segments N_segments-1, N_segments (Innervated by last muscle)
        activation_last_d = self.m.DorsalMuscleOutput(N_MUSCLES) / 2.0
        activation_last_v = self.m.VentralMuscleOutput(N_MUSCLES) / 2.0
        self.b.SetDorsalSegmentActivation(N_SEGMENTS - 1, activation_last_d)
        self.b.SetVentralSegmentActivation(N_SEGMENTS - 1, activation_last_v)
        self.b.SetDorsalSegmentActivation(N_SEGMENTS, activation_last_d)
        self.b.SetVentralSegmentActivation(N_SEGMENTS, activation_last_v)

        # --- Time and Step Count ---
        self.t += StepSize
        self.step_count += 1

    def CoMx(self) -> float:
        """Calculates the X coordinate of the Center of Mass."""
        temp = 0.0
        if N_RODS == 0: # Avoid division by zero
            return 0.0
        # Loop from 1 to N_RODS (inclusive), matching C++ loop
        for i in range(1, N_RODS + 1):
            temp += self.b.X(i) # Call the X() method of the WormBody instance
        return temp / N_RODS

    def CoMy(self) -> float:
        """Calculates the Y coordinate of the Center of Mass."""
        temp = 0.0
        if N_RODS == 0: # Avoid division by zero
            return 0.0
        # Loop from 1 to N_RODS (inclusive), matching C++ loop
        for i in range(1, N_RODS + 1):
             temp += self.b.Y(i) # Call the Y() method of the WormBody instance
        return temp / N_RODS

    def _calculate_angle_diff(self, a1, a2):
        """ Helper to calculate angle difference, handling wrap-around at +/- PI."""
        # Normalize angles to be within [-pi, pi] might simplify, but follow C++ logic:
        if a1 > math.pi / 2 and a2 < -math.pi / 2:
            a = (a1 - 2 * math.pi) - a2
        elif a1 < -math.pi / 2 and a2 > math.pi / 2:
            a = a1 - (a2 - 2 * math.pi)
        else:
            a = a1 - a2
        return a

    def Curvature(self, c: TVector[float]):
        """Calculates curvature along the body and stores it in vector c."""
        # Expects c to be pre-allocated with the correct size (N_SEGMENTS-3)//2 + 1 ? -> 23 if N_SEGMENTS=48
        k = c.LowerBound() # Start writing at the lower bound of c
        # Loop through segments where curvature can be calculated (needs 2 points before and 2 after)
        for i in range(3, N_SEGMENTS - 1, 2): # C++ loop: i=3, 5, ..., N_segments-3(or-2)
            dx1 = self.b.X(i) - self.b.X(i - 2)
            dy1 = self.b.Y(i) - self.b.Y(i - 2)
            dx2 = self.b.X(i + 2) - self.b.X(i)
            dy2 = self.b.Y(i + 2) - self.b.Y(i)

            a1 = math.atan2(dy1, dx1)
            a2 = math.atan2(dy2, dx2)

            a = self._calculate_angle_diff(a1, a2)

            # Calculate distance between the endpoints of the two segments used
            seg_dx = self.b.X(i + 2) - self.b.X(i - 2)
            seg_dy = self.b.Y(i + 2) - self.b.Y(i - 2)
            seg = math.sqrt(seg_dx**2 + seg_dy**2)

            if seg > 1e-9: # Avoid division by zero
                curvature_val = (2 * math.sin(a) / seg) / 1000.0 # Matches C++ scaling
            else:
                curvature_val = 0.0

            if k <= c.UpperBound():
                c[k] = curvature_val
            else:
                print("Warning: Curvature vector overflow", file=sys.stderr)
                break # Stop writing if vector is too small
            k += 1

    def AngleCurvature(self, c: TVector[float]):
        """Calculates angle differences along the body and stores it in vector c."""
        k = c.LowerBound()
        for i in range(3, N_SEGMENTS - 1, 2):
            dx1 = self.b.X(i) - self.b.X(i - 2)
            dy1 = self.b.Y(i) - self.b.Y(i - 2)
            dx2 = self.b.X(i + 2) - self.b.X(i)
            dy2 = self.b.Y(i + 2) - self.b.Y(i)

            a1 = math.atan2(dy1, dx1)
            a2 = math.atan2(dy2, dx2)

            a = self._calculate_angle_diff(a1, a2)

            if k <= c.UpperBound():
                c[k] = a
            else:
                print("Warning: AngleCurvature vector overflow", file=sys.stderr)
                break
            k += 1

    def Orientation(self) -> float:
        """Calculates the overall orientation from tail to head."""
        # Use HEAD and TAIL constants which depend on N_SEGMENTS
        # Assumes X, Y indices correspond to rod numbers (1 to N_rods)
        # Need to map HEAD/TAIL segments to rod indices (often HEAD=1, TAIL=N_rods)
        head_rod_idx = HEAD # Typically rod 1 is head
        tail_rod_idx = N_RODS # Typically last rod is tail
        dx = self.b.X(head_rod_idx) - self.b.X(tail_rod_idx)
        dy = self.b.Y(head_rod_idx) - self.b.Y(tail_rod_idx)
        return math.atan2(dy, dx)

    # --- Dump Methods ---

    def DumpBodyState(self, ofs: TextIO, skips: int):
        """Dumps body state (X, Y, Phi) to file stream if skips steps passed."""
        if self.step_count % skips == 0:
            ofs.write(f"{self.t:.4f}") # Use more precision if needed

            for i in range(1, N_RODS + 1):
                ofs.write(f" {self.b.X(i):.8f} {self.b.Y(i):.8f} {self.b.Phi(i):.8f}")
            ofs.write("\n")

    def DumpActState(self, ofs: TextIO, skips: int):
        """Dumps SR, Neuron, and Muscle activations to file stream."""
        if self.step_count % skips == 0:
            ofs.write(f"{self.t:.4f}")
            # Stretch receptors
            for i in range(1, self.sr.NSR + 1):
                 ofs.write(f" {self.sr.A_D_sr(i):.8f} {self.sr.A_V_sr(i):.8f}"
                           f" {self.sr.B_D_sr(i):.8f} {self.sr.B_V_sr(i):.8f}")
            # Ventral Cord Motor Neurons Outputs
            for i in range(1, N_UNITS + 1):
                for j in range(1, N_NEURONSPERUNIT + 1):
                    ofs.write(f" {self.n.NeuronOutput(nn(j, i)):.8f}")
            # Muscles Outputs
            for i in range(1, N_MUSCLES + 1):
                 ofs.write(f" {self.m.DorsalMuscleOutput(i):.8f}"
                           f" {self.m.VentralMuscleOutput(i):.8f}")
            ofs.write("\n")

    def DumpCurvature(self, ofs: TextIO, skips: int):
        """Dumps calculated curvature to file stream."""
        if self.step_count % skips == 0:
            ofs.write(f"{self.t:.4f}")
            # Re-calculate curvature for dumping (same logic as Curvature method)
            for i in range(3, N_SEGMENTS - 1, 2):
                dx1 = self.b.X(i) - self.b.X(i - 2)
                dy1 = self.b.Y(i) - self.b.Y(i - 2)
                dx2 = self.b.X(i + 2) - self.b.X(i)
                dy2 = self.b.Y(i + 2) - self.b.Y(i)

                a1 = math.atan2(dy1, dx1)
                a2 = math.atan2(dy2, dx2)
                a = self._calculate_angle_diff(a1, a2)

                seg_dx = self.b.X(i + 2) - self.b.X(i - 2)
                seg_dy = self.b.Y(i + 2) - self.b.Y(i - 2)
                seg = math.sqrt(seg_dx**2 + seg_dy**2)

                if seg > 1e-9:
                    curvature_val = (2 * math.sin(a) / seg) / 1000.0
                else:
                    curvature_val = 0.0
                ofs.write(f" {curvature_val:.8f}")
            ofs.write("\n")


    def DumpVoltage(self, ofs: TextIO, skips: int):
        """Dumps Neuron States (voltages/internal states) to file stream."""
        if self.step_count % skips == 0:
            ofs.write(f"{self.t:.4f}")
            # Ventral Cord Motor Neurons States
            for i in range(1, N_UNITS + 1):
                for j in range(1, N_NEURONSPERUNIT + 1):
                    ofs.write(f" {self.n.NeuronState(nn(j, i)):.8f}")
            ofs.write("\n")

    def DumpParams(self, ofs: TextIO):
        """Writes key simulation parameters (weights, biases, gains) to file."""
        ofs.write("--- Worm Simulation Parameters ---\n")
        # Assuming accessors exist in placeholder/actual NervousSystem
        # Using nn(X, 1) to get params for the first unit as representative
        ofs.write("Time-constants (Unit 1):\n")
        ofs.write(f"  DA: {self.n.NeuronTimeConstant(nn(DA, 1))}\n")
        ofs.write(f"  DB: {self.n.NeuronTimeConstant(nn(DB, 1))}\n")
        ofs.write(f"  DD: {self.n.NeuronTimeConstant(nn(DD, 1))}\n")
        ofs.write(f"  VD: {self.n.NeuronTimeConstant(nn(VD, 1))}\n")
        ofs.write(f"  VA: {self.n.NeuronTimeConstant(nn(VA, 1))}\n")
        ofs.write(f"  VB: {self.n.NeuronTimeConstant(nn(VB, 1))}\n\n")

        ofs.write("Biases (Unit 1):\n")
        ofs.write(f"  DA/VA: {self.n.NeuronBias(nn(DA, 1))}\n") # DA/VA share bias
        ofs.write(f"  DB/VB: {self.n.NeuronBias(nn(DB, 1))}\n") # DB/VB share bias
        ofs.write(f"  DD/VD: {self.n.NeuronBias(nn(DD, 1))}\n") # DD/VD share bias
        ofs.write("\n")

        ofs.write("Self connections (Unit 1):\n")
        ofs.write(f"  DA/VA->DA/VA: {self.n.ChemicalSynapseWeight(nn(DA, 1), nn(DA, 1))}\n") # Shared
        ofs.write(f"  DB/VB->DB/VB: {self.n.ChemicalSynapseWeight(nn(DB, 1), nn(DB, 1))}\n") # Shared
        ofs.write(f"  DD/VD->DD/VD: {self.n.ChemicalSynapseWeight(nn(DD, 1), nn(DD, 1))}\n") # Shared
        ofs.write("\n")

        ofs.write("Interneuron properties:\n")
        ofs.write(f"  AVA active state: {self.AVA_act}\n")
        ofs.write(f"  AVB active state: {self.AVB_act}\n")
        ofs.write(f"  AVA inactive state: {self.AVA_inact}\n")
        ofs.write(f"  AVB inactive state: {self.AVB_inact}\n\n")

        # C++ DumpParams showed VA->VD and VB->VD which look like typos based on constructor?
        # Assuming C++ constructor logic is correct for connectivity:
        ofs.write("Chemical Connections (Unit 1):\n")
        ofs.write(f"  DA->VD: {self.n.ChemicalSynapseWeight(nn(DA, 1), nn(VD, 1))}\n")
        ofs.write(f"  VA->DD: {self.n.ChemicalSynapseWeight(nn(VA, 1), nn(DD, 1))}\n")
        ofs.write(f"  VB->DD: {self.n.ChemicalSynapseWeight(nn(VB, 1), nn(DD, 1))}\n")
        ofs.write(f"  DB->VD: {self.n.ChemicalSynapseWeight(nn(DB, 1), nn(VD, 1))}\n")
        ofs.write(f"  VD->VA: {self.n.ChemicalSynapseWeight(nn(VD, 1), nn(VA, 1))}\n")
        ofs.write(f"  DD->DA: {self.n.ChemicalSynapseWeight(nn(DD, 1), nn(DA, 1))}\n")
        ofs.write(f"  VD->VB: {self.n.ChemicalSynapseWeight(nn(VD, 1), nn(VB, 1))}\n")
        ofs.write(f"  DD->DB: {self.n.ChemicalSynapseWeight(nn(DD, 1), nn(DB, 1))}\n")
        ofs.write("\n")

        ofs.write("Gap Junctions (Unit 1 to 2):\n")
        ofs.write(f"  DD-DD+1: {self.n.ElectricalSynapseWeight(nn(DD, 1), nn(DD, 2))}\n")
        ofs.write(f"  VD-VD+1: {self.n.ElectricalSynapseWeight(nn(VD, 1), nn(VD, 2))}\n")
        ofs.write(f"  VB-VB+1: {self.n.ElectricalSynapseWeight(nn(VB, 1), nn(VB, 2))}\n")
        ofs.write(f"  DB-DB+1: {self.n.ElectricalSynapseWeight(nn(DB, 1), nn(DB, 2))}\n\n")

        ofs.write("Stretch Receptors Gains:\n")
        ofs.write(f"  A-class SR Gain: {self.sr.SR_A_gain}\n")
        ofs.write(f"  B-class SR Gain: {self.sr.SR_B_gain}\n\n")

        ofs.write("NMJ weights:\n")
        ofs.write(f"  DA/VA: {self.NMJ_DA}\n") # Shared weight
        ofs.write(f"  DB/VB: {self.NMJ_DB}\n") # Shared weight
        ofs.write(f"  DD: {self.NMJ_DD}\n")
        ofs.write(f"  VD: {self.NMJ_VD}\n\n")
        ofs.write("-------------------------------\n")