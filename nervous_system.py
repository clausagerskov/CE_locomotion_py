import math
import random
import sys
from typing import Optional, TextIO

# Assuming vector_matrix.py is in the same directory or accessible
from vector_matrix import TVector, TMatrix

# Helper class to mimic the C++ weightentry struct
class WeightEntry:
    def __init__(self, from_neuron: int = 0, weight_val: float = 0.0):
        self.from_neuron = from_neuron
        self.weight = weight_val

    def __str__(self):
        return f"(from={self.from_neuron}, w={self.weight:.4f})"

    def __repr__(self):
        return f"WeightEntry(from_neuron={self.from_neuron}, weight_val={self.weight})"


# Sigmoid and Inverse Sigmoid functions
def sigmoid(x: float) -> float:
    """Standard logistic sigmoid function."""
    try:
        return 1.0 / (1.0 + math.exp(-x))
    except OverflowError:
        # Handle large negative x gracefully
        return 0.0 if x < 0 else 1.0

def inverse_sigmoid(y: float) -> float:
    """Inverse of the sigmoid function (logit)."""
    # Ensure input is within the valid range (0, 1) exclusive
    if y <= 0.0 or y >= 1.0:
        # Handle edge cases: return large magnitude values or raise error
        # Returning large values might be more robust in simulation if outputs hit bounds
        epsilon = 1e-15 # Small epsilon
        y = max(epsilon, min(y, 1.0 - epsilon))
        # raise ValueError(f"Input to inverse_sigmoid must be in (0, 1), got {y}")
    try:
        return math.log(y / (1.0 - y))
    except ValueError: # Should be caught by the initial check, but as fallback
         raise ValueError(f"Math domain error in inverse_sigmoid for y={y}")


class NervousSystem:
    """
    Python port of the NervousSystem class (CTRNN implementation).
    Models the continuous-time dynamics of a neural network.
    """

    def __init__(self, size: int = 0, maxchemconns: int = -1, maxelecconns: int = -1):
        """
        Initializes the NervousSystem.
        Args:
            size: Number of neurons in the circuit.
            maxchemconns: Max chemical connections per neuron (-1 for size).
            maxelecconns: Max electrical connections per neuron (-1 for maxchemconns).
        """
        # Initialize attributes to default/empty states
        self.size: int = 0
        self.maxchemconns: int = 0
        self.maxelecconns: int = 0

        self.states: TVector[float] = TVector()
        self.paststates: TVector[float] = TVector() # State at previous time step (for gap junctions)
        self.outputs: TVector[float] = TVector()
        self.biases: TVector[float] = TVector()
        self.gains: TVector[float] = TVector()
        self.taus: TVector[float] = TVector()
        self.Rtaus: TVector[float] = TVector() # Reciprocal of taus (1/tau)
        self.externalinputs: TVector[float] = TVector()

        # Connection counts and weights (sparse storage)
        self.NumChemicalConns: TVector[int] = TVector()
        self.chemicalweights: TMatrix[WeightEntry] = TMatrix() # Stores WeightEntry objects

        self.NumElectricalConns: TVector[int] = TVector()
        self.electricalweights: TMatrix[WeightEntry] = TMatrix() # Stores WeightEntry objects

        # Temp vectors used in RK4 (if implemented) - allocate even if only Euler is used now
        self.TempStates: TVector[float] = TVector()
        self.TempOutputs: TVector[float] = TVector()
        self.k1: TVector[float] = TVector()
        self.k2: TVector[float] = TVector()
        self.k3: TVector[float] = TVector()
        self.k4: TVector[float] = TVector()

        # Set size and allocate vectors
        self.SetCircuitSize(size, maxchemconns, maxelecconns)

    def SetCircuitSize(self, newsize: int, newmaxchemconns: int, newmaxelecconns: int):
        """Resizes the circuit and associated vectors/matrices."""
        if newsize < 0: raise ValueError("Circuit size cannot be negative")
        self.size = newsize

        # Determine max connections
        if newmaxchemconns == -1:
            self.maxchemconns = self.size
        else:
            self.maxchemconns = min(newmaxchemconns, self.size)
        if self.maxchemconns < 0: self.maxchemconns = 0 # Ensure non-negative

        if newmaxelecconns == -1:
            self.maxelecconns = self.maxchemconns # Max elec defaults to max chem
        else:
            self.maxelecconns = min(newmaxelecconns, self.maxchemconns)
        if self.maxelecconns < 0: self.maxelecconns = 0 # Ensure non-negative


        # Resize vectors (using 1-based indexing)
        lb = 1
        ub = self.size
        self.states.SetBounds(lb, ub)
        self.states.FillContents(0.0)
        self.paststates.SetBounds(lb, ub)
        self.paststates.FillContents(0.0)
        self.outputs.SetBounds(lb, ub)
        self.outputs.FillContents(0.0)
        self.biases.SetBounds(lb, ub)
        self.biases.FillContents(0.0)
        self.gains.SetBounds(lb, ub)
        self.gains.FillContents(1.0) # Default gain = 1
        self.taus.SetBounds(lb, ub)
        self.taus.FillContents(1.0) # Default tau = 1
        self.Rtaus.SetBounds(lb, ub)
        self.Rtaus.FillContents(1.0) # Default Rtau = 1/1 = 1
        self.externalinputs.SetBounds(lb, ub)
        self.externalinputs.FillContents(0.0)

        # Connection counts
        self.NumChemicalConns.SetBounds(lb, ub)
        self.NumChemicalConns.FillContents(0)
        self.NumElectricalConns.SetBounds(lb, ub)
        self.NumElectricalConns.FillContents(0)

        # Weight matrices (sparse representation)
        # Rows: receiving neuron (1 to size)
        # Columns: connection index (1 to maxchemconns/maxelecconns)
        self.chemicalweights.SetBounds(lb, ub, 1, self.maxchemconns)
        self.electricalweights.SetBounds(lb, ub, 1, self.maxelecconns)
        # Initialize TMatrix entries with default WeightEntry objects
        default_entry = WeightEntry()
        if self.chemicalweights.RowSize() > 0 and self.chemicalweights.ColumnSize() > 0:
            self.chemicalweights.FillContents(default_entry)
        if self.electricalweights.RowSize() > 0 and self.electricalweights.ColumnSize() > 0:
            self.electricalweights.FillContents(default_entry)


        # Temp vectors for potential RK4 integration
        self.TempStates.SetBounds(lb, ub)
        self.TempOutputs.SetBounds(lb, ub)
        self.k1.SetBounds(lb, ub)
        self.k2.SetBounds(lb, ub)
        self.k3.SetBounds(lb, ub)
        self.k4.SetBounds(lb, ub)

    # --- Accessors ---
    def CircuitSize(self) -> int:
        return self.size

    def NeuronState(self, i: int) -> float:
        """Gets the internal state (e.g., voltage) of neuron i (1-based)."""
        if 1 <= i <= self.size:
            return self.states[i]
        else:
            raise IndexError(f"Neuron index {i} out of bounds [1, {self.size}]")

    def SetNeuronState(self, i: int, value: float):
        """Sets the internal state of neuron i and updates its output accordingly."""
        if 1 <= i <= self.size:
            self.states[i] = value
            # Update output based on new state
            self.outputs[i] = sigmoid(self.gains[i] * (self.states[i] + self.biases[i]))
        else:
            raise IndexError(f"Neuron index {i} out of bounds [1, {self.size}]")

    def NeuronOutput(self, i: int) -> float:
        """Gets the output (activation) of neuron i (1-based)."""
        if 1 <= i <= self.size:
            return self.outputs[i]
        else:
            # Optionally return 0 or raise error for out-of-bounds access
            # Returning 0 might be safer if code calling this expects it
            # print(f"Warning: NeuronOutput index {i} out of bounds [1, {self.size}]")
            return 0.0
            # raise IndexError(f"Neuron index {i} out of bounds [1, {self.size}]")

    def SetNeuronOutput(self, i: int, value: float):
        """Sets the output of neuron i and updates its internal state accordingly."""
        if 1 <= i <= self.size:
            # Clip output to avoid issues with inverse_sigmoid at 0 or 1
            epsilon = 1e-15
            clipped_value = max(epsilon, min(value, 1.0 - epsilon))
            # if value <= 0.0 or value >= 1.0: # Optional warning
            #    print(f"Warning: SetNeuronOutput called with value {value} for neuron {i}. Clipping to ({epsilon}, {1-epsilon}).")

            self.outputs[i] = clipped_value # Store clipped value
            # Update state based on new output
            gain_i = self.gains[i]
            if gain_i == 0:
                # Avoid division by zero; state is ill-defined. Set to 0?
                self.states[i] = 0.0
                 # print(f"Warning: Gain is zero for neuron {i}. Cannot accurately set state from output.")
            else:
                self.states[i] = inverse_sigmoid(self.outputs[i]) / gain_i - self.biases[i]
        else:
            raise IndexError(f"Neuron index {i} out of bounds [1, {self.size}]")

    def NeuronBias(self, i: int) -> float:
        if 1 <= i <= self.size: return self.biases[i]
        else: raise IndexError(f"Neuron index {i} out of bounds [1, {self.size}]")

    def SetNeuronBias(self, i: int, value: float):
        if 1 <= i <= self.size: self.biases[i] = value
        else: raise IndexError(f"Neuron index {i} out of bounds [1, {self.size}]")

    def NeuronGain(self, i: int) -> float:
        if 1 <= i <= self.size: return self.gains[i]
        else: raise IndexError(f"Neuron index {i} out of bounds [1, {self.size}]")

    def SetNeuronGain(self, i: int, value: float):
        if 1 <= i <= self.size: self.gains[i] = value
        else: raise IndexError(f"Neuron index {i} out of bounds [1, {self.size}]")

    def NeuronTimeConstant(self, i: int) -> float:
        if 1 <= i <= self.size: return self.taus[i]
        else: raise IndexError(f"Neuron index {i} out of bounds [1, {self.size}]")

    def SetNeuronTimeConstant(self, i: int, value: float):
        if 1 <= i <= self.size:
            if value <= 0: raise ValueError(f"Time constant must be positive, got {value}")
            self.taus[i] = value
            self.Rtaus[i] = 1.0 / value # Update reciprocal tau
        else: raise IndexError(f"Neuron index {i} out of bounds [1, {self.size}]")

    def NeuronExternalInput(self, i: int) -> float:
        if 1 <= i <= self.size: return self.externalinputs[i]
        else: raise IndexError(f"Neuron index {i} out of bounds [1, {self.size}]")

    def SetNeuronExternalInput(self, i: int, value: float):
        if 1 <= i <= self.size: self.externalinputs[i] = value
        else: raise IndexError(f"Neuron index {i} out of bounds [1, {self.size}]")

    def ChemicalSynapseWeight(self, from_neuron: int, to_neuron: int) -> float:
        """Gets the weight of the chemical synapse from->to (returns 0 if no connection)."""
        if not (1 <= to_neuron <= self.size): return 0.0 # Target neuron out of bounds

        num_conns = self.NumChemicalConns[to_neuron]
        for i in range(1, num_conns + 1): # Check existing connections for this neuron
            entry = self.chemicalweights[to_neuron][i]
            if entry.from_neuron == from_neuron:
                return entry.weight
        return 0.0 # Connection not found

    def SetChemicalSynapseWeight(self, from_neuron: int, to_neuron: int, value: float):
        """Sets the weight of the chemical synapse from->to."""
        if not (1 <= from_neuron <= self.size and 1 <= to_neuron <= self.size):
            raise IndexError("Neuron index out of bounds for chemical synapse.")

        # Check if the connection already exists
        num_conns = self.NumChemicalConns[to_neuron]
        for i in range(1, num_conns + 1):
            entry = self.chemicalweights[to_neuron][i]
            if entry.from_neuron == from_neuron:
                entry.weight = value # Update existing entry
                return

        # Connection doesn't exist, add a new one if space allows
        if num_conns >= self.maxchemconns:
            raise RuntimeError(f"Maximum chemical synapses ({self.maxchemconns}) exceeded for neuron {to_neuron}")

        # Add new entry
        new_conn_idx = num_conns + 1
        self.chemicalweights[to_neuron][new_conn_idx] = WeightEntry(from_neuron, value)
        self.NumChemicalConns[to_neuron] += 1 # Increment count

    def ElectricalSynapseWeight(self, from_neuron: int, to_neuron: int) -> float:
        """Gets the weight of the electrical synapse between from<->to (returns 0 if no connection)."""
        # Electrical synapses are symmetric, query based on 'to_neuron' list
        if not (1 <= to_neuron <= self.size): return 0.0

        num_conns = self.NumElectricalConns[to_neuron]
        for i in range(1, num_conns + 1):
            entry = self.electricalweights[to_neuron][i]
            if entry.from_neuron == from_neuron:
                return entry.weight
        return 0.0

    def _InternalSetElectricalSynapseWeight(self, from_neuron: int, to_neuron: int, value: float):
        """Internal helper: Sets one direction of the electrical synapse."""
        if not (1 <= from_neuron <= self.size and 1 <= to_neuron <= self.size):
            raise IndexError("Neuron index out of bounds for electrical synapse.")

        # Check if the connection already exists in 'to_neuron's list
        num_conns = self.NumElectricalConns[to_neuron]
        for i in range(1, num_conns + 1):
            entry = self.electricalweights[to_neuron][i]
            if entry.from_neuron == from_neuron:
                entry.weight = value # Update existing entry
                return

        # Connection doesn't exist, add a new one if space allows
        if num_conns >= self.maxelecconns:
            raise RuntimeError(f"Maximum electrical synapses ({self.maxelecconns}) exceeded for neuron {to_neuron}")

        # Add new entry
        new_conn_idx = num_conns + 1
        self.electricalweights[to_neuron][new_conn_idx] = WeightEntry(from_neuron, value)
        self.NumElectricalConns[to_neuron] += 1

    def SetElectricalSynapseWeight(self, n1: int, n2: int, value: float):
        """Sets the weight of the electrical synapse (gap junction) between n1 and n2."""
        if value < 0:
            raise ValueError(f"Electrical synapse weight between {n1} and {n2} cannot be negative: {value}")

        # Set the connection in both directions for the sparse storage
        self._InternalSetElectricalSynapseWeight(n1, n2, value)
        self._InternalSetElectricalSynapseWeight(n2, n1, value)


    # --- Control Methods ---
    def RandomizeCircuitState(self, lb: float, ub: float, rs: Optional[random.Random] = None):
        """Randomizes the internal states of all neurons."""
        rand_gen = rs if rs is not None else random # Use provided generator or default
        for i in range(1, self.size + 1):
            self.SetNeuronState(i, rand_gen.uniform(lb, ub)) # Updates output too

    def RandomizeCircuitOutput(self, lb: float, ub: float, rs: Optional[random.Random] = None):
        """Randomizes the outputs of all neurons."""
        rand_gen = rs if rs is not None else random
        for i in range(1, self.size + 1):
            self.SetNeuronOutput(i, rand_gen.uniform(lb, ub)) # Updates state too

    def EulerStep(self, stepsize: float):
        """
        Integrates the circuit dynamics one step using Euler's method.
        dState/dt = (1/Tau) * (SummedInputs - State)
        SummedInputs = Bias + ExtInput + Sum(ChemWeight * Output_From) + Sum(ElecWeight * (State_From - State_This))
        """
        if stepsize <= 0:
             raise ValueError("Stepsize must be positive for EulerStep")

        # 1. Store current states to use for gap junction calculation (past states)
        # Note: A simple copy might suffice if TVector copy is deep enough, or use loop/numpy slice
        for i in range(1, self.size + 1):
            self.paststates[i] = self.states[i]
        # self.paststates = self.states.copy() # If using numpy or TVector has deep copy

        # 2. Calculate state change for all neurons
        delta_states = TVector(1, self.size) # Temporary storage for state changes
        delta_states.FillContents(0.0)

        for i in range(1, self.size + 1): # For each neuron 'i' (receiving neuron)
            # Calculate total input
            total_input = self.externalinputs[i]

            # Input from chemical synapses
            num_chem_conns = self.NumChemicalConns[i]
            for j in range(1, num_chem_conns + 1):
                entry = self.chemicalweights[i][j]
                # Use NeuronOutput to get the output of the pre-synaptic neuron
                presynaptic_output = self.NeuronOutput(entry.from_neuron)
                total_input += entry.weight * presynaptic_output

            # Input from electrical synapses (gap junctions)
            num_elec_conns = self.NumElectricalConns[i]
            for j in range(1, num_elec_conns + 1):
                entry = self.electricalweights[i][j]
                # Use PAST state of the connected neuron and PAST state of this neuron
                state_from = self.paststates[entry.from_neuron]
                state_this = self.paststates[i]
                total_input += entry.weight * (state_from - state_this)

            # Calculate change in state using Euler's method formula
            # dState/dt = Rtau * (total_input + bias - state) ??? NO, C++ code is:
            # dState/dt = Rtau * (total_input_WITHOUT_bias - state)
            # Bias is added INSIDE the sigmoid later. Let's follow the C++ EulerStep.
            # The total_input calculated above already includes chemical & electrical effects.
            # The CTRNN equation is: tau * dS/dt = -S + Bias + Sum(Weights * Output) + Sum(Gap * (Sj-Si)) + ExtInput
            # So: dS/dt = Rtau * (-S + Bias + Sum(Chem) + Sum(Gap) + ExtInput)
            # It seems the C++ code sums ExtInput, ChemInput, GapInput into 'input'
            # Then calculates: dS/dt = Rtau * (input - state)
            # This implies 'input' must implicitly contain the bias term for the equation to match.
            # Re-checking C++: No, bias is added in sigmoid step.
            # Let's stick to the C++ code's Euler step implementation:
            # Rate of change = Rtau * (sum_of_inputs - current_state)
            # Where sum_of_inputs = external + chemical_input + electrical_input
            # Note: Bias is NOT included in the state update differential equation here.
            dydt = self.Rtaus[i] * (total_input - self.states[i])
            delta_states[i] = stepsize * dydt


        # 3. Update states (simultaneously, using calculated changes)
        for i in range(1, self.size + 1):
            self.states[i] += delta_states[i]

        # 4. Update outputs based on new states and biases
        for i in range(1, self.size + 1):
            self.outputs[i] = sigmoid(self.gains[i] * (self.states[i] + self.biases[i]))


    # RK4Step method omitted as it was commented out in the C++ source provided
    # def RK4Step(self, stepsize): pass