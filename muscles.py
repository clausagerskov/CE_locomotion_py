import math
# Assuming vector_matrix.py is in the same directory or accessible
from vector_matrix import TVector, TMatrix

class Muscles:
    """
    Python port of the Muscles class.
    Models muscle activation dynamics using first-order Euler integration.
    Manages dorsal and ventral muscle states separately.
    """

    def __init__(self, nmuscles: int = 24, t_muscle: float = 0.1):
        """
        Initializes the Muscles instance.
        Args:
            nmuscles: Number of muscle pairs (dorsal/ventral) along the body.
            t_muscle: Time constant for muscle activation dynamics.
        """
        # Initialize attributes with defaults first
        self.Nmuscles: int = 0
        self.T_muscle: float = 0.0
        self.V_muscle: TMatrix[float] = TMatrix() # Activation state [muscle_idx][1=Dorsal, 2=Ventral]
        self.V_input: TMatrix[float] = TMatrix()  # Input signal [muscle_idx][1=Dorsal, 2=Ventral]

        # Call SetMuscleParams to properly initialize based on arguments
        self.SetMuscleParams(nmuscles, t_muscle)

    def SetMuscleParams(self, nmuscles: int, t_muscle: float):
        """
        Sets the parameters and initializes/resizes the internal matrices.
        """
        if nmuscles <= 0:
            raise ValueError("Number of muscles must be positive.")
        if t_muscle <= 0:
            raise ValueError("Muscle time constant must be positive.")

        self.Nmuscles = nmuscles
        self.T_muscle = t_muscle

        # Initialize/Resize TMatrix instances
        # Rows: 1 to Nmuscles (muscle index)
        # Columns: 1 (Dorsal), 2 (Ventral)
        self.V_muscle.SetBounds(1, self.Nmuscles, 1, 2)
        self.V_input.SetBounds(1, self.Nmuscles, 1, 2)

        # Optional: Fill with default values (0.0) - done in InitializeMuscleState
        self.InitializeMuscleState()

    def InitializeMuscleState(self):
        """Initializes muscle activation and input to zero."""
        # Check if matrices have been allocated (have non-zero size)
        if self.V_muscle.RowSize() > 0 and self.V_muscle.ColumnSize() > 0:
            self.V_muscle.FillContents(0.0)
        if self.V_input.RowSize() > 0 and self.V_input.ColumnSize() > 0:
            self.V_input.FillContents(0.0)

    def SetDorsalMuscleInput(self, muscle: int, input_val: float):
        """Sets the input signal for a specific dorsal muscle (1-based index)."""
        if 1 <= muscle <= self.Nmuscles:
            self.V_input[muscle][1] = input_val # Access row 'muscle', column 1 (Dorsal)
        else:
            print(f"Warning: SetDorsalMuscleInput index {muscle} out of bounds [1, {self.Nmuscles}]")

    def SetVentralMuscleInput(self, muscle: int, input_val: float):
        """Sets the input signal for a specific ventral muscle (1-based index)."""
        if 1 <= muscle <= self.Nmuscles:
            self.V_input[muscle][2] = input_val # Access row 'muscle', column 2 (Ventral)
        else:
            print(f"Warning: SetVentralMuscleInput index {muscle} out of bounds [1, {self.Nmuscles}]")

    def SetDorsalMuscleActivation(self, muscle: int, activation: float):
        """Directly sets the activation state for a specific dorsal muscle (1-based index)."""
        if 1 <= muscle <= self.Nmuscles:
            self.V_muscle[muscle][1] = activation # Access row 'muscle', column 1 (Dorsal)
        else:
            print(f"Warning: SetDorsalMuscleActivation index {muscle} out of bounds [1, {self.Nmuscles}]")

    def SetVentralMuscleActivation(self, muscle: int, activation: float):
        """Directly sets the activation state for a specific ventral muscle (1-based index)."""
        if 1 <= muscle <= self.Nmuscles:
            self.V_muscle[muscle][2] = activation # Access row 'muscle', column 2 (Ventral)
        else:
            print(f"Warning: SetVentralMuscleActivation index {muscle} out of bounds [1, {self.Nmuscles}]")

    def DorsalMuscleOutput(self, muscle: int) -> float:
        """Gets the current activation output of a specific dorsal muscle (1-based index)."""
        if 1 <= muscle <= self.Nmuscles:
            return self.V_muscle[muscle][1] # Access row 'muscle', column 1 (Dorsal)
        else:
            print(f"Warning: DorsalMuscleOutput index {muscle} out of bounds [1, {self.Nmuscles}]")
            return 0.0 # Default value

    def VentralMuscleOutput(self, muscle: int) -> float:
        """Gets the current activation output of a specific ventral muscle (1-based index)."""
        if 1 <= muscle <= self.Nmuscles:
            return self.V_muscle[muscle][2] # Access row 'muscle', column 2 (Ventral)
        else:
            print(f"Warning: VentralMuscleOutput index {muscle} out of bounds [1, {self.Nmuscles}]")
            return 0.0 # Default value

    def EulerStep(self, StepSize: float):
        """
        Performs one first-order Euler integration step for all muscle activations.
        V_muscle += StepSize * (V_input - V_muscle) / T_muscle
        """
        if self.T_muscle <= 0:
            print("Warning: T_muscle is zero or negative in Muscles.EulerStep. Skipping update.")
            return

        # Loop through muscles (1-based index)
        for i in range(1, self.Nmuscles + 1):
            # Dorsal muscle update (column 1)
            current_act_d = self.V_muscle[i][1]
            input_d = self.V_input[i][1]
            delta_act_d = StepSize * ((input_d - current_act_d) / self.T_muscle)
            self.V_muscle[i][1] = current_act_d + delta_act_d

            # Ventral muscle update (column 2)
            current_act_v = self.V_muscle[i][2]
            input_v = self.V_input[i][2]
            delta_act_v = StepSize * ((input_v - current_act_v) / self.T_muscle)
            self.V_muscle[i][2] = current_act_v + delta_act_v

            # Optional: Clip activation to a valid range (e.g., [0, 1] or [-1, 1] depending on model)
            # self.V_muscle[i][1] = max(0.0, min(self.V_muscle[i][1], 1.0))
            # self.V_muscle[i][2] = max(0.0, min(self.V_muscle[i][2], 1.0))