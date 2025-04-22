import math
import numpy as np
import sys

# --- Configuration Flags ---
BBC_STRICT = False # Set to True to match BBC code exactly where NOTES indicate differences

# --- Settable Constants (from C++) ---
MEDIUM = 0.0                        # Normalized medium drag coefficient (0=water, 1=agar)
L_WORM = 1.0e-3                     # Length of worm in m
N_SEGMENTS = 50                     # Number of segments
R_MIN = 40.0e-6                     # Minor radius of prolate ellipse body in m
C_AGAR_PAR_TOTAL = 3.2e-3           # Total tangential drag coefficient for agar in kg/s
C_AGAR_PERP_TOTAL = 128e-3          # Total rod normal drag coefficient in agar in kg/s
C_WATER_PAR_TOTAL = 3.3e-6          # Total rod tangential drag coefficient for water in kg/s
C_WATER_PERP_TOTAL = 5.2e-6         # Total rod normal drag coefficient for water in kg/s
KAPPA_L = (10.0e-3 * N_SEGMENTS) / 24 # Lateral spring constant in kg/s^2 (Force = k*x, units seem off in C++ comment?) - Assuming kg/s^2 based on usage
KAPPA_D = 350 * KAPPA_L             # Diagonal spring constant in kg/s^2
KAPPA_M0 = 0 * KAPPA_L             # Baseline active muscle spring constant in kg/s^2 default 20
BETA_L = 0.025 * KAPPA_L            # Lateral passive damping constant (Force = b*v) -> kg/s
BETA_D = 0.01 * KAPPA_D             # Diagonal passive damping constant in kg/s
BETA_M0 = 100 * BETA_L              # Baseline active damping constant in kg/s
DELTA_M = 0.65                      # Rest muscle length scaling constant

# --- Derived Constants ---
N_RODS = N_SEGMENTS + 1             # Number of rods (mass points)
N_STATES = 3 * N_RODS               # Total number of states (x, y, phi per rod)
L_SEG = L_WORM / N_SEGMENTS         # Length of an individual segment in m
D_MIN = 2 * R_MIN                   # Minor diameter of prolate ellipse body in m

# Handle BBC_STRICT for drag calculation (NOTE 1)
if BBC_STRICT:
    DENOMINATOR = (2 * N_SEGMENTS + 1)
else:
    DENOMINATOR = (2 * (N_SEGMENTS + 1)) # Fix potential C macro bug expansion

if DENOMINATOR == 0: raise ValueError("N_SEGMENTS leads to zero denominator")

C_AGAR_PAR = C_AGAR_PAR_TOTAL / DENOMINATOR
C_AGAR_PERP = C_AGAR_PERP_TOTAL / DENOMINATOR
C_WATER_PAR = C_WATER_PAR_TOTAL / DENOMINATOR
C_WATER_PERP = C_WATER_PERP_TOTAL / DENOMINATOR

C_PAR = (C_AGAR_PAR - C_WATER_PAR) * MEDIUM + C_WATER_PAR    # Per rod tangential drag coefficient in kg/s
C_PERP = (C_AGAR_PERP - C_WATER_PERP) * MEDIUM + C_WATER_PERP # Per rod normal drag coefficient in kg/s

# --- Global Arrays (Initialized by InitializeBodyConstants) ---
# These store per-rod or per-segment geometric constants
R = np.zeros(N_RODS)                  # Rod radii in m
L_D0 = np.zeros(N_SEGMENTS)           # Rest length of each diagonal element in m
L_L0 = np.zeros(N_SEGMENTS)           # Rest length of each lateral element in m
L_min = np.zeros(N_SEGMENTS)          # Minimal length of each lateral muscle in m
L_L0_minus_L_min = np.zeros(N_SEGMENTS) # Precomputed difference

# --- Global Initialization Function ---
_body_constants_initialized = False

def InitializeBodyConstants():
    """
    Initializes global arrays R, L_D0, L_L0, L_min based on constants.
    Must be called once before creating any WormBody instances.
    """
    global R, L_D0, L_L0, L_min, L_L0_minus_L_min, _body_constants_initialized
    if _body_constants_initialized:
        # print("Warning: Body constants already initialized.")
        return

    NS2 = N_SEGMENTS / 2.0
    LS2 = L_SEG * L_SEG
    for i in range(N_RODS):
        # Ensure argument for acos is within [-1, 1] due to potential float issues
        acos_arg = max(-1.0, min(1.0, (i - NS2) / (NS2 + 0.2)))
        # Calculate radius based on prolate ellipse shape (from BBC code logic)
        R[i] = R_MIN * abs(math.sin(math.acos(acos_arg)))

    for i in range(N_SEGMENTS):
        r_diff = R[i] - R[i+1]
        L_L0[i] = math.sqrt(LS2 + r_diff * r_diff) # Lateral rest length

        r_sum = R[i] + R[i+1]
        L_D0[i] = math.sqrt(LS2 + r_sum * r_sum) # Diagonal rest length

        if D_MIN == 0: raise ValueError("D_MIN cannot be zero")
        L_min[i] = L_L0[i] * (1 - DELTA_M * r_sum / D_MIN) # Minimal muscle length

        L_L0_minus_L_min[i] = L_L0[i] - L_min[i] # Precompute difference

    _body_constants_initialized = True
    print("Body constants initialized.")


# --- The WormBody Class ---
class WormBody:
    def __init__(self):
        if not _body_constants_initialized:
            raise RuntimeError("WormBody constants not initialized. Call InitializeBodyConstants() first.")

        # --- Instance variables ---
        self.t: float = 0.0

        # State vectors (position, velocity/derivative, residual)
        self.Z = np.zeros(N_STATES)       # State: [x0,y0,phi0, x1,y1,phi1, ...]
        self.dZ = np.zeros(N_STATES)      # State derivative (velocity)
        self.Residuals = np.zeros(N_STATES)# Residual vector F(Z, dZ) for DAE

        # Muscle activations (per segment)
        self.A_D_M = np.zeros(N_SEGMENTS) # Dorsal muscle activation
        self.A_V_M = np.zeros(N_SEGMENTS) # Ventral muscle activation

        # Intermediate kinematic arrays (per rod or segment)
        self.sinPhi = np.zeros(N_RODS)
        self.cosPhi = np.zeros(N_RODS)

        self.L_D_D = np.zeros(N_SEGMENTS) # Length of Dorsal Diagonal elements
        self.dL_D_D = np.zeros(N_SEGMENTS)# Rate of change of L_D_D
        self.L_V_D = np.zeros(N_SEGMENTS) # Length of Ventral Diagonal elements
        self.dL_V_D = np.zeros(N_SEGMENTS)# Rate of change of L_V_D
        self.uD_D_x = np.zeros(N_SEGMENTS)# Unit vector x-comp for Dorsal Diagonal
        self.uD_D_y = np.zeros(N_SEGMENTS)# Unit vector y-comp for Dorsal Diagonal
        self.uD_V_x = np.zeros(N_SEGMENTS)# Unit vector x-comp for Ventral Diagonal
        self.uD_V_y = np.zeros(N_SEGMENTS)# Unit vector y-comp for Ventral Diagonal

        self.L_D_L = np.zeros(N_SEGMENTS) # Length of Dorsal Lateral elements
        self.dL_D_L = np.zeros(N_SEGMENTS)# Rate of change of L_D_L
        self.L_V_L = np.zeros(N_SEGMENTS) # Length of Ventral Lateral elements
        self.dL_V_L = np.zeros(N_SEGMENTS)# Rate of change of L_V_L
        self.uL_D_x = np.zeros(N_SEGMENTS)# Unit vector x-comp for Dorsal Lateral
        self.uL_D_y = np.zeros(N_SEGMENTS)# Unit vector y-comp for Dorsal Lateral
        self.uL_V_x = np.zeros(N_SEGMENTS)# Unit vector x-comp for Ventral Lateral
        self.uL_V_y = np.zeros(N_SEGMENTS)# Unit vector y-comp for Ventral Lateral

        # Intermediate force arrays (per rod)
        self.f_D_x = np.zeros(N_RODS)     # Total force on Dorsal endpoint (x)
        self.f_D_y = np.zeros(N_RODS)     # Total force on Dorsal endpoint (y)
        self.f_V_x = np.zeros(N_RODS)     # Total force on Ventral endpoint (x)
        self.f_V_y = np.zeros(N_RODS)     # Total force on Ventral endpoint (y)

        # Temporary arrays used within calculations (avoids re-allocation)
        self._p_D_x = np.zeros(N_RODS)
        self._p_D_y = np.zeros(N_RODS)
        self._p_V_x = np.zeros(N_RODS)
        self._p_V_y = np.zeros(N_RODS)
        self._dp_D_x = np.zeros(N_RODS)
        self._dp_D_y = np.zeros(N_RODS)
        self._dp_V_x = np.zeros(N_RODS)
        self._dp_V_y = np.zeros(N_RODS)

        self._f_D_D_x = np.zeros(N_SEGMENTS)
        self._f_D_D_y = np.zeros(N_SEGMENTS)
        self._f_V_D_x = np.zeros(N_SEGMENTS)
        self._f_V_D_y = np.zeros(N_SEGMENTS)
        self._f_D_L_x = np.zeros(N_SEGMENTS)
        self._f_D_L_y = np.zeros(N_SEGMENTS)
        self._f_V_L_x = np.zeros(N_SEGMENTS)
        self._f_V_L_y = np.zeros(N_SEGMENTS)
        self._f_D_M_x = np.zeros(N_SEGMENTS)
        self._f_D_M_y = np.zeros(N_SEGMENTS)
        self._f_V_M_x = np.zeros(N_SEGMENTS)
        self._f_V_M_y = np.zeros(N_SEGMENTS)

    # --- Accessors (using 1-based indexing for external interface) ---
    def time(self) -> float:
        return self.t

    def X(self, i: int) -> float:
        """Returns X coordinate of rod i (1-based)."""
        if 1 <= i <= N_RODS:
            return self.Z[3 * (i - 1)]
        else:
            raise IndexError(f"X(i): Rod index {i} out of bounds [1, {N_RODS}]")

    def Y(self, i: int) -> float:
        """Returns Y coordinate of rod i (1-based)."""
        if 1 <= i <= N_RODS:
            return self.Z[3 * (i - 1) + 1]
        else:
            raise IndexError(f"Y(i): Rod index {i} out of bounds [1, {N_RODS}]")

    def Phi(self, i: int) -> float:
        """Returns angle Phi of rod i (1-based)."""
        if 1 <= i <= N_RODS:
            return self.Z[3 * (i - 1) + 2]
        else:
            raise IndexError(f"Phi(i): Rod index {i} out of bounds [1, {N_RODS}]")

    def SetDorsalSegmentActivation(self, i: int, a: float):
        """Sets dorsal muscle activation for segment i (1-based). Clips to [0, 1]."""
        if 1 <= i <= N_SEGMENTS:
            self.A_D_M[i - 1] = max(0.0, min(a, 1.0))
        else:
            raise IndexError(f"SetDorsalSegmentActivation: Segment index {i} out of bounds [1, {N_SEGMENTS}]")

    def SetVentralSegmentActivation(self, i: int, a: float):
        """Sets ventral muscle activation for segment i (1-based). Clips to [0, 1]."""
        if 1 <= i <= N_SEGMENTS:
            self.A_V_M[i - 1] = max(0.0, min(a, 1.0))
        else:
            raise IndexError(f"SetVentralSegmentActivation: Segment index {i} out of bounds [1, {N_SEGMENTS}]")

    def DorsalSegmentLength(self, i: int) -> float:
        """Gets current dorsal lateral element length for segment i (1-based)."""
        if 1 <= i <= N_SEGMENTS:
            return self.L_D_L[i - 1]
        else:
            raise IndexError(f"DorsalSegmentLength: Segment index {i} out of bounds [1, {N_SEGMENTS}]")

    def VentralSegmentLength(self, i: int) -> float:
        """Gets current ventral lateral element length for segment i (1-based)."""
        if 1 <= i <= N_SEGMENTS:
            return self.L_V_L[i - 1]
        else:
            raise IndexError(f"VentralSegmentLength: Segment index {i} out of bounds [1, {N_SEGMENTS}]")

    def RestingLength(self, i: int) -> float:
        """Gets resting lateral element length for segment i (1-based)."""
        if 1 <= i <= N_SEGMENTS:
            # Access the global L_L0 array (0-based)
            return L_L0[i - 1]
        else:
            raise IndexError(f"RestingLength: Segment index {i} out of bounds [1, {N_SEGMENTS}]")

    # --- Control Methods ---
    def InitializeBodyState(self):
        """
        Initializes the state of the body to a straight line along X-axis.
        Sets initial dZ to zero, as required by the integrator.
        """
        self.t = 0.0
        for i in range(N_RODS): # 0-based loop for internal arrays
            i3 = 3 * i
            self.Z[i3] = i * L_SEG   # x coordinate
            self.Z[i3 + 1] = 0.0     # y coordinate
            self.Z[i3 + 2] = math.pi / 2.0 # phi angle (rods perpendicular to body axis)
            self.dZ[i3] = 0.0
            self.dZ[i3 + 1] = 0.0
            self.dZ[i3 + 2] = 0.0

        # Reset muscle activations
        self.A_D_M.fill(0.0)
        self.A_V_M.fill(0.0)

        # Update kinematics based on initial Z (dZ is 0)
        # Need to ensure dZ is 0 *before* calling kinematics if it uses dZ
        self.dZ.fill(0.0)
        self._UpdateKinematics()


    def StepBody(self, h: float):
        """Advances the body simulation by one time step h."""
        self.A_D_M = np.zeros(self.A_D_M.size)  # debug
        self.A_V_M = np.zeros(self.A_V_M.size)  # debug
        self._SemiImplicitBackwardEulerDAEStep(h)

    # --- Internal Calculation Methods ---
    def _F(self, start_rod: int = 0, end_rod: int = N_RODS):
        """Helper to update kinematics, forces, and residuals for a range of rods."""
        self._UpdateKinematics(start_rod, end_rod)
        self._UpdateForces(start_rod, end_rod)
        self._UpdateResiduals(start_rod, end_rod)

    def _UpdateKinematics(self, start_rod: int = 0, end_rod: int = N_RODS):
        """
        Updates kinematic variables (lengths, velocities, unit vectors)
        based on the current state Z and its derivative dZ.
        Operates on the range [start_rod, end_rod) for rods (0-based).
        """
        # Determine the range of segments affected
        start_seg = max(0, start_rod - 1)
        end_seg = min(N_SEGMENTS, end_rod) # Need kinematics up to segment end_rod-1

        # Use temporary arrays stored as instance variables
        p_D_x, p_D_y = self._p_D_x, self._p_D_y
        p_V_x, p_V_y = self._p_V_x, self._p_V_y
        dp_D_x, dp_D_y = self._dp_D_x, self._dp_D_y
        dp_V_x, dp_V_y = self._dp_V_x, self._dp_V_y

        # 1. Update rod endpoint positions (p) and velocities (dp)
        # Need data for rods from start_seg to end_seg (inclusive)
        rod_loop_start = start_seg
        rod_loop_end = min(N_RODS, end_seg + 1)

        for i in range(rod_loop_start, rod_loop_end):
            i3 = 3 * i
            x, y, phi = self.Z[i3], self.Z[i3 + 1], self.Z[i3 + 2]
            dx, dy, dphi = self.dZ[i3], self.dZ[i3 + 1], self.dZ[i3 + 2]

            # Pre-calculate sin/cos
            self.sinPhi[i] = math.sin(phi)
            self.cosPhi[i] = math.cos(phi)

            # Rod endpoints (Dorsal/Ventral)
            r_cosPhi = R[i] * self.cosPhi[i]
            r_sinPhi = R[i] * self.sinPhi[i]
            p_D_x[i] = x + r_cosPhi
            p_D_y[i] = y + r_sinPhi
            p_V_x[i] = x - r_cosPhi
            p_V_y[i] = y - r_sinPhi

            # Endpoint velocities
            arm = R[i] * dphi
            dxtemp = -arm * self.sinPhi[i] # Component due to rotation
            dytemp = arm * self.cosPhi[i]  # Component due to rotation
            dp_D_x[i] = dx + dxtemp
            dp_D_y[i] = dy + dytemp
            dp_V_x[i] = dx - dxtemp
            dp_V_y[i] = dy - dytemp

        # 2. Update element lengths (L), unit vectors (u), and length rates (dL)
        # Loop over affected segments [start_seg, end_seg)
        for i in range(start_seg, end_seg):
            i1 = i + 1 # Index of the next rod

            # --- Diagonal elements (Dorsal rod i to Ventral rod i+1, and vice versa) ---
            # D_D: Dorsal(i) -> Ventral(i+1)
            D_x = p_V_x[i1] - p_D_x[i]
            D_y = p_V_y[i1] - p_D_y[i]
            self.L_D_D[i] = math.sqrt(D_x * D_x + D_y * D_y)
            if self.L_D_D[i] < 1e-12: # Avoid division by zero
                self.uD_D_x[i] = 1.0; self.uD_D_y[i] = 0.0 # Arbitrary unit vector
            else:
                inv_L_D_D = 1.0 / self.L_D_D[i]
                self.uD_D_x[i] = D_x * inv_L_D_D
                self.uD_D_y[i] = D_y * inv_L_D_D
            # Rate of change (velocity projection)
            self.dL_D_D[i] = ((dp_V_x[i1] - dp_D_x[i]) * self.uD_D_x[i] +
                              (dp_V_y[i1] - dp_D_y[i]) * self.uD_D_y[i])

            # V_D: Ventral(i) -> Dorsal(i+1)
            V_x = p_D_x[i1] - p_V_x[i]
            V_y = p_D_y[i1] - p_V_y[i]
            self.L_V_D[i] = math.sqrt(V_x * V_x + V_y * V_y)
            if self.L_V_D[i] < 1e-12:
                 self.uD_V_x[i] = 1.0; self.uD_V_y[i] = 0.0
            else:
                inv_L_V_D = 1.0 / self.L_V_D[i]
                self.uD_V_x[i] = V_x * inv_L_V_D
                self.uD_V_y[i] = V_y * inv_L_V_D
            self.dL_V_D[i] = ((dp_D_x[i1] - dp_V_x[i]) * self.uD_V_x[i] +
                              (dp_D_y[i1] - dp_V_y[i]) * self.uD_V_y[i])

            # --- Lateral elements (Dorsal rod i to Dorsal rod i+1, etc.) ---
            # D_L: Dorsal(i) -> Dorsal(i+1)
            D_x = p_D_x[i1] - p_D_x[i]
            D_y = p_D_y[i1] - p_D_y[i]
            self.L_D_L[i] = math.sqrt(D_x * D_x + D_y * D_y)
            if self.L_D_L[i] < 1e-12:
                 self.uL_D_x[i] = 1.0; self.uL_D_y[i] = 0.0
            else:
                inv_L_D_L = 1.0 / self.L_D_L[i]
                self.uL_D_x[i] = D_x * inv_L_D_L
                self.uL_D_y[i] = D_y * inv_L_D_L
            self.dL_D_L[i] = ((dp_D_x[i1] - dp_D_x[i]) * self.uL_D_x[i] +
                              (dp_D_y[i1] - dp_D_y[i]) * self.uL_D_y[i])

            # V_L: Ventral(i) -> Ventral(i+1)
            V_x = p_V_x[i1] - p_V_x[i]
            V_y = p_V_y[i1] - p_V_y[i]
            self.L_V_L[i] = math.sqrt(V_x * V_x + V_y * V_y)
            if self.L_V_L[i] < 1e-12:
                 self.uL_V_x[i] = 1.0; self.uL_V_y[i] = 0.0
            else:
                inv_L_V_L = 1.0 / self.L_V_L[i]
                self.uL_V_x[i] = V_x * inv_L_V_L
                self.uL_V_y[i] = V_y * inv_L_V_L
            self.dL_V_L[i] = ((dp_V_x[i1] - dp_V_x[i]) * self.uL_V_x[i] +
                              (dp_V_y[i1] - dp_V_y[i]) * self.uL_V_y[i])


    def _UpdateForces(self, start_rod: int = 0, end_rod: int = N_RODS):
        """
        Updates the internal forces based on kinematics and muscle activation.
        Operates on the range [start_rod, end_rod) for rods (0-based).
        """
        # Determine the range of segments affected
        start_seg = max(0, start_rod - 1)
        end_seg = min(N_SEGMENTS, end_rod) # Need forces from segment start_rod-1 to end_rod-1

        # Use temporary arrays stored as instance variables
        f_D_D_x, f_D_D_y = self._f_D_D_x, self._f_D_D_y
        f_V_D_x, f_V_D_y = self._f_V_D_x, self._f_V_D_y
        f_D_L_x, f_D_L_y = self._f_D_L_x, self._f_D_L_y
        f_V_L_x, f_V_L_y = self._f_V_L_x, self._f_V_L_y
        f_D_M_x, f_D_M_y = self._f_D_M_x, self._f_D_M_y
        f_V_M_x, f_V_M_y = self._f_V_M_x, self._f_V_M_y


        # 1. Update individual element forces (springs and dampers)
        # Loop over affected segments [start_seg, end_seg)
        for i in range(start_seg, end_seg):
            # --- Diagonal passive forces (spring + damper) ---
            # D_D element: Dorsal(i) -> Ventral(i+1)
            force_mag = KAPPA_D * (L_D0[i] - self.L_D_D[i]) - BETA_D * self.dL_D_D[i]
            f_D_D_x[i] = force_mag * self.uD_D_x[i]
            f_D_D_y[i] = force_mag * self.uD_D_y[i]
            # V_D element: Ventral(i) -> Dorsal(i+1)
            force_mag = KAPPA_D * (L_D0[i] - self.L_V_D[i]) - BETA_D * self.dL_V_D[i]
            f_V_D_x[i] = force_mag * self.uD_V_x[i]
            f_V_D_y[i] = force_mag * self.uD_V_y[i]

            # --- Lateral passive forces (nonlinear spring + damper) ---
            # D_L element: Dorsal(i) -> Dorsal(i+1)
            stretch = L_L0[i] - self.L_D_L[i]
            # NOTE 2: Nonlinear term for compression (matches BBC code)
            if stretch < 0:
                stretch4 = stretch * stretch * stretch * stretch
                stretch += 16.0 * stretch4 # Adds significant stiffness in compression
            force_mag = KAPPA_L * stretch - BETA_L * self.dL_D_L[i]
            f_D_L_x[i] = force_mag * self.uL_D_x[i]
            f_D_L_y[i] = force_mag * self.uL_D_y[i]
            # V_L element: Ventral(i) -> Ventral(i+1)
            stretch = L_L0[i] - self.L_V_L[i]
            if stretch < 0:
                stretch4 = stretch * stretch * stretch * stretch
                stretch += 16.0 * stretch4
            force_mag = KAPPA_L * stretch - BETA_L * self.dL_V_L[i]
            f_V_L_x[i] = force_mag * self.uL_V_x[i]
            f_V_L_y[i] = force_mag * self.uL_V_y[i]

            # --- Lateral active muscle forces (spring + damper, modulated by activation) ---
            # D_M element (alongside D_L): Dorsal(i) -> Dorsal(i+1)
            act_d = self.A_D_M[i]
            target_len_d = L_L0[i] - act_d * L_L0_minus_L_min[i] # Target length decreases with activation
            force_mag = (KAPPA_M0 * act_d * (target_len_d - self.L_D_L[i]) -
                         BETA_M0 * act_d * self.dL_D_L[i])
            f_D_M_x[i] = force_mag * self.uL_D_x[i]
            f_D_M_y[i] = force_mag * self.uL_D_y[i]
            # V_M element (alongside V_L): Ventral(i) -> Ventral(i+1)
            act_v = self.A_V_M[i]
            target_len_v = L_L0[i] - act_v * L_L0_minus_L_min[i]
            force_mag = (KAPPA_M0 * act_v * (target_len_v - self.L_V_L[i]) -
                         BETA_M0 * act_v * self.dL_V_L[i])
            f_V_M_x[i] = force_mag * self.uL_V_x[i]
            f_V_M_y[i] = force_mag * self.uL_V_y[i]


        # 2. Update total forces on rod endpoints by summing element forces
        # Loop over affected rods [start_rod, end_rod)
        # Need to handle boundary conditions carefully
        for i in range(start_rod, end_rod):
            f_d_x_i = 0.0; f_d_y_i = 0.0
            f_v_x_i = 0.0; f_v_y_i = 0.0

            # Forces from element i (connecting rod i to rod i+1)
            if i < N_SEGMENTS: # Check if element i exists
                f_d_x_i += -f_D_D_x[i] - f_D_L_x[i] - f_D_M_x[i] # Force on Dorsal(i) from element i
                f_d_y_i += -f_D_D_y[i] - f_D_L_y[i] - f_D_M_y[i]
                f_v_x_i += -f_V_D_x[i] - f_V_L_x[i] - f_V_M_x[i] # Force on Ventral(i) from element i
                f_v_y_i += -f_V_D_y[i] - f_V_L_y[i] - f_V_M_y[i]

            # Forces from element i-1 (connecting rod i-1 to rod i)
            if i > 0: # Check if element i-1 exists
                i_minus_1 = i - 1
                f_d_x_i += f_V_D_x[i_minus_1] + f_D_L_x[i_minus_1] + f_D_M_x[i_minus_1] # Force on Dorsal(i) from element i-1
                f_d_y_i += f_V_D_y[i_minus_1] + f_D_L_y[i_minus_1] + f_D_M_y[i_minus_1]
                f_v_x_i += f_D_D_x[i_minus_1] + f_V_L_x[i_minus_1] + f_V_M_x[i_minus_1] # Force on Ventral(i) from element i-1
                f_v_y_i += f_D_D_y[i_minus_1] + f_V_L_y[i_minus_1] + f_V_M_y[i_minus_1]

            self.f_D_x[i] = f_d_x_i
            self.f_D_y[i] = f_d_y_i
            self.f_V_x[i] = f_v_x_i
            self.f_V_y[i] = f_v_y_i


    def _UpdateResiduals(self, start_rod: int = 0, end_rod: int = N_RODS):
        """
        Updates the DAE residuals F(Z, dZ) based on current forces and state derivatives.
        Residual = dZ - calculated_velocity. Should be zero if Z, dZ satisfy the equations.
        Operates on the range [start_rod, end_rod) for rods (0-based).
        """
        if C_PERP == 0 or C_PAR == 0:
             raise ValueError("Drag coefficients C_PERP or C_PAR are zero.")

        inv_C_PERP = 1.0 / C_PERP
        inv_C_PAR = 1.0 / C_PAR

        for i in range(start_rod, end_rod):
            i3 = 3 * i

            # Decompose forces on Dorsal/Ventral endpoints into Parallel/Perpendicular
            # to the rod's orientation (phi).
            # Parallel is along axis from Ventral to Dorsal endpoint.
            # Perpendicular is 90 degrees counter-clockwise from parallel.
            # Note: BBC paper/code definitions seem slightly different, careful comparison needed.
            # Let's follow the C++ code's decomposition:
            # Par = Fx*sin(Phi) - Fy*cos(Phi)
            # Perp = Fx*cos(Phi) + Fy*sin(Phi)

            sinP = self.sinPhi[i]
            cosP = self.cosPhi[i]

            f_D_par = self.f_D_x[i] * sinP - self.f_D_y[i] * cosP
            f_D_perp = self.f_D_x[i] * cosP + self.f_D_y[i] * sinP
            f_V_par = self.f_V_x[i] * sinP - self.f_V_y[i] * cosP
            f_V_perp = self.f_V_x[i] * cosP + self.f_V_y[i] * sinP

            # Even and Odd components of parallel forces
            f_even_par = (f_V_par + f_D_par) / 2.0 # Related to translation
            f_odd_par = (f_V_par - f_D_par) / 2.0  # Related to rotation

            # Calculate velocities based on forces and drag coefficients
            # Velocity of Center of Mass (CoM) in Perp/Par frame
            V_CoM_perp = (f_D_perp + f_V_perp) * inv_C_PERP
            V_CoM_par = 2.0 * f_even_par * inv_C_PAR # Factor of 2 from BBC

            # Convert CoM velocity back to Lab frame (X, Y)
            calc_dx = V_CoM_par * sinP + V_CoM_perp * cosP
            calc_dy = -V_CoM_par * cosP + V_CoM_perp * sinP

            # Calculate angular velocity (dphi/dt)
            # NOTE 3: Discrepancy between BBC code and paper. Using paper/thesis version unless BBC_STRICT.
            if R[i] == 0: raise ValueError(f"Rod radius R[{i}] is zero.")
            if BBC_STRICT:
                # Formula from BBC code (units seem problematic: Force/Drag / Length?)
                if math.pi * 2.0 * R[i] == 0: raise ValueError("Denominator zero in BBC_STRICT dphi calc")
                calc_dphi = (f_odd_par * inv_C_PAR) / (math.pi * 2.0 * R[i])
            else:
                # Formula from BBC paper/thesis (units seem more plausible: Force / (Radius*Drag) -> Torque / (Moment*DragCoeff) -> AngVel)
                calc_dphi = (2.0 * f_odd_par) / (R[i] * C_PAR) # Effectively inv_C_PAR included

            # Calculate Residuals: F(Z, dZ) = dZ - calculated_velocity
            self.Residuals[i3] = self.dZ[i3] - calc_dx
            self.Residuals[i3 + 1] = self.dZ[i3 + 1] - calc_dy
            self.Residuals[i3 + 2] = self.dZ[i3 + 2] - calc_dphi

    def _SemiImplicitBackwardEulerDAEStep(self, h: float):
        """
        Performs one step using Semi-Implicit Backward Euler for DAEs.
        Solves M * deltaZ = F(Z_n, dZ_n) where M = -[(1/h)*dF/dZ' + dF/dZ].
        Assumes Z, dZ, Residuals are current at time t. Updates Z to time t+h.
        """
        if h <= 0:
            raise ValueError("Time step h must be positive")

        # Allocate Jacobian matrices (could be reused if stored on self)
        dFdZ = np.zeros((N_STATES, N_STATES))
        dFdZp = np.zeros((N_STATES, N_STATES)) # Stores dF/dZ'

        # 1. Ensure current residual F(Z_n, dZ_n) is up-to-date
        #    (Assumes kinematics is already current from previous step or init)
        self._UpdateForces()  # Depends on Z, dZ -> Kinematics
        self._UpdateResiduals() # Depends on Forces, dZ

        # Store current residual F(Z_n, dZ_n) before calculating Jacobians
        current_residuals = self.Residuals.copy()

        # 2. Compute Jacobians numerically around (Z_n, dZ_n)
        self._NumericaldFdZ(dFdZ)        # Calculates dF/dZ
        self._NumericaldFdZp(dFdZp)      # Calculates dF/dZ' (stored in dFdZp)

        # 3. Form the matrix M = -[(1/h)*dF/dZ' + dF/dZ]
        h_inv = 1.0 / h
        M = -(h_inv * dFdZp + dFdZ)

        # 4. Solve the linear system M * deltaZ = F(Z_n, dZ_n) for deltaZ
        try:
            deltaZ = np.linalg.solve(M, current_residuals)
        except np.linalg.LinAlgError:
            print("Error: Linear system solve failed (matrix M may be singular).", file=sys.stderr)
            # Handle error: Maybe take smaller step, or stop simulation
            # For now, let's just print Z and M and raise exception
            print("Current State Z:\n", self.Z)
            print("Matrix M:\n", M)
            raise RuntimeError("Singular matrix encountered in DAE solver.")

        # 5. Update the state: Z_{n+1} = Z_n + deltaZ
        self.Z += deltaZ
        # Update dZ using finite difference: dZ_{n+1} ~ deltaZ / h
        # (Or keep dZ as the value used in residual calc? The method name implies Z update only.)
        # The C++ code doesn't explicitly update dZ here, it seems dZ for the *next* step
        # is implicitly determined by the updated Z satisfying the DAE.
        # Let's recalculate dZ based on the new Z satisfying the residual F(Z_{n+1}, dZ_{n+1}) = 0
        # We can approximate dZ_{n+1} = calculated_velocity(Z_{n+1})
        # First update kinematics based on new Z, assuming dZ is approx deltaZ/h for velocity calc
        self.dZ = deltaZ / h # Update dZ estimate for next step's kinematics/forces
        self._UpdateKinematics() # Bring kinematics up-to-date with new Z and estimated new dZ

        # 6. Increment time
        self.t += h

    def _NumericaldFdZ(self, J: np.ndarray):
        """
        Numerically estimate the residual Jacobian matrix dF/dZ at the current state (Z, dZ).
        Stores result in J. Assumes Residuals vector is current.
        """
        eps = np.sqrt(np.finfo(float).eps) # Machine epsilon for float
        current_residuals = self.Residuals.copy()
        original_Z = self.Z.copy()

        for j in range(N_STATES):
            temp_Zj = original_Z[j]
            h = eps * abs(temp_Zj)
            if h == 0.0:
                h = eps

            # Perturb Z[j]
            self.Z[j] = temp_Zj + h
            h = self.Z[j] - temp_Zj # Adjust h for exact representability

            # Recalculate residuals with perturbed Z
            # Need to update kinematics, forces, residuals for the affected range
            # Optimisation: Only recalc affected range (tricky, do full recalc first)
            self._F() # Recalculate all residuals based on perturbed Z[j]

            # Compute finite difference for column j
            J[:, j] = (self.Residuals - current_residuals) / h

            # Restore original Z[j] and residuals for next iteration
            self.Z[j] = temp_Zj
            self.Residuals[:] = current_residuals # Restore full vector

        # Restore original Z completely just in case
        self.Z[:] = original_Z


    def _NumericaldFdZp(self, J: np.ndarray):
        """
        Numerically estimate the residual Jacobian matrix dF/dZ' (dZ prime)
        at the current state (Z, dZ). Stores result in J. Assumes Residuals is current.
        """
        eps = np.sqrt(np.finfo(float).eps)
        current_residuals = self.Residuals.copy()
        original_dZ = self.dZ.copy()
        # Need original Z to restore state after _F calls if _F modifies Z implicitly
        original_Z = self.Z.copy()

        for j in range(N_STATES):
            temp_dZj = original_dZ[j]
            h = eps * abs(temp_dZj)
            if h == 0.0:
                h = eps

            # Perturb dZ[j]
            self.dZ[j] = temp_dZj + h
            h = self.dZ[j] - temp_dZj

            # Recalculate residuals with perturbed dZ
            # _F updates kinematics (uses dZ), forces (uses kinematics), residuals (uses forces, dZ)
            self._F()

            # Compute finite difference for column j
            J[:, j] = (self.Residuals - current_residuals) / h

            # Restore original dZ[j] and residuals
            self.dZ[j] = temp_dZj
            self.Residuals[:] = current_residuals

        # Restore original dZ and Z completely
        self.dZ[:] = original_dZ
        self.Z[:] = original_Z

# --- End of WormBody Class ---

# Example of how to initialize constants (call this from main.py)
# if __name__ == "__main__":
#     InitializeBodyConstants()
#     # Create instance only after initializing constants
#     body = WormBody()
#     body.InitializeBodyState()
#     print("WormBody created and initialized.")
#     print("Rod 0 X:", body.X(1)) # Test 1-based access