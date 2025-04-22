import math
from vector_matrix import TVector # Assuming this is in the same directory

class StretchReceptor:
    """
    Python port of the StretchReceptor class.
    Calculates stretch receptor output based on normalized segment lengths.
    """

    def __init__(self, nSegs: int = 50, nSR: int = 10, ASRgain: float = 0.0, BSRgain: float = 0.0):
        """
        Initializes the StretchReceptor instance.
        Args:
            nSegs: Number of segments in the worm body.
            nSR: Number of stretch receptor units (typically equals N_UNITS).
            ASRgain: Gain for the A-class stretch receptors.
            BSRgain: Gain for the B-class stretch receptors.
        """
        # Initialize attributes with default values first in case SetParams isn't called immediately
        self.NSEGS: int = 0
        self.NSR: int = 0
        self.NSEGSSR: int = 6 # Number of segments that feed into one SR unit (fixed in C++)
        self.SR_A_gain: float = 0.0
        self.SR_B_gain: float = 0.0

        self.normSegLenD: TVector[float] = TVector(1, 0) # Normalized Dorsal Segment Lengths
        self.normSegLenV: TVector[float] = TVector(1, 0) # Normalized Ventral Segment Lengths

        self.A_D_sr: TVector[float] = TVector(1, 0) # Output of A-class Dorsal SRs per unit
        self.A_V_sr: TVector[float] = TVector(1, 0) # Output of A-class Ventral SRs per unit
        self.B_D_sr: TVector[float] = TVector(1, 0) # Output of B-class Dorsal SRs per unit
        self.B_V_sr: TVector[float] = TVector(1, 0) # Output of B-class Ventral SRs per unit

        # Call SetStretchReceptorParams to properly initialize based on arguments
        self.SetStretchReceptorParams(nSegs, nSR, ASRgain, BSRgain)

    def SetStretchReceptorParams(self, nSegs: int, nSR: int, ASRgain: float, BSRgain: float):
        """
        Sets the parameters and initializes/resizes the internal vectors.
        """
        self.NSEGS = nSegs          # Number of segments
        self.NSR = nSR              # Number of stretch receptors (units)
        self.NSEGSSR = 6            # Number of segments that go into a stretch receptor (fixed)
        self.SR_A_gain = ASRgain    # Stretch receptor gain for A-class
        self.SR_B_gain = BSRgain    # Stretch receptor gain for B-class

        # Initialize/Resize TVectors
        self.normSegLenD.SetBounds(1, self.NSEGS)
        self.normSegLenV.SetBounds(1, self.NSEGS)
        self.A_D_sr.SetBounds(1, self.NSR)
        self.A_V_sr.SetBounds(1, self.NSR)
        self.B_D_sr.SetBounds(1, self.NSR)
        self.B_V_sr.SetBounds(1, self.NSR)

        # Optional: Fill with default values (e.g., 0.0)
        self.normSegLenD.FillContents(0.0)
        self.normSegLenV.FillContents(0.0)
        self.A_D_sr.FillContents(0.0)
        self.A_V_sr.FillContents(0.0)
        self.B_D_sr.FillContents(0.0)
        self.B_V_sr.FillContents(0.0)


    def SetDorsalInput(self, seg: int, normlen: float):
        """Loads normalized dorsal segment length deformation information."""
        if 1 <= seg <= self.NSEGS:
            self.normSegLenD[seg] = normlen
        else:
            # Optional: Add a warning or error for out-of-bounds access
            print(f"Warning: SetDorsalInput index {seg} out of bounds [1, {self.NSEGS}]")


    def SetVentralInput(self, seg: int, normlen: float):
        """Loads normalized ventral segment length deformation information."""
        if 1 <= seg <= self.NSEGS:
            self.normSegLenV[seg] = normlen
        else:
            print(f"Warning: SetVentralInput index {seg} out of bounds [1, {self.NSEGS}]")

    def Update(self):
        """
        Updates the stretch receptor outputs based on the current segment lengths.
        This implements the specific innervation pattern described in the C++ code.
        """
        if self.NSEGSSR == 0: return # Avoid division by zero if NSEGSSR is not set properly

        # --- A-class Stretch Receptors Calculation (Matches active C++ code) ---
        # Unit 1 (Head): Receives same input as Unit 2 calculation would suggest
        d_sum = 0.0
        v_sum = 0.0
        # Segments 1 to NSEGSSR (e.g., 1 to 6)
        for j in range(1, self.NSEGSSR + 1):
             # Check bounds before accessing
            if 1 <= j <= self.NSEGS:
                d_sum += self.normSegLenD[j]
                v_sum += self.normSegLenV[j]
            else: print(f"Warning SR A1: Segment index {j} out of bounds [1, {self.NSEGS}]")

        if 1 <= self.NSR: # Check if index 1 is valid for output SR vector
            self.A_D_sr[1] = self.SR_A_gain * (d_sum / self.NSEGSSR)
            self.A_V_sr[1] = self.SR_A_gain * (v_sum / self.NSEGSSR)
        else: print("Warning: NSR < 1, cannot assign A_sr[1]")


        # Units 2 to NSR (C++ code explicitly goes to 10, assuming NSR>=10)
        # The loop structure suggests NSR is likely fixed at 10.
        max_a_unit = min(self.NSR, 10) # Ensure we don't exceed NSR bounds if NSR < 10
        for i in range(2, max_a_unit + 1):
            d_sum = 0.0
            v_sum = 0.0
            start_seg_idx_offset = (i - 2) * 4 # Calculate offset: 0 for i=2, 4 for i=3, ...
            for j in range(1, self.NSEGSSR + 1):
                seg_idx = j + start_seg_idx_offset
                # Check bounds before accessing
                if 1 <= seg_idx <= self.NSEGS:
                    d_sum += self.normSegLenD[seg_idx]
                    v_sum += self.normSegLenV[seg_idx]
                else: print(f"Warning SR A{i}: Segment index {seg_idx} out of bounds [1, {self.NSEGS}]")

            # Assign output for unit i
            if 1 <= i <= self.NSR: # Double check index i is valid for output vector
                 self.A_D_sr[i] = self.SR_A_gain * (d_sum / self.NSEGSSR)
                 self.A_V_sr[i] = self.SR_A_gain * (v_sum / self.NSEGSSR)


        # --- B-class Stretch Receptors Calculation (Matches active C++ code) ---
        # Units 1 to NSR-1 (C++ code explicitly goes 1 to 9, assuming NSR>=10)
        max_b_unit_loop1 = min(self.NSR - 1, 9)
        for i in range(1, max_b_unit_loop1 + 1):
            d_sum = 0.0
            v_sum = 0.0
            # Start sensing from segment 13 (offset 12)
            start_seg_idx_offset = 12 + (i - 1) * 4 # Offset: 12 for i=1, 16 for i=2, ...
            for j in range(1, self.NSEGSSR + 1):
                seg_idx = j + start_seg_idx_offset
                # Check bounds before accessing
                if 1 <= seg_idx <= self.NSEGS:
                    d_sum += self.normSegLenD[seg_idx]
                    v_sum += self.normSegLenV[seg_idx]
                else: print(f"Warning SR B{i}: Segment index {seg_idx} out of bounds [1, {self.NSEGS}]")

            # Assign output for unit i
            if 1 <= i <= self.NSR: # Check index i is valid for output vector
                self.B_D_sr[i] = self.SR_B_gain * (d_sum / self.NSEGSSR)
                self.B_V_sr[i] = self.SR_B_gain * (v_sum / self.NSEGSSR)


        # Unit 10 (Tail): Receives same input as Unit 9 calculation would suggest
        # Corresponds to segments starting at index 45 (j + 44)
        if self.NSR >= 10: # Only calculate if NSR is large enough
            target_b_unit = 10
            d_sum = 0.0
            v_sum = 0.0
            start_seg_idx_offset = 44 # Offset for the last group
            for j in range(1, self.NSEGSSR + 1):
                seg_idx = j + start_seg_idx_offset
                # Check bounds before accessing
                if 1 <= seg_idx <= self.NSEGS:
                    d_sum += self.normSegLenD[seg_idx]
                    v_sum += self.normSegLenV[seg_idx]
                else: print(f"Warning SR B{target_b_unit}: Segment index {seg_idx} out of bounds [1, {self.NSEGS}]")

            # Assign output for the last unit (index 10, if NSR allows)
            if 1 <= target_b_unit <= self.NSR:
                self.B_D_sr[target_b_unit] = self.SR_B_gain * (d_sum / self.NSEGSSR)
                self.B_V_sr[target_b_unit] = self.SR_B_gain * (v_sum / self.NSEGSSR)


    # --- Output Accessor Methods (match naming used in Worm.py) ---
    def A_D_sr(self, i: int) -> float:
        """Returns the output of the A-class Dorsal stretch receptor for unit i."""
        if 1 <= i <= self.NSR:
            return self.A_D_sr[i]
        else:
            print(f"Warning: A_D_sr index {i} out of bounds [1, {self.NSR}]")
            return 0.0 # Default value

    def A_V_sr(self, i: int) -> float:
        """Returns the output of the A-class Ventral stretch receptor for unit i."""
        if 1 <= i <= self.NSR:
            return self.A_V_sr[i]
        else:
            print(f"Warning: A_V_sr index {i} out of bounds [1, {self.NSR}]")
            return 0.0

    def B_D_sr(self, i: int) -> float:
        """Returns the output of the B-class Dorsal stretch receptor for unit i."""
        if 1 <= i <= self.NSR:
            return self.B_D_sr[i]
        else:
            print(f"Warning: B_D_sr index {i} out of bounds [1, {self.NSR}]")
            return 0.0

    def B_V_sr(self, i: int) -> float:
        """Returns the output of the B-class Ventral stretch receptor for unit i."""
        if 1 <= i <= self.NSR:
            return self.B_V_sr[i]
        else:
            print(f"Warning: B_V_sr index {i} out of bounds [1, {self.NSR}]")
            return 0.0