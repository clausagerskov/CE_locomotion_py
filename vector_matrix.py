import sys
import math
import pickle
from numbers import Number
from typing import TypeVar, Generic, List, Optional, Union, Sequence

# Define a generic type variable
EltType = TypeVar('EltType')

# Equivalent to C++ DEBUG macro (can be set to True for bounds checking)
DEBUG = False

class TVector(Generic[EltType]):
    """
    Python implementation of the C++ TVector class.
    Simulates arbitrary lower/upper bounds and 1-based indexing access
    by default, while using a standard 0-based Python list internally.
    """
    def __init__(self, lower_bound: int = 1, upper_bound: int = 0, data: Optional[Sequence[EltType]] = None):
        self._lb: int = 1
        self._ub: int = 0
        self._data: List[EltType] = []
        if data is not None:
            if upper_bound - lower_bound + 1 != len(data):
                raise ValueError("Data length must match bounds size")
            self._lb = lower_bound
            self._ub = upper_bound
            self._data = list(data)
        else:
            self.SetBounds(lower_bound, upper_bound)

    def Size(self) -> int:
        """Returns the number of elements in the vector."""
        return self._ub - self._lb + 1

    def SetSize(self, NewSize: int):
        """Resizes the vector, keeping the lower bound the same."""
        if NewSize < 0:
             raise ValueError(f"Invalid vector size: {NewSize}")
        self.SetBounds(self._lb, self._lb + NewSize - 1)

    def LowerBound(self) -> int:
        """Returns the lower bound index."""
        return self._lb

    def SetLowerBound(self, NewLB: int):
        """Sets the lower bound index, adjusting upper bound."""
        self.SetBounds(NewLB, self._ub)

    def UpperBound(self) -> int:
        """Returns the upper bound index."""
        return self._ub

    def SetUpperBound(self, NewUB: int):
        """Sets the upper bound index."""
        self.SetBounds(self._lb, NewUB)

    def SetBounds(self, NewLB: int, NewUB: int):
        """
        Sets the bounds of the TVector, reallocating space as necessary and
        preserving as much of the previous contents as possible.
        """
        new_len = NewUB - NewLB + 1
        if new_len < 0:
            new_len = 0 # Allow zero-length vectors

        if self._lb == NewLB and self._ub == NewUB:
            return # No change needed

        old_len = self.Size()
        old_data = self._data
        old_lb = self._lb

        self._data = [None] * new_len # Allocate new storage (Python handles type)

        # Copy overlapping data
        if old_len > 0 and new_len > 0:
            copy_start_old_idx = max(self._lb, NewLB)
            copy_end_old_idx = min(self._ub, NewUB)

            for i in range(copy_start_old_idx, copy_end_old_idx + 1):
                old_internal_idx = i - old_lb
                new_internal_idx = i - NewLB
                if 0 <= old_internal_idx < old_len:
                     self._data[new_internal_idx] = old_data[old_internal_idx]

        self._lb = NewLB
        self._ub = NewUB

    def FillContents(self, value: EltType):
        """Fills the vector with the given value."""
        for i in range(len(self._data)):
            self._data[i] = value

    def PushFront(self, value: EltType):
        """
        Pushes an element to the front.
        The whole list shifts and the last element is dropped.
        """
        if self.Size() > 0:
            self._data.insert(0, value)
            self._data.pop() # Remove the last element

    def InitializeContents(self, *args: EltType):
        """Initializes vector contents from arguments (use with caution)."""
        if len(args) != self.Size():
            raise ValueError("Number of arguments must match vector size for InitializeContents")
        for i in range(self.Size()):
            self._data[i] = args[i]

    def _check_bounds(self, index: int):
        """Internal helper for bounds checking."""
        if not (self._lb <= index <= self._ub):
            raise IndexError(f"Vector index {index} out of bounds ({self._lb} to {self._ub})")

    def __getitem__(self, index: int) -> EltType:
        """Allows access using v[index] with bounds checking."""
        if DEBUG:
            self._check_bounds(index)
        # Adjust index for 0-based internal list
        internal_index = index - self._lb
        if not (0 <= internal_index < len(self._data)):
             raise IndexError(f"Internal calculation error or Vector index {index} out of bounds ({self._lb} to {self._ub})")
        return self._data[internal_index]

    def __setitem__(self, index: int, value: EltType):
        """Allows assignment using v[index] = value with bounds checking."""
        if DEBUG:
            self._check_bounds(index)
        # Adjust index for 0-based internal list
        internal_index = index - self._lb
        if not (0 <= internal_index < len(self._data)):
             raise IndexError(f"Internal calculation error or Vector index {index} out of bounds ({self._lb} to {self._ub})")
        self._data[internal_index] = value

    def __call__(self, index: int) -> EltType:
         """Allows safe access using v(index) like C++."""
         self._check_bounds(index)
         internal_index = index - self._lb
         return self._data[internal_index]

    def __len__(self) -> int:
        """Returns the size of the vector."""
        return self.Size()

    def __eq__(self, other) -> bool:
        """Compares two TVectors for equality."""
        if not isinstance(other, TVector):
            return NotImplemented
        if self.Size() != other.Size() or self._lb != other._lb:
             # Note: C++ also checked UpperBound, but Size + LowerBound check is sufficient
            return False
        return self._data == other._data # Compare internal lists directly

    def __str__(self) -> str:
        """String representation for printing."""
        return " ".join(map(str, self._data))

    def __repr__(self) -> str:
        """Detailed representation."""
        return f"TVector(lb={self._lb}, ub={self._ub}, data={self._data})"

    def to_list(self) -> List[EltType]:
        """Returns the internal data as a standard Python list."""
        return list(self._data) # Return a copy

    def BinaryWriteVector(self, bofs): # bofs is a binary file stream opened for writing
        """Writes the vector state to a binary file stream using pickle."""
        pickle.dump(self._lb, bofs)
        pickle.dump(self._ub, bofs)
        pickle.dump(self._data, bofs) # Pickle the internal list

    def BinaryReadVector(self, bifs): # bifs is a binary file stream opened for reading
        """Reads the vector state from a binary file stream using pickle."""
        self._lb = pickle.load(bifs)
        self._ub = pickle.load(bifs)
        self._data = pickle.load(bifs) # Load the internal list
        # Basic consistency check
        if len(self._data) != max(0, self._ub - self._lb + 1):
             raise ValueError("Inconsistent data size read from binary file")

import copy # Needed for deep copy in copy constructor

# *******
# TMatrix
# *******

class TMatrix(Generic[EltType]):
    """
    Python implementation of the C++ TMatrix class.
    Uses a list of TVector objects to represent rows, allowing arbitrary
    row and column bounds similar to the C++ version.
    """
    def __init__(self,
                 row_lb: Optional[int] = None, row_ub: Optional[int] = None,
                 col_lb: Optional[int] = None, col_ub: Optional[int] = None,
                 source_matrix: Optional["TMatrix[EltType]"] = None):
        """
        Constructs a TMatrix.
        - TMatrix(): Creates an empty matrix with default bounds (1, 0, 1, 0).
        - TMatrix(rlb, rub, clb, cub): Creates a matrix with specified bounds.
        - TMatrix(source_matrix=other_matrix): Creates a deep copy of another TMatrix.
        """
        self._lb1: int = 1  # Row lower bound
        self._ub1: int = 0  # Row upper bound
        self._lb2: int = 1  # Col lower bound
        self._ub2: int = 0  # Col upper bound
        self._rows: List[TVector[EltType]] = [] # Internal storage: list of row TVectors

        if source_matrix is not None:
            # --- Copy constructor logic ---
            if not isinstance(source_matrix, TMatrix):
                raise TypeError("source_matrix must be an instance of TMatrix for copy constructor")
            # Set bounds first, which allocates rows
            self.SetBounds(source_matrix.RowLowerBound(), source_matrix.RowUpperBound(),
                           source_matrix.ColumnLowerBound(), source_matrix.ColumnUpperBound())
            # Perform deep copy of elements
            for r in range(self._lb1, self._ub1 + 1):
                for c in range(self._lb2, self._ub2 + 1):
                    # Use __getitem__ for row access and then TVector's __setitem__
                    self[r][c] = source_matrix(r, c) # Use M(r,c) for safe access

        elif row_lb is not None and row_ub is not None and \
             col_lb is not None and col_ub is not None:
            # --- Bounds constructor logic ---
            self.SetBounds(row_lb, row_ub, col_lb, col_ub)
        else:
            # --- Default constructor logic ---
            # Keep default bounds (1, 0, 1, 0), no rows allocated
            pass

    # --- Accessors ---
    def RowSize(self) -> int:
        """Returns the number of rows."""
        return max(0, self._ub1 - self._lb1 + 1)

    def SetRowSize(self, NewSize: int):
        """Resizes the number of rows, keeping row lower bound."""
        self.SetBounds(self._lb1, self._lb1 + NewSize - 1, self._lb2, self._ub2)

    def ColumnSize(self) -> int:
        """Returns the number of columns."""
        return max(0, self._ub2 - self._lb2 + 1)

    def SetColumnSize(self, NewSize: int):
        """Resizes the number of columns, keeping column lower bound."""
        self.SetBounds(self._lb1, self._ub1, self._lb2, self._lb2 + NewSize - 1)

    def SetSize(self, NewRowSize: int, NewColSize: int):
        """Resizes both rows and columns, keeping lower bounds."""
        self.SetBounds(self._lb1, self._lb1 + NewRowSize - 1,
                       self._lb2, self._lb2 + NewColSize - 1)

    def RowLowerBound(self) -> int: return self._lb1
    def SetRowLowerBound(self, newlb1: int): self.SetBounds(newlb1, self._ub1, self._lb2, self._ub2)
    def RowUpperBound(self) -> int: return self._ub1
    def SetRowUpperBound(self, newub1: int): self.SetBounds(self._lb1, newub1, self._lb2, self._ub2)
    def ColumnLowerBound(self) -> int: return self._lb2
    def SetColumnLowerBound(self, newlb2: int): self.SetBounds(self._lb1, self._ub1, newlb2, self._ub2)
    def ColumnUpperBound(self) -> int: return self._ub2
    def SetColumnUpperBound(self, newub2: int): self.SetBounds(self._lb1, self._ub1, self._lb2, newub2)

    def SetBounds(self, newlb1: int, newub1: int, newlb2: int, newub2: int):
        """
        Sets the bounds of the TMatrix. Reallocates space and does *not*
        preserve previous contents, matching C++ behavior.
        """
        # Only do it if bounds actually change
        if (newlb1 == self._lb1 and newub1 == self._ub1 and
            newlb2 == self._lb2 and newub2 == self._ub2):
            return

        # Store new bounds info
        self._lb1 = newlb1
        self._ub1 = newub1
        self._lb2 = newlb2
        self._ub2 = newub2

        num_rows = max(0, self._ub1 - self._lb1 + 1)
        num_cols = max(0, self._ub2 - self._lb2 + 1) # Size for inner TVectors

        # No negative sizes allowed (handled by max(0, ...))
        # C++ version error check:
        # if num_rows < 0 or num_cols < 0: # Should not happen with max(0,...)
        #    raise ValueError("Attempt to allocate a negative sized TMatrix")

        # Allocate new storage (list of TVectors)
        # Clear the old list first (Python garbage collection handles old TVectors)
        self._rows = []
        if num_rows > 0:
            self._rows = [TVector[EltType](self._lb2, self._ub2) for _ in range(num_rows)]


    # --- Overloaded Operators ---
    def __getitem__(self, row_index: int) -> TVector[EltType]:
        """
        Provides access to a row using M[row_index].
        Returns the TVector representing that row.
        Performs row bounds checking.
        """
        if DEBUG:
             if not (self._lb1 <= row_index <= self._ub1):
                  raise IndexError(f"Matrix row index {row_index} out of bounds ({self._lb1} to {self._ub1})")
        # Convert external row_index to internal 0-based index for the list
        internal_row_index = row_index - self._lb1
        if not (0 <= internal_row_index < len(self._rows)):
             # This might happen if bounds are set but allocation failed, or index is wrong
             raise IndexError(f"Matrix row index {row_index} maps to internal index {internal_row_index}, which is out of bounds for internal list (size {len(self._rows)}). Bounds are ({self._lb1}..{self._ub1}, {self._lb2}..{self._ub2})")
        return self._rows[internal_row_index]

    def __call__(self, row_index: int, col_index: int) -> EltType:
        """
        Provides safe element access using M(row_index, col_index).
        Performs bounds checking for both row and column.
        """
        # Get the row TVector using __getitem__, which checks row bounds
        row_vector = self[row_index]
        # Access the element in the row TVector using its () operator for column check
        return row_vector(col_index)

    # Assignment operator (Python handles by reference, use copy constructor for value copy)
    # def __assign__(self, other): # Not standard Python, use copy ctor
    #     if not isinstance(other, TMatrix):
    #         raise TypeError("Can only assign from another TMatrix")
    #     self.SetBounds(other.RowLowerBound(), other.RowUpperBound(),
    #                    other.ColumnLowerBound(), other.ColumnUpperBound())
    #     for r in range(self._lb1, self._ub1 + 1):
    #         for c in range(self._lb2, self._ub2 + 1):
    #              self[r][c] = other(r, c)
    #     return self

    def __eq__(self, other) -> bool:
        """Checks for equality with another TMatrix."""
        if not isinstance(other, TMatrix):
            return NotImplemented
        if (self.RowSize() != other.RowSize() or
            self.ColumnSize() != other.ColumnSize() or
            self.RowLowerBound() != other.RowLowerBound() or
            self.ColumnLowerBound() != other.ColumnLowerBound()):
             return False
        # Check element-wise equality
        for r in range(self._lb1, self._ub1 + 1):
             # Compare internal row TVectors directly
             if self[r] != other[r]: # Uses TVector.__eq__
                  return False
        return True

    # --- Other Stuff ---
    def FillContents(self, value: EltType):
        """Fills the entire matrix with the given value."""
        for r in range(self.RowSize()): # Iterate using 0-based internal index
            if r < len(self._rows):
                 self._rows[r].FillContents(value) # Use TVector's FillContents

    # InitializeContents not implemented due to complexity with *args in Python
    # def InitializeContents(self, v1: EltType, *args):
    #     pass

    def __str__(self) -> str:
        """String representation for printing."""
        return "\n".join(str(row) for row in self._rows)

    def __repr__(self) -> str:
        """Detailed representation."""
        return f"TMatrix(rlb={self._lb1}, rub={self._ub1}, clb={self._lb2}, cub={self._ub2}, rows={len(self._rows)})"