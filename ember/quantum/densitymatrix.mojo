##### UNTESTED ##### UNTESTED ##### UNTESTED ##### UNTESTED ##### UNTESTED ##### UNTESTED ##### 

from collections import Dict
from bit import log2_floor
from math import sqrt
from memory import UnsafePointer
from utils import BlockingSpinLock, BlockingScopedLock

from ..cplx import CMatrix, is_positive_semidefinite
from ..config import DEFAULT_TYPE, DEFAULT_TOL


struct DensityMatrix[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    Copyable, Movable, Defaultable
):
    '''A density matrix representing the state of multiple qubits.
    
    Parameters:
        type: A type for the data.
        tol: A tolerance for zero checks.
    '''

    var n_qubits: Int
    '''The number of qubits in the state.'''
    var matrix: CMatrix[Self.type]
    '''The density matrix.'''
    # See Statevector for the reasoning behind storing the lock in this manner
    var _lock: BlockingSpinLock
    '''Lock for thread-safe access to the density matrix elements.'''
    var _lock_pointer: UnsafePointer[BlockingSpinLock]
    '''Lock for thread-safe access to the density matrix elements.'''

    @always_inline
    fn __init__(out self):
        '''Initialize an empty density matrix.'''
        self.n_qubits = 0
        self.matrix = CMatrix[Self.type](rows=0, cols=0, fill_zeros=False)
        self._lock = BlockingSpinLock()
        self._lock_pointer = UnsafePointer(to=self._lock)
    
    fn __init__(out self, owned matrix: CMatrix[Self.type]) raises:
        '''Initialize a density matrix from a CMatrix.

        Args:
            matrix: The CMatrix defining the density matrix.
        '''
        if not (matrix.rows & (matrix.rows - 1) == 0):
            raise Error(
                'Matrix does not represent a multi-qubit state. '
                'Matrix dimension must be a power of 2.'
            )
        if not matrix.trace().is_close[Self.tol](1):
            raise Error('Matrix does not have trace 1.')
        if not is_positive_semidefinite[Self.type, Self.tol](matrix):
            # This also checks if the matrix is square and Hermitian
            raise Error('Matrix is not Hermitian or not positive semi-definite.')
        self.n_qubits = log2_floor(matrix.rows)
        self.matrix = matrix
        self._lock = BlockingSpinLock()
        self._lock_pointer = UnsafePointer(to=self._lock)
    
    fn __init__(out self, elem_dict: Dict[Int, ComplexScalar[Self.type]], size: Int) raises:
        '''Initialize a density matrix from a dictionary of statavector elements and a size.

        Args:
            elem_dict: A dictionary mapping statevector indices to non-zero values.
            size: The total length of the statevector.
        '''
        if not (size & (size - 1) == 0):
            raise Error(
                'Matrix does not represent a multi-qubit state. '
                'Size must be a power of 2.'
            )
        var sum_sqr: Scalar[type] = 0
        for elem in elem_dict.values():
            sum_sqr += elem[].squared_norm()
        if abs(sum_sqr - 1) >= Self.tol:
            raise Error(
                'Cannot convert statevector to density matrix. '
                'Statevector is not normalized.'
            )

        self.n_qubits = log2_floor(size)
        self._lock = BlockingSpinLock()
        self._lock_pointer = UnsafePointer(to=self._lock)
        self.matrix = CMatrix[Self.type](size, size, fill_zeros=False)
        # TODO: Check that the statevector is normalized
        with BlockingScopedLock(self._lock_pointer):
            for i in range(size):
                self.matrix.store_crd[1](i, i, elem_dict.get(i, 0).squared_norm())
                for j in range(i + 1, size):
                    var elem: ComplexScalar[Self.type] = (
                        elem_dict.get(i, 0) * elem_dict.get(j, 0).conj()
                    )
                    self.matrix.store_crd[1](i, j, elem)
                    self.matrix.store_crd[1](j, i, elem.conj())
    
    @always_inline
    fn __copyinit__(out self, existing: Self):
        '''Initialize a density matrix by copying another.

        Args:
            existing: The density matrix to copy.
        '''
        self.n_qubits = existing.n_qubits
        self.matrix = existing.matrix
        self._lock = BlockingSpinLock()
        self._lock_pointer = UnsafePointer(to=self._lock)
    
    @always_inline
    fn __moveinit__(out self, owned existing: Self):
        '''Initialize a density matrix by moving another into self.

        Args:
            existing: The density matrix to move into self.
        '''
        self.n_qubits = existing.n_qubits
        self.matrix = existing.matrix^
        self._lock = BlockingSpinLock()
        self._lock_pointer = UnsafePointer(to=self._lock)
    
    @always_inline
    fn purity(self) raises -> Scalar[Self.type]:
        '''Compute the purity of the density matrix. Pure states have purity 1, whereas maximally
        mixed states in N qubits have purity 1/2^N.

        Returns:
            The purity of self.
        '''
        with BlockingScopedLock(self._lock_pointer):
            return (self.matrix @ self.matrix).trace().re
    
    @always_inline
    fn is_pure[tol: Scalar[type] = Self.tol](self) raises -> Bool:
        '''Check if the density matrix is pure.

        Returns:
            True if self is pure, False otherwise.
        '''
        return abs(self.purity() - 1) < tol
        
    fn partial_trace(self, traced_qubits: List[Int, True]) raises -> Self:
        '''Compute the partial trace of the density matrix over specified qubits.

        Args:
            traced_qubits: List of qubit indices to trace out.

        Returns:
            New density matrix with specified qubits traced out.
        '''
        # Validate input qubits
        for idx in traced_qubits:
            if idx[] < 0 or idx[] >= self.n_qubits:
                raise Error('Invalid qubit specifier: ' + String(idx[]) + '.')

        # Calculate remaining qubits
        var remaining_qubits = List[Int, True]()
        for i in range(self.n_qubits):
            if i not in traced_qubits:
                remaining_qubits.append(i)

        var reduced_size: Int = 2 ** len(remaining_qubits)

        # Create result matrix
        var result_matrix = CMatrix[Self.type](reduced_size, reduced_size, fill_zeros=True)

        with BlockingScopedLock(self._lock_pointer):
            # Iterate over all basis states of the full system
            for row_idx in range(self.matrix.rows):
                for col_idx in range(self.matrix.cols):
                    # Convert linear indices to multi-dimensional indices
                    var row_indices: List[Int, True] = self._linear_to_multi_index(
                        row_idx, self.n_qubits
                    )
                    var col_indices: List[Int, True] = self._linear_to_multi_index(
                        col_idx, self.n_qubits
                    )

                    # Check if traced subsystems have matching indices
                    for traced_q in traced_qubits:
                        if row_indices[traced_q[]] != col_indices[traced_q[]]:
                            break
                    else:
                        # If traced subsystems match, contribute to partial trace
                        # Extract indices for remaining subsystems
                        var remaining_row_indices = List[Int, True]()
                        var remaining_col_indices = List[Int, True]()
                        for q in remaining_qubits:
                            remaining_row_indices.append(row_indices[q[]])
                            remaining_col_indices.append(col_indices[q[]])

                        # Convert back to linear indices for reduced system
                        var reduced_row: Int = self._multi_to_linear_index(
                            remaining_row_indices, len(remaining_qubits)
                        )
                        var reduced_col: Int = self._multi_to_linear_index(
                            remaining_col_indices, len(remaining_qubits)
                        )

                        # Add contribution to result matrix
                        result_matrix.store_crd[1](
                            reduced_row,
                            reduced_col,
                            (
                                result_matrix.load_crd[1](reduced_row, reduced_col) 
                                + self.matrix.load_crd[1](row_idx, col_idx)
                            ),
                        )

        return Self(result_matrix^)

    # Helper function: Convert linear index to multi-dimensional indices
    @staticmethod
    fn _linear_to_multi_index(linear_idx: Int, n_qubits: Int) -> List[Int, True]:
        '''Convert a linear index to multi-dimensional indices.

        Args:
            linear_idx: The linear index.
            n_qubits: The number of qubits in the system.

        Returns:
            The multi-index corresponding to the linear index.
        '''
        var indices = List[Int, True]()
        var remaining: Int = linear_idx

        for _ in range(n_qubits):
            indices.append(remaining & 1)  # Extract least significant bit
            remaining >>= 1  # Shift right by 1 bit

        return indices

    # Helper function: Convert multi-dimensional indices to linear index
    @staticmethod
    fn _multi_to_linear_index(indices: List[Int, True], n_qubits: Int) -> Int:
        '''Convert multi-dimensional indices to a linear index.

        Args:
            indices: The multi-index.
            n_qubits: The nujmber of qubits in the system.

        Returns:
            The linear index corresponding to the multi-index.
        '''
        var linear_idx: Int = 0

        for i in range(n_qubits):
            linear_idx += indices[i] * (2 ** i)

        return linear_idx

    # TODO: Test these two methods and benchmark to check which to keep
    fn _to_sv_data[tol: Scalar[type] = Self.tol](
        self
    ) raises -> Tuple[Dict[Int, ComplexScalar[Self.type]], Int]:
        '''Convert a pure density matrix to a statevector using its eigendecomposition.

        Parameters:
            tol: A tolerance for purity checks.

        Returns:
            The data necessary to instantiate a Statevector.
        '''
        if not self.is_pure[tol]():
            raise Error('Cannot convert to statevector. Density matrix is not pure.')
        
        var eigvals: List[ComplexScalar[Self.type], True]
        var eigvecs: CMatrix[Self.type]
        with BlockingScopedLock(self._lock_pointer):
            eigvals, eigvecs = eigs[Self.type, tol](self.matrix)

        # Find the eigenvalue closest to 1 (should be exactly 1 for pure states)
        var max_eigval: Scalar[Self.type] = 0
        var max_eigval_idx: Int = -1
        for i in range(len(eigvals)):
            if abs(eigvals[i].re) > max_eigval:
                max_eigval = abs(eigvals[i].re)
                max_eigval_idx = i
        
        # Verify the eigenvalue is approximately 1
        if abs(max_eigval - 1) >= tol:
            raise Error(
                'Max eigenvalue exceeds 1. '
                'This indicates an invalid density matrix (or a bug in the code!).'
            )
        
        var sv: CMatrix[Self.type] = eigvecs.get_column(col=max_eigval_idx)

        # Normalize the statevector (should already be normalized, but ensure it)
        var norm: Scalar[Self.type] = 0
        for i in range(sv.size):
            norm += sv.load_idx[1](i).squared_norm()
        norm = sqrt(norm)

        # Generate the dictionary of statevector elements
        var sv_elems = Dict[Int, ComplexScalar[Self.type]]()
        for i in range(sv.size):
            var elem: ComplexScalar[Self.type] = sv.load_idx[1](i) / norm
            if elem != 0:
                sv_elems[i] = elem
        
        return sv_elems^, len(sv)

    fn _to_sv_data_direct[tol: Scalar[type] = Self.tol](
        self
    ) raises -> Tuple[Dict[Int, ComplexScalar[Self.type]], Int]:
        '''Convert pure density matrix to statevector using direct extraction.
        This method is more efficient but requires the density matrix to be exactly pure.

        Parameters:
            tol: A tolerance for purity checks.

        Returns:
            The data necessary to instantiate a Statevector.
        '''
        
        if not self.is_pure[tol]():
            raise Error('Cannot convert to statevector. Density matrix is not pure.')
        
        var sv_elems = Dict[Int, ComplexScalar[Self.type]]()
        with BlockingScopedLock(self._lock_pointer):
            # Find the first non-zero diagonal element to determine phase reference
            var ref_idx: Int
            for i in range(self.matrix.rows):
                var diagonal_elem: ComplexScalar[Self.type] = self.matrix.load_crd[1](i, i)
                if diagonal_elem.squared_norm() > Self.tol:
                    ref_idx = i
                    break
            else:
                raise Error('No non-zero diagonal elements found.')
            
            # Extract statevector amplitudes
            var ref_amp: Scalar[Self.type] = sqrt(self.matrix.load_crd[1](ref_idx, ref_idx).re)
            for i in range(self.matrix.rows):
                if i == ref_idx:
                    # Use positive real amplitude for reference
                    sv_elems[i] = ComplexScalar[Self.type](ref_amp)
                    continue
                # Calculate amplitude from off-diagonal element
                if ref_amp >= Self.tol:
                    sv_elems[i] = self.matrix.load_crd[1](ref_idx, i) / ref_amp
                    continue
                # Check if this diagonal element is non-zero
                var diag_elem: ComplexScalar[Self.type] = self.matrix.load_crd[1](i, i)
                sv_elems[i] = sqrt(diag_elem.re) if abs(diag_elem.re) > Self.tol else 0
            
        # Verify the extracted statevector is normalized
        var norm_squared: Scalar[Self.type] = 0
        for amp in sv_elems.values():
            norm_squared += amp[].squared_norm()
        
        if abs(norm_squared - 1) >= Self.tol:
            raise Error('Extracted statevector is not normalized.')
        
        return sv_elems^, 2 ** self.n_qubits
