from collections import Dict
from math import sqrt
from memory import UnsafePointer
from utils import BlockingSpinLock, BlockingScopedLock

from ..cplx import ComplexScalar, CSRCMatrix, CMatrix
from ..config import DEFAULT_TYPE, DEFAULT_TOL, DEFAULT_ZERO_THRESHOD


struct Statevector[
    type: DType = DEFAULT_TYPE,
    tol: Scalar[type] = DEFAULT_TOL,
](Copyable, Movable, Representable, Stringable, Writable, Sized, Defaultable):
    '''A statevector.

    Parameters:
        type: A type for the circuit data.
        tol: A tolerance for normalization checks.
    '''

    var _size: Int
    '''The size of the statevector.'''
    var _elems: Dict[Int, ComplexScalar[Self.type]]
    '''A dictionary storing the nonzero statevector elements.'''

    # We must store both a lock and a pointer to the lock.
    # If we just do `self._lock = UnsafePointer(to=BlockingSpinLock())` as in the Mojo standard
    # library code for the BlockingSpinLock and BlockingScopedLock, the lock itself gets
    # destroyed after the __init__ and the pointer points to nothing.
    # We cannot just store the lock itself because the BlockingScopedLock constructor requires a
    # mutable reference to the lock, which means methods in which we use
    # `with BlockingScopedLock(self._lock):` must have a mutable reference to self, but this is not
    # always possible, for instance __getitem__ requires self to be borrowed.
    var _lock: BlockingSpinLock
    '''Lock for thread-safe access to the statevector elements.'''
    var _lock_pointer: UnsafePointer[BlockingSpinLock]
    '''Lock for thread-safe access to the statevector elements.'''

    @always_inline
    fn __init__(out self):
        '''Initialize an empty statevector.'''
        self._size = 0
        self._elems = Dict[Int, ComplexScalar[Self.type]]()
        self._lock = BlockingSpinLock()
        self._lock_pointer = UnsafePointer(to=self._lock)

    fn __init__(
        out self, 
        owned statevector: Dict[Int, ComplexScalar[Self.type]],
        n_elements: Int, 
        normalize: Bool = False,
    ) raises:
        '''Initialize a statevector from a dictionary of index element pairs.

        Args:
            statevector: A dictionary containing statevector elements.
            n_elements: The number of elements in the statevector.
            normalize: If True, the statevector will be automatically normalized. If False and the
                statevector is unnormalized, an error will be raised.
        '''
        self._lock = BlockingSpinLock()
        self._lock_pointer = UnsafePointer(to=self._lock)

        # Check that all indices are valid and compute the sum of squares of elements
        var sum_sqr: Scalar[Self.type] = 0
        for idx_elem in statevector.items():
            var idx: Int = idx_elem[].key
            var elem: ComplexScalar[Self.type] = idx_elem[].value
            if idx < 0 or idx >= n_elements:
                raise Error(
                    'Invalid statevector index for statevector with ' + String(n_elements)
                    + ' elements: ' + String(idx) + '.'
                )
            sum_sqr += elem.squared_norm()

        self._size = n_elements

        # Initialize the statevector, normalizing (or raising) if necessary
        self._elems = Dict[Int, ComplexScalar[Self.type]]()
        if abs(sum_sqr - 1) >= Self.tol:
            if normalize:
                var norm_factor: ComplexScalar[Self.type] = sqrt(sum_sqr)
                for idx_elem in statevector.items():
                    var element: ComplexScalar[Self.type] = idx_elem[].value / norm_factor
                    if element != 0:
                        self._elems[idx_elem[].key] = element
            else:
                raise Error('Statevector is not normalized.')
        else:
            for idx_elem in statevector.items():
                var element: ComplexScalar[Self.type] = idx_elem[].value
                if element != 0:
                    self._elems[idx_elem[].key] = element
            self._elems = statevector^
    
    fn __init__[zero_threshold: Scalar[Self.type]](
        out self, 
        owned statevector: CSRCMatrix[Self.type, zero_threshold],
        normalize: Bool = False,
        enforce_n_elements: Int = -1,
    ) raises:
        '''Initialize a statevector from a sparse matrix.

        Args:
            statevector: A sparse matrix containing statevector elements.
            normalize: If True, the statevector will be automatically normalized. If False and the
                statevector is unnormalized, an error will be raised.
            enforce_n_elements: If nonnegative, enforce that the statevector has exactly this many
                elements.
        '''
        self._lock = BlockingSpinLock()
        self._lock_pointer = UnsafePointer(to=self._lock)

        # Check that the statevector is one dimensional
        if statevector.rows != 1 and statevector.cols != 1:
            raise Error(
                'Expected 1D statevector. Received shape (' + String(statevector.rows)
                + ', ' + String(statevector.cols) + ').'
            )

        # Enforce the right number of elements were provided
        if enforce_n_elements >= 0 and statevector.size != enforce_n_elements:
            raise Error(
                'Statevector must have ' + String(enforce_n_elements) + ' elements, but '
                + String(statevector.size) + ' elements were provided.'
            )
        self._size = statevector.size

        # Make it a row vector
        if statevector.rows != 1:
            statevector = statevector.transpose()

        # Check that the statevector is normalized
        var sum_sqr: Scalar[Self.type] = 0
        for elem in statevector.v:
            sum_sqr += elem[].squared_norm()
        if abs(sum_sqr - 1) >= Self.tol:
            if normalize:
                statevector /= sqrt(sum_sqr)
            else:
                raise Error('Statevector is not normalized.')

        # Initialize the statevector
        self._elems = Dict[Int, ComplexScalar[Self.type]]()
        for i in range(statevector.cols):
            self._elems[statevector.col_idx[i]] = statevector.v[i]

    fn __init__(
        out self, 
        owned statevector: CMatrix[Self.type],
        normalize: Bool = False, 
        enforce_n_elements: Int = -1,
    ) raises:
        '''Initialize a statevector from a matrix.

        Args:
            statevector: A matrix containing statevector elements.
            normalize: If True, the statevector will be automatically normalized. If False and the
                statevector is unnormalized, an error will be raised.
            enforce_n_elements: If nonnegative, enforce that the statevector has exactly this many
                elements.
        '''
        self._lock = BlockingSpinLock()
        self._lock_pointer = UnsafePointer(to=self._lock)

        # Check that the statevector is one dimensional
        if statevector.rows != 1 and statevector.cols != 1:
            raise Error(
                'Expected 1D statevector. Received shape (' + String(statevector.rows)
                + ', ' + String(statevector.cols) + ').'
            )

        # Enforce the right number of elements were provided
        if enforce_n_elements >= 0 and statevector.size != enforce_n_elements:
            raise Error(
                'Statevector must have ' + String(enforce_n_elements) + ' elements, but '
                + String(statevector.size) + ' elements were provided.'
            )
        self._size = statevector.size

        # Make it a row vector
        if statevector.rows != 1:
            statevector = statevector.transpose()

        # Check that the statevector is normalized
        var sum_sqr: Scalar[Self.type] = 0
        for idx in range(statevector.cols):
            sum_sqr += statevector.load_idx[1](idx).squared_norm()
        if abs(sum_sqr - 1) >= Self.tol:
            if normalize:
                statevector /= sqrt(sum_sqr)
            else:
                raise Error('Statevector is not normalized.')

        # Initialize the statevector
        self._elems = Dict[Int, ComplexScalar[Self.type]]()
        for idx in range(statevector.cols):
            var element: ComplexScalar[Self.type] = statevector.load_idx[1](idx)
            if element != 0:
                self._elems[idx] = element

    fn __init__(
        out self, 
        owned statevector: List[ComplexScalar[Self.type], True],
        normalize: Bool = False,
        enforce_n_elements: Int = -1,
    ) raises:
        '''Initialize a statevector from a list of elements.

        Args:
            statevector: A list containing statevector elements.
            normalize: If True, the statevector will be automatically normalized. If False and the
                statevector is unnormalized, an error will be raised.
            enforce_n_elements: If nonnegative, enforce that the statevector has exactly this many
                elements.
        '''
        self._lock = BlockingSpinLock()
        self._lock_pointer = UnsafePointer(to=self._lock)

        # Enforce the right number of elements were provided
        if enforce_n_elements >= 0 and len(statevector) != enforce_n_elements:
            raise Error(
                'Statevector must have ' + String(enforce_n_elements) + ' elements, but '
                + String(len(statevector)) + ' elements were provided.'
            )
        self._size = len(statevector)

        # Initialize the statevector, normalizing (or raising) if necessary
        var sum_sqr: Scalar[Self.type] = 0
        for elem in statevector:
            sum_sqr += elem[].squared_norm()
        if abs(sum_sqr - 1) >= Self.tol:
            if normalize:
                self._elems = Dict[Int, ComplexScalar[Self.type]]()
                var norm_factor: ComplexScalar[Self.type] = sqrt(sum_sqr)
                for i in range(len(statevector)):
                    self._elems[i] = statevector[i] / norm_factor
            else:
                raise Error('Statevector is not normalized.')
        else:
            self._elems = Dict[Int, ComplexScalar[Self.type]]()
            for i in range(len(statevector)):
                self._elems[i] = statevector[i]

    fn __init__[tol: Scalar[Self.type]](
        out self, 
        owned statevector: Statevector[Self.type, tol],
        normalize: Bool = False,
        enforce_n_elements: Int = -1,
    ) raises:
        '''Initialize a statevector from another statevector.

        Args:
            statevector: A statevector.
            normalize: If True, the statevector will be automatically normalized. If False and the
                statevector is unnormalized, an error will be raised.
            enforce_n_elements: If nonnegative, enforce that the statevector has exactly this many
                elements.
        '''
        self._lock = BlockingSpinLock()
        self._lock_pointer = UnsafePointer(to=self._lock)

        if enforce_n_elements >= 0 and len(statevector) != enforce_n_elements:
            raise Error(
                'Statevector must have ' + String(enforce_n_elements) + ' elements, but '
                + String(len(statevector)) + ' elements were provided.'
            )

        var sum_sqr: Scalar[Self.type] = 0
        with BlockingScopedLock(statevector._lock_pointer):
            for elem in statevector._elems.values():
                sum_sqr += elem[].squared_norm()
        if abs(sum_sqr - 1) >= Self.tol:
            if normalize:
                statevector.normalize(sum_sqr)
            else:
                raise Error('Statevector is not normalized.')
        
        self._size = statevector._size
        self._elems = Dict[Int, ComplexScalar[Self.type]]()

        with BlockingScopedLock(statevector._lock_pointer):
            # TODO: Find out if we need a lock in this case
            # If we don't need a lock, we can transfer ownership of _elems instead of copying
            self._elems = statevector._elems

    fn __copyinit__(out self, existing: Self):
        '''Initialize a statevector by copying another.

        Args:
            existing: The statevector to copy.
        '''
        self._size= existing._size
        self._elems = existing._elems
        self._lock = BlockingSpinLock()
        self._lock_pointer = UnsafePointer(to=self._lock)

    fn __moveinit__(out self, owned existing: Self):
        '''Initialize a statevector by moving another into self.

        Args:
            existing: The statevector to move into self.
        '''
        self._size = existing._size
        self._elems = existing._elems
        self._lock = BlockingSpinLock()
        self._lock_pointer = UnsafePointer(to=self._lock)

    @staticmethod
    @always_inline
    fn zero(n_qubits: Int) -> Self:
        '''Initialize a statevector to the |0⟩ state.

        Args:
            n_qubits: The number of qubits in the statevector.
        '''
        var elems = Dict[Int, ComplexScalar[Self.type]]()
        elems[0] = 1

        var statevector = Self()
        statevector._size = 2 ** max(n_qubits, 0)
        statevector._elems = elems^
        return statevector^

    @always_inline
    fn __getitem__(self, idx: Int) raises -> ComplexScalar[Self.type]:
        '''Get the statevector element at index idx.

        Args:
            idx: The index of the element to get.

        Returns:
            The element at index idx.
        '''
        if idx < -self._size or idx >= self._size:
            raise Error('Invalid statevector index: ' + String(idx) + '.')

        var true_idx: Int = self._size + idx if idx < 0 else idx

        return self._get(true_idx)
    
    @always_inline
    fn _get(self, idx: Int) -> ComplexScalar[Self.type]:
        '''Get the statevector element at index idx without checking index validity. This does not
        support negative indices.

        Args:
            idx: The index of the element to get.

        Returns:
            The element at index idx.
        '''
        var elem: ComplexScalar[Self.type]
        with BlockingScopedLock(self._lock_pointer):
            return self._elems.get(idx, 0)

    @always_inline
    fn __setitem__(mut self, idx: Int, val: ComplexScalar[Self.type]) raises:
        '''Set the statevector element at index idx to val.

        Args:
            idx: The index of the element to set.
            val: The value to set.
        '''
        if idx < -self._size or idx >= self._size:
            raise Error('Invalid statevector index: ' + String(idx) + '.')

        var true_idx: Int = self._size + idx if idx < 0 else idx

        self._set(true_idx, val)

    @always_inline
    fn _set(mut self, idx: Int, val: ComplexScalar[Self.type]):
        '''Set the statevector element at index idx to val without checking index validity. This
        does not support negative indices.

        Args:
            idx: The index of the element to set.
            val: The value to set.
        '''
        with BlockingScopedLock(self._lock_pointer):
            if val != 0:
                self._elems[idx] = val
            elif idx in self._elems:
                # Default is only present to prevent the compiler from complaining about raising
                _ = self._elems.pop(idx, 0)

    @no_inline
    fn __str__(self) -> String:
        '''Return a description of the statevector as a string.

        Returns:
            A string representation of the statevector.
        '''
        return self._get_str_rep(max_lines=8, max_digits=8)

    @no_inline
    fn _get_str_rep(self, max_lines: Int, max_digits: Int) -> String:
        '''Return a description of the statevector as a string.

        Args:
            max_lines: The maximum number of statevector elements to show.
            max_digits: The maximum number of digits to show in each float.

        Returns:
            A string representation of the statevector.
        '''
        # Create a snapshot of the current state under lock
        var elems_snapshot = Dict[Int, ComplexScalar[Self.type]]()
        var n_elems_snapshot: Int
        with BlockingScopedLock(self._lock_pointer):
            n_elems_snapshot = self._size
            for idx_elem in self._elems.items():
                elems_snapshot[idx_elem[].key] = idx_elem[].value

        var str_rep: String = '['
        if n_elems_snapshot <= max_lines:
            for i in range(n_elems_snapshot):
                var elem: ComplexScalar[Self.type] = elems_snapshot.get(i, 0)
                if i > 0:
                    str_rep += ' '
                str_rep += elem._str1(elem, max_digits=max_digits) + '\n'
            str_rep = str_rep[:-1] + ']'
            return str_rep

        var bottom_lines: Int = max_lines // 2
        var top_lines: Int = max_lines - bottom_lines
        for i in range(top_lines):
            var elem: ComplexScalar[Self.type] = elems_snapshot.get(i, 0)
            if i > 0:
                str_rep += ' '
            str_rep += elem._str1(elem, max_digits=max_digits) + '\n'
        str_rep += ' ...\n'
        for i in range(n_elems_snapshot - 1, n_elems_snapshot - 1 - bottom_lines, -1):
            var elem: ComplexScalar[Self.type] = elems_snapshot.get(i, 0)
            str_rep += ' ' + elem._str1(elem, max_digits=max_digits) + '\n'
        str_rep = str_rep[:-1] + ']'
        return str_rep

    @no_inline
    fn __repr__(self) -> String:
        '''Return a string representation of the statevector.

        Returns:
            A string representation of the statevector.
        '''
        return self._get_str_rep(max_lines=64, max_digits=8)

    @no_inline
    fn write_to[W: Writer](self, mut writer: W):
        '''Write the statevector to a writer.

        Args:
            writer: The writer to write to.
        '''
        writer.write(String(self))

    @always_inline
    fn __len__(self) -> Int:
        '''Get the number of elemets in the statevector.

        Returns:
            The number of statevector elements.
        '''
        return self._size

    fn normalize(mut self, owned sum_sqr: Scalar[Self.type] = -1):
        '''Normalize the statevector in-place.

        Args:
            sum_sqr: The sum of absolute squares of the statevector elements. If negative, it will
                be computed from the statevector elements.
        '''
        with BlockingScopedLock(self._lock_pointer):
            if sum_sqr < 0:
                sum_sqr = 0
                for elem in self._elems.values():
                    sum_sqr += elem[].squared_norm()
            if abs(sum_sqr - 1) >= Self.tol:
                var norm_factor: ComplexScalar[Self.type] = sqrt(sum_sqr)
                for idx_elem in self._elems.items():
                    self._elems[idx_elem[].key] = idx_elem[].value / norm_factor

    fn _clean[zero_threshold: Scalar[Self.type] = DEFAULT_ZERO_THRESHOD](mut self):
        '''Zero elements in the statevector which have norm less than the zero threshold.

        Parameters:
            zero_threshold: The threshold below which elements are set to zero.
        '''
        with BlockingScopedLock(self._lock_pointer):
            var to_remove = List[Int, True]()
            for idx_elem in self._elems.items():
                if idx_elem[].value.norm() < zero_threshold:
                    to_remove.append(idx_elem[].key)
            for idx in to_remove:
                # Default is only present to prevent the compiler from complaining about raising
                _ = self._elems.pop(idx[], 0)
