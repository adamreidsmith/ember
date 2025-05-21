from collections import Dict
from math import sqrt

from ..cplx import ComplexScalar, CSRCMatrix, CMatrix
from ..config import DEFAULT_TYPE, DEFAULT_TOL

@value
struct Statevector[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    Copyable, Movable, Representable, Stringable, Writable, Sized, Defaultable
):
    '''A statevector.

    Parameters:
        type: A type for the circuit data.
        tol: A tolerance for normalization checks.
    '''

    var _n_elems: Int
    '''The number of elments in the statevector.'''
    var _elems: Dict[Int, ComplexScalar[Self.type]]
    '''A dictionary storing the statevector elements.'''

    @always_inline
    fn __init__(out self):
        '''Initialize an empty statevector.'''
        self._n_elems = 0
        self._elems = Dict[Int, ComplexScalar[Self.type]]()

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
        
        # Initialize the statevector, normalizing (or raising) if necessary
        if abs(sum_sqr - 1) >= Self.tol:
            if normalize:
                self._elems = Dict[Int, ComplexScalar[Self.type]]()
                var norm_factor: ComplexScalar[Self.type] = sqrt(sum_sqr)
                for idx_elem in statevector.items():
                    self._elems[idx_elem[].key] = idx_elem[].value / norm_factor
            else:
                raise Error('Statevector is not normalized.')
        else:
            self._elems = statevector^
        
        self._n_elems = n_elements
    
    fn __init__(
        out self, 
        owned statevector: CSRCMatrix[Self.type],  
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
        self._n_elems = statevector.size

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
        self._n_elems = statevector.size

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
            self._elems[idx] = statevector.load_idx[1](idx)
            
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
        # Enforce the right number of elements were provided
        if enforce_n_elements >= 0 and len(statevector) != enforce_n_elements:
            raise Error(
                'Statevector must have ' + String(enforce_n_elements) + ' elements, but ' 
                + String(len(statevector)) + ' elements were provided.'
            )
        self._n_elems = len(statevector)

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
    
    fn __init__(
        out self, 
        owned statevector: Self, 
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
        if enforce_n_elements >= 0 and len(statevector) != enforce_n_elements:
            raise Error(
                'Statevector must have ' + String(enforce_n_elements) + ' elements, but ' 
                + String(len(statevector)) + ' elements were provided.'
            )

        var sum_sqr: Scalar[Self.type] = 0
        for elem in statevector._elems.values():
            sum_sqr += elem[].squared_norm()
        if abs(sum_sqr - 1) >= Self.tol:
            if normalize:
                statevector.normalize(sum_sqr)
            else:
                raise Error('Statevector is not normalized.')
        
        self = statevector^
        
    
    @staticmethod
    @always_inline
    fn zero(n_qubits: Int) -> Self:
        '''Initialize a statevector to the |0⟩ state.

        Args:
            n_qubits: The number of qubits in the statevector.
        '''
        var elems = Dict[Int, ComplexScalar[Self.type]]()
        elems[0] = 1
        return Self(_n_elems = 2 ** max(n_qubits, 0), _elems = elems)
        
    @always_inline
    fn __getitem__(self, idx: Int) raises -> ComplexScalar[Self.type]:
        '''Get the statevector element at index idx.
        
        Args:
            idx: The index of the element to get.
        
        Returns:
            The element at index idx.
        '''
        if idx < -self._n_elems or idx >= self._n_elems:
            raise Error('Invalid statevector index: ' + String(idx) + '.')
        if idx < 0:
            return self._elems.get(self._n_elems + idx, 0)
        return self._elems.get(idx, 0)
    
    @always_inline
    fn _get(self, idx: Int) -> ComplexScalar[Self.type]:
        '''Get the statevector element at index idx without checking index validity.
        
        Args:
            idx: The index of the element to get.
        
        Returns:
            The element at index idx.
        '''
        if idx < 0:
            return self._elems.get(self._n_elems + idx, 0)
        return self._elems.get(idx, 0)
    
    @always_inline
    fn __setitem__(mut self, idx: Int, val: ComplexScalar[Self.type]) raises:
        '''Set the statevector element at index idx to val.
        
        Args:
            idx: The index of the element to set.
            val: The value to set.
        '''
        if idx < -self._n_elems or idx >= self._n_elems:
            raise Error('Invalid statevector index: ' + String(idx) + '.')
        if idx < 0:
            self._elems[self._n_elems + idx] = val
        else:
            self._elems[idx] = val

    @always_inline
    fn _set(mut self, idx: Int, val: ComplexScalar[Self.type]):
        '''Set the statevector element at index idx to val without checking index validity.
        
        Args:
            idx: The index of the element to set.
            val: The value to set.
        '''
        if idx < 0:
            self._elems[self._n_elems + idx] = val
        else:
            self._elems[idx] = val
    
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
        var str_rep: String = '['
        if self._n_elems <= max_lines:
            for i in range(self._n_elems):
                var elem: ComplexScalar[Self.type] = self._get(i)
                if i > 0:
                    str_rep += ' '
                str_rep += elem._str1(elem, max_digits=max_digits) + '\n'
            str_rep = str_rep[:-1] + ']'
            return str_rep

        var bottom_lines: Int = max_lines // 2
        var top_lines: Int = max_lines - bottom_lines
        for i in range(top_lines):
            var elem: ComplexScalar[Self.type] = self._get(i)
            if i > 0:
                str_rep += ' '
            str_rep += elem._str1(elem, max_digits=max_digits) + '\n'
        str_rep += ' ...\n'
        for i in range(self._n_elems - 1, self._n_elems - 1 - bottom_lines, -1):
            var elem: ComplexScalar[Self.type] = self._get(i)
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
        return self._n_elems
    
    fn normalize(mut self, owned sum_sqr: Scalar[Self.type] = -1):
        '''Normalize the statevector in-place.
        
        Args:
            sum_sqr: The sum of absolute squares of the statevector elements. If negative, it will
                be computed from the statevector elements.
        '''
        if sum_sqr < 0:
            sum_sqr = 0
            for elem in self._elems.values():
                sum_sqr += elem[].squared_norm()
        if abs(sum_sqr - 1) >= Self.tol:
            var norm_factor: ComplexScalar[Self.type] = sqrt(sum_sqr)
            for idx_elem in self._elems.items():
                self._elems[idx_elem[].key] = idx_elem[].value / norm_factor