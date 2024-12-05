from algorithm import parallelize, map
from math import nan, sqrt
from collections import Dict

from .complexsimd import ComplexScalar
from .cmatrix import CMatrix
from .config import DEFAULT_TOL, DEFAULT_ZERO_THRESHOD


@value
struct CSRCMatrix[type: DType, zero_threshold: Scalar[type] = DEFAULT_ZERO_THRESHOD](
    Absable,
    Formattable, 
    Sized, 
    Representable,
    StringableCollectionElement,
):
    var rows: Int
    var cols: Int
    var size: Int
    var v: List[ComplexScalar[Self.type], True]
    var col_idx: List[Int, True]
    var row_idx: List[Int, True]

    # Initialization ##################

    @always_inline
    fn __init__(inout self, rows: Int, cols: Int):
        '''Initialize a CSRCMatrix of zeros.'''
        self.rows = rows
        self.cols = cols
        self.size = rows * cols
        self.v = List[ComplexScalar[Self.type], True]()
        self.col_idx = List[Int, True]()
        self.row_idx = List[Int, True](0) * (rows + 1)

    fn __init__(
        inout self, rows: Int, cols: Int, *data: Tuple[Int, Int, ComplexScalar[Self.type]]
    ) raises:
        '''Initialize a CSRCMatrix from tuples (row index, column index, value).'''
        for t in data:
            if t[][1] >= cols or t[][1] < 0:
                raise Error(
                    'Column index ' + str(t[][1]) + ' is out of range for a matrix with ' 
                    + str(cols) + ' columns'
                )
            if t[][0] >= rows or t[][0] < 0:
                raise Error(
                    'Row index ' + str(t[][0]) + ' is out of range for a matrix with ' 
                    + str(rows) + ' rows'
                )

        self.rows = rows
        self.cols = cols
        self.size = self.rows * self.cols
        self.v = List[ComplexScalar[Self.type], True]()
        self.col_idx = List[Int, True]()
        self.row_idx = List[Int, True](0) * (rows + 1)
        
        @parameter
        fn cmp_fn(i: Int, j: Int) -> Bool:
            # Sort first by row, then by column
            var row_i: Int = data[i][0]
            var row_j: Int = data[j][0]
            if row_i == row_j:
                return data[i][1] < data[j][1]
            return row_i < row_j

        # TODO: Gotta be a better way to do this
        indices = List[Int, True]()
        for i in range(len(data)):
            if data[i][2] > Self.zero_threshold:
                indices.append(i)
        sort[cmp_fn=cmp_fn](indices)

        # Sort v and col_idx using indices
        var row_indices = List[Int, True]()
        for i in indices:
            self.v.append(data[i[]][2])
            self.col_idx.append(data[i[]][1])
            row_indices.append(data[i[]][0])

        # Remove any duplicates
        var i: Int = 1
        while i < self.v.size:
            if self.col_idx[i] == self.col_idx[i - 1] and row_indices[i] == row_indices[i - 1]:
                self.v[i - 1] += self.v.pop(i)
                _ = self.col_idx.pop(i)
                _ = row_indices.pop(i)
            else:
                i += 1

        self.row_idx = self._sorted_row_indices_to_idx(rows, row_indices)

    fn __init__(
        inout self, 
        rows: Int, 
        cols: Int, 
        row_indices: List[Int, True], 
        col_indices: List[Int, True], 
        vals: List[ComplexScalar[Self.type], True]
    ) raises:
        '''Initialize a CSRCMatrix from lists of row indices, column indices, and values.'''
        if vals.size != col_indices.size or vals.size != row_indices.size:
            raise Error(
                'Row indices, column indices, and values must be the same length; found lengths '
                + str(row_indices.size) + ', '
                + str(col_indices.size) + ', '
                + str(vals.size)
            )

        for i in range(col_indices.size):
            if col_indices[i] >= cols or col_indices[i] < 0:
                raise Error(
                    'Column index ' + str(col_indices[i]) 
                    + ' is out of range for a matrix with ' + str(cols) + ' columns'
                )
            if row_indices[i] >= rows or row_indices[i] < 0:
                raise Error(
                    'Row index ' + str(row_indices[i]) 
                    + ' is out of range for a matrix with ' + str(rows) + ' rows'
                )

        self.rows = rows
        self.cols = cols
        self.size = rows * cols
        self.v = List[ComplexScalar[Self.type], True]()
        self.col_idx = List[Int, True]()
        self.row_idx = List[Int, True](capacity=rows + 1)

        @parameter
        fn cmp_fn(i: Int, j: Int) -> Bool:
            # Sort first by row, then by column
            if row_indices[i] == row_indices[j]:
                return col_indices[i] < col_indices[j]
            return row_indices[i] < row_indices[j]

        # TODO: Gotta be a better way to do this
        indices = List[Int, True]()
        for i in range(vals.size):
            if vals[i] > Self.zero_threshold:
                indices.append(i)
        sort[cmp_fn=cmp_fn](indices)

        # Sort v and col_indices using indices
        var row_idx = List[Int, True]()
        for i in indices:
            self.v.append(vals[i[]])
            self.col_idx.append(col_indices[i[]])
            row_idx.append(row_indices[i[]])

        # Remove any duplicates
        var i: Int = 1
        while i < self.v.size:
            if self.col_idx[i] == self.col_idx[i - 1] and row_idx[i] == row_idx[i - 1]:
                self.v[i - 1] += self.v.pop(i)
                _ = self.col_idx.pop(i)
                _ = row_idx.pop(i)
            else:
                i += 1
        self.row_idx = self._sorted_row_indices_to_idx(rows, row_idx^)
    
    fn __init__(
        inout self, 
        rows: Int, 
        cols: Int, 
        sparse_dict: Dict[Int, Dict[Int, ComplexScalar[Self.type]]],
    ) raises:
        '''Initialize a CSRCMatrix from a nested dictionary representation. 
        sparse_dict is a nested dictionary of non-zero entries: {row_index: {col_index: value}}.
        '''
        self.rows = rows
        self.cols = cols
        self.size = rows * cols
        self.v = List[ComplexScalar[Self.type], True]()
        self.col_idx = List[Int, True]()
        self.row_idx = List[Int, True](0)

        # Extract and sort row indices
        var row_keys = List[Int, True]()
        for key in sparse_dict.keys():
            row_keys.append(key[])
        sort(row_keys)

        # Process each row based on sorted row keys
        var row_ptr: Int = 0
        for row in range(self.rows):
            if row_ptr < row_keys.size and row == row_keys[row_ptr]:
                row_ptr += 1

                # Extract and sort column indices
                col_keys = List[Int, True]()
                for col in sparse_dict[row].keys():
                    col_keys.append(col[])
                sort(col_keys)

                # Add non-zero elements from this row
                for col in col_keys:
                    var val: ComplexScalar[Self.type] = sparse_dict[row][col[]]
                    if val > Self.zero_threshold:
                        self.col_idx.append(col[])
                        self.v.append(val)

            # Update row index pointer
            self.row_idx.append(self.col_idx.size)
    
    fn __init__(inout self, dense: CMatrix[Self.type]):
        '''Initialize a CSRCMatrix from a dense CMatrix.'''
        self.rows = dense.rows
        self.cols = dense.cols
        self.size = dense.size
        self.v = List[ComplexScalar[Self.type], True]()
        self.col_idx = List[Int, True]()
        self.row_idx = List[Int, True](capacity=dense.rows + 1)
        self.row_idx.append(0)

        for row in range(self.rows):
            for col in range(self.cols):
                var val: ComplexScalar[Self.type] = dense.load_crd[1](row, col)
                if val > Self.zero_threshold:
                    self.v.append(val)
                    self.col_idx.append(col)
            self.row_idx.append(self.v.size)

    # Static constructors #############

    @staticmethod
    fn zeros(rows: Int, cols: Int) -> Self:
        '''Return a sparse matrix of zeros of the specified size.'''
        return Self(rows, cols)

    @staticmethod
    fn eye(rows: Int, cols: Int) raises -> Self: ###########################################################################
        '''Return a sparse identity matrix of the specified size.'''
        var result = Self(rows, cols)
        result.fill_diag(ComplexScalar[Self.type](1))
        return result

    # Conversion operations and helpers

    @always_inline
    def to_dense(self) -> CMatrix[Self.type]:
        '''Convert self to a dense CMatrix.'''
        var dense = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
        for row in range(self.rows):
            for k in range(self.row_idx[row], self.row_idx[row + 1]):
                dense[row, self.col_idx[k]] = self.v[k]
        return dense

    fn _row_idx_to_indices(self) -> List[Int, True]:
        '''Convert the list of index pointers to a list of row indices.'''
        var row_indices = List[Int, True]()
        for row in range(self.rows):
            var row_start: Int = self.row_idx[row]
            var row_end: Int = self.row_idx[row + 1]
            if row_start == row_end:
                continue
            elif row_end - row_start == 1:
                row_indices.append(row)
            else:
                row_indices.extend(List[Int, True](row) * (row_end - row_start))
        return row_indices
    
    @staticmethod
    fn _sorted_row_indices_to_idx(
        rows: Int, owned row_indices: List[Int, True]
    ) -> List[Int, True]:
        '''Convert a list of sorted row indices to a list of index pointers.'''
        var row_idx = List[Int, True](0) * (rows + 1)
        # Count the number of elements in each row
        for idx in row_indices:
            row_idx[idx[] + 1] += 1
        # Compute the cumulative sum of row_idx
        for i in range(1, rows + 1):
            row_idx[i] += row_idx[i - 1]
        return row_idx^

    # String conversion ###############

    @no_inline
    fn __str__(self) -> String:
        return (
            '<Compressed Sparse Row complex matrix with ' + str(self.n_stored())
            + ' stored elements of type ComplexScalar[' + str(Self.type) 
            + '] and shape (' + str(self.rows) + ', ' + str(self.cols) + ')>'
        )

    @no_inline
    fn __repr__(self) -> String:
        return self.__str__()
    
    @no_inline
    fn format_to(self, inout writer: Formatter):
        writer.write(str(self))
    
    # Assertions ######################

    @always_inline
    fn _assert_same_shape(self, other: Self) raises:
        if self.rows != other.rows or self.cols != other.cols:
            raise Error(
                'Incompatible matrix dimensions: (' 
                + str(self.rows) 
                + ', ' 
                + str(self.cols) 
                + ') and (' 
                + str(other.rows) 
                + ', ' 
                + str(other.cols) 
                + ')'
            )
    
    @always_inline
    fn _assert_same_shape(self, other: CMatrix[Self.type]) raises:
        if self.rows != other.rows or self.cols != other.cols:
            raise Error(
                'Incompatible matrix dimensions: (' 
                + str(self.rows) 
                + ', ' 
                + str(self.cols) 
                + ') and (' 
                + str(other.rows) 
                + ', ' 
                + str(other.cols) 
                + ')'
            )
    
    @always_inline
    fn _assert_reshape_compatible(self, new_rows: Int, new_cols: Int) raises:
        if new_rows * new_cols != self.size:
            raise Error(
                'Impatible reshape dimensions: ('
                + str(self.rows) 
                + ', ' 
                + str(self.cols) 
                + ') and (' 
                + str(new_rows) 
                + ', ' 
                + str(new_cols) 
                + ')'
            )

    # Properties ######################

    @always_inline
    fn __len__(self) -> Int:
        '''Return the total number of elements in the matrix, including the zeros.'''
        return self.size
    
    @always_inline
    fn shape(self) -> Tuple[Int, Int]:
        '''Return the tuple (rows, cols).'''
        return (self.rows, self.cols)
    
    @always_inline
    fn n_nonzero(self) -> Int:
        '''Returns the number of non-zero elements in the matrix.'''
        var nonzero: Int = self.v.size
        for v in self.v:
            if v[] <= Self.zero_threshold:
                nonzero -= 1
        return nonzero
    
    @always_inline
    fn n_stored(self) -> Int:
        '''Returns the number of elements stored in the matrix.'''
        return self.v.size

    @always_inline
    fn is_square(self) -> Bool:
        '''Return True if self is square and non-degenerate (rows, cols > 0), 
        and False otherwise.
        '''
        return self.rows == self.cols and self.rows > 0
    
    # @always_inline
    # fn is_unitary[tol: Scalar[Self.type] = DEFAULT_TOL](self) raises -> Bool:
    #     '''Return True if self is unitary, False otherwise.'''
    #     if not self.is_square():
    #         return False
    #     return (self @ self.dag() - self.eye_like()).frobenius_norm() < tol

    # @always_inline
    # fn is_hermitian[tol: Scalar[Self.type] = DEFAULT_TOL](self) raises -> Bool:
    #     '''Return True if self is Hermitian, False otherwise.'''
    #     if not self.is_square():
    #         return False
    #     return (self - self.dag()).frobenius_norm() < tol
    
    # Item access #####################

    @always_inline
    fn __getitem__(self, row: Int, col: Int) raises -> ComplexScalar[Self.type]:
        '''Get the value at the specified row and col.'''
        if row < 0 or row >= self.rows:
            raise Error('Invalid row index: ' + str(row))
        if col < 0 or col >= self.cols:
            raise Error('Invalid column index: ' + str(col))
        return self._getitem_noraise(row, col)
    
    @always_inline
    fn _getitem_linear_noraise(self, row: Int, col: Int) -> ComplexScalar[Self.type]:
        '''Get the value at the specified row and col without checking index validity.'''
        for i in range(self.row_idx[row], self.row_idx[row + 1]):
            if self.col_idx[i] == col:
                return self.v[i]
        return ComplexScalar[Self.type]()

    @always_inline
    fn _getitem_noraise(self, row: Int, col: Int) -> ComplexScalar[Self.type]:
        '''Get the value at the specified row and col without checking index validity.'''
        var left = self.row_idx[row]
        var right = self.row_idx[row + 1] - 1

        # Search for the column in the row with a binary search
        while left <= right:
            var mid: Int = (left + right) // 2
            var curr_col: Int = self.col_idx[mid]
            if curr_col == col:
                return self.v[mid]
            elif curr_col < col:
                left = mid + 1
            else:
                right = mid - 1
        return ComplexScalar[Self.type](0)
    
    @always_inline
    fn __setitem__(inout self, row: Int, col: Int, val: ComplexScalar[Self.type]) raises:
        '''Set the value at the specified row and col.'''
        if row < 0 or row >= self.rows:
            raise Error('Invalid row index: ' + str(row))
        if col < 0 or col >= self.cols:
            raise Error('Invalid column index: ' + str(col))
        self._setitem_noraise(row, col, val)
    
    fn _setitem_linear_noraise(inout self, row: Int, col: Int, val: ComplexScalar[Self.type]):
        '''Set the value at the specified row and col without checking index validity.'''
        var start: Int = self.row_idx[row]
        var end: Int = self.row_idx[row + 1]

        # Search for the column in the row
        for i in range(start, end):
            if self.col_idx[i] == col:
                if val <= Self.zero_threshold:
                    # Remove element if value is 0
                    _ = self.v.pop(i)
                    _ = self.col_idx.pop(i)
                    # Update row_idx accordingly
                    for j in range(row + 1, self.row_idx.size):
                        self.row_idx[j] -= 1
                else:
                    self.v[i] = val
                return

        # If value is zero, do nothing
        if val <= Self.zero_threshold:
            return

        # Find the correct insertion point to maintain sorted column indices
        while start < end and self.col_idx[start] < col:
            start += 1
        self.v.insert(start, val)
        self.col_idx.insert(start, col)

        # Update row_idx for subsequent rows
        for i in range(row + 1, self.row_idx.size):
            self.row_idx[i] += 1

    fn _setitem_noraise(inout self, row: Int, col: Int, val: ComplexScalar[Self.type]):
        '''Set the value at the specified row and col without checking index validity.'''
        var start: Int = self.row_idx[row]
        var end: Int = self.row_idx[row + 1]
        var left: Int = start
        var right: Int = end - 1

        # Search for the column in the row with a binary search
        while left <= right:
            var mid: Int = (left + right) // 2
            var curr_col: Int = self.col_idx[mid]
            if curr_col == col:
                if val <= Self.zero_threshold:
                    # Remove element if value is 0
                    _ = self.v.pop(mid)
                    _ = self.col_idx.pop(mid)
                    # Update row_idx accordingly
                    for j in range(row + 1, self.row_idx.size):
                        self.row_idx[j] -= 1
                else:
                    self.v[mid] = val
                return
            elif curr_col < col:
                left = mid + 1
            else:
                right = mid - 1

        # If value is zero, do nothing
        if val <= Self.zero_threshold:
            return

        # Insertion point is now at `left`
        self.v.insert(left, val)
        self.col_idx.insert(left, col)

        # Update row_idx for subsequent rows
        for i in range(row + 1, self.row_idx.size):
            self.row_idx[i] += 1

    fn __contains__(self, coords: Tuple[Int, Int]) -> Bool:
        '''Returns True if self[r, c] is stored, False otherwise.'''
        row, col = coords
        if row < 0 or row >= self.rows or col < 0 or col >= self.cols:
            return False
        
        # Binary search for the column in the row's indices
        left = self.row_idx[row]
        right = self.row_idx[row + 1]
        while left < right:
            var mid: Int = (left + right) // 2
            mid_col = self.col_idx[mid]
            if mid_col == col:
                return True
            elif mid_col < col:
                left = mid + 1
            else:
                right = mid
        return False

    @always_inline
    fn extract_row_as_sparse(self, row: Int) raises -> Self:
        '''Extract the specified row from the sparse matrix and return it as a sparse matrix.'''
        if row < 0 or row >= self.rows:
            raise Error('Invalid row index: ' + str(row))
        var row_start: Int = self.row_idx[row]
        var row_end: Int = self.row_idx[row + 1]
        return Self(
            rows=1,
            cols=self.cols,
            size=self.cols,
            v=self.v[row_start:row_end],
            col_idx=self.col_idx[row_start:row_end],
            row_idx=List[Int, True](0, row_end - row_start),
        )

    @always_inline
    fn extract_row_as_matrix(self, row: Int) raises -> CMatrix[Self.type]:
        '''Extract the specified row from the sparse matrix and return it as a CMatrix.'''
        if row < 0 or row >= self.rows:
            raise Error('Invalid row index: ' + str(row))
        var result = CMatrix[Self.type](1, self.cols, fill_zeros=True)
        var row_ptr: Int = self.row_idx[row]
        var row_end: Int = self.row_idx[row + 1]
        for c in range(self.cols):
            if row_ptr < row_end and c == self.col_idx[row_ptr]:
                result.store_idx[1](c, self.v[row_ptr])
                row_ptr += 1
        return result
    
    fn extract_column_as_sparse(self, col: Int) raises -> Self:
        '''Extract the specified column from the sparse matrix and return it as a sparse matrix.'''
        if col < 0 or col >= self.cols:
            raise Error('Invalid column index: ' + str(col))
        var result = Self(self.rows, 1)
        for row in range(self.rows):
            var row_start: Int = self.row_idx[row]
            var row_end: Int = self.row_idx[row + 1]
            for i in range(row_start, row_end):
                if self.col_idx[i] == col:
                    result.v.append(self.v[i])
                    result.col_idx.append(0)
                    break
            result.row_idx[row + 1] = result.v.size
        return result
    
    fn extract_column_as_matrix(self, col: Int) raises -> CMatrix[Self.type]:
        '''Extract the specified column from the sparse matrix and return it as a CMatrix.'''
        if col < 0 or col >= self.cols:
            raise Error('Invalid column index: ' + str(col))
        var result = CMatrix[Self.type](self.rows, 1, fill_zeros=True)
        for row in range(self.rows):
            var row_start: Int = self.row_idx[row]
            var row_end: Int = self.row_idx[row + 1]
            for i in range(row_start, row_end):
                if self.col_idx[i] == col:
                    result.store_idx[1](row, self.v[i])
                    break
        return result
    
    # Math dunders ####################

    @always_inline
    fn __neg__(self) -> Self:
        '''Defines the `-` unary negation operator. Returns -self.'''
        var result: Self = self
        for i in range(self.v.size):
            result.v[i] = -self.v[i]
        return result
    
    @always_inline
    fn __pos__(self) -> Self:
        '''Defines the `+` unary positive operator. Returns self.'''
        return self
    
    fn __add__(self, other: Self) raises -> Self:
        '''Defines the `+` add operator. Returns self + other.'''
        self._assert_same_shape(other)
        # Initialize each attribute directly to avoid populating row_idx
        var result = Self(
            rows=self.rows, 
            cols=self.cols, 
            size=self.rows * self.cols, 
            v=List[ComplexScalar[Self.type], True](), 
            col_idx=List[Int, True](), 
            row_idx=List[Int, True](0),
        )
        for row in range(self.rows):
            var self_start: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            var other_start: Int = other.row_idx[row]
            var other_end: Int = other.row_idx[row + 1]

            # Lists to store non-zero elements for this row
            var result_row_vals = List[ComplexScalar[Self.type], True]()
            var result_row_cols = List[Int, True]()

            # Pointers for both matrices
            var self_ptr: Int = self_start
            var other_ptr: Int = other_start

            # Loop until both pointers reach the end of their rows
            while self_ptr < self_end and other_ptr < other_end:
                var self_col: Int = self.col_idx[self_ptr]
                var other_col: Int = other.col_idx[other_ptr]
                if self_col < other_col:
                    # Add value from self matrix
                    result_row_vals.append(self.v[self_ptr])
                    result_row_cols.append(self_col)
                    self_ptr += 1
                elif self_col > other_col:
                    # Add value from other matrix
                    result_row_vals.append(other.v[other_ptr])
                    result_row_cols.append(other_col)
                    other_ptr += 1
                else:
                    # Add values from both matrices
                    var sum_val: ComplexScalar[Self.type] = self.v[self_ptr] + other.v[other_ptr]
                    if sum_val > Self.zero_threshold:
                        result_row_vals.append(sum_val)
                        result_row_cols.append(self_col)
                    self_ptr += 1; other_ptr += 1

            # Add remaining elements from self matrix
            while self_ptr < self_end:
                result_row_vals.append(self.v[self_ptr])
                result_row_cols.append(self.col_idx[self_ptr])
                self_ptr += 1

            # Add remaining elements from other matrix
            while other_ptr < other_end:
                result_row_vals.append(other.v[other_ptr])
                result_row_cols.append(other.col_idx[other_ptr])
                other_ptr += 1

            # Update result matrix with this row's values
            if result_row_vals:
                result.v.extend(result_row_vals)
                result.col_idx.extend(result_row_cols)
            result.row_idx.append(result.v.size)
        return result
    
    fn __add__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the `+` add operator. Returns a CMatrix of self + other.'''
        self._assert_same_shape(other)
        var result: CMatrix[Self.type] = other
        for row in range(self.rows):
            var self_start: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            for i in range(self_start, self_end):
                result.store_crd[1](
                    row, self.col_idx[i], result.load_crd[1](row, self.col_idx[i]) + self.v[i]
                )
        return result

    # # Parallelization may add more overhead than value
    # fn __add__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
    #     '''Defines the `+` add operator. Returns a CMatrix of self + other.'''
    #     self._assert_same_shape(other)
    #     var result: CMatrix[Self.type] = other
    #     @parameter
    #     fn add_row(row: Int):
    #         var self_start: Int = self.row_idx[row]
    #         var self_end: Int = self.row_idx[row + 1]
    #         for i in range(self_start, self_end):
    #             result.store_crd[1](
    #                 row, self.col_idx[i], result.load_crd[1](row, self.col_idx[i]) + self.v[i]
    #             )
    #     parallelize[add_row](self.rows, self.rows)
    #     return result

    @always_inline
    fn __add__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Defines the `+` add operator. Returns self if other == 0, raises otherwise.
        Adding a non-zero scalar to a sparse matrix is not supported.
        '''
        if other <= Self.zero_threshold:
            return self
        raise Error('Adding a non-zero scalar to a sparse matrix is not supported')
    
    @always_inline
    fn __radd__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the right `+` add operator. Returns other + self.'''
        return self + other
    
    @always_inline
    fn __radd__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Defines the right `+` add operator. Returns self if other == 0, raises otherwise.
        Adding a non-zero scalar to a sparse matrix is not supported.
        '''
        if other <= Self.zero_threshold:
            return self
        raise Error('Adding a non-zero scalar to a sparse matrix is not supported')
    
    @always_inline
    fn __iadd__(inout self, other: Self) raises:
        '''Defines the `+=` in-place add operator.'''
        self = self + other

    @always_inline
    fn __iadd__(inout self, other: CMatrix[Self.type]) raises:
        '''Raises. In-place addition of a sparse matrix and a matrix is not supported.'''
        raise Error('In-place addition of a sparse matrix and a matrix is not supported.')
    
    @always_inline
    fn __iadd__(inout self, other: ComplexScalar[Self.type]) raises:
        '''Defines the `+=` in-place add operator. Raises if other is non-zero.'''
        if other <= Self.zero_threshold:
            return
        raise Error('In-place addition of a sparse matrix and a scalar is not supported')

    fn __sub__(self, other: Self) raises -> Self:
        '''Defines the `-` subtraction operator. Returns self - other.'''
        self._assert_same_shape(other)
        # Initialize each attribute directly to avoid populating row_idx
        var result = Self(
            rows=self.rows, 
            cols=self.cols, 
            size=self.rows * self.cols, 
            v=List[ComplexScalar[Self.type], True](), 
            col_idx=List[Int, True](), 
            row_idx=List[Int, True](0),
        )
        for row in range(self.rows):
            var self_start: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            var other_start: Int = other.row_idx[row]
            var other_end: Int = other.row_idx[row + 1]

            # Lists to store non-zero elements for this row
            var result_row_vals = List[ComplexScalar[Self.type], True]()
            var result_row_cols = List[Int, True]()

            # Pointers for both matrices
            var self_ptr: Int = self_start
            var other_ptr: Int = other_start

            # Loop until both pointers reach the end of their rows
            while self_ptr < self_end and other_ptr < other_end:
                var self_col: Int = self.col_idx[self_ptr]
                var other_col: Int = other.col_idx[other_ptr]
                if self_col < other_col:
                    # Add value from self matrix
                    result_row_vals.append(self.v[self_ptr])
                    result_row_cols.append(self_col)
                    self_ptr += 1
                elif self_col > other_col:
                    # Subtract value from other matrix
                    result_row_vals.append(-other.v[other_ptr])
                    result_row_cols.append(other_col)
                    other_ptr += 1
                else:
                    # Add difference of values from both matrices
                    var diff_val: ComplexScalar[Self.type] = self.v[self_ptr] - other.v[other_ptr]
                    if diff_val > Self.zero_threshold:
                        result_row_vals.append(diff_val)
                        result_row_cols.append(self_col)
                    self_ptr += 1; other_ptr += 1

            # Add remaining elements from self matrix
            while self_ptr < self_end:
                result_row_vals.append(self.v[self_ptr])
                result_row_cols.append(self.col_idx[self_ptr])
                self_ptr += 1

            # Subtract remaining elements from other matrix
            while other_ptr < other_end:
                result_row_vals.append(-other.v[other_ptr])
                result_row_cols.append(other.col_idx[other_ptr])
                other_ptr += 1

            # Update result matrix with this row's values
            if result_row_vals:
                result.v.extend(result_row_vals)
                result.col_idx.extend(result_row_cols)
            result.row_idx.append(result.v.size)
        return result

    fn __sub__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the `-` subtraction operator. Returns a CMatrix of self - other.'''
        self._assert_same_shape(other)
        var result: CMatrix[Self.type] = -other
        for row in range(self.rows):
            var self_start: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            for i in range(self_start, self_end):
                result.store_crd[1](
                    row, self.col_idx[i], result.load_crd[1](row, self.col_idx[i]) + self.v[i]
                )
        return result

    # # Parallelization may add more overhead than value
    # fn __sub__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
    #     '''Defines the `-` subtraction operator. Returns a CMatrix of self - other.'''
    #     self._assert_same_shape(other)
    #     var result: CMatrix[Self.type] = -other
    #     @parameter
    #     fn add_row(row: Int):
    #         var self_start: Int = self.row_idx[row]
    #         var self_end: Int = self.row_idx[row + 1]
    #         for i in range(self_start, self_end):
    #             result.store_crd[1](
    #                 row, self.col_idx[i], result.load_crd[1](row, self.col_idx[i]) + self.v[i]
    #             )
    #     parallelize[add_row](self.rows, self.rows)
    #     return result

    @always_inline
    fn __sub__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Defines the `-` subtraction operator. Returns self if other == 0, raises otherwise.
        Subtracting a non-zero scalar from a sparse matrix is not supported.
        '''
        if other <= Self.zero_threshold:
            return self
        raise Error('Subtracting a non-zero scalar from a sparse matrix is not supported')
    
    @always_inline
    fn __rsub__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the right `-` subtraction operator. Returns other - self.'''
        # TODO: This could be improved
        return -self + other
    
    @always_inline
    fn __rsub__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Defines the right `-` subtraction operator. Returns -self if other == 0, raises
        otherwise. Subtracting a sparse matrix from a non-zero scalar is not supported.
        '''
        if other <= Self.zero_threshold:
            return -self
        raise Error('Subtracting a sparse matrix from a non-zero scalar is not supported.')
    
    @always_inline
    fn __isub__(inout self, other: Self) raises:
        '''Defines the `-=` in-place subtraction operator.'''
        self = self - other

    @always_inline
    fn __isub__(inout self, other: CMatrix[Self.type]) raises:
        '''Raises. In-place subtraction of a matrix from a sparse matrix is not supported.'''
        raise Error('In-place subtraction of a matrix from a sparse matrix is not supported')
    
    @always_inline
    fn __isub__(inout self, other: ComplexScalar[Self.type]) raises:
        '''Defines the `-=` in-place subtraction operator. Raises if other is non-zero.'''
        if other <= Self.zero_threshold:
            return
        raise Error('In-place subtraction of a scalar from a sparse matrix is not supported')

    fn __mul__(self, other: Self) raises -> Self:
        '''Defines the `*` product operator. Returns self * other.'''
        self._assert_same_shape(other)
        # Initialize each attribute directly to avoid populating row_idx
        var result = Self(
            rows=self.rows, 
            cols=self.cols, 
            size=self.rows * self.cols, 
            v=List[ComplexScalar[Self.type], True](), 
            col_idx=List[Int, True](), 
            row_idx=List[Int, True](0),
        )
        for row in range(self.rows):
            var self_start: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            var other_start: Int = other.row_idx[row]
            var other_end: Int = other.row_idx[row + 1]

            # Lists to store non-zero elements for this row
            var result_row_vals = List[ComplexScalar[Self.type], True]()
            var result_row_cols = List[Int, True]()

            # Pointers for both matrices
            var self_ptr: Int = self_start
            var other_ptr: Int = other_start
            while self_ptr < self_end and other_ptr < other_end:
                var self_col: Int = self.col_idx[self_ptr]
                var other_col: Int = other.col_idx[other_ptr]
                if self_col < other_col:
                    self_ptr += 1
                elif self_col > other_col:
                    other_ptr += 1
                else:
                    # Multiply values from both matrices
                    var prod_val: ComplexScalar[Self.type] = self.v[self_ptr] * other.v[other_ptr]
                    if prod_val > Self.zero_threshold:
                        result_row_vals.append(prod_val)
                        result_row_cols.append(self_col)
                    self_ptr += 1; other_ptr += 1

            # Update result matrix with this row's values
            if result_row_vals:
                result.v.extend(result_row_vals)
                result.col_idx.extend(result_row_cols)
            result.row_idx.append(result.v.size)
        return result
    
    fn __mul__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the `*` multiplication operator. Returns a CMatrix of self * other.'''
        self._assert_same_shape(other)
        var result = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
        for row in range(self.rows):
            var self_start: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            for i in range(self_start, self_end):
                result.store_crd[1](
                    row, self.col_idx[i], other.load_crd[1](row, self.col_idx[i]) * self.v[i]
                )
        return result

    # # Parallelization may add more overhead than value
    # fn __mul__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
    #     '''Defines the `*` multiplication operator. Returns a CMatrix of self * other.'''
    #     self._assert_same_shape(other)
    #     var result = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
    #     @parameter
    #     fn mul_row(row: Int):
    #         var self_start: Int = self.row_idx[row]
    #         var self_end: Int = self.row_idx[row + 1]
    #         for i in range(self_start, self_end):
    #             result.store_crd[1](
    #                 row, self.col_idx[i], other.load_crd[1](row, self.col_idx[i]) * self.v[i]
    #             )
    #     parallelize[mul_row](self.rows, self.rows)
    #     return result

    fn __mul__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `*` product operator. Returns self * other.'''
        # This method fails if other is NaN. However, following the
        # conventions of scipy.sparse, we ignore this case.
        if other == 0:
            return Self(self.rows, self.cols)
        var result = Self(self.rows, self.cols)
        for row in range(self.rows):
            for i in range(self.row_idx[row], self.row_idx[row + 1]):
                var val: ComplexScalar[Self.type] = self.v[i] * other
                if val > Self.zero_threshold:
                    result.v.append(val)
                    result.col_idx.append(self.col_idx[i])
            result.row_idx[row + 1] = result.v.size
        return result
    
    @always_inline
    fn __rmul__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the right `*` multiplication operator. Returns other * self.'''
        return self * other
    
    @always_inline
    fn __rmul__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Defines the right `*` multiplication operator. Returns other * self.'''
        return self * other

    @always_inline
    fn __imul__(inout self, other: Self) raises:
        '''Defines the `*=` in-place multiplication operator.'''
        self = self * other

    @always_inline
    fn __imul__(inout self, other: CMatrix[Self.type]) raises:
        '''Raises. In-place multiplication of a sparse matrix and a matrix is not supported.'''
        raise Error('In-place multiplication of a sparse matrix and a matrix is not supported')
    
    @always_inline
    fn __imul__(inout self, other: ComplexScalar[Self.type]):
        '''Defines the `*=` in-place multiplication operator.'''
        self = self * other

    @always_inline
    fn __truediv__(self, other: Self) raises -> Self:
        '''Raises. Dividing a sparse matrix by a sparse matrix is not supported.'''
        raise Error('Dividing a sparse matrix by a sparse matrix is not supported')

    @always_inline
    fn __truediv__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `/` divide operator. Returns self / other.'''
        return self * other.reciprocal()

    fn __truediv__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the `/` divide operator. Returns a CMatrix of self / other.'''
        self._assert_same_shape(other)
        alias _csnan = ComplexScalar[type](nan[Self.type](), nan[Self.type]())
        var result = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
        for row in range(self.rows):
            var row_ptr: Int = self.row_idx[row]
            var row_end: Int = self.row_idx[row + 1]
            for c in range(self.cols):
                if row_ptr < row_end and c == self.col_idx[row_ptr]:
                    result.store_crd[1](row, c, self.v[row_ptr] / other.load_crd[1](row, c))
                    row_ptr += 1
                elif other.load_crd[1](row, c) == 0:
                    result.store_crd[1](row, c, _csnan)
        return result

    # # Parallelization may add more overhead than value
    # fn __truediv__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
    #     '''Defines the `/` divide operator. Returns a CMatrix of self / other.'''
    #     self._assert_same_shape(other)
    #     alias _csnan = ComplexScalar[type](nan[Self.type](), nan[Self.type]())
    #     var result = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
    #     @parameter
    #     fn div_row(row: Int):
    #         var row_ptr: Int = self.row_idx[row]
    #         var row_end: Int = self.row_idx[row + 1]
    #         for c in range(self.cols):
    #             if row_ptr < row_end and c == self.col_idx[row_ptr]:
    #                 result.store_crd[1](row, c, self.v[row_ptr] / other.load_crd[1](row, c))
    #                 row_ptr += 1
    #             elif other.load_crd[1](row, c) == 0:
    #                 result.store_crd[1](row, c, _csnan)
    #     parallelize[div_row](self.rows, self.rows)
    #     return result
    
    @always_inline
    fn __rtruediv__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Raises. Dividing a matrix by a sparse matrix is not supported.'''
        raise Error('Dividing a matrix by a sparse matrix is not supported')
    
    @always_inline
    fn __rtruediv__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Raises. Dividing a scalar by a sparse matrix is not supported.'''
        raise Error('Dividing a scalar by a sparse matrix is not supported')

    @always_inline
    fn __itruediv__(inout self, other: Self) raises:
        '''Defines the `/=` in-place divide operator.'''
        self = self / other

    @always_inline
    fn __itruediv__(inout self, other: CMatrix[Self.type]) raises:
        '''Raises. In-place division of a sparse matrix and a matrix is not supported.'''
        raise Error('In-place division of a sparse matrix and a matrix is not supported')
    
    @always_inline
    fn __itruediv__(inout self, other: ComplexScalar[Self.type]):
        '''Defines the `/=` in-place divide operator.'''
        self *= other.reciprocal()
    
    @always_inline
    fn __floordiv__(self, other: Self) raises -> Self:
        '''Raises. Dividing a sparse matrix by a sparse matrix is not supported.'''
        raise Error('Dividing a sparse matrix by a sparse matrix is not supported')

    fn __floordiv__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `//` floor divide operator. Returns self // other.'''
        var result = Self(self.rows, self.cols)
        for row in range(self.rows):
            for i in range(self.row_idx[row], self.row_idx[row + 1]):
                var val: ComplexScalar[Self.type] = self.v[i] // other
                if val > Self.zero_threshold:
                    result.v.append(val)
                    result.col_idx.append(self.col_idx[i])
            result.row_idx[row + 1] = result.v.size
        return result

    fn __floordiv__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the `//` floor divide operator. Returns a CMatrix of self // other.'''
        self._assert_same_shape(other)
        alias _csnan = ComplexScalar[type](nan[Self.type](), nan[Self.type]())
        var result = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
        for row in range(self.rows):
            var row_ptr: Int = self.row_idx[row]
            var row_end: Int = self.row_idx[row + 1]
            for c in range(self.cols):
                if row_ptr < row_end and c == self.col_idx[row_ptr]:
                    result.store_crd[1](row, c, self.v[row_ptr] // other.load_crd[1](row, c))
                    row_ptr += 1
                elif other.load_crd[1](row, c) == 0:
                    result.store_crd[1](row, c, _csnan)
        return result

    # # Parallelization may add more overhead than value
    # fn __floordiv__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
    #     '''Defines the `//` floor divide operator. Returns a CMatrix of self // other.'''
    #     self._assert_same_shape(other)
    #     alias _csnan = ComplexScalar[type](nan[Self.type](), nan[Self.type]())
    #     var result = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
    #     @parameter
    #     fn div_row(row: Int):
    #         var row_ptr: Int = self.row_idx[row]
    #         var row_end: Int = self.row_idx[row + 1]
    #         for c in range(self.cols):
    #             if row_ptr < row_end and c == self.col_idx[row_ptr]:
    #                 result.store_crd[1](row, c, self.v[row_ptr] // other.load_crd[1](row, c))
    #                 row_ptr += 1
    #             elif other.load_crd[1](row, c) == 0:
    #                 result.store_crd[1](row, c, _csnan)
    #     parallelize[div_row](self.rows, self.rows)
    #     return result

    @always_inline
    fn __rfloordiv__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Raises. Dividing a matrix by a sparse matrix is not supported.'''
        raise Error('Dividing a matrix by a sparse matrix is not supported')
    
    @always_inline
    fn __rfloordiv__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Raises. Dividing a scalar by a sparse matrix is not supported.'''
        raise Error('Dividing a scalar by a sparse matrix is not supported')

    @always_inline
    fn __ifloordiv__(inout self, other: Self) raises:
        '''Defines the `//=` in-place floor divide operator.'''
        self = self // other

    @always_inline
    fn __ifloordiv__(inout self, other: CMatrix[Self.type]) raises:
        '''Raises. In-place floor division of a sparse matrix and a matrix is not supported.'''
        raise Error('In-place floor division of a sparse matrix and a matrix is not supported')
    
    @always_inline
    fn __ifloordiv__(inout self, other: ComplexScalar[Self.type]):
        '''Defines the `//=` in-place divide operator.'''
        self = self // other
    
    @always_inline
    fn __mod__(self, other: Self) raises -> Self:
        '''Raises. The modulo operation between sparse matrices is not supported.'''
        raise Error('The modulo operation between sparse matrices is not supported')
    
    fn __mod__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `%` mod operator. Returns self % other.'''
        var result = Self(self.rows, self.cols)
        for row in range(self.rows):
            for i in range(self.row_idx[row], self.row_idx[row + 1]):
                var val: ComplexScalar[Self.type] = self.v[i] % other
                if val > Self.zero_threshold:
                    result.v.append(val)
                    result.col_idx.append(self.col_idx[i])
            result.row_idx[row + 1] = result.v.size
        return result

    fn __mod__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the `%` floor divide operator. Returns a CMatrix of self % other.'''
        self._assert_same_shape(other)
        alias _csnan = ComplexScalar[type](nan[Self.type](), nan[Self.type]())
        var result = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
        for row in range(self.rows):
            var row_ptr: Int = self.row_idx[row]
            var row_end: Int = self.row_idx[row + 1]
            for c in range(self.cols):
                if row_ptr < row_end and c == self.col_idx[row_ptr]:
                    result.store_crd[1](row, c, self.v[row_ptr] % other.load_crd[1](row, c))
                    row_ptr += 1
                elif other.load_crd[1](row, c) == 0:
                    result.store_crd[1](row, c, _csnan)
        return result

    # # Parallelization may add more overhead than value
    # fn __mod__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
    #     '''Defines the `%` floor divide operator. Returns a CMatrix of self % other.'''
    #     self._assert_same_shape(other)
    #     alias _csnan = ComplexScalar[type](nan[Self.type](), nan[Self.type]())
    #     var result = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
    #     @parameter
    #     fn div_row(row: Int):
    #         var row_ptr: Int = self.row_idx[row]
    #         var row_end: Int = self.row_idx[row + 1]
    #         for c in range(self.cols):
    #             if row_ptr < row_end and c == self.col_idx[row_ptr]:
    #                 result.store_crd[1](row, c, self.v[row_ptr] % other.load_crd[1](row, c))
    #                 row_ptr += 1
    #             elif other.load_crd[1](row, c) == 0:
    #                 result.store_crd[1](row, c, _csnan)
    #     parallelize[div_row](self.rows, self.rows)
    #     return result

    @always_inline
    fn __rmod__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Raises. The modulo operation between a matrix and sparse matrix is not supported.'''
        raise Error('The modulo operation between a matrix and sparse matrix is not supported')
    
    @always_inline
    fn __rmod__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Raises. The modulo operation between a scalar and sparse matrix is not supported.'''
        raise Error('The modulo operation between a scalar and sparse matrix is not supported')

    @always_inline
    fn __imod__(inout self, other: Self) raises:
        '''Defines the `%=` in-place modulo operator.'''
        self = self % other

    @always_inline
    fn __imod__(inout self, other: CMatrix[Self.type]) raises:
        '''Raises. The in-place modulo operation of a sparse matrix and a matrix is not 
        supported.
        '''
        raise Error(
            'The in-place modulo operation of a sparse matrix and a matrix is not supported'
        )
    
    @always_inline
    fn __imod__(inout self, other: ComplexScalar[Self.type]):
        '''Defines the `//=` in-place divide operator.'''
        self = self % other

    @always_inline
    fn __divmod__(self, other: Self) raises -> Tuple[Self, Self]:
        '''Raises. The divmod operation between sparse matrices is not supported.'''
        raise Error('The divmod operation between sparse matrices is not supported')
    
    @always_inline
    fn __divmod__(self, other: ComplexScalar[Self.type]) -> Tuple[Self, Self]:
        '''Defines the divmod operator. Returns (self // other, self % other).'''
        return (self // other, self % other)

    @always_inline
    fn __divmod__(
        self, other: CMatrix[Self.type]
    ) raises -> Tuple[CMatrix[Self.type], CMatrix[Self.type]]:
        '''Defines the divmod operator. Returns a tuple of CMatrix 
        (self // other, self % other).
        '''
        return (self // other, self % other)
    
    @always_inline
    fn __rdivmod__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Raises. The divmod operation between a matrix and sparse matrix is not supported.'''
        raise Error('The divmod operation between a matrix and sparse matrix is not supported')
    
    @always_inline
    fn __rdivmod__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Raises. The divmod operation between a scalar and sparse matrix is not supported.'''
        raise Error('The divmod operation between a scalar and sparse matrix is not supported')
    
    fn _dense_vec_matmul(self, vec: CMatrix[Self.type]) -> CMatrix[Self.type]:
        '''Matrix vector multiplication of a sparse matrix with a dense vector.'''
        var result = CMatrix[type](self.rows, 1, fill_zeros=True)
        @parameter
        fn dot(row: Int):
            for i in range(self.row_idx[row], self.row_idx[row + 1]):
                result.store_idx[1](
                    row, result.load_idx[1](row) + self.v[i] * vec.load_idx[1](self.col_idx[i])
                )
        parallelize[dot](self.rows, self.rows)
        return result

    fn _dense_mat_matmul(self, mat: CMatrix[Self.type]) -> CMatrix[Self.type]:
        '''Matrix multiplication of a sparse matrix with a dense matrix.'''
        var result = CMatrix[type](self.rows, mat.cols, fill_zeros=True)
        @parameter
        fn dot(row: Int):
            for col in range(mat.cols):
                for i in range(self.row_idx[row], self.row_idx[row + 1]):
                    result.store_crd[1](
                        row,
                        col,
                        result.load_crd[1](row, col) 
                        + self.v[i] 
                        * mat.load_crd[1](self.col_idx[i], col),
                    )
        parallelize[dot](self.rows, self.rows)
        return result
    
    fn _sparse_matmul_gustavson(self, other: Self) raises -> Self:
        '''Matrix multiplication of a sparse matrix with a sparse matrix using 
        Gustavson's algorithm.
        '''
        var result_v = List[ComplexScalar[type], True]()
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)

        # Iterate through the rows of self
        for i in range(self.rows):
            workspace = Dict[Int, ComplexScalar[Self.type]]()

            # Iterate through non-zero elements in current row of self
            for j in range(self.row_idx[i], self.row_idx[i + 1]):
                var self_val: ComplexScalar[Self.type] = self.v[j]
                var self_col: Int = self.col_idx[j]

                # Iterate through non-zero elements in corresponding row of other
                for k in range(other.row_idx[self_col], other.row_idx[self_col + 1]):
                    var other_val: ComplexScalar[Self.type] = other.v[k]
                    var other_col: Int = other.col_idx[k]
                    
                    # Accumulate or create result
                    if other_col in workspace:
                        workspace[other_col] += self_val * other_val
                    else:
                        workspace[other_col] = self_val * other_val

            # Since we can't do `for col, val in sorted(workspace.items()):` like in Python,
            # we have to do something like this nonsense
            var workspace_keys = List[Int, True]()
            for k in workspace.keys():
                workspace_keys.append(k[])
            sort(workspace_keys)
            for k in workspace_keys:
                if workspace[k[]] > Self.zero_threshold:
                    result_col_idx.append(k[])
                    result_v.append(workspace[k[]])
            
            # Update row index
            result_row_idx.append(result_col_idx.size)
        return Self(
            rows=self.rows,
            cols=other.cols,
            size=self.rows * other.cols,
            v=result_v^,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )
    
    fn _sparse_matmul_hash(self, other: Self) raises -> Self:
        '''Matrix multiplication of a sparse matrix with a sparse matrix using a 
        hash-based approach.
        '''
        result = Dict[Int, Dict[Int, ComplexScalar[Self.type]]]()

        # Iterate through the non-zero elements of self
        for i in range(self.rows):
            for j in range(self.row_idx[i], self.row_idx[i + 1]):
                var self_col: Int = self.col_idx[j]

                # Multiply with the corresponding row of other
                for k in range(other.row_idx[self_col], other.row_idx[self_col + 1]):
                    var other_col: Int = other.col_idx[k]
                    # I'm assuming this is faster than always checking if the dict entries exist
                    # Best to use a defaultdict if that gets added to the standard library
                    try:
                        result[i][other_col] += self.v[j] * other.v[k]
                    except:
                        if i not in result:
                            result[i] = Dict[Int, ComplexScalar[Self.type]]()
                        if other_col not in result[i]:
                            result[i][other_col] = ComplexScalar[Self.type]()
                        result[i][other_col] += self.v[j] * other.v[k]
        return Self(self.rows, other.cols, result)

    fn _sparse_matmul_tranposed(self, other: Self) raises -> Self:
        '''Perform matrix multiplication self * other^T.'''
        if self.cols != other.cols:
            raise Error(
                'Cannot transpose-multiply a sparse matrix with ' + str(self.cols) 
                + ' columns with a sparse matrix with ' + str(other.cols) + ' columns'
            )
        var result = Self(self.rows, other.rows)
        # Element (r, c) in result is a dot product of row r from self and row c from other
        # Loop over rows of self
        for r in range(self.rows):
            var self_start: Int = self.row_idx[r]
            var self_end: Int = self.row_idx[r + 1]
            # If row r is empty, don't bother looping over other
            if self_start == self_end:
                result.row_idx[r + 1] = result.v.size
                continue
            # Loop over rows of other
            for c in range(other.rows):
                var other_start: Int = other.row_idx[c]
                var other_end: Int = other.row_idx[c + 1]
                # If row c of other is empty, skip it
                if other_start == other_end:
                    continue
                var dot = ComplexScalar[Self.type]()
                # Pointers point to the col (and value) index in row r of self and row c of other
                var self_ptr: Int = self_start
                var other_ptr: Int = other_start
                while self_ptr < self_end and other_ptr < other_end:
                    var self_col: Int = self.col_idx[self_ptr]
                    var other_col: Int = other.col_idx[other_ptr]
                    # We only need to consider the case where the column indices are the same 
                    # since products vanish otherwise
                    if self_col == other_col:
                        dot += self.v[self_ptr] * other.v[other_ptr]
                        self_ptr += 1
                        other_ptr += 1
                    elif self_col < other_col:
                        self_ptr += 1
                    else:
                        other_ptr += 1
                # Don't store values close to zero
                if dot > Self.zero_threshold:
                    result.v.append(dot)
                    result.col_idx.append(c)
            result.row_idx[r + 1] = result.v.size
        return result

    @always_inline
    fn __matmul__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Matrix multiplication of a sparse matrix with a dense matrix.'''
        if self.cols != other.rows:
            raise Error(
                'Cannot multiply a sparse matrix with ' + str(self.cols)
                + ' columns with a matrix with ' + str(other.rows) + ' rows'
            )
        if other.cols == 1:
            return self._dense_vec_matmul(other)
        return self._dense_mat_matmul(other)
    
    @always_inline
    fn __matmul__(self, other: Self) raises -> Self:
        '''Matrix multiplication of a sparse matrix with a dense matrix.'''
        if self.cols != other.rows:
            raise Error(
                'Cannot multiply a sparse matrix with ' + str(self.cols)
                + ' columns with a sparse matrix with ' + str(other.rows) + ' rows'
            )
        # TODO: Benchmark to see which of these is better
        # return self._sparse_matmul_gustavson(other)
        return self._sparse_matmul_hash(other)

    # Other math ######################

    @always_inline
    fn __abs__(self) -> Self:
        '''Applies the absolute value to each element.'''
        var result: Self = self
        for i in range(self.v.size):
            result.v[i] = abs(result.v[i])
        return result

    
    @always_inline
    fn conj(self) -> Self:
        '''Applies the conjugate to each element.'''
        var result: Self = self
        for i in range(self.v.size):
            result.v[i] = result.v[i].conj()
        return result
    
    @always_inline
    fn iconj(inout self):
        '''Conjugate the matrix in-place.'''
        for i in range(self.v.size):
            self.v[i] = self.v[i].conj()
    
    fn dag(self) -> Self:
        '''Return the conjugate-transpose of the matrix.'''
        # Count non-zero elements in each column and cumsum
        var t_row_idx = List[Int, True](0) * (self.cols + 1)
        for c in self.col_idx:
            t_row_idx[c[] + 1] += 1
        for i in range(1, self.cols + 1):
            t_row_idx[i] += t_row_idx[i - 1]
        
        # Prepare data and indices for transposed matrix
        # This data will be overwritten
        var t_v: List[ComplexScalar[Self.type], True] = self.v
        var t_col_idx: List[Int, True] = self.col_idx
        var col_current_pos = List[Int, True](0) * self.cols
        
        for row in range(self.rows):
            for i in range(self.row_idx[row], self.row_idx[row + 1]):
                var col: Int = self.col_idx[i]
                var insert_pos: Int = t_row_idx[col] + col_current_pos[col]
                t_v[insert_pos] = self.v[i].conj()
                t_col_idx[insert_pos] = row
                col_current_pos[col] += 1
        return Self(
            rows=self.cols,
            cols=self.rows,
            size=self.size,
            v=t_v^,
            col_idx=t_col_idx^,
            row_idx=t_row_idx^,
        )
    
    @always_inline
    fn dagger(self) -> Self:
        '''Return the conjugate-transpose of the matrix. Alias of dag.'''
        return self.dag()

    @always_inline
    fn sum(self) -> ComplexScalar[Self.type]:
        var s = ComplexScalar[Self.type]()
        for v in self.v:
            s += v[]
        return s

    @always_inline
    fn frobenius_norm(self) raises -> Scalar[Self.type]:
        '''Return the Frobenius norm of self.'''
        var norm = self.__abs__()
        return sqrt((norm * norm).sum().re)

    # Shape operations ################

    fn transpose(self) -> Self:
        '''Return the transpose of self.'''
        # Count non-zero elements in each column and cumsum
        var t_row_idx = List[Int, True](0) * (self.cols + 1)
        for c in self.col_idx:
            t_row_idx[c[] + 1] += 1
        for i in range(1, self.cols + 1):
            t_row_idx[i] += t_row_idx[i - 1]
        
        # Prepare data and indices for transposed matrix
        # This data will be overwritten
        var t_v: List[ComplexScalar[Self.type], True] = self.v
        var t_col_idx: List[Int, True] = self.col_idx
        var col_current_pos = List[Int, True](0) * self.cols
        
        for row in range(self.rows):
            for i in range(self.row_idx[row], self.row_idx[row + 1]):
                var col: Int = self.col_idx[i]
                var insert_pos: Int = t_row_idx[col] + col_current_pos[col]
                t_v[insert_pos] = self.v[i]
                t_col_idx[insert_pos] = row
                col_current_pos[col] += 1
        return Self(
            rows=self.cols,
            cols=self.rows,
            size=self.size,
            v=t_v^,
            col_idx=t_col_idx^,
            row_idx=t_row_idx^,
        )

    fn reshape(self, new_rows: Int, new_cols: Int) raises -> Self:
        '''Return a reshaped matrix.'''
        self._assert_reshape_compatible(new_rows, new_cols)
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)
        var prev_row: Int = 0
        for r in range(self.rows):
            for i in range(self.row_idx[r], self.row_idx[r + 1]):
                # Compute element's index in row major format
                var idx = r * self.cols + self.col_idx[i]
                var new_row = idx // new_cols
                while prev_row < new_row:
                    result_row_idx.append(result_col_idx.size)
                    prev_row += 1
                result_col_idx.append(idx % new_cols)
        while result_row_idx.size < new_rows + 1:
            result_row_idx.append(result_col_idx.size)
        return Self(
            rows=new_rows, 
            cols=new_cols, 
            size=self.size, 
            v=self.v, 
            col_idx=result_col_idx^, 
            row_idx=result_row_idx^,
        )

    # Fill operations #################

    @always_inline
    fn fill_zero(inout self):
        '''Fill the sparse matrix with zeros in-place.'''
        self.v = List[ComplexScalar[Self.type], True]()
        self.col_idx = List[Int, True]()
        self.row_idx = List[Int, True](0) * (self.rows + 1)

    @always_inline
    fn zeros_like(self) -> Self:
        '''Return a matrix of zeros with the same shape as self.'''
        return Self(rows=self.rows, cols=self.cols)
    
    fn inset(
        self, 
        owned matrix: CMatrix[Self.type], 
        row: Int, 
        col: Int, 
        row_stride: Int = 1, 
        col_stride: Int = 1
    ) raises -> Self:
        '''Inset a matrix into a sparse matrix at position (row, col).'''
        var end_row: Int = row + (matrix.rows - 1) * row_stride 
        if row < 0 or end_row > self.rows - 1:
            raise Error('Cannot inset matrix: row dimension exceeds sparse matrix rows')
        var end_col: Int = col + (matrix.cols - 1) * col_stride 
        if col < 0 or end_col > self.cols - 1:
            raise Error('Cannot inset matrix: column dimension exceeds sparse matrix columns')

        if self.v.size == 0:
            return self._inset_into_zero(matrix^, row, col, row_stride, col_stride)
        if matrix.cols == 1:
            return self._inset_column(matrix^, row, col, row_stride)
        return self._inset(matrix^, row, end_row, col, end_col, row_stride, col_stride)
    
    fn _inset_into_zero(
        self, 
        owned matrix: CMatrix[Self.type], 
        row: Int, 
        col: Int, 
        row_stride: Int = 1, 
        col_stride: Int = 1
    ) -> Self:
        '''Inset a matrix into a sparse zero matrix at position (row, col).'''
        var result_v = List[ComplexScalar[Self.type], True]()
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0) * (row + 1)
        for r in range(matrix.rows):
            for c in range(matrix.cols):
                var val: ComplexScalar[Self.type] = matrix.load_crd[1](r, c)
                if val > Self.zero_threshold:
                    result_v.append(val)
                    result_col_idx.append(col + c * col_stride)
            if r and row_stride > 1:
                result_row_idx.extend(List[Int, True](result_row_idx[-1]) * (row_stride - 1))
            result_row_idx.append(result_v.size)
        var n_remaining_rows = self.rows - 1 - row - (matrix.rows - 1) * row_stride
        if n_remaining_rows > 0:
            result_row_idx.extend(List[Int, True](result_row_idx[-1]) * n_remaining_rows)
        return Self(
            rows=self.rows,
            cols=self.cols,
            size=self.size,
            v=result_v^,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )
    
    fn _inset_column(self,
        owned matrix: CMatrix[Self.type], 
        row: Int, 
        col: Int, 
        row_stride: Int = 1, 
    ) -> Self:
        '''Inset a column matrix into a sparse zero matrix at position (row, col).'''
        # If matrix is a column vector, we have to perform a search for the column index 
        # on each row anyway, so just do it the simple way with the setter
        var result: Self = self
        for r in range(matrix.rows):
            var new_val: ComplexScalar[Self.type] = matrix.load_idx[1](r)
            var result_row_idx: Int = row + r * row_stride
            if new_val == 0 and result.__contains__((result_row_idx, col)):
                result._setitem_noraise(result_row_idx, col, 0)
            elif new_val != 0:
                result._setitem_noraise(result_row_idx, col, new_val)
        return result
        
    fn _inset(
        self, 
        owned matrix: CMatrix[Self.type], 
        row: Int, 
        end_row: Int, 
        col: Int, 
        end_col: Int, 
        row_stride: Int = 1, 
        col_stride: Int = 1,
    ) -> Self:
        '''Inset a matrix into a sparse zero matrix at position (row, col).'''
        var result_v = List[ComplexScalar[Self.type], True]()
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)
        # Loop over rows of sparse matrix
        for r in range(self.rows):
            var self_start: Int = self.row_idx[r]
            var self_end: Int = self.row_idx[r + 1]
            if r < row or r > end_row or (r - row) % row_stride != 0:
                # Outside the inset
                result_v.extend(self.v[self_start:self_end])
                result_col_idx.extend(self.col_idx[self_start:self_end])
            else:
                var self_ptr: Int = self_start
                for c in range(self.cols):
                    var val_in_loc: Bool = self_ptr < self_end and c == self.col_idx[self_ptr]
                    var outside_inset: Bool = c < col or c > end_col or (c - col) % col_stride != 0
                    if val_in_loc and outside_inset:
                        result_v.append(self.v[self_ptr])
                        result_col_idx.append(c)
                    elif not outside_inset:
                        var val: ComplexScalar[Self.type] = matrix.load_crd[1](
                            (r - row) // row_stride, (c - col) // col_stride
                        )
                        if val > Self.zero_threshold:
                            result_v.append(val)
                            result_col_idx.append(c)
                    if val_in_loc:
                        self_ptr += 1
            result_row_idx.append(result_v.size)
        return Self(
            rows=self.rows,
            cols=self.cols,
            size=self.size,
            v=result_v^,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )

    fn fill_diag(inout self, val: ComplexScalar[Self.type], offset: Int = 0):
        '''Fill the diagonal at index offset with val in-place.'''
        var n_diag_elements: Int = (
            min(self.rows, self.cols - offset)
            if offset >= 0
            else min(self.rows + offset, self.cols)
        )
        var start_row: Int
        var start_col: Int
        if offset >= 0:
            start_row = 0
            start_col = offset
        else:
            start_row = -offset
            start_col = 0
        for i in range(n_diag_elements):
            # Setter performs a binary search each row
            self._setitem_noraise(start_row + i, start_col + i, val)

    @always_inline
    fn eye_like(self) -> Self:
        '''Return an identity matrix with the same shape as self.'''
        var result: Self = Self(rows=self.rows, cols=self.cols)
        result.fill_diag(ComplexScalar[Self.type](1))
        return result

    # Comparison operators ############

    # TODO: Make it better
    fn matrix_equals(self, other: Self) -> Bool:
        '''Returns True if self is the same shape as other and all elements 
        are equal, False otherwise.
        '''
        if self.rows != other.rows or self.cols != other.cols:
            return False
        if self.v.size != other.v.size:
            return False
        for i in range(self.v.size):
            if self.v[i] != other.v[i] or self.col_idx[i] != other.col_idx[i]:
                return False
        for i in range(self.row_idx.size):
            if self.row_idx[i] != other.row_idx[i]:
                return False
        return True
    
    fn is_close[tol: Scalar[Self.type] = DEFAULT_TOL](self, other: Self) -> Bool:
        '''Returns True if self is the same shape as other and corresponding elements 
        are within tol of each other, False otherwise.
        '''
        if self.rows != other.rows or self.cols != other.cols:
            return False
        for row in range(self.rows):
            var self_ptr: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            var other_ptr: Int = other.row_idx[row]
            var other_end: Int = other.row_idx[row + 1]

            while self_ptr < self_end and other_ptr < other_end:
                var self_col: Int = self.col_idx[self_ptr]
                var other_col: Int = other.col_idx[other_ptr]

                if self_col == other_col:
                    if not self.v[self_ptr].is_close[tol](other.v[other_ptr]):
                        return False
                    self_ptr += 1
                    other_ptr += 1
                elif self_col < other_col:
                    if not self.v[self_ptr].is_close[tol](0):
                        return False
                    self_ptr += 1
                else:
                    if not other.v[other_ptr].is_close[tol](0):
                        return False
                    other_ptr += 1

            while self_ptr < self_end:
                if not self.v[self_ptr].is_close[tol](0):
                    return False
                self_ptr += 1
            while other_ptr < other_end:
                if not other.v[other_ptr].is_close[tol](0):
                    return False
                other_ptr += 1
        return True

    fn is_close[tol: Scalar[Self.type] = DEFAULT_TOL](
        self, other: ComplexScalar[Self.type]
    ) -> Bool:
        '''Returns True all elements of self are within tol of other, False otherwise.'''
        if self.v.size < self.size and not other.is_close[tol](0):
            return False
        for v in self.v:
            if not other.is_close[tol](v[]):
                return False
        return True

    fn __eq__(self, other: Self) raises -> Self:
        '''Returns a CSRCMatrix with ones in the positions in which self 
        is equal to other and zeros elsewhere.
        '''
        self._assert_same_shape(other)
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)
        for row in range(self.rows):
            var self_ptr: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            var other_ptr: Int = other.row_idx[row]
            var other_end: Int = other.row_idx[row + 1]
            var current_col: Int = 0
            # Track whether we've processed all columns for this row
            while current_col < self.cols:
                # Get current column for each matrix (or max column if done)
                var self_col: Int = self.cols if self_ptr >= self_end else self.col_idx[self_ptr]
                var other_col: Int = self.cols 
                    if other_ptr >= other_end 
                    else other.col_idx[other_ptr]
                if current_col < self_col and current_col < other_col:
                    # Both matrices have zero at this column
                    result_col_idx.append(current_col)
                elif current_col == self_col and current_col == other_col:
                    # Both matrices have non-zero values at this column
                    if self.v[self_ptr] == other.v[other_ptr]:
                        result_col_idx.append(current_col)
                    self_ptr += 1
                    other_ptr += 1
                elif current_col == self_col:
                    # Only self has a non-zero value at this column
                    self_ptr += 1
                else:
                    # Only other has a non-zero value at this column
                    other_ptr += 1
                # Move to the next column
                current_col += 1
            result_row_idx.append(result_col_idx.size)
        return Self(
            rows=self.rows,
            cols=self.cols,
            size=self.size,
            v=List[ComplexScalar[Self.type], True](1) * result_col_idx.size,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )

    fn __eq__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Returns a CSRCMatrix with ones in the positions in which self 
        is equal to other and zeros elsewhere.
        '''
        var other_is_zero: Bool = other == 0
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)
        for row in range(self.rows):
            var self_ptr: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            for col in range(self.cols):
                var self_col = self.cols if self_ptr >= self_end else self.col_idx[self_ptr]
                if col < self_col:
                    # Self is zero at this column
                    if other_is_zero:
                        result_col_idx.append(col)
                else:
                    # Self is non-zero at this column
                    if other == self.v[self_ptr]:
                        result_col_idx.append(col)
                    self_ptr += 1
            result_row_idx.append(result_col_idx.size)
        return Self(
            rows=self.rows, 
            cols=self.cols, 
            size=self.size, 
            v=List[ComplexScalar[Self.type], True](1) * result_col_idx.size,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )

    fn __ne__(self, other: Self) raises -> Self:
        '''Returns a CSRCMatrix with ones in the positions in which self 
        is not equal to other and zeros elsewhere.
        '''
        self._assert_same_shape(other)
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)
        for row in range(self.rows):
            var self_ptr: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            var other_ptr: Int = other.row_idx[row]
            var other_end: Int = other.row_idx[row + 1]
            var current_col: Int = 0
            # Track whether we've processed all columns for this row
            while current_col < self.cols:
                # Get current column for each matrix (or max column if done)
                var self_col: Int = self.cols if self_ptr >= self_end else self.col_idx[self_ptr]
                var other_col: Int = self.cols 
                    if other_ptr >= other_end 
                    else other.col_idx[other_ptr]
                if current_col == self_col and current_col == other_col:
                    # Both matrices have non-zero values at this column
                    if self.v[self_ptr] != other.v[other_ptr]:
                        result_col_idx.append(current_col)
                    self_ptr += 1
                    other_ptr += 1
                elif current_col == self_col:
                    # Only self has a non-zero value at this column
                    result_col_idx.append(current_col)
                    self_ptr += 1
                elif current_col == other_col:
                    # Only other has a non-zero value at this column
                    result_col_idx.append(current_col)
                    other_ptr += 1
                # Move to the next column
                current_col += 1
            result_row_idx.append(result_col_idx.size)
        return Self(
            rows=self.rows, 
            cols=self.cols, 
            size=self.size, 
            v=List[ComplexScalar[Self.type], True](1) * result_col_idx.size,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )

    fn __ne__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Returns a CSRCMatrix with ones in the positions in which self 
        is not equal to other and zeros elsewhere.
        '''
        var other_is_non_zero: Bool = other != 0
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)
        for row in range(self.rows):
            var self_ptr: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            for col in range(self.cols):
                var self_col = self.cols if self_ptr >= self_end else self.col_idx[self_ptr]
                if col < self_col:
                    # Self is zero at this column
                    if other_is_non_zero:
                        result_col_idx.append(col)
                else:
                    # Self is non-zero at this column
                    if other != self.v[self_ptr]:
                        result_col_idx.append(col)
                    self_ptr += 1
            result_row_idx.append(result_col_idx.size)
        return Self(
            rows=self.rows,
            cols=self.cols, 
            size=self.size, 
            v=List[ComplexScalar[Self.type], True](1) * result_col_idx.size,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )

    fn __gt__(self, other: Self) raises -> Self:
        '''Returns a CSRCMatrix with ones in the positions in which self 
        is greater than other and zeros elsewhere.
        '''
        self._assert_same_shape(other)
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)
        for row in range(self.rows):
            var self_ptr: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            var other_ptr: Int = other.row_idx[row]
            var other_end: Int = other.row_idx[row + 1]
            var current_col: Int = 0
            # Track whether we've processed all columns for this row
            while current_col < self.cols:
                # Get current column for each matrix (or max column if done)
                var self_col: Int = self.cols if self_ptr >= self_end else self.col_idx[self_ptr]
                var other_col: Int = self.cols 
                    if other_ptr >= other_end 
                    else other.col_idx[other_ptr]
                if current_col == self_col and current_col == other_col:
                    # Both matrices have non-zero values at this column
                    if self.v[self_ptr] > other.v[other_ptr]:
                        result_col_idx.append(current_col)
                    self_ptr += 1
                    other_ptr += 1
                elif current_col == self_col:
                    # Only self has a non-zero value at this column
                    result_col_idx.append(current_col)
                    self_ptr += 1
                elif current_col == other_col:
                    # Only other has a non-zero value at this column
                    other_ptr += 1
                # Move to the next column
                current_col += 1
            result_row_idx.append(result_col_idx.size)
        return Self(
            rows=self.rows,
            cols=self.cols,
            size=self.size,
            v=List[ComplexScalar[Self.type], True](1) * result_col_idx.size,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )
    
    fn __gt__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Returns a CSRCMatrix with ones in the positions in which self 
        is greater than other and zeros elsewhere.
        '''
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)
        for row in range(self.rows):
            var self_ptr: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            for col in range(self.cols):
                var self_col = self.cols if self_ptr >= self_end else self.col_idx[self_ptr]
                if col == self_col:
                    # Self is non-zero at this column
                    if self.v[self_ptr] > other:
                        result_col_idx.append(col)
                    self_ptr += 1
            result_row_idx.append(result_col_idx.size)
        return Self(
            rows=self.rows,
            cols=self.cols, 
            size=self.size, 
            v=List[ComplexScalar[Self.type], True](1) * result_col_idx.size,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )
    
    fn __ge__(self, other: Self) raises -> Self:
        '''Returns a CSRCMatrix with ones in the positions in which self 
        is greater than or equal to other and zeros elsewhere.
        '''
        self._assert_same_shape(other)
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)
        for row in range(self.rows):
            var self_ptr: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            var other_ptr: Int = other.row_idx[row]
            var other_end: Int = other.row_idx[row + 1]
            var current_col: Int = 0
            # Track whether we've processed all columns for this row
            while current_col < self.cols:
                # Get current column for each matrix (or max column if done)
                var self_col: Int = self.cols if self_ptr >= self_end else self.col_idx[self_ptr]
                var other_col: Int = self.cols 
                    if other_ptr >= other_end 
                    else other.col_idx[other_ptr]
                if current_col < self_col and current_col < other_col:
                    # Both matrices have zero at this column
                    result_col_idx.append(current_col)
                elif current_col == self_col and current_col == other_col:
                    # Both matrices have non-zero values at this column
                    if self.v[self_ptr] >= other.v[other_ptr]:
                        result_col_idx.append(current_col)
                    self_ptr += 1
                    other_ptr += 1
                elif current_col == self_col:
                    # Only self has a non-zero value at this column
                    result_col_idx.append(current_col)
                    self_ptr += 1
                else:
                    # Only other has a non-zero value at this column
                    other_ptr += 1
                # Move to the next column
                current_col += 1
            result_row_idx.append(result_col_idx.size)
        return Self(
            rows=self.rows,
            cols=self.cols,
            size=self.size,
            v=List[ComplexScalar[Self.type], True](1) * result_col_idx.size,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )

    fn __ge__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Returns a CSRCMatrix with ones in the positions in which self 
        is greater than or equal to other and zeros elsewhere.
        '''
        var other_is_zero: Bool = other == 0
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)
        for row in range(self.rows):
            var self_ptr: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            for col in range(self.cols):
                var self_col = self.cols if self_ptr >= self_end else self.col_idx[self_ptr]
                if col < self_col:
                    # Self is zero at this column
                    if other_is_zero:
                        result_col_idx.append(col)
                else:
                    # Self is non-zero at this column
                    if self.v[self_ptr] >= other:
                        result_col_idx.append(col)
                    self_ptr += 1
            result_row_idx.append(result_col_idx.size)
        return Self(
            rows=self.rows, 
            cols=self.cols, 
            size=self.size, 
            v=List[ComplexScalar[Self.type], True](1) * result_col_idx.size,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )
    
    fn __lt__(self, other: Self) raises -> Self:
        '''Returns a CSRCMatrix with ones in the positions in which self 
        is less than other and zeros elsewhere.
        '''
        self._assert_same_shape(other)
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)
        for row in range(self.rows):
            var self_ptr: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            var other_ptr: Int = other.row_idx[row]
            var other_end: Int = other.row_idx[row + 1]
            var current_col: Int = 0
            # Track whether we've processed all columns for this row
            while current_col < self.cols:
                # Get current column for each matrix (or max column if done)
                var self_col: Int = self.cols if self_ptr >= self_end else self.col_idx[self_ptr]
                var other_col: Int = self.cols 
                    if other_ptr >= other_end 
                    else other.col_idx[other_ptr]
                if current_col == self_col and current_col == other_col:
                    # Both matrices have non-zero values at this column
                    if self.v[self_ptr] < other.v[other_ptr]:
                        result_col_idx.append(current_col)
                    self_ptr += 1
                    other_ptr += 1
                elif current_col == self_col:
                    # Only self has a non-zero value at this column
                    self_ptr += 1
                elif current_col == other_col:
                    # Only other has a non-zero value at this column
                    result_col_idx.append(current_col)
                    other_ptr += 1
                # Move to the next column
                current_col += 1
            result_row_idx.append(result_col_idx.size)
        return Self(
            rows=self.rows,
            cols=self.cols,
            size=self.size,
            v=List[ComplexScalar[Self.type], True](1) * result_col_idx.size,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )

    fn __lt__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Returns a CSRCMatrix with ones in the positions in which self 
        is less than other and zeros elsewhere.
        '''
        var other_is_non_zero: Bool = other != 0
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)
        for row in range(self.rows):
            var self_ptr: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            for col in range(self.cols):
                var self_col = self.cols if self_ptr >= self_end else self.col_idx[self_ptr]
                if col < self_col:
                    # Self is zero at this column
                    if other_is_non_zero:
                        result_col_idx.append(col)
                else:
                    # Self is non-zero at this column
                    if self.v[self_ptr] < other:
                        result_col_idx.append(col)
                    self_ptr += 1
            result_row_idx.append(result_col_idx.size)
        return Self(
            rows=self.rows, 
            cols=self.cols, 
            size=self.size, 
            v=List[ComplexScalar[Self.type], True](1) * result_col_idx.size,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )
    
    fn __le__(self, other: Self) raises -> Self:
        '''Returns a CSRCMatrix with ones in the positions in which self 
        is less than or equal to other and zeros elsewhere.
        '''
        self._assert_same_shape(other)
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)
        for row in range(self.rows):
            var self_ptr: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            var other_ptr: Int = other.row_idx[row]
            var other_end: Int = other.row_idx[row + 1]
            var current_col: Int = 0
            # Track whether we've processed all columns for this row
            while current_col < self.cols:
                # Get current column for each matrix (or max column if done)
                var self_col: Int = self.cols if self_ptr >= self_end else self.col_idx[self_ptr]
                var other_col: Int = self.cols 
                    if other_ptr >= other_end 
                    else other.col_idx[other_ptr]
                if current_col < self_col and current_col < other_col:
                    # Both matrices have zero at this column
                    result_col_idx.append(current_col)
                elif current_col == self_col and current_col == other_col:
                    # Both matrices have non-zero values at this column
                    if self.v[self_ptr] <= other.v[other_ptr]:
                        result_col_idx.append(current_col)
                    self_ptr += 1
                    other_ptr += 1
                elif current_col == self_col:
                    # Only self has a non-zero value at this column
                    self_ptr += 1
                else:
                    # Only other has a non-zero value at this column
                    result_col_idx.append(current_col)
                    other_ptr += 1
                # Move to the next column
                current_col += 1
            result_row_idx.append(result_col_idx.size)
        return Self(
            rows=self.rows,
            cols=self.cols,
            size=self.size,
            v=List[ComplexScalar[Self.type], True](1) * result_col_idx.size,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )

    fn __le__(self, other: ComplexScalar[Self.type]) raises -> Self:
        '''Returns a CSRCMatrix with ones in the positions in which self 
        is less than or equal to other and zeros elsewhere.
        '''
        var result_col_idx = List[Int, True]()
        var result_row_idx = List[Int, True](0)
        for row in range(self.rows):
            var self_ptr: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            for col in range(self.cols):
                var self_col = self.cols if self_ptr >= self_end else self.col_idx[self_ptr]
                if col < self_col:
                    # Self is zero at this column
                    result_col_idx.append(col)
                else:
                    # Self is non-zero at this column
                    if self.v[self_ptr] <= other:
                        result_col_idx.append(col)
                    self_ptr += 1
            result_row_idx.append(result_col_idx.size)
        return Self(
            rows=self.rows, 
            cols=self.cols, 
            size=self.size, 
            v=List[ComplexScalar[Self.type], True](1) * result_col_idx.size,
            col_idx=result_col_idx^,
            row_idx=result_row_idx^,
        )
