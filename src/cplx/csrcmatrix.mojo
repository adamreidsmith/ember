from algorithm import parallelize
from math import nan
from collections import Dict

from .complexsimd import ComplexScalar
from .cmatrix import CMatrix


@value
struct CSRCMatrix[type: DType, zero_threshold: Scalar[type] = 1e-15]:
# (
#     Absable,
#     Formattable, 
#     Sized, 
#     Representable,
#     StringableCollectionElement,
# ):
    # alias _tol_default: Scalar[Self.type] = 1e-12
    # Values less than or equal to Self.zero_threshold are considered zero

    var rows: Int
    var cols: Int
    var size: Int
    var v: List[ComplexScalar[Self.type], True]
    var col_idx: List[Int, True]
    var row_idx: List[Int, True]

    alias _csnan = ComplexScalar[type](nan[Self.type](), nan[Self.type]())

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

    fn __init__(inout self, rows: Int, cols: Int, *data: Tuple[Int, Int, ComplexScalar[Self.type]]) raises:
        '''Initialize a CSRCMatrix from tuples (row index, column index, value).'''
        for t in data:
            if t[][1] >= cols or t[][1] < 0:
                raise Error('Column index ' + str(t[][1]) + ' is out of range for a matrix with ' + str(cols) + ' columns')
            if t[][0] >= rows or t[][0] < 0:
                raise Error('Row index ' + str(t[][0]) + ' is out of range for a matrix with ' + str(rows) + ' rows')

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

        # Sort v and col_idx using indices and count the number of elements in each row
        var row_indices = List[Int, True]()
        for i in indices:
            self.v.append(data[i[]][2])
            self.col_idx.append(data[i[]][1])
            row_indices.append(data[i[]][0])

        # Remove any duplicates
        var i: Int = 1
        while i < len(self.v):
            if self.col_idx[i] == self.col_idx[i - 1] and row_indices[i] == row_indices[i - 1]:
                self.v[i - 1] += self.v.pop(i)
                _ = self.col_idx.pop(i)
                _ = row_indices.pop(i)
            else:
                i += 1
        
        # Count the number of elements in each row
        for idx in row_indices:
            self.row_idx[idx[] + 1] += 1

        # Compute the cumulative sum of self.row_idx
        for i in range(1, rows + 1):
            self.row_idx[i] += self.row_idx[i - 1]

    fn __init__(
        inout self, 
        rows: Int, 
        cols: Int, 
        row_idx: List[Int, True], 
        col_idx: List[Int, True], 
        val: List[ComplexScalar[Self.type], True]
    ) raises:
        '''Initialize a CSRCMatrix from lists of row indices, column indices, and values.'''
        if len(val) != len(col_idx) or len(val) != len(row_idx):
            raise Error(
                'Row indices, column indices, and values must be the same length; found lengths ' 
                + str(len(row_idx)) + ', '
                + str(len(col_idx)) + ', '
                + str(len(val))
            )

        for i in range(len(col_idx)):
            if col_idx[i] >= cols or col_idx[i] < 0:
                raise Error('Column index ' + str(col_idx[i]) + ' is out of range for a matrix with ' + str(cols) + ' columns')
            if row_idx[i] >= rows or row_idx[i] < 0:
                raise Error('Row index ' + str(row_idx[i]) + ' is out of range for a matrix with ' + str(rows) + ' rows')

        self.rows = rows
        self.cols = cols
        self.size = rows * cols
        self.v = List[ComplexScalar[Self.type], True]()
        self.col_idx = List[Int, True]()
        self.row_idx = List[Int, True](0) * (rows + 1)

        @parameter
        fn cmp_fn(i: Int, j: Int) -> Bool:
            # Sort first by row, then by column
            if row_idx[i] == row_idx[j]:
                return col_idx[i] < col_idx[j]
            return row_idx[i] < row_idx[j]

        # TODO: Gotta be a better way to do this
        indices = List[Int, True]()
        for i in range(len(val)):
            if val[i] > Self.zero_threshold:
                indices.append(i)
        sort[cmp_fn=cmp_fn](indices)

        # Sort v and col_idx using indices
        var row_indices = List[Int, True]()
        for i in indices:
            self.v.append(val[i[]])
            self.col_idx.append(col_idx[i[]])
            row_indices.append(row_idx[i[]])

        # Remove any duplicates
        var i: Int = 1
        while i < len(self.v):
            if self.col_idx[i] == self.col_idx[i - 1] and row_indices[i] == row_indices[i - 1]:
                self.v[i - 1] += self.v.pop(i)
                _ = self.col_idx.pop(i)
                _ = row_indices.pop(i)
            else:
                i += 1
        
        # Count the number of elements in each row
        for idx in row_indices:
            self.row_idx[idx[] + 1] += 1

        # Compute the cumulative sum of self.row_idx
        for i in range(1, rows + 1):
            self.row_idx[i] += self.row_idx[i - 1]
    
    fn __init__(
        inout self, rows: Int, cols: Int, sparse_dict: Dict[Int, Dict[Int, ComplexScalar[Self.type]]]
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
            if row_ptr < len(row_keys) and row == row_keys[row_ptr]:
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
            self.row_idx.append(len(self.col_idx))
    
    fn __init__(inout self, dense: CMatrix[Self.type]):
        '''Initialize a CSRCMatrix from a dense CMatrix.'''
        self.rows = dense.rows
        self.cols = dense.cols
        self.size = dense.size
        self.v = List[ComplexScalar[Self.type], True]()
        self.col_idx = List[Int, True]()
        self.row_idx = List[Int, True](0)

        for row in range(self.rows):
            for col in range(self.cols):
                var val: ComplexScalar[Self.type] = dense.load_crd[1](row, col)
                if val > Self.zero_threshold:
                    self.v.append(val)
                    self.col_idx.append(col)
            self.row_idx.append(len(self.v))

    
    # # These should be the same as what's created by @value
    # @always_inline
    # fn __copyinit__(inout self, existing: Self):
    #     self.rows = existing.rows
    #     self.cols = existing.cols
    #     self.size = existing.size
    #     self.v = existing.v
    #     self.col_idx = existing.col_idx
    #     self.row_idx = existing.row_idx
    
    # @always_inline
    # fn __moveinit__(inout self, owned existing: Self):
    #     self.rows = existing.rows
    #     self.cols = existing.cols
    #     self.size = existing.size
    #     self.v = existing.v^
    #     self.col_idx = existing.col_idx^
    #     self.row_idx = existing.row_idx^

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
    fn _assert_matmul_compatible(self, other: Self) raises:
        if self.cols != other.rows:
            raise Error(
                'Cannot multiply matrices with shapes (' 
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
    fn _assert_matmul_compatible(self, other: CMatrix[Self.type]) raises:
        if self.cols != other.rows:
            raise Error(
                'Cannot multiply matrices with shapes (' 
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
        var nonzero: Int = len(self.v)
        for v in self.v:
            if v[] <= Self.zero_threshold:
                nonzero -= 1
        return nonzero
    
    @always_inline
    fn n_stored(self) -> Int:
        '''Returns the number of elements stored in the matrix.'''
        return len(self.v)

    @always_inline
    fn is_square(self) -> Bool:
        '''Return True if self is square and non-degenerate (rows, cols > 0), and False otherwise.'''
        return self.rows == self.cols and self.rows > 0
    
    # @always_inline
    # fn is_unitary[tol: Scalar[Self.type] = Self._tol_default](self) raises -> Bool:
    #     '''Return True if self is unitary, False otherwise.'''
    #     if not self.is_square():
    #         return False
    #     return (self @ self.dag() - self.eye_like()).frobenius_norm() < tol

    # @always_inline
    # fn is_hermitian[tol: Scalar[Self.type] = Self._tol_default](self) raises -> Bool:
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

        # Search for the column in the row
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
                    for j in range(row + 1, len(self.row_idx)):
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
        for i in range(row + 1, len(self.row_idx)):
            self.row_idx[i] += 1

    fn _setitem_noraise(inout self, row: Int, col: Int, val: ComplexScalar[Self.type]):
        '''Set the value at the specified row and col without checking index validity.'''
        var start: Int = self.row_idx[row]
        var end: Int = self.row_idx[row + 1]
        var left: Int = start
        var right: Int = end - 1

        # Search for the column in the row
        while left <= right:
            var mid: Int = (left + right) // 2
            var curr_col: Int = self.col_idx[mid]
            if curr_col == col:
                if val <= Self.zero_threshold:
                    # Remove element if value is 0
                    _ = self.v.pop(mid)
                    _ = self.col_idx.pop(mid)
                    # Update row_idx accordingly
                    for j in range(row + 1, len(self.row_idx)):
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
        for i in range(row + 1, len(self.row_idx)):
            self.row_idx[i] += 1

    @always_inline
    fn __contains__(self, coords: Tuple[Int, Int]) raises -> Bool:
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
    
    @always_inline
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
            result.row_idx[row + 1] = len(result.v)
        return result
    
    @always_inline
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

    # Helper functions ################

    @always_inline
    def to_dense(self) -> CMatrix[Self.type]:
        '''Convert self to a dense CMatrix.'''
        var dense = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
        for row in range(self.rows):
            for k in range(self.row_idx[row], self.row_idx[row + 1]):
                dense[row, self.col_idx[k]] = self.v[k]
        return dense
    
    # Math dunders ####################

    @always_inline
    fn __neg__(self) -> Self:
        '''Defines the `-` unary negation operator. Returns -self.'''
        var result: Self = self
        for i in range(len(self.v)):
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
            result.row_idx.append(len(result.v))
        return result
    
    @always_inline
    fn __add__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the `+` add operator. Returns a CMatrix of self + other.'''
        self._assert_same_shape(other)
        var result: CMatrix[Self.type] = other
        for row in range(self.rows):
            var self_start: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            for i in range(self_start, self_end):
                result.store_crd[1](row, self.col_idx[i], result.load_crd[1](row, self.col_idx[i]) + self.v[i])
        return result

    # # Parallelization may add more overhead than value
    # @always_inline
    # fn __add__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
    #     '''Defines the `+` add operator. Returns a CMatrix of self + other.'''
    #     self._assert_same_shape(other)
    #     var result: CMatrix[Self.type] = other
    #     @parameter
    #     fn add_row(row: Int):
    #         var self_start: Int = self.row_idx[row]
    #         var self_end: Int = self.row_idx[row + 1]
    #         for i in range(self_start, self_end):
    #             result.store_crd[1](row, self.col_idx[i], result.load_crd[1](row, self.col_idx[i]) + self.v[i])
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
            result.row_idx.append(len(result.v))
        return result

    @always_inline
    fn __sub__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the `-` subtraction operator. Returns a CMatrix of self - other.'''
        self._assert_same_shape(other)
        var result: CMatrix[Self.type] = -other
        for row in range(self.rows):
            var self_start: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            for i in range(self_start, self_end):
                result.store_crd[1](row, self.col_idx[i], result.load_crd[1](row, self.col_idx[i]) + self.v[i])
        return result

    # # Parallelization may add more overhead than value
    # @always_inline
    # fn __sub__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
    #     '''Defines the `-` subtraction operator. Returns a CMatrix of self - other.'''
    #     self._assert_same_shape(other)
    #     var result: CMatrix[Self.type] = -other
    #     @parameter
    #     fn add_row(row: Int):
    #         var self_start: Int = self.row_idx[row]
    #         var self_end: Int = self.row_idx[row + 1]
    #         for i in range(self_start, self_end):
    #             result.store_crd[1](row, self.col_idx[i], result.load_crd[1](row, self.col_idx[i]) + self.v[i])
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
        '''Defines the right `-` subtraction operator. Returns -self if other == 0, raises otherwise.
        Subtracting a sparse matrix from a non-zero scalar is not supported.
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
            result.row_idx.append(len(result.v))
        return result
    
    @always_inline
    fn __mul__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the `*` multiplication operator. Returns a CMatrix of self * other.'''
        self._assert_same_shape(other)
        var result = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
        for row in range(self.rows):
            var self_start: Int = self.row_idx[row]
            var self_end: Int = self.row_idx[row + 1]
            for i in range(self_start, self_end):
                result.store_crd[1](row, self.col_idx[i], other.load_crd[1](row, self.col_idx[i]) * self.v[i])
        return result

    # # Parallelization may add more overhead than value
    # @always_inline
    # fn __mul__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
    #     '''Defines the `*` multiplication operator. Returns a CMatrix of self * other.'''
    #     self._assert_same_shape(other)
    #     var result = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
    #     @parameter
    #     fn mul_row(row: Int):
    #         var self_start: Int = self.row_idx[row]
    #         var self_end: Int = self.row_idx[row + 1]
    #         for i in range(self_start, self_end):
    #             result.store_crd[1](row, self.col_idx[i], other.load_crd[1](row, self.col_idx[i]) * self.v[i])
    #     parallelize[mul_row](self.rows, self.rows)
    #     return result

    @always_inline
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
            result.row_idx[row + 1] = len(result.v)
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

    @always_inline
    fn __truediv__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the `/` divide operator. Returns a CMatrix of self / other.'''
        self._assert_same_shape(other)
        var result = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
        for row in range(self.rows):
            var row_ptr: Int = self.row_idx[row]
            var row_end: Int = self.row_idx[row + 1]
            for c in range(self.cols):
                if row_ptr < row_end and c == self.col_idx[row_ptr]:
                    result.store_crd[1](row, c, self.v[row_ptr] / other.load_crd[1](row, c))
                    row_ptr += 1
                elif other.load_crd[1](row, c) == 0:
                    result.store_crd[1](row, c, Self._csnan)
        return result

    # # Parallelization may add more overhead than value
    # @always_inline
    # fn __truediv__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
    #     '''Defines the `/` divide operator. Returns a CMatrix of self / other.'''
    #     self._assert_same_shape(other)
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
    #                 result.store_crd[1](row, c, Self._csnan)
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

    @always_inline
    fn __floordiv__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `//` floor divide operator. Returns self // other.'''
        var result = Self(self.rows, self.cols)
        for row in range(self.rows):
            for i in range(self.row_idx[row], self.row_idx[row + 1]):
                var val: ComplexScalar[Self.type] = self.v[i] // other
                if val > Self.zero_threshold:
                    result.v.append(val)
                    result.col_idx.append(self.col_idx[i])
            result.row_idx[row + 1] = len(result.v)
        return result

    @always_inline
    fn __floordiv__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the `//` floor divide operator. Returns a CMatrix of self // other.'''
        self._assert_same_shape(other)
        var result = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
        for row in range(self.rows):
            var row_ptr: Int = self.row_idx[row]
            var row_end: Int = self.row_idx[row + 1]
            for c in range(self.cols):
                if row_ptr < row_end and c == self.col_idx[row_ptr]:
                    result.store_crd[1](row, c, self.v[row_ptr] // other.load_crd[1](row, c))
                    row_ptr += 1
                elif other.load_crd[1](row, c) == 0:
                    result.store_crd[1](row, c, Self._csnan)
        return result

    # # Parallelization may add more overhead than value
    # @always_inline
    # fn __floordiv__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
    #     '''Defines the `//` floor divide operator. Returns a CMatrix of self // other.'''
    #     self._assert_same_shape(other)
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
    #                 result.store_crd[1](row, c, Self._csnan)
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
    
    @always_inline
    fn __mod__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `%` mod operator. Returns self % other.'''
        var result = Self(self.rows, self.cols)
        for row in range(self.rows):
            for i in range(self.row_idx[row], self.row_idx[row + 1]):
                var val: ComplexScalar[Self.type] = self.v[i] % other
                if val > Self.zero_threshold:
                    result.v.append(val)
                    result.col_idx.append(self.col_idx[i])
            result.row_idx[row + 1] = len(result.v)
        return result

    @always_inline
    fn __mod__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
        '''Defines the `%` floor divide operator. Returns a CMatrix of self % other.'''
        self._assert_same_shape(other)
        var result = CMatrix[Self.type](self.rows, self.cols, fill_zeros=True)
        for row in range(self.rows):
            var row_ptr: Int = self.row_idx[row]
            var row_end: Int = self.row_idx[row + 1]
            for c in range(self.cols):
                if row_ptr < row_end and c == self.col_idx[row_ptr]:
                    result.store_crd[1](row, c, self.v[row_ptr] % other.load_crd[1](row, c))
                    row_ptr += 1
                elif other.load_crd[1](row, c) == 0:
                    result.store_crd[1](row, c, Self._csnan)
        return result

    # # Parallelization may add more overhead than value
    # @always_inline
    # fn __mod__(self, other: CMatrix[Self.type]) raises -> CMatrix[Self.type]:
    #     '''Defines the `%` floor divide operator. Returns a CMatrix of self % other.'''
    #     self._assert_same_shape(other)
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
    #                 result.store_crd[1](row, c, Self._csnan)
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
        '''Raises. The in-place modulo operation of a sparse matrix and a matrix is not supported.'''
        raise Error('The in-place modulo operation of a sparse matrix and a matrix is not supported')
    
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
    fn __divmod__(self, other: CMatrix[Self.type]) raises -> Tuple[CMatrix[Self.type], CMatrix[Self.type]]:
        '''Defines the divmod operator. Returns a tuple of CMatrix (self // other, self % other).'''
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
                result.store_idx[1](row, result.load_idx[1](row) + self.v[i] * vec.load_idx[1](self.col_idx[i]))
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
                        result.load_crd[1](row, col) + self.v[i] * mat.load_crd[1](self.col_idx[i], col),
                    )
        parallelize[dot](self.rows, self.rows)
        return result
    
    fn _sparse_matmul_gustavson(self, other: Self) raises -> Self:
        '''Matrix multiplication of a sparse matrix with a sparse matrix using Gustavson's algorithm.'''
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
            result_row_idx.append(len(result_col_idx))
        return Self(
            rows=self.rows,
            cols=other.cols,
            size=self.rows * other.cols,
            v=result_v,
            col_idx=result_col_idx,
            row_idx=result_row_idx,
        )
    
    fn _sparse_matmul_hash(self, other: Self) raises -> Self:
        '''Matrix multiplication of a sparse matrix with a sparse matrix using a hash-based approach.'''
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
                result.row_idx[r + 1] = len(result.v)
                continue
            # Loop over rows of other
            for c in range(other.rows):
                var other_start: Int = other.row_idx[c]
                var other_end: Int = other.row_idx[c + 1]
                # If row c of other is empty, skip it
                if other_start == other_end:
                    continue
                var dot = ComplexScalar[Self.type]()
                # Pointers point to the column (and value) index in row r of self and row c of other
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
            result.row_idx[r + 1] = len(result.v)
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

    # # Other math ######################

    # @always_inline
    # fn __abs__(self) -> Self:
    #     '''Returns a matrix with the absolute value applied to each element.'''
    #     if self._is_col_dominant:
    #         @parameter
    #         fn abs_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
    #             return self.load_crd[simd_width](r, c).__abs__()
    #         return self._parallelize_vectorize_op[abs_r]()
    #     else:
    #         @parameter
    #         fn abs_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
    #             return self.strided_load_idx[simd_width](r * self.cols + c, self.cols).__abs__()
    #         return self._parallelize_vectorize_op[abs_c]()
    
    # @always_inline
    # fn conj(self) -> Self:
    #     '''Return the conjugate of the matrix.'''
    #     var result: Self = self
    #     @parameter
    #     fn conj_simd[simd_width: Int](idx: Int):
    #         result.im.store[width=simd_width](idx, -self.im.load[width=simd_width](idx))
    #     vectorize[conj_simd, simdwidthof[Self.type]()](self.size)
    #     return result
    
    # @always_inline
    # fn iconj(self):
    #     '''Conjugate the matrix in-place.'''
    #     @parameter
    #     fn conj_simd[simd_width: Int](idx: Int):
    #         self.im.store[width=simd_width](idx, -self.im.load[width=simd_width](idx))
    #     vectorize[conj_simd, simdwidthof[Self.type]()](self.size)
    
    # @always_inline
    # fn dag(self) -> Self:
    #     '''Return the conjugate-transpose of the matrix.'''
    #     var result = Self(rows=self.cols, cols=self.rows, fill_zeros=False)
    #     @parameter
    #     fn transpose_row(r: Int):
    #         @parameter
    #         fn transpose_col[simd_width: Int](c: Int):
    #             result.strided_store_idx[simd_width](
    #                 c * result.cols + r, 
    #                 result.cols, 
    #                 self.load_crd[simd_width](r, c).conj(),
    #             )
    #         vectorize[transpose_col, simdwidthof[Self.type]()](self.cols)
    #     parallelize[transpose_row](self.rows, self.rows)
    #     return result
    
    # # alias dagger = Self.dag
    # @always_inline
    # fn dagger(self) -> Self:
    #     '''Return the conjugate-transpose of the matrix. Alias of dag.'''
    #     return self.dag()

    # # # TODO: Make it better
    # # # This can most definitely be done better with something like
    # # # buffer.Buffer or algorithm.reduce
    # @always_inline
    # fn sum(self) -> ComplexScalar[Self.type]:
    #     if self._is_col_dominant:
    #         var row_sums = Self(rows=self.rows, cols=1, fill_zeros=False)
    #         @parameter
    #         fn sum_rows(r: Int):
    #             var row_sum = ComplexScalar[Self.type]()
    #             for c in range(self.cols):
    #                 row_sum += self[r, c]
    #             row_sums[r] = row_sum
    #         parallelize[sum_rows](self.rows, self.rows)
    #         var total = ComplexScalar[Self.type]()
    #         for r in range(self.rows):
    #             total += row_sums[r]
    #         return total 
    #     else:
    #         var col_sums = Self(rows=1, cols=self.cols, fill_zeros=False)
    #         @parameter
    #         fn sum_cols(c: Int):
    #             var col_sum = ComplexScalar[Self.type]()
    #             for r in range(self.rows):
    #                 col_sum += self[r, c]
    #             col_sums[c] = col_sum
    #         parallelize[sum_cols](self.cols, self.cols)
    #         var total = ComplexScalar[Self.type]()
    #         for c in range(self.cols):
    #             total += col_sums[c]
    #         return total 
    
    # @always_inline
    # fn echelon[tol: Scalar[Self.type] = Self._tol_default](self) -> Self:
    #     '''Return the row echelon form of self.'''
    #     var A = self
    #     var h: Int = 0
    #     var k: Int = 0
    #     while h < self.rows and k < self.cols:
    #         var i_max: Scalar[Self.type] = A[h, k].norm()
    #         var i_argmax: Int = h
    #         for i in range(h + 1, self.rows):
    #             var i_norm: Scalar[Self.type] = A[i, k].norm()
    #             if i_norm > i_max:
    #                 i_max = i_norm
    #                 i_argmax = i
    #         if A[i_argmax, k].norm() < tol:
    #             k += 1
    #         else:
    #             for c in range(self.cols):
    #                 var p: ComplexScalar[Self.type] = A[h, c]
    #                 A[h, c] = A[i_argmax, c]
    #                 A[i_argmax, c] = p
    #             for i in range(h + 1, self.rows):
    #                 var f: ComplexScalar[type] = A[i, k] / A[h, k]
    #                 A[i, k] = ComplexScalar[type](0)
    #                 for j in range(k + 1, self.cols):
    #                     A[i, j] -= A[h, j] * f
    #             h += 1; k += 1
    #     return A

    # @always_inline
    # fn det[tol: Scalar[type] = Self._tol_default](self) raises -> ComplexScalar[Self.type]:
    #     '''Return the determinant of self.'''
    #     if self.rows != self.cols:
    #         raise Error('Only square matrices have determinants')
    #     var echelon: Self = self.echelon()
    #     var d = ComplexScalar[Self.type](1)
    #     for i in range(self.rows):
    #         d *= echelon[i, i]
    #     return d
    
    # @always_inline
    # fn determinant[tol: Scalar[type] = Self._tol_default](self) raises -> ComplexScalar[Self.type]:
    #     return self.det[tol]()
    
    # @always_inline
    # fn inv[tol: SIMD[Self.type, 1] = Self._tol_default](self) raises -> Self:
    #     '''Return the inverse of a square matrix.'''
    #     alias zero = ComplexScalar[Self.type](0)
    #     alias one = ComplexScalar[Self.type](1)

    #     if self.rows != self.cols:
    #         raise Error('Only square matrices are invertible')
        
    #     var augmented = Self(rows=self.rows, cols=2 * self.cols, fill_zeros=True)
    #     for r in range(self.rows):
    #         memcpy(augmented.re + r * augmented.cols, self.re + r * self.cols, self.cols)
    #         memcpy(augmented.im + r * augmented.cols, self.im + r * self.cols, self.cols)
    #         augmented[r, self.cols + r] = one
        
    #     for i in range(self.rows):
    #         if augmented[i, i].norm() < tol:
    #             for j in range(i + 1, self.rows):
    #                 if augmented[j, i].norm() >= tol:
    #                     var i_start: Int = i * augmented.cols
    #                     var j_start: Int = j * augmented.cols
    #                     var row_i = Self(rows=1, cols=augmented.cols, fill_zeros=False)
    #                     memcpy(row_i.re, augmented.re + i_start, augmented.cols)
    #                     memcpy(row_i.im, augmented.im + i_start, augmented.cols)
    #                     memcpy(augmented.re + i_start, augmented.re + j_start, augmented.cols)
    #                     memcpy(augmented.im + i_start, augmented.im + j_start, augmented.cols)
    #                     memcpy(augmented.re + j_start, row_i.re, augmented.cols)
    #                     memcpy(augmented.im + j_start, row_i.im, augmented.cols)
    #                     break
    #             else:
    #                 raise Error('Matrix is not invertible')

    #         var pivot: ComplexScalar[Self.type] = augmented[i, i]
    #         for j in range(augmented.cols):
    #             augmented[i, j] /= pivot

    #         for j in range(self.rows):
    #             if j != i:
    #                 var factor: ComplexScalar[Self.type] = augmented[j, i]
    #                 for k in range(augmented.cols):
    #                     augmented[j, k] -= factor * augmented[i, k]
        
    #     var result = Self(rows=self.rows, cols=self.cols, fill_zeros=False)
    #     for i in range(self.rows):
    #         memcpy(result.re + i * self.cols, augmented.re + i * augmented.cols + self.cols, self.cols)
    #         memcpy(result.im + i * self.cols, augmented.im + i * augmented.cols + self.cols, self.cols)
    #     return result

    # @always_inline
    # fn inverse(self) raises -> Self:
    #     '''Return the inverse of a square matrix. Alias of inv.'''
    #     return self.inv()

    # @always_inline
    # fn frobenius_norm(self) raises -> Scalar[Self.type]:
    #     '''Return the Frobenius norm of self.'''
    #     var norm = self.__abs__()
    #     return sqrt((norm * norm).sum().re)

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
            var row_start: Int = self.row_idx[row]
            var row_end: Int = self.row_idx[row + 1]
            for i in range(row_start, row_end):
                var col: Int = self.col_idx[i]
                var val: ComplexScalar[Self.type] = self.v[i]
                var insert_pos: Int = t_row_idx[col] + col_current_pos[col]
                t_v[insert_pos] = val
                t_col_idx[insert_pos] = row
                col_current_pos[col] += 1
        
        return Self(
            rows=self.cols,
            cols=self.rows,
            size=self.size,
            v=t_v,
            col_idx=t_col_idx,
            row_idx=t_row_idx
        )


    # @always_inline
    # fn reshape(self, new_rows: Int, new_cols: Int) raises -> Self:
    #     '''Return a reshaped matrix.'''
    #     self._assert_reshape_compatible(new_rows, new_cols)
    #     var result= Self(rows=new_rows, cols=new_cols, fill_zeros=False)
    #     memcpy(result.re, self.re, self.size)
    #     memcpy(result.im, self.im, self.size)
    #     return result
    
    # @always_inline
    # fn ireshape(inout self, new_rows: Int, new_cols: Int) raises:
    #     '''Reshape self, in-place.'''
    #     self._assert_reshape_compatible(new_rows, new_cols)
    #     self.rows = new_rows
    #     self.cols = new_cols
    #     self._is_col_dominant = self.cols >= self.rows
    
    # @always_inline
    # fn flatten_to_row(self) -> Self:
    #     '''Return a flattened row matrix.'''
    #     var result = Self(rows=1, cols=self.size, fill_zeros=False)
    #     memcpy(result.re, self.re, self.size)
    #     memcpy(result.im, self.im, self.size)
    #     return result

    # @always_inline
    # fn iflatten_to_row(inout self):
    #     '''Flatten self to a row, in-place.'''
    #     self.rows = 1
    #     self.cols = self.size
    #     self._is_col_dominant = True
    
    # @always_inline
    # fn flatten_to_column(self) -> Self:
    #     '''Return a flattened column matrix.'''
    #     var result = Self(rows=self.size, cols=1, fill_zeros=False)
    #     memcpy(result.re, self.re, self.size)
    #     memcpy(result.im, self.im, self.size)
    #     return result

    # @always_inline
    # fn iflatten_to_column(inout self):
    #     '''Flatten self to a row, in-place.'''
    #     self.rows = self.size
    #     self.cols = 1
    #     self._is_col_dominant = False
    
    # @always_inline
    # fn transpose(self) -> Self:
    #     '''Return the transpose of the matrix.'''
    #     var result = Self(rows=self.cols, cols=self.rows, fill_zeros=False)
    #     @parameter
    #     fn transpose_row(r: Int):
    #         for c in range(self.cols):
    #             result[c, r] = self[r, c]
    #     parallelize[transpose_row](self.rows, self.rows)
    #     return result

    # # TODO: Make an efficient in-place transpose

    # # Fill operations #################

    # @always_inline
    # fn fill_zero(self):
    #     '''Fill a matrix with zeros in-place.'''
    #     memset_zero(self.re.address, self.size)
    #     memset_zero(self.im.address, self.size)
    
    # @always_inline
    # fn fill(self, val: ComplexScalar[Self.type]):
    #     '''Fill a matrix with val in-place.'''
    #     @parameter
    #     fn fill_row(r: Int):
    #         @parameter
    #         fn fill_col[simd_width: Int](c: Int):
    #             self.store_crd[simd_width](r, c, ComplexSIMD[Self.type, simd_width](val.re, val.im))
    #         vectorize[fill_col, simdwidthof[Self.type]()](self.cols)
    #     parallelize[fill_row](self.rows, self.rows)
    
    # @always_inline
    # fn fill_one(self):
    #     '''Fill a matrix with ones in-place.'''
    #     self.fill(ComplexScalar[Self.type](1, 0))
    
    # @always_inline
    # fn fill_i(self):
    #     '''Fill a matrix with ones in-place.'''
    #     self.fill(ComplexScalar[Self.type](0, 1))

    # @always_inline
    # fn fill_diag(self, val: ComplexScalar[Self.type], offset: Int = 0):
    #     '''Fill the diagonal at index offset of a matrix with val in-place.'''
    #     var stride = self.cols + 1
    #     var n_diag_elements: Int = (
    #         min(self.rows, self.cols - offset) 
    #         if offset >= 0 
    #         else min(self.rows + offset, self.cols)
    #     )
    #     var idx: Int
    #     @parameter
    #     fn fill_diag_simd_rp[simd_width: Int](p: Int):
    #         idx = p * stride + offset if offset >= 0 else p * stride - offset * self.cols
    #         self.strided_store_idx[simd_width](
    #             idx, stride, ComplexSIMD[Self.type, simd_width](val.re, val.im)
    #         )
    #     vectorize[fill_diag_simd_rp, simdwidthof[Self.type]()](n_diag_elements)

    # @always_inline
    # fn fill_range(self):
    #     memset_zero(self.im.address, self.size)
    #     for idx in range(self.size):
    #         self.re.store[width=1](idx, SIMD[type, 1](idx))
    
    # @always_inline
    # fn fill_range(self, start: Int):
    #     memset_zero(self.im.address, self.size)
    #     for idx in range(self.size):
    #         self.re.store[width=1](idx, SIMD[type, 1](idx + start))
    
    # @always_inline
    # fn fill_range(self, start: Int, step: Int):
    #     memset_zero(self.im.address, self.size)
    #     for idx in range(self.size):
    #         self.re.store[width=1](idx, SIMD[type, 1](step * idx + start))

    # @always_inline
    # fn zeros_like(self) -> Self:
    #     '''Return a matrix of zeros with the same shape as self.'''
    #     return Self(rows=self.rows, cols=self.cols, fill_zeros=True)
    
    # @always_inline
    # fn ones_like(self) -> Self:
    #     '''Return a matrix of ones with the same shape as self.'''
    #     var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=False)
    #     result.fill_one()
    #     return result
    
    # @always_inline
    # fn i_like(self) -> Self:
    #     '''Return a matrix of i's with the same shape as self.'''
    #     var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=False)
    #     result.fill_i()
    #     return result
    
    # @always_inline
    # fn eye_like(self) -> Self:
    #     '''Return an identity matrix with the same shape as self.'''
    #     var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=True)
    #     result.fill_diag(ComplexScalar[Self.type](1, 0))
    #     return result

    # # Comparison operators ############

    # # TODO: Make it better
    # @always_inline
    # fn matrix_equals(self, other: Self) -> Bool:
    #     '''Returns True if self is the same shape as other and all elements 
    #     are equal, False otherwise.
    #     '''
    #     if self.rows != other.rows or self.cols != other.cols:
    #         return False
    #     for r in range(self.rows):
    #         for c in range(self.cols):
    #             if self[r, c] != other[r, c]:
    #                 return False
    #     return True
    
    # # TODO: Make it better
    # @always_inline
    # fn is_close[tol: Scalar[Self.type] = Self._tol_default](self, other: Self) -> Bool:
    #     '''Returns True if self is the same shape as other and corresponding elements 
    #     are within tol of each other, False otherwise.
    #     '''
    #     if self.rows != other.rows or self.cols != other.cols:
    #         return False
    #     for r in range(self.rows):
    #         for c in range(self.cols):
    #             if not self[r, c].is_close(other[r, c], tol):
    #                 return False
    #     return True

    # # TODO: Make it better
    # @always_inline
    # fn is_close[tol: Scalar[Self.type] = Self._tol_default](self, other: ComplexScalar[Self.type]) -> Bool:
    #     '''Returns True all elements of self are within tol of other, False otherwise.'''
    #     for r in range(self.rows):
    #         for c in range(self.cols):
    #             if not self[r, c].is_close(other, tol):
    #                 return False
    #     return True

    # @always_inline
    # fn __eq__(self, other: Self) raises -> Self:
    #     '''Returns a CMatrix with ones in the positions in which self 
    #     is equal to other and zeros elsewhere. The returned matrix has
    #     the order of self.
    #     '''
    #     self._assert_same_shape(other)
    #     alias one = SIMD[Self.type, 1](1)
    #     var result: Self = Self(self.rows, self.cols, fill_zeros=True)
    #     if self._is_col_dominant:
    #         @parameter
    #         fn row_eq(r: Int):
    #             for c in range(self.cols):
    #                 if self[r, c] == other[r, c]:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[row_eq](self.rows, self.rows)
    #     else:
    #         @parameter
    #         fn col_eq(c: Int):
    #             for r in range(self.rows):
    #                 if self[r, c] == other[r, c]:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[col_eq](self.cols, self.cols)
    #     return result
    
    # @always_inline
    # fn __eq__(self, other: ComplexScalar[Self.type]) -> Self:
    #     '''Returns a CMatrix with ones in the positions in which self 
    #     is equal to other and zeros elsewhere.
    #     '''
    #     alias one = SIMD[Self.type, 1](1)
    #     var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=True)
    #     if self._is_col_dominant:
    #         @parameter
    #         fn row_eq(r: Int):
    #             for c in range(self.cols):
    #                 if self[r, c] == other:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[row_eq](self.rows, self.rows)
    #     else:
    #         @parameter
    #         fn col_eq(c: Int):
    #             for r in range(self.rows):
    #                 if self[r, c] == other:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[col_eq](self.cols, self.cols)
    #     return result
    
    # @always_inline
    # fn __ne__(self, other: Self) raises -> Self:
    #     '''Returns a CMatrix with ones in the positions in which self 
    #     is not equal to other and zeros elsewhere. The returned matrix has
    #     the order of self.
    #     '''
    #     self._assert_same_shape(other)
    #     alias one = SIMD[Self.type, 1](1)
    #     var result: Self = Self(self.rows, self.cols, fill_zeros=True)
    #     if self._is_col_dominant:
    #         @parameter
    #         fn row_eq(r: Int):
    #             for c in range(self.cols):
    #                 if self[r, c] != other[r, c]:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[row_eq](self.rows, self.rows)
    #     else:
    #         @parameter
    #         fn col_eq(c: Int):
    #             for r in range(self.rows):
    #                 if self[r, c] != other[r, c]:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[col_eq](self.cols, self.cols)
    #     return result
    
    # @always_inline
    # fn __ne__(self, other: ComplexScalar[Self.type]) -> Self:
    #     '''Returns a CMatrix with ones in the positions in which self 
    #     is not equal to other and zeros elsewhere.
    #     '''
    #     alias one = SIMD[Self.type, 1](1)
    #     var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=True)
    #     if self._is_col_dominant:
    #         @parameter
    #         fn row_eq(r: Int):
    #             for c in range(self.cols):
    #                 if self[r, c] != other:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[row_eq](self.rows, self.rows)
    #     else:
    #         @parameter
    #         fn col_eq(c: Int):
    #             for r in range(self.rows):
    #                 if self[r, c] != other:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[col_eq](self.cols, self.cols)
    #     return result
    
    # @always_inline
    # fn __gt__(self, other: Self) raises -> Self:
    #     '''Returns a CMatrix with ones in the positions in which self 
    #     is greater than other and zeros elsewhere. The returned matrix has
    #     the order of self.
    #     '''
    #     self._assert_same_shape(other)
    #     alias one = SIMD[Self.type, 1](1)
    #     var result: Self = Self(self.rows, self.cols, fill_zeros=True)
    #     if self._is_col_dominant:
    #         @parameter
    #         fn row_eq(r: Int):
    #             for c in range(self.cols):
    #                 if self[r, c] > other[r, c]:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[row_eq](self.rows, self.rows)
    #     else:
    #         @parameter
    #         fn col_eq(c: Int):
    #             for r in range(self.rows):
    #                 if self[r, c] > other[r, c]:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[col_eq](self.cols, self.cols)
    #     return result
    
    # @always_inline
    # fn __gt__(self, other: ComplexScalar[Self.type]) -> Self:
    #     '''Returns a CMatrix with ones in the positions in which self 
    #     is greater than other and zeros elsewhere.
    #     '''
    #     alias one = SIMD[Self.type, 1](1)
    #     var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=True)
    #     if self._is_col_dominant:
    #         @parameter
    #         fn row_eq(r: Int):
    #             for c in range(self.cols):
    #                 if self[r, c] > other:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[row_eq](self.rows, self.rows)
    #     else:
    #         @parameter
    #         fn col_eq(c: Int):
    #             for r in range(self.rows):
    #                 if self[r, c] > other:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[col_eq](self.cols, self.cols)
    #     return result
    
    # @always_inline
    # fn __ge__(self, other: Self) raises -> Self:
    #     '''Returns a CMatrix with ones in the positions in which self 
    #     is greater than or equal to other and zeros elsewhere. The returned matrix has
    #     the order of self.
    #     '''
    #     self._assert_same_shape(other)
    #     alias one = SIMD[Self.type, 1](1)
    #     var result: Self = Self(self.rows, self.cols, fill_zeros=True)
    #     if self._is_col_dominant:
    #         @parameter
    #         fn row_eq(r: Int):
    #             for c in range(self.cols):
    #                 if self[r, c] >= other[r, c]:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[row_eq](self.rows, self.rows)
    #     else:
    #         @parameter
    #         fn col_eq(c: Int):
    #             for r in range(self.rows):
    #                 if self[r, c] >= other[r, c]:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[col_eq](self.cols, self.cols)
    #     return result
    
    # @always_inline
    # fn __ge__(self, other: ComplexScalar[Self.type]) -> Self:
    #     '''Returns a CMatrix with ones in the positions in which self 
    #     is greater than or equal to other and zeros elsewhere.
    #     '''
    #     alias one = SIMD[Self.type, 1](1)
    #     var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=True)
    #     if self._is_col_dominant:
    #         @parameter
    #         fn row_eq(r: Int):
    #             for c in range(self.cols):
    #                 if self[r, c] >= other:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[row_eq](self.rows, self.rows)
    #     else:
    #         @parameter
    #         fn col_eq(c: Int):
    #             for r in range(self.rows):
    #                 if self[r, c] >= other:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[col_eq](self.cols, self.cols)
    #     return result
    
    # @always_inline
    # fn __lt__(self, other: Self) raises -> Self:
    #     '''Returns a CMatrix with ones in the positions in which self 
    #     is less than other and zeros elsewhere. The returned matrix has
    #     the order of self.
    #     '''
    #     self._assert_same_shape(other)
    #     alias one = SIMD[Self.type, 1](1)
    #     var result: Self = Self(self.rows, self.cols, fill_zeros=True)
    #     if self._is_col_dominant:
    #         @parameter
    #         fn row_eq(r: Int):
    #             for c in range(self.cols):
    #                 if self[r, c] < other[r, c]:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[row_eq](self.rows, self.rows)
    #     else:
    #         @parameter
    #         fn col_eq(c: Int):
    #             for r in range(self.rows):
    #                 if self[r, c] < other[r, c]:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[col_eq](self.cols, self.cols)
    #     return result
    
    # @always_inline
    # fn __lt__(self, other: ComplexScalar[Self.type]) -> Self:
    #     '''Returns a CMatrix with ones in the positions in which self 
    #     is less than other and zeros elsewhere.
    #     '''
    #     alias one = SIMD[Self.type, 1](1)
    #     var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=True)
    #     if self._is_col_dominant:
    #         @parameter
    #         fn row_eq(r: Int):
    #             for c in range(self.cols):
    #                 if self[r, c] < other:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[row_eq](self.rows, self.rows)
    #     else:
    #         @parameter
    #         fn col_eq(c: Int):
    #             for r in range(self.rows):
    #                 if self[r, c] < other:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[col_eq](self.cols, self.cols)
    #     return result
    
    # @always_inline
    # fn __le__(self, other: Self) raises -> Self:
    #     '''Returns a CMatrix with ones in the positions in which self 
    #     is less than or equal to other and zeros elsewhere. The returned matrix has
    #     the order of self.
    #     '''
    #     self._assert_same_shape(other)
    #     alias one = SIMD[Self.type, 1](1)
    #     var result: Self = Self(self.rows, self.cols, fill_zeros=True)
    #     if self._is_col_dominant:
    #         @parameter
    #         fn row_eq(r: Int):
    #             for c in range(self.cols):
    #                 if self[r, c] <= other[r, c]:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[row_eq](self.rows, self.rows)
    #     else:
    #         @parameter
    #         fn col_eq(c: Int):
    #             for r in range(self.rows):
    #                 if self[r, c] <= other[r, c]:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[col_eq](self.cols, self.cols)
    #     return result
    
    # @always_inline
    # fn __le__(self, other: ComplexScalar[Self.type]) -> Self:
    #     '''Returns a CMatrix with ones in the positions in which self 
    #     is less than or equal to other and zeros elsewhere.
    #     '''
    #     alias one = SIMD[Self.type, 1](1)
    #     var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=True)
    #     if self._is_col_dominant:
    #         @parameter
    #         fn row_eq(r: Int):
    #             for c in range(self.cols):
    #                 if self[r, c] <= other:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[row_eq](self.rows, self.rows)
    #     else:
    #         @parameter
    #         fn col_eq(c: Int):
    #             for r in range(self.rows):
    #                 if self[r, c] <= other:
    #                     result.re.store[width=1](r * self.cols + c, one)
    #         parallelize[col_eq](self.cols, self.cols)
    #     return result