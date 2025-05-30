from random import randn_float64
from math import sqrt
from memory import memset_zero, memcpy, UnsafePointer
from algorithm import parallelize, vectorize
from os import Atomic
from sys import simdwidthof

from .complexsimd import ComplexScalar, ComplexSIMD
from ..config import DEFAULT_TYPE, DEFAULT_TOL


@value
struct CMatrix[type: DType = DEFAULT_TYPE](
    Absable,
    Sized,
    Representable,
    CollectionElement,
    Stringable,
    Writable,
):
    '''Complex matrices with fast, efficient operations.
    
    Parameters:
        type: A type for the matrix data.
    '''

    var re: UnsafePointer[Scalar[Self.type]]
    '''The real components of the matrix elements.'''
    var im: UnsafePointer[Scalar[Self.type]]
    '''The imaginary components of the matrix elements.'''
    var rows: Int
    '''The number of rows in the matrix.'''
    var cols: Int
    '''The number of columns in the matrix.'''
    var size: Int
    '''The total number of elements in the matrix.'''

    # Initialization ##################

    fn __init__[__: None = None](out self, rows: Int, cols: Int, fill_zeros: Bool = True):
        '''Initialize a zero matrix with a defined size.
        
        Args:
            rows: The number of rows in the matrix.
            cols: The number of columns in the matrix.
            fill_zeros: Whether or not to fill the matrix with zeros. If False, memory is allocated 
                but remains uninitialized.
        '''
        self.rows = rows
        self.cols = cols
        self.size = self.rows * self.cols
        self.re = UnsafePointer[Scalar[Self.type]].alloc(self.size)
        self.im = UnsafePointer[Scalar[Self.type]].alloc(self.size)
        if fill_zeros:
            memset_zero(self.re.address, self.size)
            memset_zero(self.im.address, self.size)

    fn __init__(out self, rows: Int, cols: Int, *data: ComplexScalar[Self.type]) raises:
        '''Initialize a matrix with a defined size and data.
        
        Args:
            rows: The number of rows in the matrix.
            cols: The number of columns in the matrix.
            data: The data to fill the matrix. The matrix is filled in row major order, with any
                remaining unfilled elements set to zero.
        '''
        if len(data) > rows * cols:
            raise Error('More elements provided than fit in the matrix')
        self.rows = rows
        self.cols = cols
        self.size = self.rows * self.cols
        self.re = UnsafePointer[Scalar[Self.type]].alloc(self.size)
        self.im = UnsafePointer[Scalar[Self.type]].alloc(self.size)
        for idx in range(len(data)):
            self.store_idx[1](idx, data[idx])
        for idx in range(len(data), self.size):
            self.store_idx[1](idx, 0)
    
    fn __init__(out self, rows: Int, cols: Int, data: List[ComplexScalar[Self.type], True]) raises:
        '''Initialize a matrix with a defined size and data.
        
        Args:
            rows: The number of rows in the matrix.
            cols: The number of columns in the matrix.
            data: The data to fill the matrix. The matrix is filled in row major order, with any
                remaining unfilled elements set to zero.
        '''
        if len(data) > rows * cols:
            raise Error('More elements provided than fit in the matrix')
        self.rows = rows
        self.cols = cols
        self.size = self.rows * self.cols
        self.re = UnsafePointer[Scalar[Self.type]].alloc(self.size)
        self.im = UnsafePointer[Scalar[Self.type]].alloc(self.size)
        for idx in range(len(data)):
            self.store_idx[1](idx, data[idx])
        for idx in range(len(data), self.size):
            self.store_idx[1](idx, 0)

    fn __init__(out self, data: List[ComplexScalar[Self.type], True]):
        '''Initialize a row vector with specified data.

        Args:
            data: The data to fill the row vector.
        '''
        self.rows = 1
        self.cols = len(data)
        self.size = self.rows * self.cols
        self.re = UnsafePointer[Scalar[Self.type]].alloc(self.size)
        self.im = UnsafePointer[Scalar[Self.type]].alloc(self.size)
        for idx in range(len(data)):
            self.store_idx[1](idx, data[idx])

    fn __init__(out self, data: List[List[ComplexScalar[Self.type], True]]) raises:
        '''Initialize a matrix from a list of lists. Each sublist must have the same length.

        Args:
            data: A list of lists defining the matrix values.
        '''
        self.rows = len(data)
        self.cols = len(data[0])
        for row in data[1:]:
            if len(row[]) != self.cols:
                raise Error('All sub-list of `data` must be the same length')
        self.size = self.rows * self.cols
        self.re = UnsafePointer[Scalar[Self.type]].alloc(self.size)
        self.im = UnsafePointer[Scalar[Self.type]].alloc(self.size)
        for r in range(self.rows):
            for c in range(self.cols):
                self.store_crd[1](r, c, data[r][c])

    fn __copyinit__(out self, existing: Self):
        '''Initialize a matrix by copying another matrix.

        Args:
            existing: The matrix to copy into self.
        '''
        self.rows = existing.rows
        self.cols = existing.cols
        self.size = existing.size
        self.re = UnsafePointer[Scalar[Self.type]].alloc(self.size)
        self.im = UnsafePointer[Scalar[Self.type]].alloc(self.size)
        # memcpy(self.re, existing.re, self.size)
        # memcpy(self.im, existing.im, self.size)
        @parameter
        fn cpy_data[simd_width: Int](idx: Int):
            self.store_idx[simd_width](idx, existing.load_idx[simd_width](idx))
        vectorize[cpy_data, simdwidthof[Self.type]()](self.size)

    fn __moveinit__(out self, owned existing: Self):
        '''Initialize a matrix by moving data from another matrix.

        Args:
            existing: The matrix to move into self.
        '''
        self.rows = existing.rows
        self.cols = existing.cols
        self.size = existing.size
        self.re = existing.re
        self.im = existing.im

    fn __del__(owned self):
        '''Delete the matrix and free memory.'''
        # # This loop may not be necessary since scalar types are trivial
        # for idx in range(self.size):
        #     (self.re + idx).destroy_pointee()
        #     (self.im + idx).destroy_pointee()
        if self.re:
            self.re.free()
        if self.im:
            self.im.free()
    
    @always_inline
    fn _should_parallelize(self) -> Bool:
        '''Determine whether or not we should parallelize operations in a row-major format.

        Returns:
            Whether or not we should parallelize.
        '''
        # Never parallelize small matrices or matrices with few rows
        if self.size < 1000:
            return False
        if self.rows <= 2:
            return False
        
        # For medium matrices, only parallelize if there are many rows
        if self.size < 40000:
            return self.rows >= 8
        
        # For big matrices parallelize in most cases
        return self.rows >= 4
    
    # Static constructors #############

    @staticmethod
    @always_inline
    fn eye(rows: Int, cols: Int) -> Self:
        '''Initialize an identity matrix with the specified dimensions.
        
        Args:
            rows: The number of rows in the matrix.
            cols: The number of columns in the matrix.
        
        Returns:
            A matrix with ones on the main diagonal and zeros elsewhere.
        '''
        var result = Self(rows, cols, fill_zeros=True)
        result.fill_diag(1)
        return result

    @staticmethod
    @always_inline
    fn zeros(rows: Int, cols: Int) -> Self:
        '''Initialize an matrix of zeros with the specified dimensions.
        
        Args:
            rows: The number of rows in the matrix.
            cols: The number of columns in the matrix.
        
        Returns:
            A matrix fileld with zeros.
        '''
        return Self(rows, cols, fill_zeros=True)

    @staticmethod
    @always_inline
    fn ones(rows: Int, cols: Int) -> Self:
        '''Initialize an matrix of ones with the specified dimensions.
        
        Args:
            rows: The number of rows in the matrix.
            cols: The number of columns in the matrix.
        
        Returns:
            A matrix filled with ones.
        '''
        var result = Self(rows, cols, fill_zeros=False)
        result.fill_one()
        return result

    @staticmethod
    @always_inline
    fn i(rows: Int, cols: Int) -> Self:
        '''Initialize an matrix of i's with the specified dimensions.
        
        Args:
            rows: The number of rows in the matrix.
            cols: The number of columns in the matrix.
        
        Returns:
            A matrix filled with the imaginary unit i.
        '''
        var result = Self(rows, cols, fill_zeros=False)
        result.fill_i()
        return result

    @staticmethod
    @always_inline
    fn arange(rows: Int, cols: Int) -> Self:
        '''Initialize a matrix with values from range(self.size) in row-major order.
        
        Args:
            rows: The number of rows in the matrix.
            cols: The number of columns in the matrix.
        
        Returns:
            A matrix with values from range(self.size) in row-major order.
        '''
        var result = Self(rows, cols, fill_zeros=False)
        result.fill_range()
        return result

    @staticmethod
    @always_inline
    fn arange(rows: Int, cols: Int, start: Int) -> Self:
        '''Initialize a matrix with values from range(start, start + self.size) in row-major order.
        
        Args:
            rows: The number of rows in the matrix.
            cols: The number of columns in the matrix.
            start: The start value for the range.
        
        Returns:
            A matrix with values from range(start, start + self.size) in row-major order.
        '''
        var result = Self(rows, cols, fill_zeros=False)
        result.fill_range(start)
        return result

    @staticmethod
    @always_inline
    fn arange(rows: Int, cols: Int, start: Int, step: Int) -> Self:
        '''Initialize a matrix with values from range(start, start + self.size * step, step)
        in row-major order.
        
        Args:
            rows: The number of rows in the matrix.
            cols: The number of columns in the matrix.
            start: The start value for the range.
            step: The step size for the range.
        
        Returns:
            A matrix with values from range(start, start + self.size * step, step)
            in row-major order.
        
        '''
        var result = Self(rows, cols, fill_zeros=False)
        result.fill_range(start, step)
        return result
    
    @staticmethod
    fn random(rows: Int, cols: Int, mean: Scalar[type] = 0, stdev: Scalar[type] = 1) -> Self:
        '''Initialize a random matrix. Real and imaginary parts of all elements are chosen
        from a normal distribution with the specified mean and standard deviation.

        Args:
            rows: The number of rows in the matrix.
            cols: The number of columns in the matrix.
            mean: The mean of the normal distribution.
            stdev: The standard deviation of the normal distribution.

        Returns:
            A random matrix.
        '''
        var result = Self(rows=rows, cols=cols, fill_zeros=False)
        for idx in range(result.size):
            var val = ComplexScalar[Self.type](
                (
                    randn_float64(Float64(mean), Float64(stdev)), 
                    randn_float64(Float64(mean), Float64(stdev)),
                )
            )
            result.store_idx[1](idx, val)
        return result^


    # Assertions ######################

    @always_inline
    fn _assert_same_shape(self, other: Self) raises:
        '''Raise if other is not the same shape as self.'''
        if self.size == 0 and other.size == 0:
            return
        if self.rows != other.rows or self.cols != other.cols:
            raise Error(
                'Incompatible matrix dimensions: ('
                + String(self.rows) + ', ' + String(self.cols) + ') and ('
                + String(other.rows) + ', ' + String(other.cols) + ')'
            )

    @always_inline
    fn _assert_matmul_compatible(self, other: Self) raises:
        '''Raise if the matrix product self * other is not possible due to incompatible matrix
        dimensions.
        '''
        if self.cols != other.rows:
            raise Error(
                'Cannot multiply matrices with shapes ('
                + String(self.rows) + ', ' + String(self.cols) + ') and ('
                + String(other.rows) + ', ' + String(other.cols) + ')'
            )

    @always_inline
    fn _assert_reshape_compatible(self, new_rows: Int, new_cols: Int) raises:
        '''Raise if self cannot be reshaped to a matrix with new_rows rows and new_cols columns.'''
        if new_rows * new_cols != self.size:
            raise Error(
                'Impatible reshape dimensions: ('
                + String(self.rows) + ', ' + String(self.cols) + ') and ('
                + String(new_rows) + ', ' + String(new_cols) + ')'
            )

    # Properties ######################

    @always_inline
    fn __len__(self) -> Int:
        '''Return the total number of elements in the matrix.
        
        Returns:
            The total number of elements in the matrix.
        '''
        return self.size

    @always_inline
    fn shape(self) -> Tuple[Int, Int]:
        '''Return the tuple (rows, cols).
        
        Returns:
            The shape of the matrix as a tuple.    
        '''
        return (self.rows, self.cols)

    @always_inline
    fn is_square(self) -> Bool:
        '''Check if the matrix is square.
        
        Return:
            True if self is square and non-degenerate, and False otherwise.
        '''
        return self.rows == self.cols and self.rows > 0

    @always_inline
    fn is_unitary[tol: Scalar[Self.type] = DEFAULT_TOL](self) raises -> Bool:
        '''Check if the matrix is unitary within the specified tolerance. Specifically, check if
        the matrix is square and the Frobenius norm of `self * self.dagger() - I` is less than tol.

        Parameters:
            tol: The tolerance for the unitarity check.
        
        Returns:
            True if self is unitary, False otherwise.
        '''
        if not self.is_square():
            return False
        return (self @ self.dag() - self.eye_like()).frobenius_norm() < tol

    @always_inline
    fn is_hermitian[tol: Scalar[Self.type] = DEFAULT_TOL](self) raises -> Bool:
        '''Check if the matrix is Hermitian within the specified tolerance. Specifically, check if
        the matrix is square and the Frobenius norm of `self - self.dagger()` is less than tol.

        Parameters:
            tol: The tolerance for the Hermitian check.
        
        Returns:
            True if self is Hermitian, False otherwise.
        '''
        if not self.is_square():
            return False
        return (self - self.dag()).frobenius_norm() < tol

    # String conversion ###############

    @no_inline
    fn __str__(self) -> String:
        '''Convert the matrix to a string.
        
        Returns:
            A string representation of the matrix.
        '''
        return self._get_str_rep(max_lines=6)

    @no_inline
    fn _get_str_rep(self, max_lines: Int) -> String:
        '''Get a string representation of the matrix, showing at most max_lines of data.'''
        alias dots: String = '...'
        if self.size == 0:
            return '[]'
        var mat = List[List[String, True], True]()
        var hide_rows: Bool = self.rows > max_lines and max_lines > 0
        var hide_cols: Bool = self.cols > max_lines and max_lines > 0
        for r in range(max_lines if hide_rows else self.rows):
            if r == max_lines // 2 and hide_rows:
                mat.append(List[String, True](dots))
            mat.append(List[String, True]())
            var r_idx: Int = r if r < max_lines // 2 or not hide_rows else self.rows + r - max_lines
            for c in range(max_lines if hide_cols else self.cols):
                if c == max_lines // 2 and hide_cols:
                    mat[-1].append(dots)
                var c_idx: Int = c if c < max_lines // 2 or not hide_cols else self.cols + c - max_lines
                mat[-1].append(String(self.load_crd[1](r_idx, c_idx)))

        var max_col_widths = List[Int, True]()
        for c in range(len(mat[0])):
            var max_w: Int = 0
            for r in range(len(mat)):
                if mat[r][0] == dots:
                    continue
                max_w = max(max_w, len(mat[r][c]))
            max_col_widths.append(max_w)

        var out: String = '['
        for r in range(len(mat)):
            var line: String
            if mat[r][0] == dots:
                line = ' ...\n'
            else:
                line = ' [' if r > 0 else '['
                for c in range(len(mat[0])):
                    line += mat[r][c].ljust(max_col_widths[c])
                    if c < len(mat[0]) - 1:
                        line += ', '
                line += ']\n'
            out += line
        out = out[:-1]
        out += ']\n'
        out += ' DType: ' + String(Self.type) + ' | Size: (' + String(self.rows) + ', ' + String(self.cols) + ')'
        return out

    @no_inline
    fn __repr__(self) -> String:
        '''Convert the matrix to a string.
        
        Returns:
            A string representation of the matrix.
        '''
        return self._get_str_rep(max_lines=0)

    @no_inline
    fn write_to[W: Writer](self, mut writer: W):
        '''Write the matrix to a writer.
        
        Args:
            writer: The writer to write to.
        '''
        writer.write(String(self))

    # Item access #####################

    @always_inline
    fn __getitem__(self, row: Int, col: Int) raises -> ComplexScalar[Self.type]:
        '''Get the value at the specified row and col.
        
        Args:
            row: The row index.
            col: The column index.
        
        Returns:
            The value at position (row, col).
        '''
        if row < 0 or row >= self.rows:
            raise Error('Invalid row index: ' + String(row))
        if col < 0 or col >= self.cols:
            raise Error('Invalid column index: ' + String(col))
        return self.load_idx[1](row * self.cols + col)

    @always_inline
    fn __getitem__(self, idx: Int) raises -> ComplexScalar[Self.type]:
        '''Get the value at the specified index. Matrices are indexed from index 0 at postion
        (0, 0) to index size - 1 at position (rows - 1, cols - 1) in row-major order.

        Args:
            idx: The index of the value to get.
        
        Returns:
            The value at index idx.
        '''
        if idx < 0 or idx >= self.size:
            raise Error('Invalid index: '+ String(idx))
        return self.load_idx[1](idx)

    @always_inline
    fn __setitem__(self, row: Int, col: Int, val: ComplexScalar[Self.type]) raises:
        '''Set the value at the specified row and col.
        
        Args:
            row: The row index.
            col: The column index.
            val: The value to set.
        '''
        if row < 0 or row >= self.rows:
            raise Error('Invalid row index: ' + String(row))
        if col < 0 or col >= self.cols:
            raise Error('Invalid column index: ' + String(col))
        self.store_crd[1](row, col, val)

    @always_inline
    fn __setitem__(self, idx: Int, val: ComplexScalar[Self.type]) raises:
        '''Set the value at the specified index. Matrices are indexed from index 0 at postion 
        (0, 0) to index size - 1 at position (rows - 1, cols - 1) in row-major order.

        Args:
            idx: The index of the item to set.
            val: The value to set.
        '''
        if idx < 0 or idx >= self.size:
            raise Error('Invalid index: '+ String(idx))
        self.store_idx[1](idx, val)

    @always_inline
    fn load_idx[width: Int](self, idx: Int) -> ComplexSIMD[Self.type, width]:
        '''Load width lanes of the matrix at index idx as a complex SIMD vector.
        
        Parameters:
            width: The SIMD width to load.
        
        Args:
            idx: The index at which to load the SIMD vector.
        
        Returns:
            The SIMD vector.
        '''
        return ComplexSIMD(
            self.re.load[width=width](idx),
            self.im.load[width=width](idx),
        )

    @always_inline
    fn store_idx[width: Int](self, idx: Int, val: ComplexSIMD[Self.type, width]):
        '''Store width lanes of a complex SIMD vector in the matrix at index idx.
        
        Parameters:
            width: The SIMD width.
        
        Args:
            idx: The index at which to store the SIMD vector.
            val: The complex SIMD vector to store.
        '''
        self.re.store(idx, val.re)
        self.im.store(idx, val.im)

    @always_inline
    fn load_crd[width: Int](self, row: Int, col: Int) -> ComplexSIMD[Self.type, width]:
        '''Load width lanes of the ComplexSIMD at position (row, col).'''
        return self.load_idx[width](row * self.cols + col)

    @always_inline
    fn store_crd[width: Int](self, row: Int, col: Int, val: ComplexSIMD[Self.type, width]):
        '''Store val in the ComplexSIMD at position (row, col).
        
        Parameters:
            width: The SIMD width.
        
        Args:
            row: The row index.
            col: The column index.
            val: The complex SIMD vector to store.
        '''
        self.store_idx[width](row * self.cols + col, val)

    @always_inline
    fn strided_load_idx[width: Int](self, idx: Int, stride: Int) -> ComplexSIMD[Self.type, width]:
        '''Load width lanes starting at index idx spaced by stride.
        
        Parameters:
            width: The SIMD width.
        
        Args:
            idx: The index to load.
            stride: The stride between loads.
        
        Returns:
            The complex SIMD vector.
        '''
        return ComplexSIMD(
            (self.re + idx).strided_load[width=width](stride),
            (self.im + idx).strided_load[width=width](stride),
        )

    @always_inline
    fn strided_store_idx[width: Int](self, idx: Int, stride: Int, val: ComplexSIMD[Self.type, width]):
        '''Store width lanes starting at index idx spaced by stride.
        
        Parameters:
            width: The SIMD width.
        
        Args:
            idx: The index at which to store the SIMD vector.
            stride: The stride between stores.
            val: The complex SIMD vector to store.
        '''
        (self.re + idx).strided_store[width=width](val.re, stride)
        (self.im + idx).strided_store[width=width](val.im, stride)
    
    fn get_column(
        self, col: Int, row_start: Int = 0, row_end: Int = -1
    ) raises -> CMatrix[Self.type]:
        '''Get the specified column between the specified row indices.

        Args:
            col: The index of the column to return.
            row_start: The row start index. Default is 0.
            row_end: The row end index. If negative, it is taken to be self.rows. Default is -1.
        
        Returns:
            The specified column values.
        '''
        if row_start < 0 or row_start >= self.rows:
            raise Error('Invalid row start index: ' + String(row_start))
        if row_end > self.rows:
            raise Error('Invalid row end index: ' + String(row_end))
        if col < 0 or col >= self.cols:
            raise Error('Invalid column index: ' + String(col))
        if row_end >= 0 and row_start >= row_end:
            return CMatrix[type](rows=0, cols=0)
        return self._get_column(col, row_start, row_end)
    
    fn _get_column(
        self, col: Int, owned row_start: Int = 0, owned row_end: Int = -1
    ) -> CMatrix[Self.type]:
        '''Get the specified column between the specified row indices without checking index
        validity.

        Args:
            col: The index of the column to return.
            row_start: The row start index. Default is 0.
            row_end: The row end index. If negative, it is taken to be self.rows. Default is -1.
        
        Returns:
            The specified column values.
        '''
        if row_end < 0:
            row_end = self.rows
        
        var result = CMatrix[type](rows=row_end - row_start, cols=1, fill_zeros=False)
        for i in range(result.rows):
            result.store_idx[1](i, self.load_crd[1](row_start + i, col))
        return result
    
    fn get_row(
        self, row: Int, col_start: Int = 0, col_end: Int = -1
    ) raises -> CMatrix[Self.type]:
        '''Get the specified row between the specified column indices.

        Args:
            row: The index of the row to return.
            col_start: The column start index. Default is 0.
            col_end: The column end index. If negative, it is taken to be self.cols. Default is -1.
        
        Returns:
            The specified row values.
        '''
        if col_start < 0 or col_start >= self.cols:
            raise Error('Invalid column start index: ' + String(col_start))
        if col_end > self.cols:
            raise Error('Invalid column end index: ' + String(col_end))
        if row < 0 or row >= self.rows:
            raise Error('Invalid row index: ' + String(row))
        if col_end >= 0 and col_start >= col_end:
            return CMatrix[type](rows=0, cols=0)
        return self._get_row(row, col_start, col_end)
    
    fn _get_row(
        self, row: Int, owned col_start: Int = 0, owned col_end: Int = -1
    ) -> CMatrix[Self.type]:
        '''Get the specified column between the specified row indices without checking index
        validity.

        Args:
            row: The index of the row to return.
            col_start: The column start index. Default is 0.
            col_end: The column end index. If negative, it is taken to be self.cols. Default is -1.

        Returns:
            The specified row values.
        '''
        if col_end < 0:
            col_end = self.cols
        
        var result = CMatrix[type](rows=1, cols=col_end - col_start, fill_zeros=False)
        for i in range(result.cols):
            result.store_idx[1](i, self.load_crd[1](row, col_start + i))
        return result
    
    fn get_block(
        self, row_start: Int = 0, row_end: Int = -1, col_start: Int = 0, col_end: Int = -1
    ) raises -> CMatrix[Self.type]:
        '''Get the specified block of the matrix.

        Args:
            row_start: The row start index. Default is 0.
            row_end: The row end index. If negative, it is taken to be self.rows, Default is -1.
            col_start: The column start index. Default is 0.
            col_end: The column end index. If negative, it is taken to be self.cols, Default is -1.
        
        Returns:
            The block with rows between row_start and row_end and columns between col_start 
            and col_end.
        '''
        if row_start < 0 or row_start >= self.rows:
            raise Error('Invalid row start index: ' + String(row_start))
        if col_start < 0 or col_start >= self.cols:
            raise Error('Invalid column start index: ' + String(col_start))
        if row_end > self.rows:
            raise Error('Invalid row end index: ' + String(row_end))
        if col_end > self.cols:
            raise Error('Invalid column end index: ' + String(col_end))
        if (row_end >= 0 and row_start >= row_end) or (col_end >= 0 and col_start >= col_end):
            return CMatrix[Self.type](rows=0, cols=0)
        return self._get_block(row_start, row_end, col_start, col_end)
    
    fn _get_block(
        self,
        owned row_start: Int = 0,
        owned row_end: Int = -1,
        owned col_start: Int = 0,
        owned col_end: Int = -1,
    ) -> CMatrix[Self.type]:
        '''Get the specified block of the matrix without cehcking index validity.

        Args:
            row_start: The row start index. Default is 0.
            row_end: The row end index. If negative, it is taken to be self.rows, Default is -1.
            col_start: The column start index. Default is 0.
            col_end: The column end index. If negative, it is taken to be self.cols, Default is -1.
        
        Returns:
            The block with rows between row_start and row_end and columns between col_start 
            and col_end.
        '''
        if row_end < 0:
            row_end = self.rows
        if col_end < 0:
            col_end = self.cols
        
        var result = CMatrix[Self.type](
            rows=row_end - row_start, cols=col_end - col_start, fill_zeros=False
        )
        for r in range(result.rows):
            for c in range(result.cols):
                result.store_crd[1](r, c, self.load_crd[1](row_start + r, col_start + c))
        return result
    
    fn set_block(self, row: Int, col: Int, block: CMatrix[Self.type]) raises:
        '''Set the specified block in self with the supplied block.

        Args:
            row: The row start index.
            col: The column start index.
            block: The block of values to set.
        '''
        if block.rows <= 0 or block.cols <= 0:
            return
        if row < 0 or row >= self.rows:
            raise Error('Invalid row start index: ' + String(row))
        if col < 0 or col >= self.cols:
            raise Error('Invalid column start index: ' + String(col))
        if row + block.rows > self.rows or col + block.cols > self.cols:
            raise Error('Block overflows matrix at specified position.')
        self._set_block(row, col, block)
    
    fn _set_block(self, row: Int, col: Int, block: CMatrix[Self.type]):
        '''Set the specified block in self with the supplied block without checking index or block
        validity.

        Args:
            row: The row start index.
            col: The column start index.
            block: The block of values to set.
        '''
        for r in range(block.rows):
            for c in range(block.cols):
                self.store_crd[1](row + r, col + c, block.load_crd[1](r, c))
    
    # Helper functions ################

    fn _parallelize_vectorize_op[
        func: fn[width: Int](r: Int, c: Int) capturing -> ComplexSIMD[Self.type, width],
        vectorize_only: Bool = False,
    ](self) -> Self:
        '''A helper function for parallelizing over rows and vectorizing over columns.
        
        Parameters:
            func: The function to apply at each coordinate.
            vectorize_only: If True, only vectorization will be applied.

        Returns:
            A CMatrix in which entry (r, c) holds the value func(r, c).
        '''
        var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=False)
        if vectorize_only or not self._should_parallelize():
            for r in range(self.rows):
                @parameter
                fn op_col0[simd_width: Int](c: Int):
                    result.store_crd[simd_width](r, c, func[simd_width](r, c))
                vectorize[op_col0, simdwidthof[Self.type]()](self.cols)
        else:
            @parameter
            fn op_row(r: Int):
                @parameter
                fn op_col1[simd_width: Int](c: Int):
                    result.store_crd[simd_width](r, c, func[simd_width](r, c))
                vectorize[op_col1, simdwidthof[Self.type]()](self.cols)
            parallelize[op_row](self.rows)
        return result

    fn _parallelize_vectorize_op_inplace[
        func: fn[width: Int](r: Int, c: Int) capturing -> ComplexSIMD[Self.type, width],
        vectorize_only: Bool = False,
    ](self):
        '''A helper function for parallelizing over rows and vectorizing over columns.
        Stores the value func(r, c) at each position (r, c) in the matrix.
        
        Parameters:
            func: The function to apply at each coordinate.
            vectorize_only: If True, only vectorization will be applied.
        '''
        if vectorize_only or not self._should_parallelize():
            for r in range(self.rows):
                @parameter
                fn op_col0[simd_width: Int](c: Int):
                    self.store_crd[simd_width](r, c, func[simd_width](r, c))
                vectorize[op_col0, simdwidthof[Self.type]()](self.cols)
        else:
            @parameter
            fn op_row(r: Int):
                @parameter
                fn op_col1[simd_width: Int](c: Int):
                    self.store_crd[simd_width](r, c, func[simd_width](r, c))
                vectorize[op_col1, simdwidthof[Self.type]()](self.cols)
            parallelize[op_row](self.rows)

    # Math dunders ####################

    fn __neg__(self) -> Self:
        '''Defines the `-` unary negation operator.
        
        Returns:
            The negative of self.
        '''
        @parameter
        fn neg_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return -self.load_crd[simd_width](r, c)
        return self._parallelize_vectorize_op[neg_r, True]()  # Small operation, vectorize only

    @always_inline
    fn __pos__(self) -> Self:
        '''Defines the `+` unary positive operator.
        
        Returns:
            The positive of self.
        '''
        return self

    fn __add__(self, other: Self) raises -> Self:
        '''Defines the `+` add operator.
        
        Args:
            other: The matrix to add to self.
        
        Returns:
            The sum of self and other.
        '''
        self._assert_same_shape(other)
        @parameter
        fn add_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) + other.load_crd[simd_width](r, c)
        return self._parallelize_vectorize_op[add_r]()

    fn __add__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `+` add operator.
        
        Args:
            other: The number to add to self.
        
        Returns:
            The sum of self and other.
        '''
        @parameter
        fn add_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) + other
        return self._parallelize_vectorize_op[add_r]()

    fn __sub__(self, other: Self) raises -> Self:
        '''Defines the `-` subtraction operator.
        
        Args:
            other: The matrix to subtract from self.
        
        Returns:
            The difference of self and other.
        '''
        self._assert_same_shape(other)
        @parameter
        fn sub_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) - other.load_crd[simd_width](r, c)
        return self._parallelize_vectorize_op[sub_r]()

    fn __sub__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `-` subtraction operator.
        
        Args:
            other: The number to subtract from self.
        
        Returns:
            The difference of self and other.
        '''
        @parameter
        fn sub_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) - other
        return self._parallelize_vectorize_op[sub_r]()

    fn __mul__(self, other: Self) raises -> Self:
        '''Defines the `*` product operator.
        
        Args:
            other: The matrix to multiply with self.
        
        Returns:
            The elementwise product of self and other.
        '''
        self._assert_same_shape(other)
        @parameter
        fn mul_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) * other.load_crd[simd_width](r, c)
        return self._parallelize_vectorize_op[mul_r]()

    fn __mul__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `*` product operator.
        
        Args:
            other: The number to multiply with self.
        
        Returns:
            The elementwise product of self and other.
        '''
        @parameter
        fn mul_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) * other
        return self._parallelize_vectorize_op[mul_r]()

    fn __truediv__(self, other: Self) raises -> Self:
        '''Defines the `/` divide operator.
        
        Args:
            other: The matrix to divide self by.
        
        Returns:
            The elementwise quotient of self with other.
        '''
        self._assert_same_shape(other)
        @parameter
        fn div_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) / other.load_crd[simd_width](r, c)
        return self._parallelize_vectorize_op[div_r]()

    @always_inline
    fn __truediv__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `/` divide operator.
        
        Args:
            other: The number to divide self by.
        
        Returns:
            The elementwise quotient of self with other.
        '''
        return self.__mul__(other.reciprocal())

    fn __floordiv__(self, other: Self) raises -> Self:
        '''Defines the `//` floor divide operator.
        
        Args:
            other: The matrix to floor divide self by.
        
        Returns:
            The elementwise floor division of self with other.
        '''
        self._assert_same_shape(other)
        @parameter
        fn fdiv_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) // other.load_crd[simd_width](r, c)
        return self._parallelize_vectorize_op[fdiv_r]()

    fn __floordiv__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `//` floor divide operator.
        
        Args:
            other: The number to floor divide self by.
        
        Returns:
            The elementwise floor division of self with other.
        '''
        @parameter
        fn fdiv_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) // other
        return self._parallelize_vectorize_op[fdiv_r]()

    fn __mod__(self, other: Self) raises -> Self:
        '''Defines the `%` mod operator.
        
        Args:
            other: The matrix to modulo with self.

        Returns:
            The elementwise modulo of self with other.
        '''
        self._assert_same_shape(other)
        @parameter
        fn mod_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) % other.load_crd[simd_width](r, c)
        return self._parallelize_vectorize_op[mod_r]()

    fn __mod__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `%` mod operator.
        
        Args:
            other: The number to modulo with self.
        
        Returns:
            The elementwise modulo of self with other.
        
        '''
        @parameter
        fn mod_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) % other
        return self._parallelize_vectorize_op[mod_r]()

    @always_inline
    fn __divmod__(self, other: Self) raises -> Tuple[Self, Self]:
        '''Defines the divmod operator.
        
        Args:
            other: The matrix to divmod with self.
        
        Returns:
            The elementwise divmod of self with other.
        '''
        return (self // other, self % other)

    @always_inline
    fn __divmod__(self, other: ComplexScalar[Self.type]) -> Tuple[Self, Self]:
        '''Defines the divmod operator.
        
        Args:
            other: The matrix to divmod with self.
        
        Returns:
            The elementwise divmod of self with other.
        '''
        return (self // other, self % other)

    fn __matmul__(self, other: Self) raises -> Self:
        '''Defines the `@` matrix multiply operator.
        
        Args:
            other: The matrix to multiply with self.
        
        Returns:
            The matrix product of self with other.
        '''
        self._assert_matmul_compatible(other)
        var result = Self(rows=self.rows, cols=other.cols, fill_zeros=True)
        @parameter
        fn calc_row(r: Int):
            for k in range(self.cols):
                @parameter
                fn dot[simd_width: Int](c: Int):
                    result.store_crd[simd_width](
                        r,
                        c,
                        result.load_crd[simd_width](r, c)
                        + self.load_crd[1](r, k)
                        * other.load_crd[simd_width](k, c),
                    )
                vectorize[dot, simdwidthof[Self.type]()](result.cols)
        # Always parallelize matmul
        parallelize[calc_row](result.rows)
        return result

    # In-place math dunders ###########

    fn __iadd__(self, other: Self) raises:
        '''Defines the `+=` in-place add operator.
        
        Args:
            other: The matrix to add to self.
        '''
        self._assert_same_shape(other)
        @parameter
        fn add_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) + other.load_crd[simd_width](r, c)
        self._parallelize_vectorize_op_inplace[add_r]()

    fn __iadd__(self, other: ComplexScalar[Self.type]):
        '''Defines the `+=` in-place add operator.
        
        Args:
            other: The number to add to self.
        '''
        @parameter
        fn add_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) + other
        self._parallelize_vectorize_op_inplace[add_r]()

    fn __isub__(self, other: Self) raises:
        '''Defines the `-=` in-place subtraction operator.
        
        Args:
            other: The matrix to subtract from self.
        '''
        self._assert_same_shape(other)
        @parameter
        fn sub_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) - other.load_crd[simd_width](r, c)
        self._parallelize_vectorize_op_inplace[sub_r]()

    fn __isub__(self, other: ComplexScalar[Self.type]):
        '''Defines the `-=` in-place subtraction operator.
        
        Args:
            other: The number to subtract from self.
        '''
        @parameter
        fn sub_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) - other
        self._parallelize_vectorize_op_inplace[sub_r]()

    fn __imul__(self, other: Self) raises:
        '''Defines the `*=` in-place product operator.
        
        Args:
            other: The matrix to multiply elementwise with self.
        '''
        self._assert_same_shape(other)
        @parameter
        fn mul_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) * other.load_crd[simd_width](r, c)
        self._parallelize_vectorize_op_inplace[mul_r]()

    fn __imul__(self, other: ComplexScalar[Self.type]):
        '''Defines the `*=` in-place product operator.
        
        Args:
            other: The number to multiply elementwise with self.
        '''
        @parameter
        fn mul_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) * other
        self._parallelize_vectorize_op_inplace[mul_r]()

    fn __itruediv__(self, other: Self) raises:
        '''Defines the `/=` in-place divide operator.
        
        Args:
            other: The matrix to divide self by elementwise.
        '''
        self._assert_same_shape(other)
        @parameter
        fn div_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) / other.load_crd[simd_width](r, c)
        self._parallelize_vectorize_op_inplace[div_r]()

    @always_inline
    fn __itruediv__(self, other: ComplexScalar[Self.type]):
        '''Defines the `/=` in-place divide operator.
        
        Args:
            other: The number to divide self by elementwise.
        '''
        self.__imul__(other.reciprocal())

    fn __ifloordiv__(self, other: Self) raises:
        '''Defines the `//=` in-place floor divide operator.
        
        Args:
            other: The matrix to floor divide self by elementwise.
        '''
        self._assert_same_shape(other)
        @parameter
        fn fdiv_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) // other.load_crd[simd_width](r, c)
        self._parallelize_vectorize_op_inplace[fdiv_r]()

    fn __ifloordiv__(self, other: ComplexScalar[Self.type]):
        '''Defines the `//=` in-place floor divide operator.
        
        Args:
            other: The number to divide self by elementwise.
        '''
        @parameter
        fn fdiv_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) // other
        self._parallelize_vectorize_op_inplace[fdiv_r]()

    fn __imod__(self, other: Self) raises:
        '''Defines the `%=` in-place mod operator.
        
        Args:
            other: The matrix to modulo self by elementwise.
        '''
        self._assert_same_shape(other)
        @parameter
        fn mod_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) % other.load_crd[simd_width](r, c)
        self._parallelize_vectorize_op_inplace[mod_r]()

    fn __imod__(self, other: ComplexScalar[Self.type]):
        '''Defines the `%=` in-place mod operator.
        
        Args:
            other: The number to modulo self by elementwise.
        '''
        @parameter
        fn mod_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c) % other
        self._parallelize_vectorize_op_inplace[mod_r]()

    # TODO: Make it better
    @always_inline
    fn __imatmul__(mut self, other: Self) raises:
        '''Defines the `@=` in-place matrix multiply operator. Only valid for square matrices.
        
        Args:
            other: The matrix to matrix multiply with self.
        '''
        self = self @ other
        # self._assert_matmul_compatible(other)
        # var result = Self(rows=self.rows, cols=other.cols, fill_zeros=True)
        # @parameter
        # fn calc_row(r: Int):
        #     for k in range(self.cols):
        #         @parameter
        #         fn dot[simd_width: Int](c: Int):
        #             result.store_crd[simd_width](
        #                 r,
        #                 c,
        #                 result.load_crd[simd_width](r, c)
        #                 + self.load_crd[1](r, k)
        #                 * other.load_crd[simd_width](k, c),
        #             )
        #         vectorize[dot, simdwidthof[Self.type]()](result.cols)
        # parallelize[calc_row](result.rows)
        # self = result^

    # Right math dunders ##############

    @always_inline
    fn __radd__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `+` right add operator.
        
        Args:
            other: The number to add to self.

        Returns:
            The sum of other and self.
        '''
        return self + other

    fn __rsub__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `-` right subtraction operator.
        
        Args:
            other: The number to subtract from self.

        Returns:
            The difference of other and self.
        '''
        @parameter
        fn sub_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return other - self.load_crd[simd_width](r, c)
        return self._parallelize_vectorize_op[sub_r]()

    @always_inline
    fn __rmul__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `*` right product operator.
        
        Args:
            other: The number to multiply with self.

        Returns:
            The product of other and self.
        '''
        return self * other

    fn __rtruediv__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `/` right divide operator.
        
        Args:
            other: The number to divide by self.

        Returns:
            The quotient of other and self.
        '''
        @parameter
        fn div_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return other / self.load_crd[simd_width](r, c)
        return self._parallelize_vectorize_op[div_r]()

    fn __rfloordiv__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `//` right floor divide operator.
        
        Args:
            other: The number to floow divide by self.

        Returns:
            The floor division of other and self.
        '''
        @parameter
        fn fdiv_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return other // self.load_crd[simd_width](r, c)
        return self._parallelize_vectorize_op[fdiv_r]()

    fn __rmod__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `%` right mod operator.
        
        Args:
            other: The number to modulo by self.

        Returns:
            The modulo of other and self.
        '''
        @parameter
        fn mod_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return other % self.load_crd[simd_width](r, c)
        return self._parallelize_vectorize_op[mod_r]()

    @always_inline
    fn __rdivmod__(self, other: ComplexScalar[Self.type]) -> Tuple[Self, Self]:
        '''Defines the right divmod operator.
        
        Args:
            other: The number to divmod with self.

        Returns:
            The divmod of other and self.
        '''
        return (other // self, other % self)

    # Other math ######################

    fn __abs__(self) -> Self:
        '''Applies the absolute value to each element.

        Returns:
            The absolute value of self.
        '''
        @parameter
        fn abs_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return self.load_crd[simd_width](r, c).__abs__()
        return self._parallelize_vectorize_op[abs_r]()

    fn conj(self) -> Self:
        '''The conjugate of the matrix.
        
        Returns:
            The matrix with the conjugate applied to each element.
        '''
        var result: Self = self
        @parameter
        fn conj_simd[simd_width: Int](idx: Int):
            result.im.store(idx, -self.im.load[width=simd_width](idx))
        vectorize[conj_simd, simdwidthof[Self.type]()](self.size)
        return result

    fn iconj(self):
        '''Conjugate the matrix in-place.'''
        @parameter
        fn conj_simd[simd_width: Int](idx: Int):
            self.im.store(idx, -self.im.load[width=simd_width](idx))
        vectorize[conj_simd, simdwidthof[Self.type]()](self.size)

    fn dag(self) -> Self:
        '''The conjugate-transpose of the matrix.
        
        Returns:
            The conjugate-transpose of self.
        '''
        var result = Self(rows=self.cols, cols=self.rows, fill_zeros=False)
        @parameter
        fn transpose_row(r: Int):
            @parameter
            fn transpose_col[simd_width: Int](c: Int):
                result.strided_store_idx[simd_width](
                    c * result.cols + r,
                    result.cols,
                    self.load_crd[simd_width](r, c).conj(),
                )
            vectorize[transpose_col, simdwidthof[Self.type]()](self.cols)
        parallelize[transpose_row](self.rows)
        return result

    @always_inline
    fn dagger(self) -> Self:
        '''The conjugate-transpose of the matrix. Alias of dag().
        
        Returns:
            The conjugate-transpose of self.
        '''
        return self.dag()

    @always_inline
    fn sum(self) raises -> ComplexScalar[Self.type]:
        '''Sum all elements in the matrix.

        Returns:
            The sum of all elements in the matrix.
        '''
        var row_sums = Self(rows=self.rows, cols=1, fill_zeros=False)
        @parameter
        fn sum_rows(r: Int):
            var row_sum = ComplexScalar[Self.type]()
            for c in range(self.cols):
                row_sum += self.load_crd[1](r, c)
            row_sums.store_idx[1](r, row_sum)
        parallelize[sum_rows](self.rows)
        var total = ComplexScalar[Self.type]()
        for r in range(self.rows):
            total += row_sums.load_idx[1](r)
        return total

    fn trace(self) -> ComplexScalar[Self.type]:
        '''Compute the trace of a matrix. This sums the elements along the main diagonal,
        regardless of the shape of the matrix.

        Returns:
            The trace of the matrix.
        '''
        var trace = ComplexScalar[Self.type].zero()
        for i in range(min(self.rows, self.cols)):
            trace += self.load_crd[1](i, i)
        return trace

    fn echelon[tol: Scalar[Self.type] = DEFAULT_TOL](self) -> Self:
        '''Compute the row echelon form of a matrix.
        
        Parameters:
            tol: The tolerance for zero checks.
        
        Returns:
            The row echelon form of self.
        '''
        var A = self
        var h: Int = 0
        var k: Int = 0
        while h < self.rows and k < self.cols:
            var i_max: Scalar[Self.type] = A.load_crd[1](h, k).norm()
            var i_argmax: Int = h
            for i in range(h + 1, self.rows):
                var i_norm: Scalar[Self.type] = A.load_crd[1](i, k).norm()
                if i_norm > i_max:
                    i_max = i_norm
                    i_argmax = i
            if A.load_crd[1](i_argmax, k).norm() < tol:
                k += 1
            else:
                for c in range(self.cols):
                    var p: ComplexScalar[Self.type] = A.load_crd[1](h, c)
                    A.store_crd[1](h, c, A.load_crd[1](i_argmax, c))
                    A.store_crd[1](i_argmax, c, p)
                for i in range(h + 1, self.rows):
                    var f: ComplexScalar[Self.type] = A.load_crd[1](i, k) / A.load_crd[1](h, k)
                    A.store_crd[1](i, k, 0)
                    for j in range(k + 1, self.cols):
                        A.store_crd[1](i, j, A.load_crd[1](i, j) - A.load_crd[1](h, j) * f)
                h += 1; k += 1
        return A

    fn det[tol: Scalar[Self.type] = DEFAULT_TOL](self) raises -> ComplexScalar[Self.type]:
        '''Compute the determinant of a matrix.
        
        Parameters:
            tol: The tolerance for zero checks when computing row echelon form.
        
        Returns:
            The determinant of self.
        '''
        if self.rows != self.cols:
            raise Error('Only square matrices have determinants')
        var echelon: Self = self.echelon()
        var d = ComplexScalar[Self.type](1)
        for i in range(self.rows):
            d *= echelon.load_crd[1](i, i)
        return d

    fn determinant[tol: Scalar[Self.type] = DEFAULT_TOL](self) raises -> ComplexScalar[Self.type]:
        '''Compute the determinant of a matrix. Alias of det().
        
        Parameters:
            tol: The tolerance for zero checks when computing row echelon form.
        
        Returns:
            The determinant of self.
        '''
        return self.det[tol]()

    fn inv[tol: Scalar[Self.type] = DEFAULT_TOL](self) raises -> Self:
        '''Compute the inverse of a square matrix.
        
        Parameters:
            tol: The tolerance for zero checks.

        Returns:
            The inverse of self.
        '''
        alias zero = ComplexScalar[Self.type](0)
        alias one = ComplexScalar[Self.type](1)

        if self.rows != self.cols:
            raise Error('Only square matrices are invertible')

        var augmented = Self(rows=self.rows, cols=2 * self.cols, fill_zeros=True)
        for r in range(self.rows):
            memcpy(augmented.re + r * augmented.cols, self.re + r * self.cols, self.cols)
            memcpy(augmented.im + r * augmented.cols, self.im + r * self.cols, self.cols)
            augmented.store_crd[1](r, self.cols + r, one)

        for i in range(self.rows):
            if augmented[i, i].norm() < tol:
                for j in range(i + 1, self.rows):
                    if augmented[j, i].norm() >= tol:
                        var i_start: Int = i * augmented.cols
                        var j_start: Int = j * augmented.cols
                        var row_i = Self(rows=1, cols=augmented.cols, fill_zeros=False)
                        memcpy(row_i.re, augmented.re + i_start, augmented.cols)
                        memcpy(row_i.im, augmented.im + i_start, augmented.cols)
                        memcpy(augmented.re + i_start, augmented.re + j_start, augmented.cols)
                        memcpy(augmented.im + i_start, augmented.im + j_start, augmented.cols)
                        memcpy(augmented.re + j_start, row_i.re, augmented.cols)
                        memcpy(augmented.im + j_start, row_i.im, augmented.cols)
                        break
                else:
                    raise Error('Matrix is not invertible')

            var pivot: ComplexScalar[Self.type] = augmented.load_crd[1](i, i)
            for j in range(augmented.cols):
                augmented.store_crd[1](i, j, augmented.load_crd[1](i, j) / pivot)

            for j in range(self.rows):
                if j != i:
                    var factor: ComplexScalar[Self.type] = augmented.load_crd[1](j, i)
                    for k in range(augmented.cols):
                        augmented.store_crd[1](
                            j,
                            k,
                            augmented.load_crd[1](j, k) - augmented.load_crd[1](i, k) * factor
                        )

        var result = Self(rows=self.rows, cols=self.cols, fill_zeros=False)
        for i in range(self.rows):
            memcpy(result.re + i * self.cols, augmented.re + i * augmented.cols + self.cols, self.cols)
            memcpy(result.im + i * self.cols, augmented.im + i * augmented.cols + self.cols, self.cols)
        return result

    @always_inline
    fn inverse[tol: Scalar[Self.type] = DEFAULT_TOL](self) raises -> Self:
        '''Compute the inverse of a square matrix. Alias of inv().
        
        Parameters:
            tol: The tolerance for zero checks.

        Returns:
            The inverse of self.
        '''
        return self.inv()

    fn frobenius_norm(self) -> Scalar[Self.type]:
        '''Compute the Frobenius norm of self.
        
        Returns:
            The Frobenius norm of self.
        '''
        var total = Atomic[Self.type](0)
        @parameter
        fn sum_row_squared_norms(r: Int):
            var row_sum: Scalar[Self.type] = 0
            for c in range(self.cols):
                row_sum += self.load_crd[1](r, c).squared_norm()
            _ = total.fetch_add(row_sum)
        parallelize[sum_row_squared_norms](self.rows)
        return sqrt(total.load())

    # Shape operations ################

    @always_inline
    fn reshape(self, new_rows: Int, new_cols: Int) raises -> Self:
        '''Reshape the matrix.
        
        Args:
            new_rows: The number of rows in the new matrix.
            new_cols: The number of columns in the new matrix.
        
        Returns:
            The reshaped matrix.
        '''
        self._assert_reshape_compatible(new_rows, new_cols)
        var result = Self(rows=new_rows, cols=new_cols, fill_zeros=False)
        memcpy(result.re, self.re, self.size)
        memcpy(result.im, self.im, self.size)
        return result

    @always_inline
    fn ireshape(mut self, new_rows: Int, new_cols: Int) raises:
        '''Reshape the matrix, in-place.
        
        Args:
            new_rows: The number of rows in the new matrix.
            new_cols: The number of columns in the new matrix.
        '''
        self._assert_reshape_compatible(new_rows, new_cols)
        self.rows = new_rows
        self.cols = new_cols

    @always_inline
    fn flatten_to_row(self) -> Self:
        '''Flatten the matrix to a row vector.
        
        Returns:
            The flattened matrix.
        '''
        var result = Self(rows=1, cols=self.size, fill_zeros=False)
        memcpy(result.re, self.re, self.size)
        memcpy(result.im, self.im, self.size)
        return result

    @always_inline
    fn iflatten_to_row(mut self):
        '''Flatten the matrix to a row vector, in-place.'''
        self.rows = 1
        self.cols = self.size

    @always_inline
    fn flatten_to_column(self) -> Self:
        '''Flatten the matrix to a column vector.
        
        Returns:
            The flattened matrix.
        '''
        var result = Self(rows=self.size, cols=1, fill_zeros=False)
        memcpy(result.re, self.re, self.size)
        memcpy(result.im, self.im, self.size)
        return result

    @always_inline
    fn iflatten_to_column(mut self):
        '''Flatten the matrix to a column vector, in-place.'''
        self.rows = self.size
        self.cols = 1

    fn transpose(self) -> Self:
        '''Transpose the matrix.
        
        Returns:
            The transpose of self.
        '''
        var result = Self(rows=self.cols, cols=self.rows, fill_zeros=False)
        @parameter
        fn transpose_row(r: Int):
            for c in range(self.cols):
                result.store_crd[1](c, r, self.load_crd[1](r, c))
        parallelize[transpose_row](self.rows)
        return result

    fn itranspose(mut self) raises:
        '''Compute the transpose of the matrix in-place.'''
        if self.size == 0:
            # For empty matrices, just swap dimensions
            self.rows, self.cols = self.cols, self.rows
            return

        if self.rows == self.cols:
            for r in range(self.rows):
                # The diagonal elements (r == c) remain unchanged.
                for c in range(r + 1, self.cols):
                    var idx1: Int = r * self.cols + c  # Index of self[r, c]
                    var idx2: Int = c * self.cols + r  # Index of self[c, r]

                    var temp_re_val: Scalar[Self.type] = self.re.load[width=1](idx1)
                    self.re.store(idx1, self.re.load[width=1](idx2))
                    self.re.store(idx2, temp_re_val)
                    var temp_im_val: Scalar[Self.type] = self.im.load[width=1](idx1)
                    self.im.store(idx1, self.im.load[width=1](idx2))
                    self.im.store(idx2, temp_im_val)
            return

        # Non-square matrix transpose requires temporary storage for re-arrangement
        var original_rows: Int = self.rows
        var original_cols: Int = self.cols
        var new_rows: Int = original_cols  # Transposed rows
        var new_cols: Int = original_rows  # Transposed columns

        # Allocate temporary buffers for the transposed data
        var temp_re_ptr = UnsafePointer[Scalar[Self.type]].alloc(self.size)
        var temp_im_ptr = UnsafePointer[Scalar[Self.type]].alloc(self.size)
        # # Initialize temporary buffers to zero to ensure clean state
        # memset_zero(temp_re_ptr, self.size)
        # memset_zero(temp_im_ptr, self.size)

        @parameter
        fn process_original_row(r_orig: Int):
            # r_orig is the row index in the original matrix
            @parameter
            fn process_original_col_block[simd_width: Int](c_orig_block_start: Int):
                # c_orig_block_start is the starting column index for a SIMD block in the
                # original matrix.

                # Load a row segment from self
                # self[r_orig, c_orig_block_start ... c_orig_block_start + simd_width - 1]
                var val_simd: ComplexSIMD[Self.type, simd_width] = self.load_crd[simd_width](
                    r_orig, c_orig_block_start
                )
                
                # Store this SIMD vector as a column segment in the temporary buffers
                var target_idx_start: Int = c_orig_block_start * new_cols + r_orig
                
                # The stride to move to the next element in the same column of the transposed
                # view is new_cols
                (temp_re_ptr + target_idx_start).strided_store[width=simd_width](
                    val_simd.re, new_cols
                )
                (temp_im_ptr + target_idx_start).strided_store[width=simd_width](
                    val_simd.im, new_cols
                )
            vectorize[process_original_col_block, simdwidthof[Self.type]()](original_cols)
        parallelize[process_original_row](original_rows)

        # Free old memory
        self.re.free()
        self.im.free()

        self.re = temp_re_ptr
        self.im = temp_im_ptr

        self.rows = new_rows
        self.cols = new_cols

    # Fill operations #################

    @always_inline
    fn fill_zero(self):
        '''Fill a matrix with zeros in-place.'''
        memset_zero(self.re.address, self.size)
        memset_zero(self.im.address, self.size)

    fn fill(self, val: ComplexScalar[Self.type]):
        '''Fill a matrix with val in-place.
        
        Args:
            val: The value to fill the amtrix with.
        '''
        @parameter
        fn fill_col[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
            return ComplexSIMD[Self.type, simd_width](val.re, val.im)
        self._parallelize_vectorize_op_inplace[fill_col, True]()  # Small operation, vectorize only

    @always_inline
    fn fill_one(self):
        '''Fill a matrix with ones in-place.'''
        self.fill(1)

    @always_inline
    fn fill_i(self):
        '''Fill a matrix with the imaginary unit i in-place.'''
        self.fill(ComplexScalar[Self.type]().i())

    fn fill_diag(self, val: ComplexScalar[Self.type], offset: Int = 0):
        '''Fill the diagonal at index offset of a matrix with val in-place.
        
        Args:
            val: The value to fill along the diagonal.
            offset: The offset of the diagonal to fill.
        '''
        var stride = self.cols + 1
        var n_diag_elements: Int = (
            min(self.rows, self.cols - offset)
            if offset >= 0
            else min(self.rows + offset, self.cols)
        )
        @parameter
        fn fill_diag_simd_rp[simd_width: Int](p: Int):
            var idx: Int = p * stride + offset if offset >= 0 else p * stride - offset * self.cols
            self.strided_store_idx[simd_width](
                idx, stride, ComplexSIMD[Self.type, simd_width](val.re, val.im)
            )
        vectorize[fill_diag_simd_rp, simdwidthof[Self.type]()](n_diag_elements)

    fn fill_range(self):
        '''Fill the matrix with values from range(self.size) in row-major order, in-place.'''
        memset_zero(self.im.address, self.size)
        for idx in range(self.size):
            self.re.store(idx, SIMD[Self.type, 1](idx))

    fn fill_range(self, start: Int):
        '''Fill the matrix with values from range(start, start + self.size) in row-major order, 
        in-place.
        
        Args:
            start: The start value for the range.
        '''
        memset_zero(self.im.address, self.size)
        for idx in range(self.size):
            self.re.store(idx, SIMD[Self.type, 1](idx + start))

    fn fill_range(self, start: Int, step: Int):
        '''Fill the matrix with values from range(start, start + self.size * step, step)
        in row-major order, in-place.

        Args:
            start: The start value for the range.
            step: The step size for the range.
        '''
        memset_zero(self.im.address, self.size)
        for idx in range(self.size):
            self.re.store(idx, SIMD[Self.type, 1](step * idx + start))

    @always_inline
    fn range_like(self) -> Self:
        '''Return a matrix with values from range(self.size) in row-major order.
        
        Returns:
            A range matrix with the same shape as self.
        '''
        var result = Self(self.rows, self.cols, fill_zeros=False)
        result.fill_range()
        return result

    @always_inline
    fn range_like(self, start: Int) -> Self:
        '''Return a matrix with values from range(start, start + self.size) in row-major order.
        
        Args:
            start: The start value for the range.
        
        Returns:
            A range matrix with the same chape as self.
        '''
        var result = Self(self.rows, self.cols, fill_zeros=False)
        result.fill_range(start)
        return result
    
    @always_inline
    fn range_like(self, start: Int, step: Int) -> Self:
        '''Return a matrix with values from range(start, start + self.size * step, step)
        in row-major order.
        
        Args:
            start: The start value for the range.
            step: The step size for the range.
        
        Returns:
            A range matrix with the same chape as self.
        '''
        var result = Self(self.rows, self.cols, fill_zeros=False)
        result.fill_range(start, step)
        return result

    @always_inline
    fn zeros_like(self) -> Self:
        '''Return a matrix of zeros with the same shape as self.
        
        Returns:
            A matrix of zeros.
        '''
        return Self(rows=self.rows, cols=self.cols, fill_zeros=True)

    @always_inline
    fn ones_like(self) -> Self:
        '''Return a matrix of ones with the same shape as self.
        
        Returns:
            A matrix of ones.
        '''
        var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=False)
        result.fill_one()
        return result

    @always_inline
    fn i_like(self) -> Self:
        '''Return a matrix of i's with the same shape as self.
        
        Returns:
            A matrix of i's.
        '''
        var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=False)
        result.fill_i()
        return result

    @always_inline
    fn eye_like(self) -> Self:
        '''Return an identity matrix with the same shape as self.
        
        Returns:
            An identity matrix.
        '''
        var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=True)
        result.fill_diag(1)
        return result

    # Comparison operators ############

    fn matrix_equals(self, other: Self) -> Bool:
        '''Check two matrices for equality.

        Args:
            other: The matrix to compare with self.
        
        Returns:
            True if self is the same shape as other and all elements are equal, False otherwise.
        '''
        if self.rows != other.rows or self.cols != other.cols:
            return False
        for idx in range(self.size):
            if self.load_idx[1](idx) != other.load_idx[1](idx):
                return False
        return True

    fn is_close[tol: Scalar[Self.type] = DEFAULT_TOL](self, other: Self) -> Bool:
        '''Check if all elements of two matrices are within tol of eachother.

        Parameters:
            tol: The tolerance for checking closeness.
        
        Args:
            other: The matrix to compare with self.
        
        Returns:
            True if self is the same shape as other and all elements are within tol of eachother,
            False otherwise.
        '''
        if self.rows != other.rows or self.cols != other.cols:
            return False
        for idx in range(self.size):
            if not self.load_idx[1](idx).is_close[tol](other.load_idx[1](idx)):
                return False
        return True

    fn is_close[tol: Scalar[Self.type] = DEFAULT_TOL](self, other: ComplexScalar[Self.type]) -> Bool:
        '''Check if all elements of the matrix are within tol of other.

        Parameters:
            tol: The tolerance for checking closeness.
        
        Args:
            other: The number to compare with self.
        
        Returns:
            True if all elements of self are within tol of other, False otherwise.
        '''
        for idx in range(self.size):
            if not self.load_idx[1](idx).is_close[tol](other):
                return False
        return True

    fn __eq__(self, other: Self) raises -> Self:
        '''Find positions where two matrices are equal.
        
        Args:
            other: The matrix to compare with self.
        
        Returns:
            A matrix with ones in the positions in which self is equal to other and zeros
            elsewhere.
        '''
        self._assert_same_shape(other)
        alias one = SIMD[Self.type, 1](1)
        var result: Self = Self(self.rows, self.cols, fill_zeros=True)
        if self._should_parallelize():
            @parameter
            fn row_eq(r: Int):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) == other.load_crd[1](r, c):
                        result.re.store(r * self.cols + c, one)
            parallelize[row_eq](self.rows)
        else:
            for r in range(self.rows):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) == other.load_crd[1](r, c):
                        result.re.store(r * self.cols + c, one)
        return result

    fn __eq__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Find positions where self is equal to other.
        
        Args:
            other: The number to compare with self.
        
        Returns:
            A matrix with ones in the positions in which self is equal to other and zeros
            elsewhere.
        '''
        alias one = SIMD[Self.type, 1](1)
        var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=True)
        if self._should_parallelize():
            @parameter
            fn row_eq(r: Int):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) == other:
                        result.re.store(r * self.cols + c, one)
            parallelize[row_eq](self.rows)
        else:
            for r in range(self.rows):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) == other:
                        result.re.store(r * self.cols + c, one)
        return result

    fn __ne__(self, other: Self) raises -> Self:
        '''Find positions where two matrices are not equal.
        
        Args:
            other: The matrix to compare with self.
        
        Returns
            A matrix with ones in the positions in which self is not equal to other and zeros
            elsewhere.
        '''
        self._assert_same_shape(other)
        alias one = SIMD[Self.type, 1](1)
        var result: Self = Self(self.rows, self.cols, fill_zeros=True)
        if self._should_parallelize():
            @parameter
            fn row_eq(r: Int):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) != other.load_crd[1](r, c):
                        result.re.store(r * self.cols + c, one)
            parallelize[row_eq](self.rows)
        else:
            for r in range(self.rows):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) != other.load_crd[1](r, c):
                        result.re.store(r * self.cols + c, one)
        return result

    fn __ne__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Find positions where self is not equal to other.
        
        Args:
            other: The number to compare with self.
        
        Returns:
            A matrix with ones in the positions in which self is not equal to other and zeros
            elsewhere.
        '''
        alias one = SIMD[Self.type, 1](1)
        var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=True)
        if self._should_parallelize():
            @parameter
            fn row_eq(r: Int):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) != other:
                        result.re.store(r * self.cols + c, one)
            parallelize[row_eq](self.rows)
        else:
            for r in range(self.rows):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) != other:
                        result.re.store(r * self.cols + c, one)
        return result

    fn __gt__(self, other: Self) raises -> Self:
        '''Find positions where self is greater than other.
        
        Args:
            other: The matrix to compare with self.
        
        Returns:
            A matrix with ones in the positions in which self is greater than other and zeros
            elsewhere.
        '''
        self._assert_same_shape(other)
        alias one = SIMD[Self.type, 1](1)
        var result: Self = Self(self.rows, self.cols, fill_zeros=True)
        if self._should_parallelize():
            @parameter
            fn row_eq(r: Int):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) > other.load_crd[1](r, c):
                        result.re.store(r * self.cols + c, one)
            parallelize[row_eq](self.rows)
        else:
            for r in range(self.rows):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) > other.load_crd[1](r, c):
                        result.re.store(r * self.cols + c, one)
        return result

    fn __gt__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Find positions where self is greater than other.
        
        Args:
            other: The number to compare with self.
        
        Returns:
            A matrix with ones in the positions in which self is greater than other and zeros
            elsewhere.
        '''
        alias one = SIMD[Self.type, 1](1)
        var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=True)
        if self._should_parallelize():
            @parameter
            fn row_eq(r: Int):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) > other:
                        result.re.store(r * self.cols + c, one)
            parallelize[row_eq](self.rows)
        else:
            for r in range(self.rows):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) > other:
                        result.re.store(r * self.cols + c, one)
        return result

    fn __ge__(self, other: Self) raises -> Self:
        '''Find positions where self is greater than or equal to other.
        
        Args:
            other: The matrix to compare with self.
        
        Returns:
            A matrix with ones in the positions in which self is greater than or equal to other
            and zeros elsewhere.
        '''
        self._assert_same_shape(other)
        alias one = SIMD[Self.type, 1](1)
        var result: Self = Self(self.rows, self.cols, fill_zeros=True)
        if self._should_parallelize():
            @parameter
            fn row_eq(r: Int):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) >= other.load_crd[1](r, c):
                        result.re.store(r * self.cols + c, one)
            parallelize[row_eq](self.rows)
        else:
            for r in range(self.rows):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) >= other.load_crd[1](r, c):
                        result.re.store(r * self.cols + c, one)
        return result

    fn __ge__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Find positions where self is greater than or equal to other.
        
        Args:
            other: The number to compare with self.
        
        Returns:
            A matrix with ones in the positions in which self is greater than or equal to other
            and zeros elsewhere.
        '''
        alias one = SIMD[Self.type, 1](1)
        var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=True)
        if self._should_parallelize():
            @parameter
            fn row_eq(r: Int):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) >= other:
                        result.re.store(r * self.cols + c, one)
            parallelize[row_eq](self.rows)
        else:
            for r in range(self.rows):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) >= other:
                        result.re.store(r * self.cols + c, one)
        return result

    fn __lt__(self, other: Self) raises -> Self:
        '''Find positions where self is less than other.
        
        Args:
            other: The matrix to compare with self.
        
        Returns:
            A matrix with ones in the positions in which self is less than other and zeros
            elsewhere.
        '''
        self._assert_same_shape(other)
        alias one = SIMD[Self.type, 1](1)
        var result: Self = Self(self.rows, self.cols, fill_zeros=True)
        if self._should_parallelize():
            @parameter
            fn row_eq(r: Int):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) < other.load_crd[1](r, c):
                        result.re.store(r * self.cols + c, one)
            parallelize[row_eq](self.rows)
        else:
            for r in range(self.rows):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) < other.load_crd[1](r, c):
                        result.re.store(r * self.cols + c, one)
        return result

    fn __lt__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Find positions where self is less than other.
        
        Args:
            other: The number to compare with self.
        
        Returns:
            A matrix with ones in the positions in which self is less than other and zeros
            elsewhere.
        '''
        alias one = SIMD[Self.type, 1](1)
        var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=True)
        if self._should_parallelize():
            @parameter
            fn row_eq(r: Int):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) < other:
                        result.re.store(r * self.cols + c, one)
            parallelize[row_eq](self.rows)
        else:
            for r in range(self.rows):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) < other:
                        result.re.store(r * self.cols + c, one)
        return result

    fn __le__(self, other: Self) raises -> Self:
        '''Find positions where self is less than or equal to other.
        
        Args:
            other: The matrix to compare with self.
        
        Returns:
            A matrix with ones in the positions in which self is less than or equal to other
            and zeros elsewhere.
        '''
        self._assert_same_shape(other)
        alias one = SIMD[Self.type, 1](1)
        var result: Self = Self(self.rows, self.cols, fill_zeros=True)
        if self._should_parallelize():
            @parameter
            fn row_eq(r: Int):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) <= other.load_crd[1](r, c):
                        result.re.store(r * self.cols + c, one)
            parallelize[row_eq](self.rows)
        else:
            for r in range(self.rows):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) <= other.load_crd[1](r, c):
                        result.re.store(r * self.cols + c, one)
        return result

    fn __le__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Find positions where self is less than or equal to other.
        
        Args:
            other: The number to compare with self.
        
        Returns:
            A matrix with ones in the positions in which self is less than or equal to other
            and zeros elsewhere.
        '''
        alias one = SIMD[Self.type, 1](1)
        var result: Self = Self(rows=self.rows, cols=self.cols, fill_zeros=True)
        if self._should_parallelize():
            @parameter
            fn row_eq(r: Int):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) <= other:
                        result.re.store(r * self.cols + c, one)
            parallelize[row_eq](self.rows)
        else:
            for r in range(self.rows):
                for c in range(self.cols):
                    if self.load_crd[1](r, c) <= other:
                        result.re.store(r * self.cols + c, one)
        return result

    # Conversion ######################

    fn to_list(self) -> List[ComplexScalar[Self.type], True]:
        '''Flatten the matrix and convert it to a list of complex values.

        Returns:
            A list of matrix values.
        '''
        var result = List[ComplexScalar[Self.type], True](capacity=self.size)
        for idx in range(self.size):
            result.append(self.load_idx[1](idx))
        return result
    
    fn to_2d_list(self) -> List[List[ComplexScalar[Self.type], True]]:
        '''Return the matrix as a list of lists of complex values.

        Returns:
            A list of lists of matrix values.
        '''
        var result = List[List[ComplexScalar[Self.type], True]](capacity=self.rows)
        for r in range(self.rows):
            var row = List[ComplexScalar[Self.type], True](capacity=self.cols)
            for c in range(self.cols):
                row.append(self.load_crd[1](r, c))
            result.append(row)
        return result
