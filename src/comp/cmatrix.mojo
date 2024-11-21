from memory import memset_zero, memcpy
from algorithm import parallelize, vectorize
from sys import simdwidthof

from .complexsimd import ComplexScalar, ComplexSIMD

@value
struct CMatrix[type: DType, rows: Int, cols: Int](
    Absable,
    Formattable, 
    Sized, 
    RepresentableCollectionElement,
    StringableCollectionElement,
    Defaultable,
):    
    alias size: Int = Self.rows * Self.cols
    alias _is_col_dominant: Bool = Self.cols >= Self.rows

    var re: UnsafePointer[Scalar[Self.type]]
    var im: UnsafePointer[Scalar[Self.type]]

    # Initialization ##################

    @always_inline
    fn __init__(inout self):
        self.__init__(fill_zeros=True)

    @always_inline
    fn __init__(inout self, fill_zeros: Bool = True):
        self.re = UnsafePointer[Scalar[Self.type]].alloc(Self.size)
        self.im = UnsafePointer[Scalar[Self.type]].alloc(Self.size)
        if fill_zeros:
            memset_zero(self.re.address, Self.size)
            memset_zero(self.im.address, Self.size)
    
    @always_inline
    fn __init__(inout self, *data: ComplexScalar[Self.type]):
        self.re = UnsafePointer[Scalar[Self.type]].alloc(Self.size)
        self.im = UnsafePointer[Scalar[Self.type]].alloc(Self.size)
        for idx in range(len(data)):
            self[idx] = data[idx]
    
    @always_inline
    fn __init__(inout self, data: List[ComplexScalar[Self.type], True]):
        self.re = UnsafePointer[Scalar[Self.type]].alloc(Self.size)
        self.im = UnsafePointer[Scalar[Self.type]].alloc(Self.size)
        for idx in range(len(data)):
            self[idx] = data[idx]

    @always_inline
    fn __copyinit__(inout self, existing: Self):
        self.re = UnsafePointer[Scalar[Self.type]].alloc(Self.size)
        self.im = UnsafePointer[Scalar[Self.type]].alloc(Self.size)
        memcpy(self.re, existing.re, Self.size)
        memcpy(self.im, existing.im, Self.size)
    
    @always_inline
    fn __moveinit__(inout self, owned existing: Self):
        self.re = existing.re
        self.im = existing.im
    
    @always_inline
    fn __del__(owned self):
        for idx in range(Self.size):
            (self.re + idx).destroy_pointee()
            (self.im + idx).destroy_pointee()
        self.re.free()
        self.im.free()
    
    # Properties ######################

    @always_inline
    fn __len__(self) -> Int:
        return Self.size
    
    @always_inline
    fn shape(self) -> Tuple[Int, Int]:
        return (Self.rows, Self.cols)

    # String conversion ###############

    @no_inline
    fn __str__(self) -> String:
        return self._get_str_rep(max_lines=6)
    
    @no_inline
    fn _get_str_rep(self, max_lines: Int) -> String:
        alias dots: String = '...'
        if Self.size == 0:
            return '[]'
        var mat = List[List[String, True], True]()
        var hide_rows: Bool = Self.rows > max_lines and max_lines > 0
        var hide_cols: Bool = Self.cols > max_lines and max_lines > 0
        for r in range(max_lines if hide_rows else Self.rows):
            if r == max_lines // 2 and hide_rows:
                mat.append(List[String, True](dots))
            mat.append(List[String, True]())
            var r_idx: Int = r if r < max_lines // 2 or not hide_rows else Self.rows + r - max_lines
            for c in range(max_lines if hide_cols else Self.cols):
                if c == max_lines // 2 and hide_cols:
                    mat[-1].append(dots)
                var c_idx: Int = c if c < max_lines // 2 or not hide_cols else Self.cols + c - max_lines
                mat[-1].append(str(self[r_idx, c_idx]))
        
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
        out += ' DType: ' + str(Self.type) + ' | Size: (' + str(Self.rows) + ', ' + str(Self.cols) + ')'
        return out

    @no_inline
    fn __repr__(self) -> String:
        return self._get_str_rep(max_lines=0)
    
    @no_inline
    fn format_to(self, inout writer: Formatter):
        writer.write(str(self))
    
    # Item access #####################

    @always_inline
    fn __getitem__(self, row: Int, col: Int) -> ComplexScalar[Self.type]:
        '''Get the value at the specified row and col.'''
        return self.load_idx[1](row * Self.cols + col)
    
    @always_inline
    fn __getitem__(self, idx: Int) -> ComplexScalar[Self.type]:
        '''Get the value at the specified index. 
        Matrices are indexed from index 0 at postion (0, 0) to index size - 1 
        at position (rows - 1, cols - 1) in row-major order.
        '''
        return self.load_idx[1](idx)
    
    @always_inline
    fn __setitem__(self, row: Int, col: Int, val: ComplexScalar[Self.type]):
        '''Set the value at the specified row and col.'''
        self.store_crd[1](row, col, val)
    
    @always_inline
    fn __setitem__(self, idx: Int, val: ComplexScalar[Self.type]):
        '''Set the value at the specified index. 
        Matrices are indexed from index 0 at postion (0, 0) to index size - 1 
        at position (rows - 1, cols - 1) in row-major order.
        '''
        self.store_idx[1](idx, val)

    @always_inline
    fn load_idx[width: Int](self, idx: Int) -> ComplexSIMD[Self.type, width]:
        '''Load width lanes of the ComplexSIMD at index idx.'''
        return ComplexSIMD(
            self.re.load[width=width](idx), 
            self.im.load[width=width](idx),
        )

    @always_inline
    fn store_idx[width: Int](self, idx: Int, val: ComplexSIMD[Self.type, width]):
        '''Store val in the ComplexSIMD at index idx.'''
        self.re.store[width=width](idx, val.re)
        self.im.store[width=width](idx, val.im)

    @always_inline
    fn load_crd[width: Int](self, row: Int, col: Int) -> ComplexSIMD[Self.type, width]:
        '''Load width lanes of the ComplexSIMD at position (row, col).'''
        return self.load_idx[width](row * Self.cols + col)

    @always_inline
    fn store_crd[width: Int](self, row: Int, col: Int, val: ComplexSIMD[Self.type, width]):
        '''Store val in the ComplexSIMD at position (row, col).'''
        self.store_idx[width](row * Self.cols + col, val)
    
    @always_inline
    fn strided_load_idx[width: Int](self, idx: Int, stride: Int) -> ComplexSIMD[Self.type, width]:
        '''Load width lanes starting at index idx spaced by stride.'''
        return ComplexSIMD(
            (self.re + idx).strided_load[width=width](stride),
            (self.im + idx).strided_load[width=width](stride),
        )
    
    @always_inline
    fn strided_store_idx[width: Int](self, idx: Int, stride: Int, val: ComplexSIMD[Self.type, width]):
        '''Store width lanes starting at index idx spaced by stride.'''
        (self.re + idx).strided_store[width=width](val.re, stride)
        (self.im + idx).strided_store[width=width](val.im, stride)

    # Helper functions ################

    @always_inline
    fn _parallelize_vectorize_op[
        func: fn[width: Int](r: Int, c: Int) capturing -> ComplexSIMD[Self.type, width]
    ](self) -> Self:
        '''Returns a CMatrix in which position (r, c) holds the value func(r, c).'''
        var result: Self = Self(fill_zeros=False)
        if Self._is_col_dominant:
            @parameter
            fn op_row(r: Int):
                @parameter
                fn op_col[simd_width: Int](c: Int):
                    result.store_crd[simd_width](r, c, func[simd_width](r, c))
                vectorize[op_col, simdwidthof[Self.type]()](Self.cols)
            parallelize[op_row](Self.rows, Self.rows)
            return result
        else:
            @parameter
            fn op_col(c: Int):
                @parameter
                fn op_row[simd_width: Int](r: Int):
                    result.strided_store_idx[simd_width](r * Self.cols + c, Self.cols, func[simd_width](r, c))
                vectorize[op_row, simdwidthof[Self.type]()](Self.rows)
            parallelize[op_col](Self.cols, Self.cols)
            return result
    
    @always_inline
    fn _parallelize_vectorize_op_inplace[
        func: fn[width: Int](r: Int, c: Int) capturing -> ComplexSIMD[Self.type, width]
    ](self):  
        '''Stores the value func(r, c) at each position (r, c) in self.'''
        if Self._is_col_dominant:
            @parameter
            fn op_row(r: Int):
                @parameter
                fn op_col[simd_width: Int](c: Int):
                    self.store_crd[simd_width](r, c, func[simd_width](r, c))
                vectorize[op_col, simdwidthof[Self.type]()](Self.cols)
            parallelize[op_row](Self.rows, Self.rows)
        else:
            @parameter
            fn op_col(c: Int):
                @parameter
                fn op_row[simd_width: Int](r: Int):
                    self.strided_store_idx[simd_width](r * Self.cols + c, Self.cols, func[simd_width](r, c))
                vectorize[op_row, simdwidthof[Self.type]()](Self.rows)
            parallelize[op_col](Self.cols, Self.cols)
    
    # Math dunders ####################

    @always_inline
    fn __neg__(self) -> Self:
        '''Defines the `-` unary negation operator. Returns -self.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn neg_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return -self.load_crd[simd_width](r, c)
            return self._parallelize_vectorize_op[neg_r]()
        else:
            @parameter
            fn neg_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return -self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols)
            return self._parallelize_vectorize_op[neg_c]()
    
    @always_inline
    fn __pos__(self) -> Self:
        '''Defines the `+` unary positive operator. Returns self.'''
        return self
    
    @always_inline
    fn __add__(self, other: Self) -> Self:
        '''Defines the `+` add operator. Returns self + other.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn add_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) + other.load_crd[simd_width](r, c)
            return self._parallelize_vectorize_op[add_r]()
        else:
            @parameter
            fn add_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                var idx: Int = r * Self.cols + c
                return (
                    self.strided_load_idx[simd_width](idx, Self.cols)
                    + other.strided_load_idx[simd_width](idx, Self.cols)
                )
            return self._parallelize_vectorize_op[add_c]()
    
    @always_inline
    fn __add__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `+` add operator. Returns self + other.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn add_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) + other
            return self._parallelize_vectorize_op[add_r]()
        else:
            @parameter
            fn add_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols) + other
            return self._parallelize_vectorize_op[add_c]()
    
    @always_inline
    fn __sub__(self, other: Self) -> Self:
        '''Defines the `-` subtraction operator. Returns self - other.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn sub_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) - other.load_crd[simd_width](r, c)
            return self._parallelize_vectorize_op[sub_r]()
        else:
            @parameter
            fn sub_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                var idx: Int = r * Self.cols + c
                return (
                    self.strided_load_idx[simd_width](idx, Self.cols)
                    - other.strided_load_idx[simd_width](idx, Self.cols)
                )
            return self._parallelize_vectorize_op[sub_c]()
    
    @always_inline
    fn __sub__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `-` subtraction operator. Returns self - other.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn sub_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) - other
            return self._parallelize_vectorize_op[sub_r]()
        else:
            @parameter
            fn sub_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols) - other
            return self._parallelize_vectorize_op[sub_c]()
    
    @always_inline
    fn __mul__(self, other: Self) -> Self:
        '''Defines the `*` product operator. Returns self * other.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn mul_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) * other.load_crd[simd_width](r, c)
            return self._parallelize_vectorize_op[mul_r]()
        else:
            @parameter
            fn mul_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                var idx: Int = r * Self.cols + c
                return (
                    self.strided_load_idx[simd_width](idx, Self.cols)
                    * other.strided_load_idx[simd_width](idx, Self.cols)
                )
            return self._parallelize_vectorize_op[mul_c]()
    
    @always_inline
    fn __mul__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `*` product operator. Returns self * other.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn mul_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                # print(r, c, self.load_crd[simd_width](r, c), other, self.load_crd[simd_width](r, c) * other)
                return self.load_crd[simd_width](r, c) * other
            return self._parallelize_vectorize_op[mul_r]()
        else:
            @parameter
            fn mul_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols) * other
            return self._parallelize_vectorize_op[mul_c]()

    @always_inline
    fn __truediv__(self, other: Self) -> Self:
        '''Defines the `/` divide operator. Returns self / other.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn div_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) / other.load_crd[simd_width](r, c)
            return self._parallelize_vectorize_op[div_r]()
        else:
            @parameter
            fn div_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                var idx: Int = r * Self.cols + c
                return (
                    self.strided_load_idx[simd_width](idx, Self.cols)
                    / other.strided_load_idx[simd_width](idx, Self.cols)
                )
            return self._parallelize_vectorize_op[div_c]()

    @always_inline
    fn __truediv__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `/` divide operator. Returns self / other.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn div_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) / other
            return self._parallelize_vectorize_op[div_r]()
        else:
            @parameter
            fn div_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols) / other
            return self._parallelize_vectorize_op[div_c]()
    
    @always_inline
    fn __floordiv__(self, other: Self) -> Self:
        '''Defines the `//` floor divide operator. Returns self // other.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn fdiv_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) // other.load_crd[simd_width](r, c)
            return self._parallelize_vectorize_op[fdiv_r]()
        else:
            @parameter
            fn fdiv_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                var idx: Int = r * Self.cols + c
                return (
                    self.strided_load_idx[simd_width](idx, Self.cols)
                    // other.strided_load_idx[simd_width](idx, Self.cols)
                )
            return self._parallelize_vectorize_op[fdiv_c]()

    @always_inline
    fn __floordiv__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `//` floor divide operator. Returns self // other.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn fdiv_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) // other
            return self._parallelize_vectorize_op[fdiv_r]()
        else:
            @parameter
            fn fdiv_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols) // other
            return self._parallelize_vectorize_op[fdiv_c]()
    
    @always_inline
    fn __mod__(self, other: Self) -> Self:
        '''Defines the `%` mod operator. Returns self % other.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn mod_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) % other.load_crd[simd_width](r, c)
            return self._parallelize_vectorize_op[mod_r]()
        else:
            @parameter
            fn mod_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                var idx: Int = r * Self.cols + c
                return (
                    self.strided_load_idx[simd_width](idx, Self.cols)
                    % other.strided_load_idx[simd_width](idx, Self.cols)
                )
            return self._parallelize_vectorize_op[mod_c]()
    
    @always_inline
    fn __mod__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `%` mod operator. Returns self % other.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn mod_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) % other
            return self._parallelize_vectorize_op[mod_r]()
        else:
            @parameter
            fn mod_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols) % other
            return self._parallelize_vectorize_op[mod_c]()
    
    @always_inline
    fn __divmod__(self, other: Self) -> Tuple[Self, Self]:
        '''Defines the divmod operator. Returns (self // other, self % other).'''
        return (self // other, self % other)
    
    @always_inline
    fn __divmod__(self, other: ComplexScalar[Self.type]) -> Tuple[Self, Self]:
        '''Defines the divmod operator. Returns (self // other, self % other).'''
        return (self // other, self % other)
    
    @always_inline
    fn __matmul__(
        self, other: CMatrix[Self.type, rows=Self.cols],
    ) -> CMatrix[Self.type, Self.rows, other.cols]:
        '''Defines the `@` matrix multiple operator. Returns self @ other.'''
        var result = CMatrix[Self.type, Self.rows, other.cols](fill_zeros=True)
        @parameter
        fn calc_row(r: Int):
            for k in range(Self.cols):
                @parameter
                fn dot[simd_width: Int](c: Int):
                    result.store_crd[simd_width](
                        r, 
                        c, 
                        result.load_crd[simd_width](r, c) + self[r, k] * other.load_crd[simd_width](k, c),
                    )
                vectorize[dot, simdwidthof[Self.type]()](result.cols)
        parallelize[calc_row](result.rows, result.rows)
        return result
    
    # In-place math dunders ###########

    @always_inline
    fn __iadd__(self, other: Self):
        '''Defines the `+=` in-place add operator.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn add_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) + other.load_crd[simd_width](r, c)
            self._parallelize_vectorize_op_inplace[add_r]()
        else:
            @parameter
            fn add_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                var idx: Int = r * Self.cols + c
                return (
                    self.strided_load_idx[simd_width](idx, Self.cols) 
                    + other.strided_load_idx[simd_width](idx, Self.cols)
                )
            self._parallelize_vectorize_op_inplace[add_c]()
    
    @always_inline
    fn __iadd__(self, other: ComplexScalar[Self.type]):
        '''Defines the `+=` in-place add operator.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn add_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) + other
            self._parallelize_vectorize_op_inplace[add_r]()
        else:
            @parameter
            fn add_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols) + other
            self._parallelize_vectorize_op_inplace[add_c]()

    @always_inline
    fn __isub__(self, other: Self):
        '''Defines the `-=` in-place subtraction operator.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn sub_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) - other.load_crd[simd_width](r, c)
            self._parallelize_vectorize_op_inplace[sub_r]()
        else:
            @parameter
            fn sub_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                var idx: Int = r * Self.cols + c
                return (
                    self.strided_load_idx[simd_width](idx, Self.cols) 
                    - other.strided_load_idx[simd_width](idx, Self.cols)
                )
            self._parallelize_vectorize_op_inplace[sub_c]()
    
    @always_inline
    fn __isub__(self, other: ComplexScalar[Self.type]):
        '''Defines the `-=` in-place subtraction operator.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn sub_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) - other
            self._parallelize_vectorize_op_inplace[sub_r]()
        else:
            @parameter
            fn sub_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols) - other
            self._parallelize_vectorize_op_inplace[sub_c]()

    @always_inline
    fn __imul__(self, other: Self):
        '''Defines the `*=` in-place product operator.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn mul_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) * other.load_crd[simd_width](r, c)
            self._parallelize_vectorize_op_inplace[mul_r]()
        else:
            @parameter
            fn mul_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                var idx: Int = r * Self.cols + c
                return (
                    self.strided_load_idx[simd_width](idx, Self.cols) 
                    * other.strided_load_idx[simd_width](idx, Self.cols)
                )
            self._parallelize_vectorize_op_inplace[mul_c]()
    
    @always_inline
    fn __imul__(self, other: ComplexScalar[Self.type]):
        '''Defines the `*=` in-place product operator.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn mul_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) * other
            self._parallelize_vectorize_op_inplace[mul_r]()
        else:
            @parameter
            fn mul_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols) * other
            self._parallelize_vectorize_op_inplace[mul_c]()

    @always_inline
    fn __itruediv__(self, other: Self):
        '''Defines the `/=` in-place divide operator.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn div_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) / other.load_crd[simd_width](r, c)
            self._parallelize_vectorize_op_inplace[div_r]()
        else:
            @parameter
            fn div_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                var idx: Int = r * Self.cols + c
                return (
                    self.strided_load_idx[simd_width](idx, Self.cols) 
                    / other.strided_load_idx[simd_width](idx, Self.cols)
                )
            self._parallelize_vectorize_op_inplace[div_c]()
    
    @always_inline
    fn __itruediv__(self, other: ComplexScalar[Self.type]):
        '''Defines the `/=` in-place divide operator.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn div_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) / other
            self._parallelize_vectorize_op_inplace[div_r]()
        else:
            @parameter
            fn div_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols) / other
            self._parallelize_vectorize_op_inplace[div_c]()

    @always_inline
    fn __ifloordiv__(self, other: Self):
        '''Defines the `//=` in-place floor divide operator.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn fdiv_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) // other.load_crd[simd_width](r, c)
            self._parallelize_vectorize_op_inplace[fdiv_r]()
        else:
            @parameter
            fn fdiv_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                var idx: Int = r * Self.cols + c
                return (
                    self.strided_load_idx[simd_width](idx, Self.cols) 
                    // other.strided_load_idx[simd_width](idx, Self.cols)
                )
            self._parallelize_vectorize_op_inplace[fdiv_c]()
    
    @always_inline
    fn __ifloordiv__(self, other: ComplexScalar[Self.type]):
        '''Defines the `//=` in-place floor divide operator.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn fdiv_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) // other
            self._parallelize_vectorize_op_inplace[fdiv_r]()
        else:
            @parameter
            fn fdiv_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols) // other
            self._parallelize_vectorize_op_inplace[fdiv_c]()

    @always_inline
    fn __imod__(self, other: Self):
        '''Defines the `%=` in-place mod operator.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn mod_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) % other.load_crd[simd_width](r, c)
            self._parallelize_vectorize_op_inplace[mod_r]()
        else:
            @parameter
            fn mod_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                var idx: Int = r * Self.cols + c
                return (
                    self.strided_load_idx[simd_width](idx, Self.cols) 
                    % other.strided_load_idx[simd_width](idx, Self.cols)
                )
            self._parallelize_vectorize_op_inplace[mod_c]()

    @always_inline
    fn __imod__(self, other: ComplexScalar[Self.type]):
        '''Defines the `%=` in-place mod operator.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn mod_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c) % other
            self._parallelize_vectorize_op_inplace[mod_r]()
        else:
            @parameter
            fn mod_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols) % other
            self._parallelize_vectorize_op_inplace[mod_c]()
    
    # TODO: Make it better
    @always_inline
    fn __imatmul__(self, other: Self) raises:
        '''Defines the `@=` in-place matrix multiply operator. Only valid for square matrices.'''
        if Self.rows != Self.cols:
            raise Error('In-place matrix multiplication is only available for square matrices')
        var result = Self(fill_zeros=True)
        @parameter
        fn calc_row(r: Int):
            for k in range(Self.cols):
                @parameter
                fn dot[simd_width: Int](c: Int):
                    result.store_crd[simd_width](
                        r, 
                        c, 
                        result.load_crd[simd_width](r, c) + self[r, k] * other.load_crd[simd_width](k, c),
                    )
                vectorize[dot, simdwidthof[Self.type]()](result.cols)
        parallelize[calc_row](result.rows, result.rows)
        # # Produces segmentation fault...
        # @parameter
        # fn cpy[simd_width: Int](idx: Int):
        #     self.store_idx[simd_width](idx, result.load_idx[simd_width](idx))
        # vectorize[cpy, simdwidthof[Self.type]()](Self.size)

        # # Copies memory incorrectly...
        # memcpy(self.re, result.re, Self.size)
        # memcpy(self.im, result.im, Self.size)

        # So for now, we do this...
        for idx in range(Self.size):
            self[idx] = result[idx]
    
    # Right math dunders ##############

    @always_inline
    fn __radd__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `+` right add operator. Returns self + other.'''
        return self + other
    
    @always_inline
    fn __rsub__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `-` right subtraction operator. Returns other - self.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn sub_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return other - self.load_crd[simd_width](r, c)
            return self._parallelize_vectorize_op[sub_r]()
        else:
            @parameter
            fn sub_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return other - self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols)
            return self._parallelize_vectorize_op[sub_c]()

    @always_inline
    fn __rmul__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `*` right product operator. Returns self * other.'''
        return self * other

    @always_inline
    fn __rtruediv__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `/` right divide operator. Returns other / self.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn div_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return other / self.load_crd[simd_width](r, c)
            return self._parallelize_vectorize_op[div_r]()
        else:
            @parameter
            fn div_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return other / self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols)
            return self._parallelize_vectorize_op[div_c]()

    @always_inline
    fn __rfloordiv__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `//` right floor divide operator. Returns other // self.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn fdiv_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return other // self.load_crd[simd_width](r, c)
            return self._parallelize_vectorize_op[fdiv_r]()
        else:
            @parameter
            fn fdiv_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return other // self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols)
            return self._parallelize_vectorize_op[fdiv_c]()

    @always_inline
    fn __rmod__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Defines the `%` right mod operator. Returns other % self.'''
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn mod_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return other % self.load_crd[simd_width](r, c)
            return self._parallelize_vectorize_op[mod_r]()
        else:
            @parameter
            fn mod_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return other % self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols)
            return self._parallelize_vectorize_op[mod_c]()

    @always_inline
    fn __rdivmod__(self, other: ComplexScalar[Self.type]) -> Tuple[Self, Self]:
        '''Defines the right divmod operator. Returns (other // self, other % self).'''
        return (other // self, other % self)

    # Other math ######################

    @always_inline
    fn __abs__(self) -> Self:
        @parameter
        if Self._is_col_dominant:
            @parameter
            fn abs_r[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.load_crd[simd_width](r, c).__abs__()
            return self._parallelize_vectorize_op[abs_r]()
        else:
            @parameter
            fn abs_c[simd_width: Int](r: Int, c: Int) -> ComplexSIMD[Self.type, simd_width]:
                return self.strided_load_idx[simd_width](r * Self.cols + c, Self.cols).__abs__()
            return self._parallelize_vectorize_op[abs_c]()
    
    @always_inline
    fn conj(self) -> Self:
        '''Return the conjugate of the matrix.'''
        var result: Self = self
        @parameter
        fn conj_simd[simd_width: Int](idx: Int):
            result.im.store[width=simd_width](idx, -self.im.load[width=simd_width](idx))
        vectorize[conj_simd, simdwidthof[Self.type]()](Self.size)
        return result
    
    @always_inline
    fn iconj(self):
        '''Conjugate the matrix in-place.'''
        @parameter
        fn conj_simd[simd_width: Int](idx: Int):
            self.im.store[width=simd_width](idx, -self.im.load[width=simd_width](idx))
        vectorize[conj_simd, simdwidthof[Self.type]()](Self.size)
    
    @always_inline
    fn dag(self) -> CMatrix[Self.type, Self.cols, Self.rows]:
        '''Return the conjugate-transpose of the matrix.'''
        var result = CMatrix[Self.type, Self.cols, Self.rows](fill_zeros=False)
        @parameter
        fn transpose_row(r: Int):
            @parameter
            fn transpose_col[simd_width: Int](c: Int):
                result.strided_store_idx[simd_width](
                    c * result.cols + r, 
                    result.cols, 
                    self.load_crd[simd_width](r, c).conj(),
                )
            vectorize[transpose_col, simdwidthof[Self.type]()](Self.cols)
        parallelize[transpose_row](Self.rows, Self.rows)
        return result
    
    # # TODO: Make it better
    # # This can most definitely be done better with something like
    # # buffer.Buffer or algorithm.reduce
    @always_inline
    fn sum(self) -> ComplexScalar[Self.type]:
        if self._is_col_dominant:
            var row_sums = CMatrix[Self.type, Self.rows, 1](fill_zeros=False)
            @parameter
            fn sum_rows(r: Int):
                var row_sum = ComplexScalar[Self.type]()
                for c in range(Self.cols):
                    row_sum += self[r, c]
                row_sums[r] = row_sum
            parallelize[sum_rows](Self.rows, Self.rows)
            var total = ComplexScalar[Self.type]()
            for r in range(Self.rows):
                total += row_sums[r]
            return total 
        else:
            var col_sums = CMatrix[Self.type, 1, cols](fill_zeros=False)
            @parameter
            fn sum_cols(c: Int):
                var col_sum = ComplexScalar[Self.type]()
                for r in range(Self.rows):
                    col_sum += self[r, c]
                col_sums[c] = col_sum
            parallelize[sum_cols](Self.cols, Self.cols)
            var total = ComplexScalar[Self.type]()
            for c in range(Self.cols):
                total += col_sums[c]
            return total 

    # TODO: Implement determinant method
    
    @always_inline
    fn inv(self, tol: SIMD[Self.type, 1] = 1e-12) raises -> Self:
        '''Return the inverse of a square matrix.'''
        alias zero = ComplexScalar[Self.type](0)
        alias one = ComplexScalar[Self.type](1)

        if Self.rows != Self.cols:
            raise Error('Only square matrices are invertible')
        
        var augmented = CMatrix[Self.type, Self.rows, 2 * Self.cols](fill_zeros=True)
        for r in range(Self.rows):
            memcpy(augmented.re + r * augmented.cols, self.re + r * Self.cols, Self.cols)
            memcpy(augmented.im + r * augmented.cols, self.im + r * Self.cols, Self.cols)
            augmented[r, Self.cols + r] = one
        
        for i in range(Self.rows):
            if augmented[i, i].norm() < tol:
                for j in range(i + 1, Self.rows):
                    if augmented[j, i].norm() >= tol:
                        var i_start: Int = i * augmented.cols
                        var j_start: Int = j * augmented.cols
                        var row_i = CMatrix[Self.type, 1, augmented.cols](fill_zeros=False)
                        memcpy(row_i.re, augmented.re + i_start, augmented.cols)
                        memcpy(row_i.im, augmented.im + i_start, augmented.cols)
                        memcpy(augmented.re + i_start, augmented.re + j_start, augmented.cols)
                        memcpy(augmented.im + i_start, augmented.im + j_start, augmented.cols)
                        memcpy(augmented.re + j_start, row_i.re, augmented.cols)
                        memcpy(augmented.im + j_start, row_i.im, augmented.cols)
                        break
                else:
                    raise Error('Matrix is not invertible')

            var pivot: ComplexScalar[Self.type] = augmented[i, i]
            for j in range(augmented.cols):
                augmented[i, j] /= pivot

            for j in range(Self.rows):
                if j != i:
                    var factor: ComplexScalar[Self.type] = augmented[j, i]
                    for k in range(augmented.cols):
                        augmented[j, k] -= factor * augmented[i, k]
        
        var result = CMatrix[Self.type, Self.rows, Self.cols](fill_zeros=False)
        for i in range(Self.rows):
            memcpy(result.re + i * Self.cols, augmented.re + i * augmented.cols + Self.cols, Self.cols)
            memcpy(result.im + i * Self.cols, augmented.im + i * augmented.cols + Self.cols, Self.cols)
        return result

    @always_inline
    fn inverse(self) raises -> Self:
        '''Return the inverse of a square matrix. Alias of inv.'''
        return self.inv()

    # Shape operations ################

    @always_inline
    fn reshape[new_rows: Int, new_cols: Int](self) raises -> CMatrix[Self.type, new_rows, new_cols]:
        '''Return a reshaped matrix.'''
        if new_rows * new_cols != Self.size:
            raise Error('Incompatible reshape dimension')
        var result: CMatrix[Self.type, new_rows, new_cols] = CMatrix[Self.type, new_rows, new_cols](fill_zeros=False)
        memcpy(result.re, self.re, Self.size)
        memcpy(result.im, self.im, Self.size)
        return result
    
    @always_inline
    fn flatten_to_row(self) -> CMatrix[Self.type, 1, Self.size]:
        '''Return a flattened row matrix.'''
        var result: CMatrix[Self.type, 1, Self.size] = CMatrix[Self.type, 1, Self.size](fill_zeros=False)
        memcpy(result.re, self.re, Self.size)
        memcpy(result.im, self.im, Self.size)
        return result
    
    @always_inline
    fn flatten_to_column(self) -> CMatrix[Self.type, Self.size, 1]:
        '''Return a flattened column matrix.'''
        var result: CMatrix[Self.type, Self.size, 1] = CMatrix[Self.type, Self.size, 1](fill_zeros=False)
        memcpy(result.re, self.re, Self.size)
        memcpy(result.im, self.im, Self.size)
        return result
    
    @always_inline
    fn transpose(self) -> CMatrix[Self.type, Self.cols, Self.rows]:
        '''Return the transpose of the matrix.'''
        var result = CMatrix[Self.type, Self.cols, Self.rows](fill_zeros=False)
        @parameter
        fn transpose_row(r: Int):
            for c in range(Self.cols):
                result[c, r] = self[r, c]
        parallelize[transpose_row](Self.rows, Self.rows)
        return result
    
    # # Not actually faster...
    # @always_inline
    # fn transpose(self) -> CMatrix[Self.type, Self.cols, Self.rows]:
    #     '''Return the transpose of the matrix.'''
    #     var result:CMatrix[Self.type, Self.cols, Self.rows] = CMatrix[Self.type, Self.cols, Self.rows](fill_zeros=False)
    #     @parameter
    #     fn transpose_row(r: Int):
    #         @parameter
    #         fn transpose_col[simd_width: Int](c: Int):
    #             result.strided_store_idx[simd_width](
    #                 c * result.cols + r, 
    #                 result.cols, 
    #                 self.load_crd[simd_width](r, c),
    #             )
    #         vectorize[transpose_col, simdwidthof[Self.type]()](Self.cols)
    #     parallelize[transpose_row](Self.rows, Self.rows)
    #     return result

    # Fill operations #################

    @always_inline
    fn fill_zero(self):
        '''Fill a matrix with zeros in-place.'''
        memset_zero(self.re.address, Self.size)
        memset_zero(self.im.address, Self.size)
    
    @always_inline
    fn fill(self, val: ComplexScalar[Self.type]):
        '''Fill a matrix with val in-place.'''
        @parameter
        fn fill_row(r: Int):
            @parameter
            fn fill_col[simd_width: Int](c: Int):
                self.store_crd[simd_width](r, c, ComplexSIMD[Self.type, simd_width](val.re, val.im))
            vectorize[fill_col, simdwidthof[Self.type]()](Self.cols)
        parallelize[fill_row](Self.rows, Self.rows)
    
    @always_inline
    fn fill_one(self):
        '''Fill a matrix with ones in-place.'''
        self.fill(ComplexScalar[Self.type](1, 0))
    
    @always_inline
    fn fill_i(self):
        '''Fill a matrix with ones in-place.'''
        self.fill(ComplexScalar[Self.type](0, 1))

    @always_inline
    fn fill_diag(self, val: ComplexScalar[Self.type], offset: Int = 0):
        '''Fill the diagonal at index offset of a matrix with val in-place.'''
        var stride = Self.cols + 1
        var n_diag_elements: Int = (
            min(Self.rows, Self.cols - offset) 
            if offset >= 0 
            else min(Self.rows + offset, Self.cols)
        )
        var idx: Int
        @parameter
        fn fill_diag_simd_rp[simd_width: Int](p: Int):
            idx = p * stride + offset if offset >= 0 else p * stride - offset * Self.cols
            self.strided_store_idx[simd_width](
                idx, stride, ComplexSIMD[Self.type, simd_width](val.re, val.im)
            )
        vectorize[fill_diag_simd_rp, simdwidthof[Self.type]()](n_diag_elements)

    @always_inline
    fn fill_range(self):
        memset_zero(self.im.address, Self.size)
        for idx in range(Self.size):
            self.re.store[width=1](idx, SIMD[type, 1](idx))
    
    @always_inline
    fn fill_range(self, start: Int):
        memset_zero(self.im.address, Self.size)
        for idx in range(Self.size):
            self.re.store[width=1](idx, SIMD[type, 1](idx + start))
    
    @always_inline
    fn fill_range(self, start: Int, step: Int):
        memset_zero(self.im.address, Self.size)
        for idx in range(Self.size):
            self.re.store[width=1](idx, SIMD[type, 1](step * idx + start))

    @always_inline
    fn zeros_like(self) -> Self:
        '''Return a matrix of zeros with the same shape as self.'''
        return Self(fill_zeros=True)
    
    @always_inline
    fn ones_like(self) -> Self:
        '''Return a matrix of ones with the same shape as self.'''
        var result: Self = Self(fill_zeros=False)
        result.fill_one()
        return result
    
    @always_inline
    fn i_like(self) -> Self:
        '''Return a matrix of i's with the same shape as self.'''
        var result: Self = Self(fill_zeros=False)
        result.fill_i()
        return result
    
    @always_inline
    fn eye_like(self) -> Self:
        '''Return an identity matrix with the same shape as self.'''
        var result: Self = Self(fill_zeros=True)
        result.fill_diag(ComplexScalar[Self.type](1, 0))
        return result

    # Comparison operators ############

    # TODO: Make it better
    @always_inline
    fn matrix_equals(self, other: Self) -> Bool:
        '''Returns True if self is the same shape as other and all elements 
        are equal, False otherwise.
        '''
        for r in range(Self.rows):
            for c in range(Self.cols):
                if self[r, c] != other[r, c]:
                    return False
        return True
    
    # TODO: Make it better
    @always_inline
    fn is_close(self, other: Self, tol: Scalar[Self.type] = 1e-10) -> Bool:
        '''Returns True if self is the same shape as other and all elements 
        are equal, False otherwise.
        '''
        for r in range(Self.rows):
            for c in range(Self.cols):
                if not self[r, c].is_close(other[r, c], tol):
                    return False
        return True

    # TODO: Make it better
    @always_inline
    fn is_close(self, other: ComplexScalar[Self.type], tol: Scalar[Self.type] = 1e-10) -> Bool:
        '''Returns True if self is the same shape as other and all elements 
        are equal, False otherwise.
        '''
        for r in range(Self.rows):
            for c in range(Self.cols):
                if not self[r, c].is_close(other, tol):
                    return False
        return True

    @always_inline
    fn __eq__(self, other: Self) -> Self:
        '''Returns a CMatrix with ones in the positions in which self 
        is equal to other and zeros elsewhere. The returned matrix has
        the order of self.
        '''
        alias one = SIMD[Self.type, 1](1)
        var result: Self = Self(fill_zeros=True)
        @parameter
        fn row_eq(r: Int):
            for c in range(Self.cols):
                if self[r, c] == other[r, c]:
                    result.re.store[width=1](r * Self.cols + c, one)
        parallelize[row_eq](Self.rows, Self.rows)
        return result
    
    @always_inline
    fn __eq__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Returns a CMatrix with ones in the positions in which self 
        is equal to other and zeros elsewhere.
        '''
        var result: Self = Self(fill_zeros=True)
        @parameter
        fn row_eq(r: Int):
            for c in range(Self.cols):
                if self[r, c] == other:
                    result.re.store[width=1](r * Self.cols + c, SIMD[Self.type, 1](1))
        parallelize[row_eq](Self.rows, Self.rows)
        return result
    
    @always_inline
    fn __ne__(self, other: Self) -> Self:
        '''Returns a CMatrix with ones in the positions in which self 
        is not equal to other and zeros elsewhere. The returned matrix has
        the order of self.
        '''
        var result: Self = Self(fill_zeros=True)
        @parameter
        fn row_eq(r: Int):
            for c in range(Self.cols):
                if self[r, c] != other[r, c]:
                    result.re.store[width=1](r * Self.cols + c, SIMD[Self.type, 1](1))
        parallelize[row_eq](Self.rows, Self.rows)
        return result
    
    @always_inline
    fn __ne__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Returns a CMatrix with ones in the positions in which self 
        is not equal to other and zeros elsewhere.
        '''
        var result: Self = Self(fill_zeros=True)
        @parameter
        fn row_eq(r: Int):
            for c in range(Self.cols):
                if self[r, c] != other:
                    result.re.store[width=1](r * Self.cols + c, SIMD[Self.type, 1](1))
        parallelize[row_eq](Self.rows, Self.rows)
        return result
    
    @always_inline
    fn __gt__(self, other: Self) -> Self:
        '''Returns a CMatrix with ones in the positions in which self 
        is greater than other and zeros elsewhere. The returned matrix has
        the order of self.
        '''
        var result: Self = Self(fill_zeros=True)
        @parameter
        fn row_eq(r: Int):
            for c in range(Self.cols):
                if self[r, c] > other[r, c]:
                    result.re.store[width=1](r * Self.cols + c, SIMD[Self.type, 1](1))
        parallelize[row_eq](Self.rows, Self.rows)
        return result
    
    @always_inline
    fn __gt__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Returns a CMatrix with ones in the positions in which self 
        is greater than other and zeros elsewhere.
        '''
        var result: Self = Self(fill_zeros=True)
        @parameter
        fn row_eq(r: Int):
            for c in range(Self.cols):
                if self[r, c] > other:
                    result.re.store[width=1](r * Self.cols + c, SIMD[Self.type, 1](1))
        parallelize[row_eq](Self.rows, Self.rows)
        return result
    
    @always_inline
    fn __ge__(self, other: Self) -> Self:
        '''Returns a CMatrix with ones in the positions in which self 
        is greater than or equal to other and zeros elsewhere. The returned matrix has
        the order of self.
        '''
        var result: Self = Self(fill_zeros=True)
        @parameter
        fn row_eq(r: Int):
            for c in range(Self.cols):
                if self[r, c] >= other[r, c]:
                    result.re.store[width=1](r * Self.cols + c, SIMD[Self.type, 1](1))
        parallelize[row_eq](Self.rows, Self.rows)
        return result
    
    @always_inline
    fn __ge__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Returns a CMatrix with ones in the positions in which self 
        is greater than or equal to other and zeros elsewhere.
        '''
        var result: Self = Self(fill_zeros=True)
        @parameter
        fn row_eq(r: Int):
            for c in range(Self.cols):
                if self[r, c] >= other:
                    result.re.store[width=1](r * Self.cols + c, SIMD[Self.type, 1](1))
        parallelize[row_eq](Self.rows, Self.rows)
        return result
    
    @always_inline
    fn __lt__(self, other: Self) -> Self:
        '''Returns a CMatrix with ones in the positions in which self 
        is less than other and zeros elsewhere. The returned matrix has
        the order of self.
        '''
        var result: Self = Self(fill_zeros=True)
        @parameter
        fn row_eq(r: Int):
            for c in range(Self.cols):
                if self[r, c] < other[r, c]:
                    result.re.store[width=1](r * Self.cols + c, SIMD[Self.type, 1](1))
        parallelize[row_eq](Self.rows, Self.rows)
        return result
    
    @always_inline
    fn __lt__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Returns a CMatrix with ones in the positions in which self 
        is less than other and zeros elsewhere.
        '''
        var result: Self = Self(fill_zeros=True)
        @parameter
        fn row_eq(r: Int):
            for c in range(Self.cols):
                if self[r, c] < other:
                    result.re.store[width=1](r * Self.cols + c, SIMD[Self.type, 1](1))
        parallelize[row_eq](Self.rows, Self.rows)
        return result
    
    @always_inline
    fn __le__(self, other: Self) -> Self:
        '''Returns a CMatrix with ones in the positions in which self 
        is less than or equal to other and zeros elsewhere. The returned matrix has
        the order of self.
        '''
        var result: Self = Self(fill_zeros=True)
        @parameter
        fn row_eq(r: Int):
            for c in range(Self.cols):
                if self[r, c] <= other[r, c]:
                    result.re.store[width=1](r * Self.cols + c, SIMD[Self.type, 1](1))
        parallelize[row_eq](Self.rows, Self.rows)
        return result
    
    @always_inline
    fn __le__(self, other: ComplexScalar[Self.type]) -> Self:
        '''Returns a CMatrix with ones in the positions in which self 
        is less than or equal to other and zeros elsewhere.
        '''
        var result: Self = Self(fill_zeros=True)
        @parameter
        fn row_eq(r: Int):
            for c in range(Self.cols):
                if self[r, c] <= other:
                    result.re.store[width=1](r * Self.cols + c, SIMD[Self.type, 1](1))
        parallelize[row_eq](Self.rows, Self.rows)
        return result
