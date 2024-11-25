from math import inf as INF, ceil, log
from collections import Dict, List
from sys import simdwidthof
from algorithm import parallelize, vectorize
from memory import memcpy

from .complexsimd import ComplexScalar
from .sizedcmatrix import SizedCMatrix


fn kron_sequential[type: DType, a_rows: Int, a_cols: Int, b_rows: Int, b_cols: Int](
    a: SizedCMatrix[type, a_rows, a_cols],
    b: SizedCMatrix[type, b_rows, b_cols],
) -> SizedCMatrix[type, a_rows * b_rows, a_cols * b_cols]:
    '''Implements the Kronecker product using the naive algorithm.'''
    var result = SizedCMatrix[type, a_rows * b_rows, a_cols * b_cols](fill_zeros=False)
    var block: SizedCMatrix[type, b.rows, b.cols]
    for r_a in range(a_rows):
        for c_a in range(a_cols):
            block = a[r_a, c_a] * b
            for r_b in range(b_rows):
                for c_b in range(b_cols):
                    result[r_a * b_rows + r_b, c_a * b_cols + c_b] = block[r_b, c_b]
    return result


fn kron[type: DType, a_rows: Int, a_cols: Int, b_rows: Int, b_cols: Int](
    a: SizedCMatrix[type, a_rows, a_cols],
    b: SizedCMatrix[type, b_rows, b_cols],
) -> SizedCMatrix[type, a_rows * b_rows, a_cols * b_cols]:
    '''Implements the Kronecker product of a with b.'''
    @parameter
    if a._is_col_dominant:
        return _kron_par_a_cols(a, b)
    else:
        return _kron_par_a_rows(a, b)


fn _kron_par_a_rows[type: DType, a_rows: Int, a_cols: Int, b_rows: Int, b_cols: Int](
    a: SizedCMatrix[type, a_rows, a_cols],
    b: SizedCMatrix[type, b_rows, b_cols],
) -> SizedCMatrix[type, a_rows * b_rows, a_cols * b_cols]:
    '''Implements the Kronecker product of a with b parallelized over the rows or a.'''
    var result = SizedCMatrix[type, a_rows * b_rows, a_cols * b_cols](fill_zeros=False)
    @parameter
    fn par_row_a(r_a: Int):
        for c_a in range(a_cols):
            var a_elem: ComplexScalar[type] = a[r_a, c_a]
            @parameter
            fn vec_col_b[simd_width: Int](c_b: Int):
                for r_b in range(b_rows):
                    result.store_crd[simd_width](
                        r_a * b_rows + r_b, c_a * b_cols + c_b,
                        a_elem * b.load_crd[simd_width](r_b, c_b),
                    )
            vectorize[vec_col_b, simdwidthof[type]()](b_cols)
    parallelize[par_row_a](a_rows, a_rows)
    return result


fn _kron_par_a_cols[type: DType, a_rows: Int, a_cols: Int, b_rows: Int, b_cols: Int](
    a: SizedCMatrix[type, a_rows, a_cols],
    b: SizedCMatrix[type, b_rows, b_cols],
) -> SizedCMatrix[type, a_rows * b_rows, a_cols * b_cols]:
    '''Implements the Kronecker product of a with b parallelized over the columns or a.'''
    var result = SizedCMatrix[type, a_rows * b_rows, a_cols * b_cols](fill_zeros=False)
    @parameter
    fn par_col_a(c_a: Int):
        for r_a in range(a_rows):
            var a_elem: ComplexScalar[type] = a[r_a, c_a]
            @parameter
            fn vec_col_b[simd_width: Int](c_b: Int):
                for r_b in range(b_rows):
                    result.store_crd[simd_width](
                        r_a * b_rows + r_b, c_a * b_cols + c_b,
                        a_elem * b.load_crd[simd_width](r_b, c_b),
                    )
            vectorize[vec_col_b, simdwidthof[type]()](b_cols)
    parallelize[par_col_a](a_cols, a_cols)
    return result


fn swap_rows[type: DType, rows: Int, cols: Int](
    A: SizedCMatrix[type, rows, cols], 
    r1: Int, 
    r2: Int
) raises -> SizedCMatrix[type, rows, cols]:
    '''Swap rows r1 and r2 in A and return the result.'''
    if r1 == r2:
        return A
    if not 0 <= r1 < rows or not 0 <= r2 < rows:
        raise Error('Invalid row index')
    var result = A
    # # memcpy seems unreliable atm
    # memcpy(result.re + r1 * result.cols, A.re + r2 * A.cols, A.cols)
    # memcpy(result.im + r1 * result.cols, A.im + r2 * A.cols, A.cols)
    # memcpy(result.re + r2 * result.cols, A.re + r1 * A.cols, A.cols)
    # memcpy(result.im + r2 * result.cols, A.im + r1 * A.cols, A.cols)
    for c in range(cols):
        result[r1, c] = A[r2, c]
        result[r2, c] = A[r1, c]
    return result 


fn swap_rows_inplace[type: DType, rows: Int, cols: Int](
    A: SizedCMatrix[type, rows, cols], 
    r1: Int, 
    r2: Int
) raises:
    '''Swap rows r1 and r2 in A in-place.'''
    if r1 == r2:
        return
    if not 0 <= r1 < rows or not 0 <= r2 < rows:
        raise Error('Invalid row index')
    # # This doesn't work...
    # var r1_start: Int = r1 * cols
    # var r2_start: Int = r2 * cols
    # var row_r1 = SizedCMatrix[type, 1, cols](fill_zeros=False)
    # memcpy(row_r1.re, A.re + r1_start, cols)
    # memcpy(row_r1.im, A.im + r1_start, cols)
    # # print(A)  # But somehow adding this print statement makes it work??
    # memcpy(A.re + r1_start, A.re + r2_start, cols)
    # memcpy(A.im + r1_start, A.im + r2_start, cols)
    # memcpy(A.re + r2_start, row_r1.re, cols)
    # memcpy(A.im + r2_start, row_r1.im, cols)
    for c in range(cols):
        var p: ComplexScalar[type] = A[r1, c]
        A[r1, c] = A[r2, c]
        A[r2, c] = p


fn swap_cols[type: DType, rows: Int, cols: Int](
    A: SizedCMatrix[type, rows, cols], 
    c1: Int, 
    c2: Int
) raises -> SizedCMatrix[type, rows, cols]:
    '''Swap columns c1 and c2 in A and return the result.'''
    if c1 == c2:
        return A
    if not 0 <= c1 < cols or not 0 <= c2 < cols:
        raise Error('Invalid col index')
    var result = A
    for r in range(rows):
        result[r, c1] = A[r, c2]
        result[r, c2] = A[r, c1]
    return result 


fn swap_cols_inplace[type: DType, rows: Int, cols: Int](
    A: SizedCMatrix[type, rows, cols], 
    c1: Int, 
    c2: Int
) raises:
    '''Swap columns c1 and c2 in A in-place.'''
    if c1 == c2:
        return
    if not 0 <= c1 < cols or not 0 <= c2 < cols:
        raise Error('Invalid col index')
    var p: ComplexScalar[type]
    for r in range(rows):
        p = A[r, c1]
        A[r, c1] = A[r, c2]
        A[r, c2] = p


fn swap_vals[type: DType, rows: Int, cols: Int](
    A: SizedCMatrix[type, rows, cols], 
    r1: Int, 
    c1: Int, 
    r2: Int, 
    c2: Int
) raises -> SizedCMatrix[type, rows, cols]:
    '''Swap values at (r1, c1) and (r2, c2) in A and return the result.'''
    if r1 == r2 and c1 == c2:
        return A
    if not 0 <= r1 < rows or not 0 <= r2 < rows:
        raise Error('Invalid row index')
    if not 0 <= c1 < cols or not 0 <= c2 < cols:
        raise Error('Invalid col index')
    var result = A
    result[r1, c1] = A[r2, c2]
    result[r2, c2] = A[r1, c1]
    return result 


fn swap_vals_inplace[type: DType, rows: Int, cols: Int](
    inout A: SizedCMatrix[type, rows, cols], 
    r1: Int, 
    c1: Int, 
    r2: Int, 
    c2: Int
) raises:
    '''Swap values at (r1, c1) and (r2, c2) in-place.'''
    if r1 == r2 and c1 == c2:
        return
    if not 0 <= r1 < rows or not 0 <= r2 < rows:
        raise Error('Invalid row index')
    if not 0 <= c1 < cols or not 0 <= c2 < cols:
        raise Error('Invalid col index')
    var p1: ComplexScalar[type] = A[r1, c1]
    A[r1, c1] = A[r2, c2]
    A[r2, c2] = p1 


fn augmented_ref[type: DType, rows: Int, A_cols: Int, B_cols: Int, tol: Scalar[type] = 1e-15](
    A: SizedCMatrix[type, rows, A_cols], B: SizedCMatrix[type, rows, B_cols]
) raises -> SizedCMatrix[type, rows, A_cols + B_cols]:
    '''Computes the row echelon form of the augmented matrix [A|B].'''
    var Aaug: SizedCMatrix[type, rows, A_cols + B_cols] = hstack(A, B)
    var h: Int = 0  # initialize pivot row
    var k: Int = 0  # initialize pivot column

    while h < rows and k < A_cols:
        # /* Find the k-th pivot: */
        var i_max: Scalar[type] = Aaug[h, k].norm()
        var i_argmax: Int = h
        for i in range(h + 1, rows):
            var i_norm: Scalar[type] = Aaug[i, k].norm()
            if i_norm > i_max:
                i_max = i_norm
                i_argmax = i
        if Aaug[i_argmax, k].norm() < tol:
            # /* No pivot in this column, pass to next column */
            k += 1
        else:
            swap_rows_inplace(Aaug, h, i_argmax)
            # /* Do for all rows below pivot: */
            for i in range(h + 1, rows):
                var f: ComplexScalar[type] = Aaug[i, k] / Aaug[h, k]
                # /* Fill with zeros the lower part of pivot column: */
                Aaug[i, k] = ComplexScalar[type](0)
                # /* Do for all remaining elements in current row: */
                for j in range(k + 1, A_cols + B_cols):
                    Aaug[i, j] -= Aaug[h, j] * f
            # /* Increase pivot row and column */
            h += 1
            k += 1
    return Aaug


# TODO: Switch to a faster algorithm
fn solve[type: DType, size: Int, B_cols: Int, tol: Scalar[type] = 1e-15](
    A: SizedCMatrix[type, size, size], B: SizedCMatrix[type, size, B_cols]
) raises -> SizedCMatrix[type, size, B_cols]:
    '''Solves a linear system of equations Ax=B via Gaussian elimination.'''
    var aug_ref: SizedCMatrix[type, size, size + B_cols] = augmented_ref[tol=tol](A, B)

    for i in range(size):
        if aug_ref[i, i].norm() < tol:
            raise Error('Matrix is singular or nearly singular')

    # Back substitution
    var X: SizedCMatrix[type, size, B_cols] = B.zeros_like()
    for c in range(B_cols):
        for i in range(size - 1, -1, -1):
            var dot = ComplexScalar[type](0)
            for j in range(i + 1, size):
                dot += aug_ref[i, j] * X[j, c]
            X[i, c] = (aug_ref[i, c + size] - dot) / aug_ref[i, i]
    return X


fn one_norm[type: DType, rows: Int, cols: Int](A: SizedCMatrix[type, rows, cols]) -> SIMD[type, 1]:
    '''Computes the 1-norm of a matrix, i.e. the maximum absolute column sum.'''
    var col_sums: SizedCMatrix[type, 1, cols] = SizedCMatrix[type, 1, cols](fill_zeros=False)
    @parameter
    fn sum_abs_cols(c: Int):
        var col_sum = ComplexScalar[type]()
        for r in range(rows):
            col_sum += abs(A[r, c])
        col_sums[c] = col_sum
    parallelize[sum_abs_cols](cols, cols)
    var max: SIMD[type, 1] = col_sums[0].re
    for c in range(1, cols):
        if (cs := col_sums[c].re) > max:
            max = cs
    return max


fn _int_matrix_power[type: DType, size: Int](owned A: SizedCMatrix[type, size, size], n: Int) -> SizedCMatrix[type, size, size]:
    '''Computes A^n for a matrix A and positive integer n.'''
    if n == 1:
        return A 
    if n % 2 == 0:
        return _int_matrix_power(A @ A, n // 2)
    else:
        return _int_matrix_power(A @ A, (n - 1) // 2) @ A


fn matrix_power[type: DType, size: Int](A: SizedCMatrix[type, size, size], n: Int) raises -> SizedCMatrix[type, size, size]:
    '''Computes A^n for a matrix A and integer n.'''
    if n == 0:
        return A.eye_like()
    var Ac = A  # copy A
    if n < 0:
        return _int_matrix_power(Ac^, -n).inv()
    return _int_matrix_power(Ac^, n)


fn mmax[type: DType, rows: Int, cols: Int](A: SizedCMatrix[type, rows, cols]) -> ComplexScalar[type]:
    '''Returns the maximum value (by modulus) of the matrix.'''
    var mx: ComplexScalar[type] = A[0]
    for idx in range(A.size):
        if (m := A[idx]) > mx:  # ComplexSIMD's are compared by norm
            mx = m
    return mx


fn mmin[type: DType, rows: Int, cols: Int](A: SizedCMatrix[type, rows, cols]) -> ComplexScalar[type]:
    '''Returns the minimum value (by modulus) of the matrix.'''
    var mx: ComplexScalar[type] = A[0]
    for idx in range(A.size):
        if (m := A[idx]) < mx:  # ComplexSIMD's are compared by norm
            mx = m
    return mx


fn hstack[type: DType, rows: Int, A_cols: Int, B_cols: Int](
    A: SizedCMatrix[type, rows, A_cols], B: SizedCMatrix[type, rows, B_cols]
) -> SizedCMatrix[type, rows, A_cols + B_cols]:
    '''Stack two matrices horizontally and return the result.'''
    var result = SizedCMatrix[type, rows, A_cols + B_cols](fill_zeros=False)
    for r in range(rows):
        # TODO: Change to memcpy for improved efficiency
        for c in range(A_cols):
            result[r, c] = A[r, c]
        for c in range(B_cols):
            result[r, c + A_cols] = B[r, c]
    return result


fn vstack[type: DType, A_rows: Int, B_rows: Int, cols: Int](
    A: SizedCMatrix[type, A_rows, cols], B: SizedCMatrix[type, B_rows, cols]
) -> SizedCMatrix[type, A_rows + B_rows, cols]:
    '''Stack two matrices vertically and return the result.'''
    var result = SizedCMatrix[type, A_rows + B_rows, cols](fill_zeros=False)
    for c in range(cols):
        # TODO: Change to memcpy for improved efficiency
        for r in range(A_rows):
            result[r, c] = A[r, c]
        for r in range(B_rows):
            result[r + A_rows, c] = B[r, c]
    return result


# Constants for expm
alias b_d = List[List[Int, True]](
    List[Int, True](120, 60, 12, 1),
    List[Int, True](30240, 15120, 3360, 420, 30, 1),
    List[Int, True](17297280, 8648640,1995840, 277200, 25200,1512, 56,1),
    List[Int, True](17643225600, 8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1),
    List[Int, True](64764752532480000, 32382376266240000, 7771770303897600, 1187353796428800, 129060195264000, 10559470521600, 670442572800, 33522128640, 1323241920, 40840800, 960960, 16380, 182, 1),
)
alias theta3 = 0.015
alias theta5 = 0.25
alias theta7 = 0.95
alias theta9 = 2.1
alias theta13 = 5.4


fn _expm_pade[type: DType, size: Int, m: Int](A: SizedCMatrix[type, size, size]) raises -> SizedCMatrix[type, size, size]:
    # alias b: List[Int, True] = b_d[(m - 3) // 2]  # This produces a segmentation fault...
    var b: List[Int, True] = b_d[(m - 3) // 2]

    var U: SizedCMatrix[type, size, size] = A.eye_like() * ComplexScalar[type](b[1])
    var V: SizedCMatrix[type, size, size] = A.eye_like() * ComplexScalar[type](b[0])
    var A2 = A @ A
    var A2n = A.eye_like()

    @parameter
    for i in range(1, m // 2 + 1):
        A2n @= A2
        U += A2n * b[2 * i + 1]
        V += A2n * b[2 * i]
    U = A @ U
    return solve(V - U, V + U)
    # return (V - U).inv() @ (V + U)


fn _expm_ss[type: DType, size: Int](A: SizedCMatrix[type, size, size], norm: SIMD[type, 1]) raises -> SizedCMatrix[type, size, size]:
    # alias b = b_d[4]  # This produces a segmentation fault...
    var b = b_d[4]
    alias inv_log_2: SIMD[type, 1] = log(2.0).cast[type]()
    
    var Ac: SizedCMatrix[type, size, size] = A
    var s: Int = max(0, int(ceil(log(norm / theta13) * inv_log_2)))
    if s > 0:
        Ac /= 2**s

    var I: SizedCMatrix[type, size, size] = Ac.eye_like()
    var A2: SizedCMatrix[type, size, size] = Ac @ Ac
    var A4: SizedCMatrix[type, size, size] = A2 @ A2
    var A6: SizedCMatrix[type, size, size] = A2 @ A4
    var U: SizedCMatrix[type, size, size] = Ac @ (A6 @ (b[13] * A6 + b[11] * A4 + b[9] * A2) + b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * I)
    var V: SizedCMatrix[type, size, size] = A6 @ (b[12] * A6 + b[10] * A4 + b[8] * A2) + b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * I
    var r13: SizedCMatrix[type, size, size] = solve(V - U, V + U)
    # var r13: SizedCMatrix[type, size, size] = (V - U).inv() @ (V + U)
    return matrix_power(r13, 2**s)


fn expm[type: DType, size: Int](A: SizedCMatrix[type, size, size]) raises -> SizedCMatrix[type, size, size]:
    '''Compute the matrix exponential of a square matrix A.
    
    Implements algorithm 10.20 from [1].

    [1] Functions of Matrices: Theory and Computation, Nicholas J. Higham, 2008.
    '''
    var norm: SIMD[type, 1] = one_norm(A)
    if norm < theta3:
        return _expm_pade[m=3](A)
    elif norm < theta5:
        return _expm_pade[m=5](A)
    elif norm < theta7:
        return _expm_pade[m=7](A)
    elif norm < theta9:
        return _expm_pade[m=9](A)
    return _expm_ss(A, norm)
    