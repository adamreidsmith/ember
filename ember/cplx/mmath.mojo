from math import inf as INF, ceil, log
from sys import simdwidthof
from algorithm import parallelize, vectorize
from memory import memcpy

from .complexsimd import ComplexScalar
from .cmatrix import CMatrix
from .csrcmatrix import CSRCMatrix
from .qr import eigvals
from ..config import DEFAULT_TOL


fn kron_sequential[type: DType](A: CMatrix[type], B: CMatrix[type]) -> CMatrix[type]:
    '''Implements the Kronecker product of A with B using the naive algorithm. By convention, 
    if A or B is degenerate, (i.e., rows=0 or cols=0) the other matrix is returned.

    Parameters:
        type: A type for the matrices.

    Args:
        A: The first matrix.
        B: The second matrix.
    
    Returns:
        The Kronecker product of A and B.
    '''
    if A.size == 0:
        return B
    if B.size == 0:
        return A
    var result = CMatrix[type](rows=A.rows * B.rows, cols=A.cols * B.cols, fill_zeros=False)
    var block: CMatrix[type]
    for r_a in range(A.rows):
        for c_a in range(A.cols):
            block = A.load_crd[1](r_a, c_a) * B
            for r_b in range(B.rows):
                for c_b in range(B.cols):
                    result.store_crd[1](
                        r_a * B.rows + r_b, c_a * B.cols + c_b, block.load_crd[1](r_b, c_b)
                    )
    return result

fn kron[type: DType](*M: CMatrix[type]) -> CMatrix[type]:
    '''Implements the Kronecker product of all supplied matrices.

    Parameters:
        type: A type for the matrices.

    Args:
        M: The matrices to take the Kronecker product of.
    
    Returns:
        The Kronecker product of all supplied matrices.
    '''
    if len(M) == 1:
        return M[0]
    var result = CMatrix[type](rows=0, cols=0, fill_zeros=False)
    for mat in M:
        result = kron(result, mat[])
    return result^

fn kron[type: DType](A: CMatrix[type], B: CMatrix[type]) -> CMatrix[type]:
    '''Implements the Kronecker product of A with B. By convention, if A or B is degenerate,
    (that is, rows=0 or cols=0) the other matrix is returned.

    Parameters:
        type: A type for the matrices.

    Args:
        A: The first matrix.
        B: The second matrix.
    
    Returns:
        The Kronecker product of A and B.
    '''
    if A.size == 0:
        return B
    if B.size == 0:
        return A
    if A.size == 1:
        if A.load_idx[1](0) == 1:
            return B
        return A.load_idx[1](0) * B
    if B.size == 1:
        if B.load_idx[1](0) == 1:
            return A
        return B.load_idx[1](0) * A
    # TODO: Benchmark to see which is better
    return _kron_par_a_rows(A, B)
    # return _kron_par_a_cols(A, B)


fn _kron_par_a_rows[type: DType](A: CMatrix[type], B: CMatrix[type]) -> CMatrix[type]:
    '''Implements the Kronecker product of A with B parallelized over the rows of A.

    Parameters:
        type: A type for the matrices.

    Args:
        A: The first matrix.
        B: The second matrix.
    
    Returns:
        The Kronecker product of A and B.
    '''
    var result = CMatrix[type](rows=A.rows * B.rows, cols=A.cols * B.cols, fill_zeros=False)
    @parameter
    fn par_row_a(r_a: Int):
        for c_a in range(A.cols):
            var a_elem: ComplexScalar[type] = A.load_crd[1](r_a, c_a)
            @parameter
            fn vec_col_b[simd_width: Int](c_b: Int):
                for r_b in range(B.rows):
                    result.store_crd[simd_width](
                        r_a * B.rows + r_b, c_a * B.cols + c_b,
                        a_elem * B.load_crd[simd_width](r_b, c_b),
                    )
            vectorize[vec_col_b, simdwidthof[type]()](B.cols)
    parallelize[par_row_a](A.rows)
    return result


fn _kron_par_a_cols[type: DType](A: CMatrix[type], B: CMatrix[type]) -> CMatrix[type]:
    '''Implements the Kronecker product of A with B parallelized over the columns of A.

    Parameters:
        type: A type for the matrices.

    Args:
        A: The first matrix.
        B: The second matrix.
    
    Returns:
        The Kronecker product of A and B.
    '''
    var result = CMatrix[type](rows=A.rows * B.rows, cols=A.cols * B.cols, fill_zeros=False)
    @parameter
    fn par_col_a(c_a: Int):
        for r_a in range(A.rows):
            var a_elem: ComplexScalar[type] = A.load_crd[1](r_a, c_a)
            @parameter
            fn vec_col_b[simd_width: Int](c_b: Int):
                for r_b in range(B.rows):
                    result.store_crd[simd_width](
                        r_a * B.rows + r_b, c_a * B.cols + c_b,
                        a_elem * B.load_crd[simd_width](r_b, c_b),
                    )
            vectorize[vec_col_b, simdwidthof[type]()](B.cols)
    parallelize[par_col_a](A.cols)
    return result


fn sparse_kron[type: DType](*M: CSRCMatrix[type]) -> CSRCMatrix[type]:
    '''Implements the Kronecker product of all supplied sparse matrices.

    Parameters:
        type: A type for the sparse matrices.

    Args:
        M: The sparse matrices to take the Kronecker product of.
    
    Returns:
        The Kronecker product of all supplied sparse matrices.
    '''
    if len(M) == 1:
        return M[0]
    var result = CSRCMatrix[type](rows=0, cols=0)
    for mat in M:
        result = sparse_kron(result, mat[])
    return result^


fn sparse_kron[type: DType](A: CSRCMatrix[type], B: CSRCMatrix[type]) -> CSRCMatrix[type]:
    '''Implements the Kronecker product of sparse matrices A with B. Note that the result
    uses the zero threshold of A.

    Parameters:
        type: A type for the sparse matrices.

    Args:
        A: The first sparse matrix.
        B: The second sparse matrix.
    
    Returns:
        The Kronecker product of A and B.
    '''
    if A.size == 0:
        return B
    if B.size == 0:
        return A
    if A.size == 1:
        if A.v:
            if A.v[0] == 1:
                return B
            return A.v[0] * B
        return B.zeros_like()
    if B.size == 1:
        if B.v:
            if B.v[0] == 1:
                return A
            return B.v[0] * A
        return A.zeros_like()
        
    var result_v = List[ComplexScalar[type], True]()
    var result_col_idx = List[Int, True]()
    var result_row_idx = List[Int, True](capacity=A.rows * B.rows + 1)
    result_row_idx.append(0)
    for r_a in range(A.rows):
        var a_start: Int = A.row_idx[r_a]
        var a_end: Int = A.row_idx[r_a + 1]
        for r_b in range(B.rows):
            var b_start: Int = B.row_idx[r_b]
            var b_end: Int = B.row_idx[r_b + 1]
            for i_a in range(a_start, a_end):
                var c_a: Int = A.col_idx[i_a]
                var v_a: ComplexScalar[type] = A.v[i_a]
                for i_b in range(b_start, b_end):
                    var prod = v_a * B.v[i_b]
                    if prod > A.zero_threshold:
                        var result_col: Int = c_a * B.cols + B.col_idx[i_b]
                        result_v.append(prod)
                        result_col_idx.append(result_col)
            result_row_idx.append(len(result_v))
    return CSRCMatrix[type](
        rows=A.rows * B.rows,
        cols=A.cols * B.cols,
        size=A.size * B.size,
        v=result_v^,
        col_idx=result_col_idx^,
        row_idx=result_row_idx^,
    )


fn swap_rows[type: DType](A: CMatrix[type], r1: Int, r2: Int) raises -> CMatrix[type]:
    '''Swap rows r1 and r2 in A and return the result.

    Parameters:
        type: A type for the matrix.

    Args:
        A: The matrix to swap rows in.
        r1: The index of the first row.
        r2: The index of the second row.
    
    Returns:
        The matrix with rows r1 and r2 swapped.
    '''
    if r1 == r2:
        return A
    if not 0 <= r1 < A.rows or not 0 <= r2 < A.rows:
        raise Error('Invalid row indices')
    var result = A
    # # memcpy seems unreliable atm
    # memcpy(result.re + r1 * result.cols, A.re + r2 * A.cols, A.cols)
    # memcpy(result.im + r1 * result.cols, A.im + r2 * A.cols, A.cols)
    # memcpy(result.re + r2 * result.cols, A.re + r1 * A.cols, A.cols)
    # memcpy(result.im + r2 * result.cols, A.im + r1 * A.cols, A.cols)
    for c in range(A.cols):
        result.store_crd[1](r1, c, A.load_crd[1](r2, c))
        result.store_crd[1](r2, c, A.load_crd[1](r1, c))
    return result


fn swap_rows_inplace[type: DType](A: CMatrix[type], r1: Int, r2: Int) raises:
    '''Swap rows r1 and r2 in A in-place.

    Parameters:
        type: A type for the matrix.

    Args:
        A: The matrix to swap rows in.
        r1: The index of the first row.
        r2: The index of the second row.
    '''
    if r1 == r2:
        return
    if not 0 <= r1 < A.rows or not 0 <= r2 < A.rows:
        raise Error('Invalid row indices')
    # # This doesn't work...
    # var r1_start: Int = r1 * cols
    # var r2_start: Int = r2 * cols
    # var row_r1 = CMatrix[type, 1, cols](fill_zeros=False)
    # memcpy(row_r1.re, A.re + r1_start, cols)
    # memcpy(row_r1.im, A.im + r1_start, cols)
    # # print(A)  # But somehow adding this print statement makes it work??
    # memcpy(A.re + r1_start, A.re + r2_start, cols)
    # memcpy(A.im + r1_start, A.im + r2_start, cols)
    # memcpy(A.re + r2_start, row_r1.re, cols)
    # memcpy(A.im + r2_start, row_r1.im, cols)
    for c in range(A.cols):
        var p: ComplexScalar[type] = A.load_crd[1](r1, c)
        A.store_crd[1](r1, c, A.load_crd[1](r2, c))
        A.store_crd[1](r2, c, p)


fn swap_cols[type: DType](A: CMatrix[type], c1: Int, c2: Int) raises -> CMatrix[type]:
    '''Swap columns c1 and c2 in A and return the result.

    Parameters:
        type: A type for the matrix.

    Args:
        A: The matrix to swap columns in.
        c1: The index of the first column.
        c2: The index of the second column.
    
    Returns:
        The matrix with columns c1 and c2 swapped.
    '''
    if c1 == c2:
        return A
    if not 0 <= c1 < A.cols or not 0 <= c2 < A.cols:
        raise Error('Invalid col indices')
    var result = A
    for r in range(A.rows):
        result.store_crd[1](r, c1, A.load_crd[1](r, c2))
        result.store_crd[1](r, c2, A.load_crd[1](r, c1))
    return result


fn swap_cols_inplace[type: DType](A: CMatrix[type], c1: Int, c2: Int) raises:
    '''Swap columns c1 and c2 in A in-place.

    Parameters:
        type: A type for the matrix.

    Args:
        A: The matrix to swap columns in.
        c1: The index of the first column.
        c2: The index of the second column.
    '''
    if c1 == c2:
        return
    if not 0 <= c1 < A.cols or not 0 <= c2 < A.cols:
        raise Error('Invalid col indices')
    for r in range(A.rows):
        var p: ComplexScalar[type] = A.load_crd[1](r, c1)
        A.store_crd[1](r, c1, A.load_crd[1](r, c2))
        A.store_crd[1](r, c2, p)


fn swap_vals[type: DType](
    A: CMatrix[type], r1: Int, c1: Int, r2: Int, c2: Int
) raises -> CMatrix[type]:
    '''Swap values at (r1, c1) and (r2, c2) in A and return the result.

    Parameters:
        type: A type for the matrix.
    
    Args:
        A: The matrix to swap values in.
        r1: The row index of the first value.
        c1: The column index of the first value.
        r2: The row index of the second value.
        c2: The column index of the second value.
    
    Returns:
        The matrix with values at (r1, c1) and (r2, c2) swapped.
    '''
    if r1 == r2 and c1 == c2:
        return A
    if not 0 <= r1 < A.rows or not 0 <= r2 < A.rows:
        raise Error('Invalid row indices')
    if not 0 <= c1 < A.cols or not 0 <= c2 < A.cols:
        raise Error('Invalid col indices')
    var result = A
    result.store_crd[1](r1, c1, A.load_crd[1](r2, c2))
    result.store_crd[1](r2, c2, A.load_crd[1](r1, c1))
    return result


fn swap_vals_inplace[type: DType](mut A: CMatrix[type], r1: Int, c1: Int, r2: Int, c2: Int) raises:
    '''Swap values at (r1, c1) and (r2, c2) in-place.

    Parameters:
        type: A type for the matrix.
    
    Args:
        A: The matrix to swap values in.
        r1: The row index of the first value.
        c1: The column index of the first value.
        r2: The row index of the second value.
        c2: The column index of the second value.
    '''
    if r1 == r2 and c1 == c2:
        return
    if not 0 <= r1 < A.rows or not 0 <= r2 < A.rows:
        raise Error('Invalid row indices')
    if not 0 <= c1 < A.cols or not 0 <= c2 < A.cols:
        raise Error('Invalid col indices')
    var p1: ComplexScalar[type] = A.load_crd[1](r1, c1)
    A.store_crd[1](r1, c1, A.load_crd[1](r2, c2))
    A.store_crd[1](r2, c2, p1)


fn augmented_ref[type: DType, tol: Scalar[type] = DEFAULT_TOL](
    A: CMatrix[type], B: CMatrix[type]
) raises -> CMatrix[type]:
    '''Computes the row echelon form of the augmented matrix [A|B].
    
    Parameters:
        type: A type for the matrices.
        tol: A tolerance for zero checks.
    
    Args:
        A: The left matrix.
        B: The right matrix.
    
    Returns:
        The row echelon form of the augmented matrix [A|B].
    '''
    var Aaug: CMatrix[type] = hstack(A, B)
    var h: Int = 0  # initialize pivot row
    var k: Int = 0  # initialize pivot column

    while h < A.rows and k < A.cols:
        # Find the k-th pivot:
        var i_max: Scalar[type] = Aaug.load_crd[1](h, k).norm()
        var i_argmax: Int = h
        for i in range(h + 1, A.rows):
            var i_norm: Scalar[type] = Aaug.load_crd[1](i, k).norm()
            if i_norm > i_max:
                i_max = i_norm
                i_argmax = i
        if Aaug.load_crd[1](i_argmax, k) < tol:
            # No pivot in this column, pass to next column
            k += 1
        else:
            swap_rows_inplace(Aaug, h, i_argmax)
            # Do for all rows below pivot:
            for i in range(h + 1, A.rows):
                var f: ComplexScalar[type] = Aaug.load_crd[1](i, k) / Aaug.load_crd[1](h, k)
                # Fill with zeros the lower part of pivot column:
                Aaug.store_crd[1](i, k, ComplexScalar[type](0))
                # Do for all remaining elements in current row:
                for j in range(k + 1, A.cols + B.cols):
                    Aaug.store_crd[1](i, j, Aaug.load_crd[1](i, j) - Aaug.load_crd[1](h, j) * f)
            # Increase pivot row and column
            h += 1; k += 1
    return Aaug


# TODO: Switch to a faster algorithm
fn solve[type: DType, tol: Scalar[type] = DEFAULT_TOL](
    A: CMatrix[type], B: CMatrix[type]
) raises -> CMatrix[type]:
    '''Solves a linear system of equations Ax=B via Gaussian elimination.
    
    Parameters:
        type: A type for the matrices.
        tol: A tolerance for zero checks.
    
    Args:
        A: The LHS matrix of the linear system.
        B: The RHS matrix of the linear system.
    
    Returns:
        The solution to the linear system.
    '''
    if A.rows != A.cols:
        raise Error('A must be a square matrix for `solve`')

    var aug_ref: CMatrix[type] = augmented_ref[type, tol](A, B)

    for i in range(A.rows):
        if aug_ref.load_crd[1](i, i) < tol:
            raise Error('Matrix is singular or nearly singular')

    # Back substitution
    var X: CMatrix[type] = B.zeros_like()
    for c in range(B.cols):
        for i in range(A.rows - 1, -1, -1):
            var dot = ComplexScalar[type](0)
            for j in range(i + 1, A.cols):
                dot += aug_ref.load_crd[1](i, j) * X.load_crd[1](j, c)
            X.store_crd[1](
                i, c, (aug_ref.load_crd[1](i, c + A.cols) - dot) / aug_ref.load_crd[1](i, i)
            )
    return X


fn one_norm[type: DType](A: CMatrix[type]) -> Scalar[type]:
    '''Computes the one-norm of a matrix, i.e. the maximum absolute column sum.
    
    Parameters:
        type: A type for the matrix.
    
    Args:
        A: The matrix to compute the one-norm of.
    
    Returns:
        The one norm of A.
    '''
    var col_sums = CMatrix[type](rows=1, cols=A.cols, fill_zeros=False)
    @parameter
    fn sum_abs_cols(c: Int):
        var col_sum = ComplexScalar[type]()
        for r in range(A.rows):
            col_sum += abs(A.load_crd[1](r, c))
        col_sums.store_idx[1](c, col_sum)
    parallelize[sum_abs_cols](A.cols, A.cols)
    var max: Scalar[type] = col_sums.load_idx[1](0).re
    for c in range(1, A.cols):
        if (cs := col_sums.load_idx[1](c).re) > max:
            max = cs
    return max


fn _int_matrix_power[type: DType](owned A: CMatrix[type], n: Int) raises -> CMatrix[type]:
    '''Computes A^n for a square matrix A and positive integer n.
    
    Parameters:
        type: A type for the matrix.
    
    Args:
        A: The matrix to compute the power of.
        n: A positive integer power.
    
    Returns:
        The matrix A to the power of n.
    '''
    if A.rows != A.cols:
        raise Error('Cannot compute a power of a non-square matrix')
    if n == 1:
        return A
    if n % 2 == 0:
        return _int_matrix_power(A @ A, n // 2)
    else:
        return _int_matrix_power(A @ A, (n - 1) // 2) @ A


fn _int_matrix_power[type: DType](owned A: CSRCMatrix[type], n: Int) raises -> CSRCMatrix[type]:
    '''Computes A^n for a square sparse matrix A and positive integer n.
    
    Parameters:
        type: A type for the sparse matrix.
    
    Args:
        A: The sparse matrix to compute the power of.
        n: A positive integer power.
    
    Returns:
        The sparse matrix A to the power of n.
    '''
    if A.rows != A.cols:
        raise Error('Cannot compute a power of a non-square matrix')
    if n == 1:
        return A
    if n % 2 == 0:
        return _int_matrix_power(A @ A, n // 2)
    else:
        return _int_matrix_power(A @ A, (n - 1) // 2) @ A


fn matrix_power[type: DType](A: CMatrix[type], n: Int) raises -> CMatrix[type]:
    '''Computes A^n for a matrix A and integer n.
    
    Parameters:
        type: A type for the matrix.
    
    Args:
        A: The matrix to compute the power of.
        n: An integer power.
    
    Returns:
        The matrix A to the power of n.
    '''
    if A.rows != A.cols:
        raise Error('Cannot compute a power of a non-square matrix')
    if n == 0:
        return A.eye_like()
    var Ac = A
    if n < 0:
        return _int_matrix_power(Ac^, -n).inv()
    return _int_matrix_power(Ac^, n)


fn matrix_power[type: DType](A: CSRCMatrix[type], n: Int) raises -> CSRCMatrix[type]:
    '''Computes A^n for a sparse matrix A and integer n.
    
    Parameters:
        type: A type for the sparse matrix.
    
    Args:
        A: The sparse matrix to compute the power of.
        n: An integer power.
    
    Returns:
        The sparse matrix A to the power of n.
    '''
    if A.rows != A.cols:
        raise Error('Cannot compute a power of a non-square matrix')
    if n < 0:
        raise Error('Negative powers of sparse matrices are not supported at this time')
    if n == 0:
        return A.eye_like()
    var Ac = A
    return _int_matrix_power(Ac^, n)


fn mmax[type: DType](A: CMatrix[type]) -> ComplexScalar[type]:
    '''Returns the maximum value (by modulus) of the matrix.
    
    Parameters:
        type: A type for the matrix.
    
    Args:
        A: The matrix to compute the max of.
    
    Returns:
        The max value of A by modulus.
    '''
    var mx: ComplexScalar[type] = A.load_idx[1](0)
    for idx in range(A.size):
        if (m := A.load_idx[1](idx)) > mx:  # ComplexSIMD's are compared by norm
            mx = m
    return mx


fn mmin[type: DType](A: CMatrix[type]) -> ComplexScalar[type]:
    '''Returns the minimum value (by modulus) of the matrix.
    
    Parameters:
        type: A type for the matrix.
    
    Args:
        A: The matrix to compute the min of.
    
    Returns:
        The min value of A by modulus.
    '''
    var mx: ComplexScalar[type] = A.load_idx[1](0)
    for idx in range(A.size):
        if (m := A.load_idx[1](idx)) < mx:  # ComplexSIMD's are compared by norm
            mx = m
    return mx


fn hstack[type: DType](A: CMatrix[type], B: CMatrix[type]) raises -> CMatrix[type]:
    '''Stack two matrices horizontally and return the result.
    
    Parameters:
        type: A type for the matrices.
    
    Args:
        A: The left matrix.
        B: The right matrix.
    
    Returns:
        The horizontally stacked matrix.
    '''
    if A.rows != B.rows:
        raise Error(
            'Invalid row dimensions for `hstack`: ' + String(A.rows) + ' and ' + String(B.rows)
        )
    var result = CMatrix[type](rows=A.rows, cols=A.cols + B.cols, fill_zeros=False)
    for r in range(A.rows):
        # TODO: Change to memcpy
        for c in range(A.cols):
            result.store_crd[1](r, c, A.load_crd[1](r, c))
        for c in range(B.cols):
            result.store_crd[1](r, c + A.cols, B.load_crd[1](r, c))
    return result


fn vstack[type: DType](A: CMatrix[type], B: CMatrix[type]) raises -> CMatrix[type]:
    '''Stack two matrices vertically and return the result.
    
    Parameters:
        type: A type for the matrices.
    
    Args:
        A: The left matrix.
        B: The right matrix.
    
    Returns:
        The vertically stacked matrix.
    '''
    if A.cols != B.cols:
        raise Error(
            'Invalid column dimensions for `vstack`: ' + String(A.cols) + ' and ' + String(B.cols)
        )
    var result = CMatrix[type](rows=A.rows + B.rows, cols=A.cols, fill_zeros=False)
    for c in range(A.cols):
        # TODO: Change to memcpy
        for r in range(A.rows):
            result.store_crd[1](r, c, A.load_crd[1](r, c))
        for r in range(B.rows):
            result.store_crd[1](r + A.rows, c, B.load_crd[1](r, c))
    return result


# Constants for expm
alias b_d = List[List[Int, True]](
    List[Int, True](120, 60, 12, 1),
    List[Int, True](30240, 15120, 3360, 420, 30, 1),
    List[Int, True](17297280, 8648640,1995840, 277200, 25200,1512, 56,1),
    List[Int, True](
        17643225600, 8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1
    ),
    List[Int, True](
        64764752532480000, 32382376266240000, 7771770303897600, 1187353796428800, 129060195264000,
        10559470521600, 670442572800, 33522128640, 1323241920, 40840800, 960960, 16380, 182, 1
    ),
)
alias theta3 = 0.015
alias theta5 = 0.25
alias theta7 = 0.95
alias theta9 = 2.1
alias theta13 = 5.4


fn _expm_pade[m: Int, type: DType](A: CMatrix[type]) raises -> CMatrix[type]:
    '''Compute the matrix exponential using Pade approximation.

    Parameters:
        m: The degree of the Pade approximation.
        type: A type for the matrix.

    Args:
        A: The matrix to exponentiate.
    
    Returns:
        The matrix exponential.
    '''
    if A.rows != A.cols:
        raise Error('Cannot exponentiate a non-square matrix')
    var b: List[Int, True] = b_d[(m - 3) // 2]

    var U: CMatrix[type] = A.eye_like() * ComplexScalar[type](b[1])
    var V: CMatrix[type] = A.eye_like() * ComplexScalar[type](b[0])
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


fn _expm_ss[type: DType](A: CMatrix[type], norm: Scalar[type]) raises -> CMatrix[type]:
    '''Compute the matrix exponential using the scaling and squaring algorithm.

    Parameters:
        type: A type for the matrix.

    Args:
        A: The matrix to exponentiate.
        norm: The one-norm of A.
    
    Returns:
        The matrix exponential.
    '''
    if A.rows != A.cols:
        raise Error('Cannot exponentiate a non-square matrix')
    # alias b = b_d[4]  # This produces a segmentation fault...
    var b = b_d[4]
    alias inv_log_2: Scalar[type] = log(2.0).cast[type]()

    var Ac: CMatrix[type] = A
    var s: Int = max(0, Int(ceil(log(norm / theta13) * inv_log_2)))
    if s > 0:
        Ac /= 2 ** s

    var I: CMatrix[type] = Ac.eye_like()
    var A2: CMatrix[type] = Ac @ Ac
    var A4: CMatrix[type] = A2 @ A2
    var A6: CMatrix[type] = A2 @ A4
    var U: CMatrix[type] = (
        Ac 
        @ (
            A6 
            @ (b[13] * A6 + b[11] * A4 + b[9] * A2) + b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * I
        )
    )
    var V: CMatrix[type] = (
        A6 
        @ (b[12] * A6 + b[10] * A4 + b[8] * A2) + b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * I
    )
    var r13: CMatrix[type] = solve(V - U, V + U)
    return matrix_power(r13, 2 ** s)


fn expm[type: DType](A: CMatrix[type]) raises -> CMatrix[type]:
    '''Compute the matrix exponential of a square matrix A.
    
    Implements algorithm 10.20 from [1].

    [1] Functions of Matrices: Theory and Computation, Nicholas J. Higham, 2008.

    Parameters:
        type: A type for the matrix.

    Args:
        A: The matrix to exponentiate.
    
    Returns:
        The matrix exponential.
    '''
    if A.rows != A.cols:
        raise Error('Cannot exponentiate a non-square matrix')
    var norm: Scalar[type] = one_norm(A)
    if norm < theta3:
        return _expm_pade[m=3](A)
    elif norm < theta5:
        return _expm_pade[m=5](A)
    elif norm < theta7:
        return _expm_pade[m=7](A)
    elif norm < theta9:
        return _expm_pade[m=9](A)
    return _expm_ss(A, norm)


fn kron_power[type: DType](owned A: CMatrix[type], n: Int) raises -> CMatrix[type]:
    '''Computes the Kronecker product of a matrix with itself n times. By convention, n=0 returns a
    degenerate matrix of size (0, 0).

    Parameters:
        type: A type for the matrix.

    Args:
        A: The matrix to compute the Kroncker power of.
        n: A nonnegative integer power.
    
    Returns:
        The Kronecker power A^n.
    '''
    if A.rows != A.cols:
        raise Error('Cannot compute a power of a non-square matrix')
    if n < 0:
        raise Error('Cannot apply an Kronecker power less than 1')
    if n == 0:
        return CMatrix[type](0, 0)
    return _kron_power(A^, n)
    

fn _kron_power[type: DType](owned A: CMatrix[type], n: Int) -> CMatrix[type]:
    '''Computes the Kronecker product of a matrix with itself n times.
    
    Parameters:
        type: A type for the matrix.

    Args:
        A: The matrix to compute the Kroncker power of.
        n: A positive integer power.
    
    Returns:
        The Kronecker power A^n.
    '''
    if n == 1:
        return A
    if n % 2 == 0:
        return _kron_power(kron(A, A), n // 2)
    else:
        return kron(_kron_power(kron(A, A), (n - 1) // 2), A)


fn kron_power[type: DType](owned A: CSRCMatrix[type], n: Int) raises -> CSRCMatrix[type]:
    '''Computes the Kronecker product of a sparse matrix with itself n times. By convention, n=0
    returns a degenerate matrix of size (0, 0).

    Parameters:
        type: A type for the sparse matrix.

    Args:
        A: The sparse matrix to compute the Kroncker power of.
        n: A nonnegative integer power.
    
    Returns:
        The Kronecker power A^n.
    '''
    if A.rows != A.cols:
        raise Error('Cannot compute a power of a non-square matrix')
    if n < 0:
        raise Error('Cannot apply negative Kronecker powers')
    if n == 0:
        return CSRCMatrix[type](0, 0)
    return _kron_power(A^, n)
    

fn _kron_power[type: DType](owned A: CSRCMatrix[type], n: Int) -> CSRCMatrix[type]:
    '''Computes the Kronecker product of a sparse matrix with itself n times. By convention, n=0
    returns a degenerate matrix of size (0, 0).

    Parameters:
        type: A type for the sparse matrix.

    Args:
        A: The sparse matrix to compute the Kroncker power of.
        n: A positive integer power.
    
    Returns:
        The Kronecker power A^n.
    '''
    if n == 1:
        return A
    if n % 2 == 0:
        return _kron_power(sparse_kron(A, A), n // 2)
    else:
        return sparse_kron(_kron_power(sparse_kron(A, A), (n - 1) // 2), A)


fn is_positive_semidefinite[type: DType, tol: Scalar[type] = DEFAULT_TOL](A: CMatrix[type]) raises -> Bool:
    '''Check if the CMatrix is positive semi-definite.

    Parameters:
        type: A type for the sparse matrix.
        tol: A tolerance for closeness checks.

    Args:
        A: The matrix to check.

    Returns:
        True if the CMatrix is positive semi-definite, False otherwise.
    '''
    if not A.is_hermitian[tol]():  # Also checks if the matrix is square
        return False
    
    var eigenvalues: List[ComplexScalar[type], True] = eigvals(A)
    for val in eigenvalues:
        if abs(val[].im) > tol:
            raise Error(
                'Non-real eigenvalue computed: ' + String(val[]) 
                + '. Matrix is not Hermitian (or bug in eigenvalue code!)'
            )
        if val[].re < -tol:
            return False
    return True
