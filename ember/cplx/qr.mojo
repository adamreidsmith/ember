from math import sqrt

from .complexsimd import ComplexScalar
from .cmatrix import CMatrix
from ..config import DEFAULT_TYPE, DEFAULT_TOL

alias sqrt2: Float64 = sqrt(Float64(2))

fn _house_gen[type: DType = DEFAULT_TYPE](
    a: CMatrix[type]
) -> Tuple[CMatrix[type], ComplexScalar[type]]:
    '''This function computes a vector u that generates a Householder reflection H = I - uu*
    satisfying Ha = -ve_1 for some complex number v. Accumulating these transformations, any
    square matrix can be reduced to upper Hessenberg form.

    Parameters:
        type: A type for the data.

    Args:
        a: A vector to generate the Householder reflection.
    
    Returns:
        A tuple of the Householder reflection generating vector u and the scalar v.
    '''
    var u: CMatrix[type] = a
    var v: ComplexScalar[type] = a.frobenius_norm()

    if v == 0:
        u.store_idx[1](0, sqrt2.cast[type]())
        return (u, v)
    
    var rho: ComplexScalar[type]
    if u.load_idx[1](0) != 0:
        rho = u.load_idx[1](0).conj()
        rho /= rho.norm()
    else:
        rho = 1

    u *= rho / v
    u.store_idx[1](0, u.load_idx[1](0) + 1)
    u /= u.load_idx[1](0).sqrt()
    v *= -rho.conj()
    return (u^, v)
    

fn _hess_reduce[type: DType = DEFAULT_TYPE](
    A: CMatrix[type]
) raises -> Tuple[CMatrix[type], CMatrix[type]]:
    '''This algorithm reduces a square matrix A of order n to upper Hessenberg form using
    Householder transformations.
    
    Parameters:
        type: A type for the data.
    
    Args:
        A: The matrix to convert to upper Hessenberg form.
    
    Returns:
        A tuple of the upper Hessenberg form of A and the associated transformations.
    '''
    var n: Int = A.rows
    var H: CMatrix[type] = A
    var Q: CMatrix[type] = A.eye_like()

    for k in range(n - 2):
        # Generate the Householder transformation that annihilates the k-th column
        # [u, H(k+1, k)] = housegen(H(k+1:n, k));
        var u: CMatrix[type]
        var Hk: ComplexScalar[type]
        u, Hk = _house_gen[type](H._get_column(col=k, row_start=k + 1, row_end=n))
        H.store_crd[1](k + 1, k, Hk)

        # Q(k+1:n, k) = u;
        Q._set_block(row=k + 1, col=k, block=u)

        # Multiply the transformation on the left
        # v = (u’)*H(k+1:n,k+1:n);
        # H(k+1:n,k+1:n) = H(k+1:n,k+1:n) - u*v;
        var H_block: CMatrix[type] = H._get_block(
            row_start=k + 1, row_end=n, col_start=k + 1, col_end=n
        )
        H._set_block(row=k + 1, col=k + 1, block=H_block - u @ (u.dag() @ H_block))
        # H(k+2:n,k) = 0;
        H._set_block(row=k + 2, col=k, block=CMatrix[type].zeros(rows=n - k - 2, cols=1))

        # Multiply the transformation on the right
        # v = H(1:n,k+1:n)*u;
        # H(1:n,k+1:n) = H(1:n,k+1:n) - v*(u’);
        H_block = H._get_block(row_start=0, row_end=n, col_start=k + 1, col_end=n)
        H._set_block(row=0, col=k + 1, block = H_block - (H_block @ u) @ u.dag())

    # Accumulate transformations on matrix Q
    var I: CMatrix[type] = A.eye_like()
    for k in range(n - 3, -1, -1):
        # u = Q(k+1:n,k);
        # v = (u’)*Q(k+1:n,k+1:n);
        # Q(k+1:n,k+1:n) = Q(k+1:n,k+1:n) - u*v;
        # Q(:,k) = I(:,k);
        var u: CMatrix[type] = Q._get_column(col=k, row_start=k + 1, row_end=n)
        var Q_block: CMatrix[type] = Q._get_block(
            row_start=k + 1, row_end=n, col_start=k + 1, col_end=n
        )
        Q._set_block(row=k + 1, col=k + 1, block=Q_block - u @ (u.dag() @ Q_block))
        Q._set_block(row=0, col=k, block=I._get_block(
            row_start=0, row_end=n, col_start=k, col_end=k + 1)
        )
        
    return (H^, Q^)


fn _backsearch[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    H: CMatrix[type], z: Int
) -> Tuple[Int, Int]:
    '''This function finds deflating rows on a complex Schur form matrix.

    Parameters:
        type: A type for the data.
        tol: A tolerance for zero checks.

    Args:
        H: Hessenberg matrix of size n x n.
        z: Index 0 < z < n.

    Returns:
        Indices i1 and i2 (i1, i2 <= z) satisfying one of the following conditions:
            1) 0 <= i1 < i2 <= z, deflate at rows i1 and i2.
            2) 0 = i1 = i2, matrix is completely deflated.
    '''
    
    var i1: Int = z
    var i2: Int = z
    var S: CMatrix[type] = H  # Create a copy to avoid modifying the original
    var normS: Scalar[type] = S.frobenius_norm()

    while i1 > 0:
        if S.load_crd[1](i1, i1 - 1).norm() < tol * normS:
            S.store_crd[1](i1, i1 - 1, 0)
            if i1 == i2:
                i2 = i1 - 1
                i1 -= 1
            else:
                return (i1, i2)
        else:
            i1 -= 1
    return (i1, i2)


fn _wilkshift[type: DType = DEFAULT_TYPE](
    a: ComplexScalar[type], b: ComplexScalar[type], c: ComplexScalar[type], d: ComplexScalar[type]
) -> ComplexScalar[type]:
    '''This function computes the Wilkinson shift of a submatrix of order 2.

    Parameters:
        type: A type for the data.

    Args:
        a: [0, 0] matrix element.
        b: [0, 1] matrix element.
        c: [1, 0] matrix element.
        d: [1, 1] matrix element.

    Returns:
        The Wilkinson shift for the supplied matrix.
    '''
    var kappa: ComplexScalar[type] = d
    var s: Scalar[type] = a.norm() + b.norm() + c.norm() + d.norm()

    if s == 0:
        return kappa
    
    var q: ComplexScalar[type] = (b / s) * (c / s)
    if q != 0:
        var p: ComplexScalar[type] = 0.5 * ((a / s) - (d / s))
        var r: ComplexScalar[type] = (p * p + q).sqrt()
        if p.re * r.re + p.im * r.im < 0:
            r = -r
        kappa -= s * (q / (p + r))
    
    return kappa


fn _rot_gen[type: DType = DEFAULT_TYPE](
    a: ComplexScalar[type], b: ComplexScalar[type]
) -> Tuple[ComplexScalar[type], ComplexScalar[type], ComplexScalar[type], ComplexScalar[type]]:
    '''This function generates a Givens rotation from elements a and b.

    Parameters:
        type: A type for the data.

    Args:
        a: The first element.
        b: The second element.

    Returns:
        The values a and b and constants c and s defining the plane rotation in a tuple.
    '''
    if b == 0:
        return (a, ComplexScalar[type](0), ComplexScalar[type](1), ComplexScalar[type](0))
    if a == 0:
        return (b, ComplexScalar[type](0), ComplexScalar[type](0), ComplexScalar[type](1))

    var mu: ComplexScalar[type] = a / a.norm()
    var tau: Scalar[type] = abs(a.re) + abs(a.im) + abs(b.re) + abs(b.im)
    var nu: Scalar[type] = tau * sqrt((a / tau).squared_norm()  + (b / tau).squared_norm())
    return (nu * mu, ComplexScalar[type](0), abs(a) / nu, mu * b.conj() / nu)


fn _rot_app[type: DType = DEFAULT_TYPE](
    c: ComplexScalar[type], s: ComplexScalar[type], x: CMatrix[type], y: CMatrix[type]
) raises -> Tuple[CMatrix[type], CMatrix[type]]:
    '''This function takes a plane rotation defined by c and s (the scalars returned by _rot_gen)
    and applies it to the vectors x and y.

    Parameters:
        type: A type for the data.

    Args:
        c: The first scalar returned by _rot_gen
        s: The second scalar returned by _rot_gen.
        x: The first vector.
        y: the second vector.

    Returns:
        The vectors x and y after application of the plane rotation.
    '''
    return (c * x + s * y, c * y - s.conj() * x)


fn _hqr[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    H: CMatrix[type], Q: CMatrix[type], maxiter: Int
) raises -> Tuple[CMatrix[type], CMatrix[type]]:
    '''This returns a unitary similar triangular matrix whose diagonals are the eigenvalues of H
    (the real Schur form).
    
    Parameters:
        type: A type for the data.
        tol: A tolerance for zero checks.

    Args:
        H: An upper Hessenberg matrix.
        Q: The transformation matrix produced by _hess_reduce.
        maxiter: The maximum number of iterations.

    Returns:
        A tuple containing the real Schur form of H and the transformation matrix Q.
    '''
    var n: Int = H.rows
    var i2: Int = n - 1
    var c_vec = CMatrix[type](1, n)
    var s_vec = CMatrix[type](1, n)

    var S: CMatrix[type] = H
    var Qn: CMatrix[type] = Q

    var it: Int = 0

    while True:
        it += 1
        if it >= maxiter:
            raise Error('Maximum number of iterations reached.')

        var old_i2: Int = i2
        var i1: Int
        # Check subdiagonal for near zeros, deflating points
        i1, i2 = _backsearch[type, tol](S, i2)

        if i2 == 0:
            # End the function if S is upper triangular
            return (S, Qn)
        
        if i2 != old_i2:
            # Set iteration number to zero if there is another deflating row
            it = 0

        # Compute Wilkinson shift
        var kappa: ComplexScalar[type] = _wilkshift[type](
            S.load_crd[1](i2 - 1, i2 - 1),
            S.load_crd[1](i2 - 1, i2),
            S.load_crd[1](i2, i2 - 1),
            S.load_crd[1](i2, i2),
        )

        # Apply shift to the element of the diagonal that is left out of the loop
        S.store_crd[1](i1, i1, S.load_crd[1](i1, i1) - kappa)

        # Loop reducing the matrix to triangular form
        for j in range(i1, i2):
            # Apply rotation so that the subdiagonal is set to zero
            var a: ComplexScalar[type]
            var b: ComplexScalar[type]
            var c: ComplexScalar[type]
            var s: ComplexScalar[type]
            a, b, c, s = _rot_gen[type](S.load_crd[1](j, j), S.load_crd[1](j + 1, j))
            S.store_crd[1](j, j, a)
            S.store_crd[1](j + 1, j, b)
            c_vec.store_idx[1](j, c)
            s_vec.store_idx[1](j, s)

            # Apply shift to diagonal
            S.store_crd[1](j + 1, j + 1, S.load_crd[1](j + 1, j + 1) - kappa)

            # Modify the involved rows
            # [S(j,j+1:n), S(j+1,j+1:n)] = rotapp(c(j), s(j), S(j,j+1:n), S(j+1,j+1:n));
            var x: CMatrix[type]
            var y: CMatrix[type]
            x, y = _rot_app[type](c, s, S._get_row(j, j + 1, n), S._get_row(j + 1, j + 1, n))
            S._set_block(j, j + 1, x)
            S._set_block(j + 1, j + 1, y)
                
        # Loop applying the back multiplication
        for k in range(i1, i2):
            # [S(1:k+1,k), S(1:k+1,k+1)] = rotapp(c(k), conj(s(k)), S(1:k+1,k), S(1:k+1,k+1));
            var x: CMatrix[type]
            var y: CMatrix[type]
            x, y = _rot_app[type](
                c_vec.load_idx[1](k),
                s_vec.load_idx[1](k).conj(),
                # S._get_column(col=k, row_start=0, row_end=k + 1),
                # S._get_column(col=k + 1, row_start=0, row_end=k + 1),
                S._get_column(col=k, row_start=0, row_end=k + 2),
                S._get_column(col=k + 1, row_start=0, row_end=k + 2),
            )
            S._set_block(row=0, col=k, block=x)
            S._set_block(row=0, col=k + 1, block=y)

            # Accumulate transformations
            # [Qn(1:n,k), Qn(1:n,k+1)] = rotapp(c(k), conj(s(k)), Qn(1:n,k), Qn(1:n,k+1));
            x, y = _rot_app[type](
                c_vec.load_idx[1](k),
                s_vec.load_idx[1](k).conj(),
                Qn._get_column(col=k, row_start=0, row_end=n),
                Qn._get_column(col=k + 1, row_start=0, row_end=n)
            )
            Qn._set_block(row=0, col=k, block=x)
            Qn._set_block(row=0, col=k + 1, block=y)

            S.store_crd[1](k, k, S.load_crd[1](k, k) + kappa)
        
        S.store_crd[1](i2, i2, S.load_crd[1](i2, i2) + kappa)

fn complex_schur[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    A: CMatrix[type]
) raises -> Tuple[CMatrix[type], CMatrix[type]]:
    '''Compute the complex Schur decomposition of a matrix. Given a square matrix A, this function
    computes an upper triangular matrix S (the Schur form of A) and a unitary matrix Q satisfying
    A = Q @ S @ Q.dag(). The main diagonal of S gives the eigenvalues of A.

    Parameters:
        type: A type for the data.
        tol: A tolerance for zero checks.

    Args:
        A: A square matrix.
    
    Returns:
        A tuple containing the Schur form of A and the unitary matrix Q.
    '''
    if not A.is_square():
        raise Error('Matrix must be square to have a Schur decomposition.')

    var H: CMatrix[type]
    var Q: CMatrix[type]
    H, Q = _hess_reduce[type](A)
    return _hqr[type, tol](H, Q, maxiter=10_000)

fn eigvals[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    A: CMatrix[type]
) raises -> List[ComplexScalar[type], True]:
    '''Computes the eigenvalus of a matrix using the complex Schur decomposition.

    Parameters:
        type: A type for the data.
        tol: A tolerance for zero checks.

    Args:
        A: A square matrix.
    
    Returns:
        A tuple containing the Schur form of A and the unitary matrix Q.
    '''
    var S: CMatrix[type]
    S, _ = complex_schur[type, tol](A)

    var eigenvalues = List[ComplexScalar[type], True](capacity=S.rows)
    for i in range(S.rows):
        eigenvalues.append(S.load_crd[1](i, i))
    return eigenvalues^
