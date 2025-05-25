from testing import assert_equal, assert_raises
from ._testing import _assert_matrix_almost_equal, _assert_matrix_equal

from ember import (
    kron,
    kron_sequential,
    sparse_kron,
    swap_cols,
    swap_cols_inplace,
    swap_rows,
    swap_rows_inplace,
    swap_vals,
    swap_vals_inplace,
    mmax,
    mmin,
    augmented_ref,
    solve,
    expm,
    one_norm,
    matrix_power,
    hstack,
    vstack,
    is_positive_semidefinite,
)
from ember import CSRCMatrix, CMatrix, ComplexScalar

alias type = DType.float64

def run_mmath_tests():
    print('Running mmath tests')
    test_swap()
    test_mmax_mmin()
    test_kron()
    test_one_norm()
    test_stack()
    test_solve()
    test_matrix_power()
    test_expm()
    test_find_smallest_eigval()
    print('All tests passed')

def test_swap():
    a = CMatrix[type].arange(40, 50)
    b = swap_cols(a, 0, 1)
    for r in range(40):
        for c in range(50):
            cs = 1 if c == 0 else 0 if c == 1 else c
            assert_equal(a[r, c].re, b[r, cs].re, 'swap_cols')
            assert_equal(a[r, c].im, b[r, cs].im, 'swap_cols')
    b = swap_cols(a, 28, 15)
    for r in range(40):
        for c in range(50):
            cs = 15 if c == 28 else 28 if c == 15 else c
            assert_equal(a[r, c].re, b[r, cs].re, 'swap_cols')
            assert_equal(a[r, c].im, b[r, cs].im, 'swap_cols')
    with assert_raises():
        _ = swap_cols(a, -1, 1)
    with assert_raises():
        _ = swap_cols(a, 1, 50)
    b = swap_rows(a, 0, 1)
    for r in range(40):
        for c in range(50):
            rs = 0 if r == 1 else 1 if r == 0 else r
            assert_equal(a[r, c].re, b[rs, c].re, 'swap_rows')
            assert_equal(a[r, c].im, b[rs, c].im, 'swap_rows')
    b = swap_rows(a, 28, 15)
    for r in range(40):
        for c in range(50):
            rs = 28 if r == 15 else 15 if r == 28 else r
            assert_equal(a[r, c].re, b[rs, c].re, 'swap_rows')
            assert_equal(a[r, c].im, b[rs, c].im, 'swap_rows')
    with assert_raises():
        _ = swap_rows(a, -1, 1)
    with assert_raises():
        _ = swap_rows(a, 1, 40)
    b = swap_vals(a, 0, 0, 5, 5)
    for r in range(40):
        for c in range(50):
            swapped_val = b[5, 5] if r == 0 and c == 0 else b[0, 0] if r == 5 and c == 5 else b[r, c]
            assert_equal(a[r, c].re, swapped_val.re, 'swap_vals')
            assert_equal(a[r, c].im, swapped_val.im, 'swap_vals')
    b = swap_vals(a, 17, 49, 3, 33)
    for r in range(40):
        for c in range(50):
            swapped_val = b[17, 49] if r == 3 and c == 33 else b[3, 33] if r == 17 and c == 49 else b[r, c]
            assert_equal(a[r, c].re, swapped_val.re, 'swap_vals')
            assert_equal(a[r, c].im, swapped_val.im, 'swap_vals')
    with assert_raises():
        _ = swap_vals(a, -1, 0, 3, 3)
    with assert_raises():
        _ = swap_vals(a, 1, 2, 3, 50)
    b = a
    swap_cols_inplace(b, 0, 1)
    for r in range(40):
        for c in range(50):
            cs = 0 if c == 1 else 1 if c == 0 else c
            assert_equal(a[r, c].re, b[r, cs].re, 'swap_cols_inplace')
            assert_equal(a[r, c].im, b[r, cs].im, 'swap_cols_inplace')
    b = a
    swap_cols_inplace(b, 28, 15)
    for r in range(40):
        for c in range(50):
            cs = 28 if c == 15 else 15 if c == 28 else c
            assert_equal(a[r, c].re, b[r, cs].re, 'swap_cols_inplace')
            assert_equal(a[r, c].im, b[r, cs].im, 'swap_cols_inplace')
    with assert_raises():
        swap_cols_inplace(b, -1, 1)
    with assert_raises():
        swap_cols_inplace(b, 1, 50)
    b = a
    swap_rows_inplace(b, 0, 1)
    for r in range(40):
        for c in range(50):
            rs = 0 if r == 1 else 1 if r == 0 else r
            assert_equal(a[r, c].re, b[rs, c].re, 'swap_rows_inplace')
            assert_equal(a[r, c].im, b[rs, c].im, 'swap_rows_inplace')
    b = a
    swap_rows_inplace(b, 28, 15)
    for r in range(40):
        for c in range(50):
            rs = 28 if r == 15 else 15 if r == 28 else r
            assert_equal(a[r, c].re, b[rs, c].re, 'swap_rows_inplace')
            assert_equal(a[r, c].im, b[rs, c].im, 'swap_rows_inplace')
    with assert_raises():
        swap_rows_inplace(b, -1, 1)
    with assert_raises():
        swap_rows_inplace(b, 1, 40)
    b = a
    swap_vals_inplace(b, 0, 0, 5, 5)
    for r in range(40):
        for c in range(50):
            swapped_val = b[5, 5] if r == 0 and c == 0 else b[0, 0] if r == 5 and c == 5 else b[r, c]
            assert_equal(a[r, c].re, swapped_val.re, 'swap_vals_inplace')
            assert_equal(a[r, c].im, swapped_val.im, 'swap_vals_inplace')
    b = a
    swap_vals_inplace(b, 17, 49, 3, 33)
    for r in range(40):
        for c in range(50):
            swapped_val = b[17, 49] if r == 3 and c == 33 else b[3, 33] if r == 17 and c == 49 else b[r, c]
            assert_equal(a[r, c].re, swapped_val.re, 'swap_vals_inplace')
            assert_equal(a[r, c].im, swapped_val.im, 'swap_vals_inplace')
    with assert_raises():
        swap_vals_inplace(a, -1, 0, 3, 3)
    with assert_raises():
        swap_vals_inplace(a, 1, 2, 3, 50)            

def test_mmax_mmin():
    a = CMatrix[type].arange(50, 50)
    mx, mn = mmax(a), mmin(a)
    assert_equal(mx.re, 50 * 50 - 1, 'mmax')
    assert_equal(mx.im, 0, 'mmax')
    assert_equal(mn.re, 0, 'mmin')
    assert_equal(mn.im, 0, 'mmin')
    a = CMatrix[type].arange(50, 50) + (17, 13)
    a *= (-12, 16)
    mx, mn = mmax(a), mmin(a)
    assert_equal(mx.re, -30400, 'mmax')
    assert_equal(mx.im, 40100, 'mmax')
    assert_equal(mn.re, -412, 'mmin')
    assert_equal(mn.im, 116, 'mmin')

def test_kron():
    a = CMatrix[type].arange(3, 3)
    b = CMatrix[type].arange(2, 2, -1, -1)
    s = CSRCMatrix[type](a)
    p = CSRCMatrix[type](b)
    k = CMatrix[type](6, 6,
        0, 0, -1, -2, -2, -4,
        0, 0, -3, -4, -6, -8,
        -3, -6, -4, -8, -5, -10,
        -9, -12, -12, -16, -15, -20,
        -6, -12, -7, -14, -8, -16,
        -18, -24, -21, -28, -24, -32,
    )
    _assert_matrix_equal(kron_sequential(a, b), k, 'kron_sequential')
    _assert_matrix_equal(kron(a, b), k, 'kron')
    _assert_matrix_equal(sparse_kron(s, p).to_dense(), k)
    a = CMatrix[type].arange(3, 2)
    b = CMatrix[type].arange(2, 2, -1, -1)
    s = CSRCMatrix[type](a)
    p = CSRCMatrix[type](b)
    k = CMatrix[type](6, 4,
        0, 0, -1, -2,
        0, 0, -3, -4,
        -2, -4, -3, -6,
        -6, -8, -9, -12,
        -4, -8, -5, -10,
        -12, -16, -15, -20,
    )
    _assert_matrix_equal(kron_sequential(a, b), k, 'kron_sequential')
    _assert_matrix_equal(kron(a, b), k, 'kron')
    _assert_matrix_equal(sparse_kron(s, p).to_dense(), k)
    c = CMatrix[type].arange(2, 5, -12, 3)
    _assert_matrix_equal(kron(a, b, c), kron(kron(a, b), c), 'kron')
    _assert_matrix_equal(kron(a, b, c), kron_sequential(kron_sequential(a, b), c), 'kron')
    _assert_matrix_equal(kron(c, b, c), kron(kron(c, b), c), 'kron')
    _assert_matrix_equal(kron(a, a, -a), kron_sequential(kron_sequential(a, a), -a), 'kron')
    asp = CSRCMatrix[type](a)
    bsp = CSRCMatrix[type](b)
    csp = CSRCMatrix[type](c)
    _assert_matrix_equal(sparse_kron(asp, bsp, csp).to_dense(), sparse_kron(sparse_kron(asp, bsp), csp).to_dense(), 'sparse_kron')
    _assert_matrix_equal(sparse_kron(asp, bsp, csp).to_dense(), kron_sequential(kron_sequential(asp.to_dense(), bsp.to_dense()), csp.to_dense()), 'sparse_kron')
    _assert_matrix_equal(sparse_kron(csp, bsp, csp).to_dense(), sparse_kron(sparse_kron(csp, bsp), csp).to_dense(), 'sparse_kron')
    _assert_matrix_equal(sparse_kron(asp, asp, -asp).to_dense(), kron_sequential(kron_sequential(asp.to_dense(), asp.to_dense()), -asp.to_dense()), 'sparse_kron')

def test_one_norm():
    a = CMatrix[type].arange(3, 3)
    assert_equal(one_norm(a), 15, 'one_norm')

def test_stack():
    a = CMatrix[type].arange(2, 4)
    b = CMatrix[type].arange(1, 4)
    s = CMatrix[type](3, 4,
        0, 1, 2, 3,
        4, 5, 6, 7,
        0, 1, 2, 3,
    )
    _assert_matrix_equal(vstack(a, b), s, 'vstack')
    a = CMatrix[type].arange(2, 4)
    b = CMatrix[type].arange(2, 3)
    s = CMatrix[type](2, 7,
        0, 1, 2, 3, 0, 1, 2,
        4, 5, 6, 7, 3, 4, 5,
    )
    _assert_matrix_equal(hstack(a, b), s, 'hstack')
    b = CMatrix[type].arange(3, 3)
    with assert_raises():
        _ = vstack(a, b)
    with assert_raises():
        _ = hstack(a, b)

def test_solve():
    a = CMatrix[type](3, 3,
        0, 1, 2,
        1, 1, 1,
        3, 0, -1,
    )
    b = CMatrix[type](3, 2,
        5, 6,
        3, 1,
        8, 8,
    )
    x = CMatrix[type](3, 2,
        5, 6.5,
        -9, -17,
        7, 11.5,
    )
    _assert_matrix_almost_equal(solve(a, b), x, 'solve')

def test_matrix_power():
    a = CMatrix[type].arange(4, 4)
    a[0, 0] = (0.0, 0.5)
    a[1, 0] = (1.0, -0.5)
    a[2, 3] = (0.97, 0.776)
    a[3, 3] = (-1.2, -1.1)
    s = CSRCMatrix(a)
    _assert_matrix_almost_equal(matrix_power(a, 4), a @ a @ a @ a, 'matrix_power', atol=1e-7)
    _assert_matrix_almost_equal(matrix_power(a, -4), (a @ a @ a @ a).inv(), 'matrix_power', atol=1e-7)
    _assert_matrix_almost_equal(matrix_power(a, 0), a.eye_like(), 'matrix_power')
    _assert_matrix_almost_equal(matrix_power(s, 4).to_dense(), a @ a @ a @ a, 'matrix_power', atol=1e-7)
    _assert_matrix_almost_equal(matrix_power(s, 0).to_dense(), a.eye_like(), 'matrix_power')

def test_expm():
    a = CMatrix[type].arange(4, 4)
    a[0, 0] = (0.0, 0.5)
    a[1, 0] = (1.0, -0.5)
    a[2, 3] = (0.97, 0.776)
    a[3, 3] = (-1.2, -1.1)
    r = CMatrix[type](4, 4, 
        (2.37727829e8, 3.66055498e7), (3.33746691e8, 5.36890868e7),
        (3.84288377e8, 6.15537626e7), (1.47880336e8, 2.84619059e7),
        (7.00694930e8, 8.80087548e7), (9.83895168e8, 1.30301443e8),
        (1.13287183e9, 1.49253690e8), (4.36337426e8, 7.14493814e7),
        (7.68077164e8, 1.54388580e8), (1.07796432e9, 2.24225171e8),
        (1.24124774e9, 2.57316196e8), (4.76943358e8, 1.14556219e8),
        (9.90871006e8, 1.08000933e8), (1.39150726e9, 1.61138164e8),
        (1.60218462e9, 1.84439887e8), (6.17420946e8, 9.07435142e7),
    )
    _assert_matrix_almost_equal(expm(a), r, 'expm', rtol=1e-8)
    _assert_matrix_almost_equal(expm(CMatrix[type].zeros(10, 10)), CMatrix[type].eye(10, 10), 'expm')

def test_find_smallest_eigval():
    # m = CMatrix.arange(4, 4, -5)
    # print(m)
    # print(find_smallest_eigval(m))
    pass