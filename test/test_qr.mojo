from math import sin, cos

from testing import assert_equal, assert_raises, assert_almost_equal, assert_true
from ._testing import _assert_matrix_almost_equal, _assert_matrix_equal

from ember import CMatrix, ComplexScalar
from ember import complex_schur, eigvals

alias type = DType.float64

def run_qr_tests():
    print('Running qr tests')
    # test_complex_schur()
    test_eigvals()
    print('All tests passed')

def test_complex_schur():
    A = CMatrix[type].arange(5, 6)
    with assert_raises(contains='must be square'):
        _ = complex_schur(A)

    test_mats = List[CMatrix[type]](
        CMatrix[type].eye(2, 2),
        CMatrix[type].eye(3, 3),
        CMatrix[type].eye(1000, 1000),
        CMatrix[type].arange(10, 10),
        CMatrix[type].arange(11, 11, -162, 3),
        CMatrix[type].arange(11, 11, -162, 3) / 27 + CMatrix[type].arange(11, 11, 48, -1) * ComplexScalar[type].i() / 13,
        CMatrix[type].arange(200, 200, 24565, -2) / 29977 + CMatrix[type].arange(200, 200, -33910, 3) * ComplexScalar[type].i() / 44581,
        CMatrix[type].random(100, 100),
        CMatrix[type].random(100, 100) * 1000,
        CMatrix[type].random(100, 100, mean=100, stdev=100),
        CMatrix[type].random(500, 500),
    )
    for A in test_mats:
        U, Q = complex_schur(A[])
        # Check if U is upper triangular
        for r in range(1, A[].rows):
            for c in range(r):
                assert_almost_equal(U[r, c].re, 0, 'complex_schur')
                assert_almost_equal(U[r, c].im, 0, 'complex_schur')
        assert_true(Q.is_unitary[tol=1e-12]())
        _assert_matrix_almost_equal(A[], Q @ U @ Q.dag(), 'complex_schur')

def test_eigvals():
    m = CMatrix[type].eye(50, 50)
    for e in eigvals(m):
        assert_almost_equal(e[].re, 1, 'eigvals')
        assert_almost_equal(e[].im, 0, 'eigvals')
    m = CMatrix[type].arange(10, 10)
    eigs = eigvals(m)
    assert_equal(len(eigs), 10, 'eigvals')
    for e in eigs:
        if not (e[].is_close[1e-12](0) or e[].is_close[1e-6](511.140380) or e[].is_close[1e-6](-16.1403801)):
            raise Error('eigvals')
    
    m = CMatrix[type].arange(10, 10) / (50 - CMatrix[type].arange(10, 10) * ComplexScalar[type].i())
    eigs = eigvals(m)
    assert_equal(len(eigs), 10, 'eigvals')
    es = List[ComplexScalar[type]](
        (4.00724214, 4.53707538),
        (-0.205403410, -0.151925426),
        (4.95299787e-05, 2.06929101e-03),
        (1.87592383e-05, -7.43379016e-06),
        (-1.16465392e-07, -1.67467462e-07),
        (-1.52420670e-09, 1.27765404e-09),
        (9.68792155e-12, 1.41532953e-11),
    )
    for e in eigs:
        for te in es:
            if e[].is_close[1e-12](0) or e[].is_close[1e-8](te[]):
                break
        else:
            raise Error('eigvals')
