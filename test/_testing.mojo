from testing import assert_almost_equal, assert_equal
from math import isnan

from ember import CMatrix, ComplexScalar

fn _assert_almost_equal[type: DType, size: Int](lhs: SIMD[type, size], rhs: SIMD[type, size], msg: String = String(''), *, atol: Float64 = 1e-10, rtol: Float64 = 0) raises:
    return assert_almost_equal[type, size](lhs, rhs, msg, atol=atol, rtol=rtol)

fn _assert_matrix_equal[type: DType](lhs: CMatrix[type], rhs: CMatrix[type], msg: String = '') raises:
    assert_equal(lhs.rows, rhs.rows, msg + ' rows')
    assert_equal(lhs.cols, rhs.cols, msg + ' cols')
    for r in range(lhs.rows):
        for c in range(lhs.cols):
            if not isnan(lhs[r, c].re) or not isnan(rhs[r, c].re):
                assert_equal(lhs[r, c].re, rhs[r, c].re, msg)
            if not isnan(lhs[r, c].im) or not isnan(rhs[r, c].im):
                assert_equal(lhs[r, c].im, rhs[r, c].im, msg)

fn _assert_matrix_equal[type: DType](lhs: CMatrix[type], rhs: ComplexScalar[type], msg: String = '') raises:
    for r in range(lhs.rows):
        for c in range(lhs.cols):
            if not isnan(lhs[r, c].re) or not isnan(rhs.re):
                assert_equal(lhs[r, c].re, rhs.re, msg)
            if not isnan(lhs[r, c].im) or not isnan(rhs.im):
                assert_equal(lhs[r, c].im, rhs.im, msg)

fn _assert_matrix_almost_equal[type: DType](lhs: CMatrix[type], rhs: CMatrix[type], msg: String = '', *, atol: Float64 = 1e-10, rtol: Float64 = 0) raises:
    assert_equal(lhs.rows, rhs.rows, msg + ' rows')
    assert_equal(lhs.cols, rhs.cols, msg + ' cols')
    for r in range(lhs.rows):
        for c in range(lhs.cols):
            if not isnan(lhs[r, c].re) or not isnan(rhs[r, c].re):
                assert_almost_equal(lhs[r, c].re, rhs[r, c].re, msg, atol=atol, rtol=rtol)
            if not isnan(lhs[r, c].im) or not isnan(rhs[r, c].im):
                assert_almost_equal(lhs[r, c].im, rhs[r, c].im, msg, atol=atol, rtol=rtol)

fn _assert_matrix_almost_equal[type: DType](lhs: CMatrix[type], rhs: ComplexScalar[type], msg: String = '', *, atol: Float64 = 1e-10, rtol: Float64 = 0) raises:
    for r in range(lhs.rows):
        for c in range(lhs.cols):
            if not isnan(lhs[r, c].re) or not isnan(rhs.re):
                assert_almost_equal(lhs[r, c].re, rhs.re, msg, atol=atol, rtol=rtol)
            if not isnan(lhs[r, c].im) or not isnan(rhs.im):
                assert_almost_equal(lhs[r, c].im, rhs.im, msg, atol=atol, rtol=rtol)