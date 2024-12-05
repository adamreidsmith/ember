from testing import assert_true, assert_false, assert_equal, assert_raises
import math

from ._testing import _assert_almost_equal, _assert_matrix_equal, _assert_matrix_almost_equal
from src import CMatrix, ComplexScalar

alias type = DType.float64

def run_cmatrix_tests():
    print('Running CMatrix tests')
    test_init()
    test_properties()
    test_item_access()
    test_arithmetic()
    test_matmul()
    test_iarithmetic()
    test_other()
    test_shape()
    test_fill()
    test_compare()
    test_static_constructors()
    print('All tests passed')

def test_init():
    m = CMatrix[type](2, 4, fill_zeros=True)
    _assert_matrix_equal(m, ComplexScalar[type](0), 'init')
    m = CMatrix[type](2, 2, ComplexScalar[type](2, 2), ComplexScalar[type](1, 3), ComplexScalar[type](-1, 4))
    assert_equal(m[0, 0].re, 2, 'init')
    assert_equal(m[0, 0].im, 2, 'init')
    assert_equal(m[0, 1].re, 1, 'init')
    assert_equal(m[0, 1].im, 3, 'init')
    assert_equal(m[1, 0].re, -1, 'init')
    assert_equal(m[1, 0].im, 4, 'init')
    assert_equal(m[1, 1].re, 0, 'init')
    assert_equal(m[1, 1].im, 0, 'init')
    m = CMatrix[type](List[ComplexScalar[type], True](ComplexScalar[type](2, 2), ComplexScalar[type](1, 3), ComplexScalar[type](-1, 4)))
    assert_equal(m.rows, 1, 'init')
    assert_equal(m.cols, 3, 'init')
    assert_equal(m[0, 0].re, 2, 'init')
    assert_equal(m[0, 0].im, 2, 'init')
    assert_equal(m[0, 1].re, 1, 'init')
    assert_equal(m[0, 1].im, 3, 'init')
    assert_equal(m[0, 2].re, -1, 'init')
    with assert_raises():
        m = CMatrix[type](List[List[ComplexScalar[type], True]](
            List[ComplexScalar[type], True](ComplexScalar[type]()),
            List[ComplexScalar[type], True](ComplexScalar[type](), ComplexScalar[type]()),
        ))
    m = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 0), ComplexScalar[type](1, 1)),
        List[ComplexScalar[type], True](ComplexScalar[type](2, 2), ComplexScalar[type](3, 3)),
    ))
    assert_equal(m.rows, 2, 'init')
    assert_equal(m.cols, 2, 'init')
    assert_equal(m[0, 0].re, 0, 'init')
    assert_equal(m[0, 0].im, 0, 'init')
    assert_equal(m[0, 1].re, 1, 'init')
    assert_equal(m[0, 1].im, 1, 'init')
    assert_equal(m[1, 0].re, 2, 'init')
    assert_equal(m[1, 0].im, 2, 'init')
    assert_equal(m[1, 1].re, 3, 'init')
    assert_equal(m[1, 1].im, 3, 'init')
    m2 = m
    _assert_matrix_equal(m, m2, 'copyinit')

def test_properties():
    m1 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 0), ComplexScalar[type](1, 1)),
        List[ComplexScalar[type], True](ComplexScalar[type](2, 2), ComplexScalar[type](3, 3)),
        List[ComplexScalar[type], True](ComplexScalar[type](4, 4), ComplexScalar[type](5, 5)),
    ))
    m2 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](math.cos(Scalar[type](3))), ComplexScalar[type](-math.sin(Scalar[type](3)))),
        List[ComplexScalar[type], True](ComplexScalar[type](math.sin(Scalar[type](3))), ComplexScalar[type](math.cos(Scalar[type](3)))),
    ))
    m3 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](1, 0), ComplexScalar[type](2, -2)),
        List[ComplexScalar[type], True](ComplexScalar[type](2, 2), ComplexScalar[type](-1, 0)),
    ))
    assert_equal(len(m1), 6, 'len')
    assert_equal(m1.shape()[0], 3, 'shape')
    assert_equal(m1.shape()[1], 2, 'shape')
    assert_false(m1.is_square(), 'is_square')
    assert_true(m2.is_square(), 'is_square')
    assert_true(m2.is_unitary(), 'is_unitary')
    assert_false(m3.is_unitary(), 'is_unitary')
    assert_true(m3.is_hermitian(), 'is_hermitian')
    assert_false(m2.is_hermitian(), 'is_hermitian')

def test_item_access():
    m1 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 0), ComplexScalar[type](0, 1)),
        List[ComplexScalar[type], True](ComplexScalar[type](1, 0), ComplexScalar[type](1, 1)),
        List[ComplexScalar[type], True](ComplexScalar[type](2, 0), ComplexScalar[type](2, 1)),
    ))
    m2 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 0), ComplexScalar[type](1, 2)),
        List[ComplexScalar[type], True](ComplexScalar[type](2, 4), ComplexScalar[type](3, 6)),
    ))
    for r in range(3):
        for c in range(2):
            assert_equal(m1[r, c].re, r, 'getitem')
            assert_equal(m1[r, c].im, c, 'getitem')
            assert_equal(m1.load_crd[1](r, c).re, r, 'getitem_noraise')
            assert_equal(m1.load_crd[1](r, c).im, c, 'getitem_noraise')
    for idx in range(4):
        assert_equal(m2[idx].re, idx, 'getitem')
        assert_equal(m2[idx].im, 2 * idx, 'getitem')
        assert_equal(m2.load_idx[1](idx).re, idx, 'getitem_noraise')
        assert_equal(m2.load_idx[1](idx).im, 2 * idx, 'getitem_noraise')
    for r in range(3):
        for c in range(2):
            m1[r, c] = ComplexScalar[type](-r, -c)
    for r in range(3):
        for c in range(2):
            assert_equal(m1[r, c].re, -r, 'setitem')
            assert_equal(m1[r, c].im, -c, 'setitem')
    for r in range(3):
        for c in range(2):
            m1.store_crd[1](r, c, ComplexScalar[type](2 * r, 2 * c))
    for r in range(3):
        for c in range(2):
            assert_equal(m1[r, c].re, 2 * r, 'setitem_noraise')
            assert_equal(m1[r, c].im, 2 * c, 'setitem_noraise')
    for idx in range(4):
        m2[idx] = ComplexScalar[type](-idx, -2 * idx)
    for idx in range(4):
        assert_equal(m2[idx].re, -idx, 'setitem')
        assert_equal(m2[idx].im, -2 * idx, 'setitem')
    for idx in range(4):
        m2.store_idx[1](idx, ComplexScalar[type](idx, 3 * idx))
    for idx in range(4):
        assert_equal(m2[idx].re, idx, 'setitem_noraise')
        assert_equal(m2[idx].im, 3 * idx, 'setitem_noraise')

def test_arithmetic():
    m1 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 0), ComplexScalar[type](1, 1)),
        List[ComplexScalar[type], True](ComplexScalar[type](2, 2), ComplexScalar[type](3, 3)),
    ))
    m2 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](1, 1), ComplexScalar[type](1, 2)),
        List[ComplexScalar[type], True](ComplexScalar[type](2, 1), ComplexScalar[type](2, 2)),
    ))
    m3 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](1, 1), ComplexScalar[type](1, 2)),
    ))
    s = ComplexScalar[type]((1.4, -2.7))
    for idx in range(4):
        assert_equal((-m1)[idx].re, -(m1[idx]).re, 'neg')
        assert_equal((-m1)[idx].im, -(m1[idx]).im, 'neg')
        assert_equal((m1 + m2)[idx].re, (m1[idx] + m2[idx]).re, 'add')
        assert_equal((m1 + m2)[idx].im, (m1[idx] + m2[idx]).im, 'add')
        assert_equal((m1 + s)[idx].re, (m1[idx] + s).re, 'add')
        assert_equal((m1 + s)[idx].im, (m1[idx] + s).im, 'add')
        assert_equal((m1 - m2)[idx].re, (m1[idx] - m2[idx]).re, 'sub')
        assert_equal((m1 - m2)[idx].im, (m1[idx] - m2[idx]).im, 'sub')
        assert_equal((m1 - s)[idx].re, (m1[idx] - s).re, 'sub')
        assert_equal((m1 - s)[idx].im, (m1[idx] - s).im, 'sub')
        assert_equal((m1 * m2)[idx].re, (m1[idx] * m2[idx]).re, 'mul')
        assert_equal((m1 * m2)[idx].im, (m1[idx] * m2[idx]).im, 'mul')
        _assert_almost_equal((m1 * s)[idx].re, (m1[idx] * s).re, 'mul')
        _assert_almost_equal((m1 * s)[idx].im, (m1[idx] * s).im, 'mul')
        assert_equal((m1 / m2)[idx].re, (m1[idx] / m2[idx]).re, 'truediv')
        assert_equal((m1 / m2)[idx].im, (m1[idx] / m2[idx]).im, 'truediv')
        _assert_almost_equal((m1 / s)[idx].re, (m1[idx] / s).re, 'truediv')
        _assert_almost_equal((m1 / s)[idx].im, (m1[idx] / s).im, 'truediv')
        assert_equal((m1 // m2)[idx].re, (m1[idx] // m2[idx]).re, 'floordiv')
        assert_equal((m1 // m2)[idx].im, (m1[idx] // m2[idx]).im, 'floordiv')
        _assert_almost_equal((m1 // s)[idx].re, (m1[idx] // s).re, 'floordiv')
        _assert_almost_equal((m1 // s)[idx].im, (m1[idx] // s).im, 'floordiv')
        assert_equal((m1 % m2)[idx].re, (m1[idx] % m2[idx]).re, 'mod')
        assert_equal((m1 % m2)[idx].im, (m1[idx] % m2[idx]).im, 'mod')
        _assert_almost_equal((m1 % s)[idx].re, (m1[idx] % s).re, 'mod')
        _assert_almost_equal((m1 % s)[idx].im, (m1[idx] % s).im, 'mod')
    with assert_raises():
        _ = m1 + m3
    with assert_raises():
        _ = m1 - m3
    with assert_raises():
        _ = m1 * m3
    with assert_raises():
        _ = m1 / m3
    with assert_raises():
        _ = m1 // m3
    with assert_raises():
        _ = m1 % m3

def test_matmul():
    m1 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 0), ComplexScalar[type](1, 1), ComplexScalar[type](0.5, 1.5)),
        List[ComplexScalar[type], True](ComplexScalar[type](2, 2), ComplexScalar[type](3, 3), ComplexScalar[type](-0.5, -1.5)),
    ))
    m2 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](1, 1), ComplexScalar[type](1, 2)),
        List[ComplexScalar[type], True](ComplexScalar[type](2, 1), ComplexScalar[type](2, 2)),
        List[ComplexScalar[type], True](ComplexScalar[type](3.3, 1.3), ComplexScalar[type](-4.4, -2.2)),
    ))
    m1m2 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0.7, 8.6), ComplexScalar[type](1.1, -3.7)),
        List[ComplexScalar[type], True](ComplexScalar[type](3.3, 7.4), ComplexScalar[type](-3.1, 25.7)),
    ))
    _assert_matrix_almost_equal(m1 @ m2, m1m2, 'matmul')
    with assert_raises():
        _ = m1 @ m1m2
    m1 @= m2
    _assert_matrix_almost_equal(m1, m1m2, 'imatmul')
    with assert_raises(contains='doesnt raise'):
        _ = m1 @ m1
        raise 'doesnt raise'
    mlarge = CMatrix[type](200, 200, fill_zeros=False)
    for idx in range(200**2):
        mlarge[idx] = ComplexScalar[type](hash(idx) % idx, hash(2 * idx) % idx)
    ilarge = mlarge.eye_like()
    _assert_matrix_almost_equal(mlarge @ ilarge, mlarge, 'matmul')
    _assert_matrix_almost_equal(ilarge @ mlarge, mlarge, 'matmul')

def test_iarithmetic():
    m1 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 0), ComplexScalar[type](1, 1)),
        List[ComplexScalar[type], True](ComplexScalar[type](2, 2), ComplexScalar[type](3, 3)),
    ))
    m2 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](1, 1), ComplexScalar[type](1, 2)),
        List[ComplexScalar[type], True](ComplexScalar[type](2, 1), ComplexScalar[type](2, 2)),
    ))
    m3 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](1, 1), ComplexScalar[type](1, 2)),
    ))
    s = ComplexScalar[type]((1.4, -2.7))
    mi = m1
    mi += m2
    _assert_matrix_equal(mi, m1 + m2, 'iadd')
    mi += s
    _assert_matrix_equal(mi, m1 + m2 + s, 'iadd')
    mi = m1
    mi -= m2
    _assert_matrix_equal(mi, m1 - m2, 'isub')
    mi -= s
    _assert_matrix_equal(mi, m1 - m2 - s, 'isub')
    mi = m1
    mi *= m2
    _assert_matrix_equal(mi, m1 * m2, 'imul')
    mi *= s
    _assert_matrix_equal(mi, m1 * m2 * s, 'imul')
    mi = m1
    mi /= m2
    _assert_matrix_equal(mi, m1 / m2, 'itruediv')
    mi /= s
    _assert_matrix_equal(mi, m1 / m2 / s, 'itruediv')
    mi = m1
    mi //= m2
    _assert_matrix_equal(mi, m1 // m2, 'ifloordiv')
    mi //= s
    _assert_matrix_equal(mi, m1 // m2 // s, 'ifloordiv')
    mi = m1
    mi %= m2
    _assert_matrix_equal(mi, m1 % m2, 'imod')
    mi %= s
    _assert_matrix_equal(mi, m1 % m2 % s, 'imod')
    with assert_raises():
        m1 += m3
    with assert_raises():
        m1 -= m3
    with assert_raises():
        m1 *= m3
    with assert_raises():
        m1 /= m3
    with assert_raises():
        m1 //= m3
    with assert_raises():
        m1 %= m3

def test_other():
    m1 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 2), ComplexScalar[type](-2, 3)),
        List[ComplexScalar[type], True](ComplexScalar[type](1, -1), ComplexScalar[type](3, 2)),
    ))
    m1abs = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](abs(ComplexScalar[type](0, 2)), abs(ComplexScalar[type](-2, 3))),
        List[ComplexScalar[type], True](abs(ComplexScalar[type](1, -1)), abs(ComplexScalar[type](3, 2))),
    ))
    _assert_matrix_almost_equal(abs(m1), m1abs, 'abs')
    m1conj = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 2).conj(), ComplexScalar[type](-2, 3).conj()),
        List[ComplexScalar[type], True](ComplexScalar[type](1, -1).conj(), ComplexScalar[type](3, 2).conj()),
    ))
    _assert_matrix_almost_equal(m1.conj(), m1conj, 'conj')
    m1conj.iconj()
    _assert_matrix_almost_equal(m1, m1conj, 'iconj')
    m1dag = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 2).conj(), ComplexScalar[type](1, -1).conj()),
        List[ComplexScalar[type], True](ComplexScalar[type](-2, 3).conj(), ComplexScalar[type](3, 2).conj()),
    ))
    _assert_matrix_almost_equal(m1.dag(), m1dag, 'dag')
    assert_equal(m1.sum().re, 2, 'sum')
    assert_equal(m1.sum().im, 6, 'sum')
    m1echelon = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 2), ComplexScalar[type](-2, 3)),
        List[ComplexScalar[type], True](ComplexScalar[type](0, 0), ComplexScalar[type](0.5, 2.5)),
    ))
    _assert_matrix_equal(m1.echelon(), m1echelon, 'echelon')
    _assert_almost_equal(m1.det().re, -5, 'det')
    _assert_almost_equal(m1.det().im, 1, 'det')
    _assert_almost_equal(m1.echelon().det().re, m1.det().re, 'det')
    _assert_almost_equal(m1.echelon().det().im, m1.det().im, 'det')
    mlarge = CMatrix[type](100, 100, fill_zeros=False)
    for idx in range(100**2):
        mlarge[idx] = ComplexScalar[type](hash(idx) % 10, hash(2 * idx) % 10)
    mlargeechelon = mlarge.echelon()
    for r in range(100):
        for c in range(100):
            if r > c:
                assert_equal(mlargeechelon[r, c].re, 0, 'echelon')
                assert_equal(mlargeechelon[r, c].im, 0, 'echelon')
    _assert_almost_equal(mlarge.det().re, mlargeechelon.det().re, 'det', rtol=1e-8)
    _assert_almost_equal(mlarge.det().im, mlargeechelon.det().im, 'det', rtol=1e-8)
    mlargeinv = mlarge.inv()
    _assert_matrix_almost_equal(mlargeinv @ mlarge, mlarge.eye_like(), 'inv')
    _assert_matrix_almost_equal(mlarge @ mlargeinv, mlarge.eye_like(), 'inv')
    _assert_almost_equal(m1.frobenius_norm(), 5.656854249492381, 'frobenius_norm')

def test_shape():
    m1 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 2), ComplexScalar[type](-2, 3)),
        List[ComplexScalar[type], True](ComplexScalar[type](1, -1), ComplexScalar[type](3, 2)),
        List[ComplexScalar[type], True](ComplexScalar[type](-4, 1), ComplexScalar[type](7, 6)),
    ))
    m1reshape = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 2), ComplexScalar[type](-2, 3), ComplexScalar[type](1, -1)),
        List[ComplexScalar[type], True](ComplexScalar[type](3, 2), ComplexScalar[type](-4, 1), ComplexScalar[type](7, 6)),
    ))
    assert_equal(m1.reshape(2, 3).rows, 2, 'reshape')
    assert_equal(m1.reshape(2, 3).cols, 3, 'reshape')
    _assert_matrix_equal(m1.reshape(2, 3), m1reshape, 'reshape')
    mi = m1
    mi.ireshape(2, 3)
    _assert_matrix_equal(mi, m1reshape, 'ireshape')
    with assert_raises():
        _ = m1.reshape(1, 1)
    with assert_raises():
        m1.ireshape(1, 1)
    assert_equal(m1.flatten_to_row().rows, 1, 'flatten_to_row')
    assert_equal(m1.flatten_to_row().cols, 6, 'flatten_to_row')
    assert_equal(m1.flatten_to_column().rows, 6, 'flatten_to_column')
    assert_equal(m1.flatten_to_column().cols, 1, 'flatten_to_column')
    mi = m1
    mi.iflatten_to_row()
    _assert_matrix_equal(mi, m1.flatten_to_row(), 'iflatten_to_row')
    mi.iflatten_to_column()
    _assert_matrix_equal(mi, m1.flatten_to_column(), 'iflatten_to_column')
    m1t = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 2), ComplexScalar[type](1, -1), ComplexScalar[type](-4, 1)),
        List[ComplexScalar[type], True](ComplexScalar[type](-2, 3), ComplexScalar[type](3, 2), ComplexScalar[type](7, 6)),
    ))
    _assert_matrix_equal(m1.transpose(), m1t, 'transpose')
    _assert_matrix_equal(m1.transpose().transpose(), m1, 'transpose')

def test_fill():
    m1 = CMatrix[type](100, 100, fill_zeros=False)
    for idx in range(100**2):
        m1[idx] = ComplexScalar[type](hash(idx) % 10, hash(2 * idx) % 10)
    s = ComplexScalar[type](19, -1345)
    m1.fill(s)
    _assert_matrix_equal(m1, s, 'fill')
    m1.fill_zero()
    _assert_matrix_equal(m1, 0, 'fill_zero')
    m1.fill_one()
    _assert_matrix_equal(m1, 1, 'fill_one')
    m1.fill_i()
    _assert_matrix_equal(m1, (0, 1), 'fill_i')
    m1.fill_zero()
    m1.fill_diag(s)
    _assert_matrix_equal(m1, m1.eye_like() * s, 'fill_diag')
    m1.fill_range()
    for idx in range(m1.size):
        assert_equal(m1[idx].re, idx, 'fill_range')
        assert_equal(m1[idx].im, 0, 'fill_range')
    m1.fill_range(12345)
    for idx in range(m1.size):
        assert_equal(m1[idx].re, idx + 12345, 'fill_range')
        assert_equal(m1[idx].im, 0, 'fill_range')
    m1.fill_range(-12345, 17)
    for idx in range(m1.size):
        assert_equal(m1[idx].re, 17 * idx - 12345, 'fill_range')
        assert_equal(m1[idx].im, 0, 'fill_range')
    _assert_matrix_equal(m1.ones_like(), 1, 'ones_like')
    _assert_matrix_equal(m1.zeros_like(), 0, 'zeros_like')
    _assert_matrix_equal(m1.i_like(), (0, 1), 'i_like')
    for r in range(m1.rows):
        for c in range(m1.cols):
            if r == c:
                assert_equal(m1.eye_like()[r, c].re, 1, 'eye_like')
                assert_equal(m1.eye_like()[r, c].im, 0, 'eye_like')
            else:
                assert_equal(m1.eye_like()[r, c].re, 0, 'eye_like')
                assert_equal(m1.eye_like()[r, c].im, 0, 'eye_like')

def test_compare():
    m1 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 0), ComplexScalar[type](1, 1)),
        List[ComplexScalar[type], True](ComplexScalar[type](2, 2), ComplexScalar[type](3, 3)),
    ))
    m2 = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](0, 0), ComplexScalar[type](1, 2)),
        List[ComplexScalar[type], True](ComplexScalar[type](2, 1), ComplexScalar[type](3, 3)),
    ))
    assert_true(m1.matrix_equals(m1), 'matrix_equals')
    assert_true(m2.matrix_equals(m2), 'matrix_equals')
    assert_false(m2.matrix_equals(m1), 'matrix_equals')
    assert_false(m1.matrix_equals(m2), 'matrix_equals')
    assert_true(m1.is_close(m1), 'is_close')
    assert_true(m2.is_close(m2), 'is_close')
    assert_false(m2.is_close(m1), 'is_close')
    assert_false(m1.is_close(m2), 'is_close')
    mi = m1
    mi[0, 0] += 0.0000001
    assert_false(m1.is_close(mi), 'is_close')
    mi = m1
    mi[0, 0] += 0.0000000000001
    assert_true(m1.is_close(mi), 'is_close')
    assert_true(m1.ones_like().is_close(1), 'is_close')
    assert_false(m1.ones_like().is_close(1.0000001), 'is_close')
    _assert_matrix_equal(m1 == m1, m1.ones_like(), 'eq')
    _assert_matrix_equal(m1 == m2, m1.eye_like(), 'eq')
    ne = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](), ComplexScalar[type](1)),
        List[ComplexScalar[type], True](ComplexScalar[type](1), ComplexScalar[type]()),
    ))
    _assert_matrix_equal(m1 != m1, m1.zeros_like(), 'ne')
    _assert_matrix_equal(m1 != m2, ne, 'ne')
    gt = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](), ComplexScalar[type]()),
        List[ComplexScalar[type], True](ComplexScalar[type](1), ComplexScalar[type]()),
    ))
    _assert_matrix_equal(m1 > m1, m1.zeros_like(), 'gt')
    _assert_matrix_equal(m1 > m2, gt, 'gt')
    ge = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](1), ComplexScalar[type]()),
        List[ComplexScalar[type], True](ComplexScalar[type](1), ComplexScalar[type](1)),
    ))
    _assert_matrix_equal(m1 >= m1, m1.ones_like(), 'ge')
    _assert_matrix_equal(m1 >= m2, ge, 'ge')

    lt = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](), ComplexScalar[type](1)),
        List[ComplexScalar[type], True](ComplexScalar[type](), ComplexScalar[type]()),
    ))
    _assert_matrix_equal(m1 < m1, m1.zeros_like(), 'lt')
    _assert_matrix_equal(m1 < m2, lt, 'lt')
    le = CMatrix[type](List[List[ComplexScalar[type], True]](
        List[ComplexScalar[type], True](ComplexScalar[type](1), ComplexScalar[type](1)),
        List[ComplexScalar[type], True](ComplexScalar[type](), ComplexScalar[type](1)),
    ))
    _assert_matrix_equal(m1 <= m1, m1.ones_like(), 'le')
    _assert_matrix_equal(m1 <= m2, le, 'le')

def test_static_constructors():
    for r in range(1, 50):
        for c in range(1, 50):
            _assert_matrix_equal(CMatrix[type].eye(r, c), CMatrix[type](r, c).eye_like(), 'eye')
            _assert_matrix_equal(CMatrix[type].zeros(r, c), CMatrix[type](r, c).zeros_like(), 'zeros')
            _assert_matrix_equal(CMatrix[type].ones(r, c), CMatrix[type](r, c).ones_like(), 'ones')
            _assert_matrix_equal(CMatrix[type].i(r, c), CMatrix[type](r, c).i_like(), 'i')
            _assert_matrix_equal(CMatrix[type].arange(r, c), CMatrix[type](r, c).range_like(), 'i')
            _assert_matrix_equal(CMatrix[type].arange(r, c), CMatrix[type](r, c).range_like(), 'i')
            for start in range(100, -100, -23):
                _assert_matrix_equal(CMatrix[type].arange(r, c, start), CMatrix[type](r, c).range_like(start), 'i')
                _assert_matrix_equal(CMatrix[type].arange(r, c, start), CMatrix[type](r, c).range_like(start), 'i')
                for step in range(-10, 10, 3):
                    _assert_matrix_equal(CMatrix[type].arange(r, c, start, step), CMatrix[type](r, c).range_like(start, step), 'i')
                    _assert_matrix_equal(CMatrix[type].arange(r, c, start, step), CMatrix[type](r, c).range_like(start, step), 'i')
