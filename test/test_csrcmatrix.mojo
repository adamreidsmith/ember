from testing import assert_true, assert_false, assert_equal, assert_raises
import math
from collections import Dict

from ._testing import _assert_matrix_equal, _assert_matrix_almost_equal
from src import CSRCMatrix, CMatrix, ComplexScalar

alias type = DType.float64

def run_csrcmatrix_tests():
    print('Running CSRCMatrix tests')
    test_init()
    test_item_access()
    test_other()
    test_arithmetic()
    test_matmul()
    test_shape()
    test_compare()
    test_fill()
    test_static_constructors()
    print('All tests passed')

def test_init():
    m1 = CSRCMatrix[type](4, 6)
    assert_equal(len(m1.v), 0, 'init')
    assert_equal(len(m1.col_idx), 0, 'init')
    assert_equal(len(m1.row_idx), 5, 'init')
    for i in range(len(m1.row_idx)):
        assert_equal(m1.row_idx[i], 0, 'init')
    v1 = List[ComplexScalar[type], True](30, 10, 20, 40, 80, 50, 60, 70)
    c1 = List[Int, True](1, 0, 1, 3, 5, 2, 3, 4)   
    r1 = List[Int, True](1, 0, 0, 1, 3, 2, 1, 2)
    m1 = CSRCMatrix[type](4, 6, r1, c1, v1)
    vt = List[ComplexScalar[type], True](10, 20, 30, 100, 50, 70, 80)
    ct = List[Int, True](0, 1, 1, 3, 2, 4, 5)   
    rt = List[Int, True](0, 2, 4, 6, 7)
    assert_equal(len(m1.row_idx), m1.rows + 1, 'init')
    assert_equal(len(m1.v), len(m1.col_idx), 'init')
    assert_equal(len(m1.v), len(vt), 'init')
    assert_equal(len(m1.col_idx), len(ct), 'init')
    assert_equal(len(m1.row_idx), len(rt), 'init')
    for i in range(len(vt)):
        assert_equal(m1.v[i].re, vt[i].re, 'init')
        assert_equal(m1.v[i].im, vt[i].im, 'init')
        assert_equal(m1.col_idx[i], ct[i], 'init')
    for i in range(len(rt)):
        assert_equal(m1.row_idx[i], rt[i], 'init')
    with assert_raises():
        m1 = CSRCMatrix[type](4, 6, r1, c1 + List[Int, True](0), v1)
    with assert_raises():
        m1 = CSRCMatrix[type](4, 6, r1, c1[:-1] + List[Int, True](-1), v1)
    with assert_raises():
        m1 = CSRCMatrix[type](4, 6, r1[:-1] + List[Int, True](10), c1, v1)
    m1 = CSRCMatrix[type](4, 6,
        (1, 1, ComplexScalar[type](30)),
        (0, 0, ComplexScalar[type](10)),
        (0, 1, ComplexScalar[type](20)),
        (1, 3, ComplexScalar[type](40)),
        (3, 5, ComplexScalar[type](80)),
        (2, 2, ComplexScalar[type](50)),
        (1, 3, ComplexScalar[type](60)),
        (2, 4, ComplexScalar[type](70)),
    )
    vt = List[ComplexScalar[type], True](10, 20, 30, 100, 50, 70, 80)
    ct = List[Int, True](0, 1, 1, 3, 2, 4, 5)   
    rt = List[Int, True](0, 2, 4, 6, 7)
    assert_equal(len(m1.row_idx), m1.rows + 1, 'init')
    assert_equal(len(m1.v), len(m1.col_idx), 'init')
    assert_equal(len(m1.v), len(vt), 'init')
    assert_equal(len(m1.col_idx), len(ct), 'init')
    assert_equal(len(m1.row_idx), len(rt), 'init')
    for i in range(len(vt)):
        assert_equal(m1.v[i].re, vt[i].re, 'init')
        assert_equal(m1.v[i].im, vt[i].im, 'init')
        assert_equal(m1.col_idx[i], ct[i], 'init')
    for i in range(len(rt)):
        assert_equal(m1.row_idx[i], rt[i], 'init')
    with assert_raises():
        m1 = CSRCMatrix[type](4, 6,
            (1, 1, ComplexScalar[type](30)),
            (0, 0, ComplexScalar[type](10)),
            (0, 1, ComplexScalar[type](20)),
            (1, 3, ComplexScalar[type](40)),
            (3, 5, ComplexScalar[type](80)),
            (2, 2, ComplexScalar[type](50)),
            (1, 3, ComplexScalar[type](60)),
            (2, 7, ComplexScalar[type](70)),
        )
    with assert_raises():
        m1 = CSRCMatrix[type](4, 6,
            (1, 1, ComplexScalar[type](30)),
            (-1, 0, ComplexScalar[type](10)),
            (0, 1, ComplexScalar[type](20)),
            (1, 3, ComplexScalar[type](40)),
            (3, 5, ComplexScalar[type](80)),
            (2, 2, ComplexScalar[type](50)),
            (1, 3, ComplexScalar[type](60)),
            (2, 4, ComplexScalar[type](70)),
        )
    v1 = List[ComplexScalar[type], True](30, 10, 20, 40, 80, 50, 60, 70)
    c1 = List[Int, True](1, 0, 1, 3, 5, 2, 3, 4)   
    r1 = List[Int, True](1, 0, 0, 1, 3, 2, 1, 2)
    m1 = CSRCMatrix[type](4, 6, r1, c1, v1)
    data = Dict[Int, Dict[Int, ComplexScalar[type]]]()
    data[0] = Dict[Int, ComplexScalar[type]]()
    data[1] = Dict[Int, ComplexScalar[type]]()
    data[2] = Dict[Int, ComplexScalar[type]]()
    data[3] = Dict[Int, ComplexScalar[type]]()
    data[0][0] = 10
    data[0][1] = 20
    data[1][1] = 30
    data[1][3] = 100
    data[2][2] = 50
    data[2][4] = 70
    data[3][5] = 80
    m1 = CSRCMatrix[type](4, 6, data)
    vt = List[ComplexScalar[type], True](10, 20, 30, 100, 50, 70, 80)
    ct = List[Int, True](0, 1, 1, 3, 2, 4, 5)   
    rt = List[Int, True](0, 2, 4, 6, 7)
    assert_equal(len(m1.v), len(vt), 'init')
    assert_equal(len(m1.col_idx), len(ct), 'init')
    assert_equal(len(m1.row_idx), len(rt), 'init')
    for i in range(len(vt)):
        assert_equal(m1.v[i].re, vt[i].re, 'init')
        assert_equal(m1.v[i].im, vt[i].im, 'init')
        assert_equal(m1.col_idx[i], ct[i], 'init')
    for i in range(len(rt)):
        assert_equal(m1.row_idx[i], rt[i], 'init')
    cm = CMatrix[type](153, 82)
    for r in range(cm.rows):
        for c in range(cm.cols):
            cm[r, c] = ComplexScalar[type]((int(hash(r + c)) % 83) - 41, (int(hash(r*c + r)) % 667) - 300)
    _assert_matrix_equal(CSRCMatrix[type](cm).to_dense(), cm, 'init')
    _assert_matrix_equal(CSRCMatrix[type](cm.zeros_like()).to_dense(), cm.zeros_like(), 'init')
    _assert_matrix_equal(CSRCMatrix[type](cm.eye_like()).to_dense(), cm.eye_like(), 'init')
    _assert_matrix_equal(CSRCMatrix[type](CMatrix[type](0, 0)).to_dense(), CMatrix[type](0, 0), 'init')

def test_item_access():
    m1 = CSRCMatrix[type](4, 6,
        (1, 1, ComplexScalar[type](30, 5)),
        (0, 0, ComplexScalar[type](10)),
        (0, 1, ComplexScalar[type](20)),
        (1, 3, ComplexScalar[type](40)),
        (3, 5, ComplexScalar[type](80)),
        (2, 2, ComplexScalar[type](50)),
        (1, 3, ComplexScalar[type](60)),
        (2, 4, ComplexScalar[type](70)),
    )
    assert_equal(m1[1, 1].re, 30, 'getitem')
    assert_equal(m1[1, 1].im, 5, 'getitem')
    assert_equal(m1[0, 0].re, 10, 'getitem')
    assert_equal(m1[0, 0].im, 0, 'getitem')
    assert_equal(m1[1, 3].re, 100, 'getitem')
    assert_equal(m1[1, 3].im, 0, 'getitem')
    assert_equal(m1[3, 2].re, 0, 'getitem')
    assert_equal(m1[3, 2].im, 0, 'getitem')
    assert_equal(m1[0, 5].re, 0, 'getitem')
    assert_equal(m1[0, 5].im, 0, 'getitem')
    with assert_raises():
        _ = m1[-1, 0]
    with assert_raises():
        _ = m1[0, 7]
    assert_equal(m1._getitem_noraise(1, 1).re, 30, 'getitem_noraise')
    assert_equal(m1._getitem_noraise(1, 1).im, 5, 'getitem_noraise')
    assert_equal(m1._getitem_noraise(0, 0).re, 10, 'getitem_noraise')
    assert_equal(m1._getitem_noraise(0, 0).im, 0, 'getitem_noraise')
    assert_equal(m1._getitem_noraise(1, 3).re, 100, 'getitem_noraise')
    assert_equal(m1._getitem_noraise(1, 3).im, 0, 'getitem_noraise')
    assert_equal(m1._getitem_noraise(3, 2).re, 0, 'getitem_noraise')
    assert_equal(m1._getitem_noraise(3, 2).im, 0, 'getitem_noraise')
    assert_equal(m1._getitem_noraise(0, 5).re, 0, 'getitem_noraise')
    assert_equal(m1._getitem_noraise(0, 5).im, 0, 'getitem_noraise')
    assert_equal(m1._getitem_linear_noraise(1, 1).re, 30, 'getitem_linear_noraise')
    assert_equal(m1._getitem_linear_noraise(1, 1).im, 5, 'getitem_linear_noraise')
    assert_equal(m1._getitem_linear_noraise(0, 0).re, 10, 'getitem_linear_noraise')
    assert_equal(m1._getitem_linear_noraise(0, 0).im, 0, 'getitem_linear_noraise')
    assert_equal(m1._getitem_linear_noraise(1, 3).re, 100, 'getitem_linear_noraise')
    assert_equal(m1._getitem_linear_noraise(1, 3).im, 0, 'getitem_linear_noraise')
    assert_equal(m1._getitem_linear_noraise(3, 2).re, 0, 'getitem_linear_noraise')
    assert_equal(m1._getitem_linear_noraise(3, 2).im, 0, 'getitem_linear_noraise')
    assert_equal(m1._getitem_linear_noraise(0, 5).re, 0, 'getitem_linear_noraise')
    assert_equal(m1._getitem_linear_noraise(0, 5).im, 0, 'getitem_linear_noraise')
    m1[0, 1] = (-90, 600)
    for r in range(4):
        for c in range(6):
            m1[r, c] = ComplexScalar[type](hash(r) % 10, hash(100 * c) % 10)
    for r in range(4):
        for c in range(6):
            assert_equal(m1[r, c].re, hash(r) % 10, 'setitem')
            assert_equal(m1[r, c].im, hash(100 * c) % 10, 'setitem')
    for r in range(4):
        for c in range(6):
            m1._setitem_noraise(r, c, ComplexScalar[type](hash(r + 5) % 10, hash(100 * c - 5) % 10))
    for r in range(4):
        for c in range(6):
            assert_equal(m1[r, c].re, hash(r + 5) % 10, 'setitem_noraise')
            assert_equal(m1[r, c].im, hash(100 * c - 5) % 10, 'setitem_noraise')
    for r in range(4):
        for c in range(6):
            m1._setitem_linear_noraise(r, c, ComplexScalar[type](hash(r - 5) % 50, hash(50 * c + 5) % 50))
    for r in range(4):
        for c in range(6):
            assert_equal(m1[r, c].re, hash(r - 5) % 50, 'setitem_linear_noraise')
            assert_equal(m1[r, c].im, hash(50 * c + 5) % 50, 'setitem_linear_noraise')
    m1 = CSRCMatrix[type](4, 6,
        (1, 1, ComplexScalar[type](30, 5)),
        (0, 0, ComplexScalar[type](10)),
        (0, 1, ComplexScalar[type](20)),
        (1, 3, ComplexScalar[type](40)),
        (3, 5, ComplexScalar[type](80)),
        (2, 2, ComplexScalar[type](50)),
        (1, 3, ComplexScalar[type](60)),
        (2, 4, ComplexScalar[type](70)),
    )
    assert_true((0, 0) in m1, 'contains')
    assert_true((1, 1) in m1, 'contains')
    assert_true((0, 1) in m1, 'contains')
    assert_true((1, 3) in m1, 'contains')
    assert_true((3, 5) in m1, 'contains')
    assert_true((2, 4) in m1, 'contains')
    assert_true((2, 2) in m1, 'contains')
    assert_false((1, 0) in m1, 'contains')
    assert_false((2, 1) in m1, 'contains')
    assert_false((3, 1) in m1, 'contains')
    assert_false((2, 3) in m1, 'contains')
    assert_false((1, 5) in m1, 'contains')
    assert_false((1, 4) in m1, 'contains')
    assert_false((0, 2) in m1, 'contains')
    assert_false((3, 0) in m1, 'contains')
    assert_false((-1, 0) in m1, 'contains')
    assert_false((4, 0) in m1, 'contains')
    assert_false((0, -1) in m1, 'contains')
    assert_false((0, 6) in m1, 'contains')
    m1[3, 5] = 0
    r0 = CSRCMatrix[type](1, m1.cols, (0, 0, ComplexScalar[type](10)), (0, 1, ComplexScalar[type](20)))
    r1 = CSRCMatrix[type](1, m1.cols, (0, 3, ComplexScalar[type](100)), (0, 1, ComplexScalar[type](30, 5)))
    r3 = CSRCMatrix[type](1, m1.cols)
    _assert_matrix_equal(m1.extract_row_as_sparse(0).to_dense(), r0.to_dense(), 'extract_row_as_sparse')
    _assert_matrix_equal(m1.extract_row_as_sparse(1).to_dense(), r1.to_dense(), 'extract_row_as_sparse')
    _assert_matrix_equal(m1.extract_row_as_sparse(3).to_dense(), r3.to_dense(), 'extract_row_as_sparse')
    _assert_matrix_equal(m1.extract_row_as_matrix(0), r0.to_dense(), 'extract_row_as_matrix')
    _assert_matrix_equal(m1.extract_row_as_matrix(1), r1.to_dense(), 'extract_row_as_matrix')
    _assert_matrix_equal(m1.extract_row_as_matrix(3), r3.to_dense(), 'extract_row_as_matrix')
    with assert_raises(contains='Invalid row'):
        _ = m1.extract_row_as_sparse(-1)
    with assert_raises(contains='Invalid row'):
        _ = m1.extract_row_as_sparse(5)
    with assert_raises(contains='Invalid row'):
        _ = m1.extract_row_as_matrix(-1)
    with assert_raises(contains='Invalid row'):
        _ = m1.extract_row_as_matrix(4)
    m1[2, 4] = 0
    m1[3, 3] = 70
    c0 = CSRCMatrix[type](m1.rows, 1, (0, 0, ComplexScalar[type](10)))
    c3 = CSRCMatrix[type](m1.rows, 1, (1, 0, ComplexScalar[type](100)), (3, 0, ComplexScalar[type](70)))
    c5 = CSRCMatrix[type](m1.rows, 1)
    _assert_matrix_equal(m1.extract_column_as_sparse(0).to_dense(), c0.to_dense(), 'extract_column_as_sparse')
    _assert_matrix_equal(m1.extract_column_as_sparse(3).to_dense(), c3.to_dense(), 'extract_column_as_sparse')
    _assert_matrix_equal(m1.extract_column_as_sparse(5).to_dense(), c5.to_dense(), 'extract_column_as_sparse')
    _assert_matrix_equal(m1.extract_column_as_matrix(0), c0.to_dense(), 'extract_column_as_matrix')
    _assert_matrix_equal(m1.extract_column_as_matrix(3), c3.to_dense(), 'extract_column_as_matrix')
    _assert_matrix_equal(m1.extract_column_as_matrix(5), c5.to_dense(), 'extract_column_as_matrix')
    with assert_raises(contains='Invalid column'):
        _ = m1.extract_column_as_sparse(-1)
    with assert_raises(contains='Invalid column'):
        _ = m1.extract_column_as_sparse(7)
    with assert_raises(contains='Invalid column'):
        _ = m1.extract_column_as_matrix(-1)
    with assert_raises(contains='Invalid column'):
        _ = m1.extract_column_as_matrix(7)

def test_other():
    v1 = List[ComplexScalar[type], True](10, 20, 30, 40, 50, 60, (70, 3), 80)
    c1 = List[Int, True](0, 1, 1, 3, 2, 3, 4, 5)   
    r1 = List[Int, True](0, 0, 1, 1, 2, 2, 2, 3)
    m1 = CSRCMatrix[type](4, 6, r1, c1, v1)
    m1d = CMatrix[type](4, 6,
        10, 20, 0, 0, 0, 0,
        0, 30, 0, 40, 0, 0,
        0, 0, 50, 60, (70, 3), 0,
        0, 0, 0, 0, 0, 80,
    )
    _assert_matrix_equal(m1.to_dense(), m1d, 'to_dense')

def test_arithmetic():
    v1 = List[ComplexScalar[type], True](10, 20, 30, 40, 50, 60, 70, 80)
    c1 = List[Int, True](0, 1, 1, 3, 2, 3, 4, 5)   
    r1 = List[Int, True](0, 0, 1, 1, 2, 2, 2, 3)
    sm1 = CSRCMatrix[type](4, 6, r1, c1, v1)
    v2 = List[ComplexScalar[type], True]((0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8))
    c2 = List[Int, True](0, 1, 4, 1, 2, 4, 4, 5)   
    r2 = List[Int, True](0, 1, 1, 3, 2, 1, 2, 3)
    sm2 = CSRCMatrix[type](4, 6, r2, c2, v2)
    for i in range(len(sm1.v)):
        assert_equal((-sm1).v[i].re, -(sm1.v[i].re), 'neg')
        assert_equal((-sm1).v[i].im, -(sm1.v[i].im), 'neg')
    for i in range(len(sm2.v)):
        assert_equal((-sm2).v[i].re, -(sm2.v[i].re), 'neg')
        assert_equal((-sm2).v[i].im, -(sm2.v[i].im), 'neg')
    cm = CMatrix[type](4, 6,
        1, 2, 3, 4, 5, 6,
        0, 0, 0, 1, 1, 1,
        0, -1, 0, -1, 0, -1,
        50, -45, 40, -35, 30, -25,
    )
    z = ComplexScalar[type]()
    s = ComplexScalar[type](14, -6)

    # sparse op sparse
    _assert_matrix_almost_equal((sm1 + sm1).to_dense(), 2 * sm1.to_dense(), 'add')
    _assert_matrix_almost_equal((sm1 + -sm1).to_dense(), sm1.to_dense().zeros_like(), 'add')
    _assert_matrix_almost_equal((sm1 + sm2).to_dense(), sm1.to_dense() + sm2.to_dense(), 'add')
    # sparse op matrix
    _assert_matrix_almost_equal(sm1 + cm, sm1.to_dense() + cm, 'add')
    _assert_matrix_almost_equal(sm2 + cm, sm2.to_dense() + cm, 'add')
    # sparse op scalar
    _assert_matrix_almost_equal((sm1 + z).to_dense(), sm1.to_dense(), 'add')
    _assert_matrix_almost_equal((sm2 + z).to_dense(), sm2.to_dense(), 'add')
    with assert_raises(contains='not supported'):
        _ = sm1 + s
    with assert_raises(contains='not supported'):
        _ = sm2 + s
    # matrix op sparse
    _assert_matrix_almost_equal(cm + sm1, cm + sm1.to_dense(), 'radd')
    _assert_matrix_almost_equal(cm + sm2, cm + sm2.to_dense(), 'radd')
    # scalar op sparse
    _assert_matrix_almost_equal((z + sm1).to_dense(), sm1.to_dense(), 'radd')
    _assert_matrix_almost_equal((z + sm2).to_dense(), sm2.to_dense(), 'radd')
    with assert_raises(contains='not supported'):
        _ = s + sm1
    with assert_raises(contains='not supported'):
        _ = s + sm2
    # sparse op= sparse
    si = sm1
    si += sm2
    _assert_matrix_almost_equal(si.to_dense(), sm1.to_dense() + sm2.to_dense(), 'iadd')
    # sparse op= matrix
    with assert_raises(contains='not supported'):
        sm1 += cm
    # sparse op= scalar
    si = sm1
    si += z
    _assert_matrix_almost_equal(si.to_dense(), sm1.to_dense(), 'iadd')
    with assert_raises(contains='not supported'):
        sm1 += s

    # sparse op sparse
    _assert_matrix_almost_equal((sm1 - sm1).to_dense(), sm1.to_dense().zeros_like(), 'sub')
    _assert_matrix_almost_equal((sm1 - sm2).to_dense(), sm1.to_dense() - sm2.to_dense(), 'sub')
    # sparse op matrix
    _assert_matrix_almost_equal(sm1 - cm, sm1.to_dense() - cm, 'sub')
    _assert_matrix_almost_equal(sm2 - cm, sm2.to_dense() - cm, 'sub')
    # sparse op scalar
    _assert_matrix_almost_equal((sm1 - z).to_dense(), sm1.to_dense(), 'sub')
    _assert_matrix_almost_equal((sm2 - z).to_dense(), sm2.to_dense(), 'sub')
    with assert_raises(contains='not supported'):
        _ = sm1 - s
    with assert_raises(contains='not supported'):
        _ = sm2 - s
    # matrix op sparse
    _assert_matrix_almost_equal(cm - sm1, cm - sm1.to_dense(), 'rsub')
    _assert_matrix_almost_equal(cm - sm2, cm - sm2.to_dense(), 'rsub')
    # scalar op sparse
    _assert_matrix_almost_equal((z - sm1).to_dense(), -sm1.to_dense(), 'rsub')
    _assert_matrix_almost_equal((z - sm2).to_dense(), -sm2.to_dense(), 'rsub')
    with assert_raises(contains='not supported'):
        _ = s - sm1
    with assert_raises(contains='not supported'):
        _ = s - sm2
    # sparse op= sparse
    si = sm1
    si -= sm2
    _assert_matrix_almost_equal(si.to_dense(), sm1.to_dense() - sm2.to_dense(), 'isub')
    # sparse op= matrix
    with assert_raises(contains='not supported'):
        sm1 -= cm
    # sparse op= scalar
    si = sm1
    si -= z
    _assert_matrix_almost_equal(si.to_dense(), sm1.to_dense(), 'isub')
    with assert_raises(contains='not supported'):
        sm1 -= s

    # sparse op sparse
    _assert_matrix_almost_equal((sm1 * sm1).to_dense(), sm1.to_dense() * sm1.to_dense(), 'mul')
    _assert_matrix_almost_equal((sm1 * sm2).to_dense(), sm1.to_dense() * sm2.to_dense(), 'mul')
    # sparse op matrix
    _assert_matrix_almost_equal(sm1 * cm, sm1.to_dense() * cm, 'mul')
    _assert_matrix_almost_equal(sm2 * cm, sm2.to_dense() * cm, 'mul')
    # sparse op scalar
    _assert_matrix_almost_equal((sm1 * z).to_dense(), sm1.to_dense().zeros_like(), 'mul')
    _assert_matrix_almost_equal((sm2 * z).to_dense(), sm2.to_dense().zeros_like(), 'mul')
    _assert_matrix_almost_equal((sm1 * s).to_dense(), sm1.to_dense() * s, 'mul')
    _assert_matrix_almost_equal((sm2 * s).to_dense(), sm2.to_dense() * s, 'mul')
    # matrix op sparse
    _assert_matrix_almost_equal(cm * sm1, cm * sm1.to_dense(), 'rmul')
    _assert_matrix_almost_equal(cm * sm2, cm * sm2.to_dense(), 'rmul')
    # scalar op sparse
    _assert_matrix_almost_equal((z * sm1).to_dense(), sm1.to_dense().zeros_like(), 'rmul')
    _assert_matrix_almost_equal((z * sm2).to_dense(), sm2.to_dense().zeros_like(), 'rmul')
    _assert_matrix_almost_equal((s * sm1).to_dense(), s * sm1.to_dense(), 'rmul')
    _assert_matrix_almost_equal((s * sm2).to_dense(), s * sm2.to_dense(), 'rmul')
    # sparse op= sparse
    si = sm1
    si *= sm2
    _assert_matrix_almost_equal(si.to_dense(), sm1.to_dense() * sm2.to_dense(), 'imul')
    # sparse op= matrix
    with assert_raises(contains='not supported'):
        sm1 *= cm
    # sparse op= scalar
    si = sm1
    si *= z
    _assert_matrix_almost_equal(si.to_dense(), sm1.to_dense().zeros_like(), 'imul')
    si = sm1
    si *= s
    _assert_matrix_almost_equal(si.to_dense(), sm1.to_dense() * s, 'imul')

    # sparse op sparse
    with assert_raises(contains='not supported'):
        _ = sm1 / sm1
    with assert_raises(contains='not supported'):
        _ = sm1 / sm2
    # sparse op matrix
    _assert_matrix_almost_equal(sm1 / cm, sm1.to_dense() / cm, 'truediv')
    _assert_matrix_almost_equal(sm2 / cm, sm2.to_dense() / cm, 'truediv')
    # sparse op scalar
    _assert_matrix_almost_equal((sm1 / s).to_dense(), sm1.to_dense() / s, 'truediv')
    _assert_matrix_almost_equal((sm2 / s).to_dense(), sm2.to_dense() / s, 'truediv')
    # matrix op sparse
    with assert_raises(contains='not supported'):
        _ = cm / sm1
    with assert_raises(contains='not supported'):
        _ = cm / sm2
    # scalar op sparse
    with assert_raises(contains='not supported'):
        _ = z / sm1
    with assert_raises(contains='not supported'):
        _ = z / sm2
    with assert_raises(contains='not supported'):
        _ = s / sm1
    with assert_raises(contains='not supported'):
        _ = s / sm2
    # sparse op= sparse
    with assert_raises(contains='not supported'):
        sm1 /= sm2
    # sparse op= matrix
    with assert_raises(contains='not supported'):
        sm1 /= cm
    # sparse op= scalar
    si = sm1
    si /= s
    _assert_matrix_almost_equal(si.to_dense(), sm1.to_dense() / s, 'itruediv')

    # sparse op sparse
    with assert_raises(contains='not supported'):
        _ = sm1 // sm1
    with assert_raises(contains='not supported'):
        _ = sm1 // sm2
    # sparse op matrix
    _assert_matrix_almost_equal(sm1 // cm, sm1.to_dense() // cm, 'floordiv')
    _assert_matrix_almost_equal(sm2 // cm, sm2.to_dense() // cm, 'floordiv')
    # sparse op scalar
    _assert_matrix_almost_equal((sm1 // s).to_dense(), sm1.to_dense() // s, 'floordiv')
    _assert_matrix_almost_equal((sm2 // s).to_dense(), sm2.to_dense() // s, 'floordiv')
    # matrix op sparse
    with assert_raises(contains='not supported'):
        _ = cm // sm1
    with assert_raises(contains='not supported'):
        _ = cm // sm2
    # scalar op sparse
    with assert_raises(contains='not supported'):
        _ = z // sm1
    with assert_raises(contains='not supported'):
        _ = z // sm2
    with assert_raises(contains='not supported'):
        _ = s // sm1
    with assert_raises(contains='not supported'):
        _ = s // sm2
    # sparse op= sparse
    with assert_raises(contains='not supported'):
        sm1 //= sm2
    # sparse op= matrix
    with assert_raises(contains='not supported'):
        sm1 //= cm
    # sparse op= scalar
    si = sm1
    si //= s
    _assert_matrix_almost_equal(si.to_dense(), sm1.to_dense() // s, 'ifloordiv')

    # sparse op sparse
    with assert_raises(contains='not supported'):
        _ = sm1 % sm1
    with assert_raises(contains='not supported'):
        _ = sm1 % sm2
    # sparse op matrix
    _assert_matrix_almost_equal(sm1 % cm, sm1.to_dense() % cm, 'mod')
    _assert_matrix_almost_equal(sm2 % cm, sm2.to_dense() % cm, 'mod')
    # sparse op scalar
    _assert_matrix_almost_equal((sm1 % s).to_dense(), sm1.to_dense() % s, 'mod')
    _assert_matrix_almost_equal((sm2 % s).to_dense(), sm2.to_dense() % s, 'mod')
    # matrix op sparse
    with assert_raises(contains='not supported'):
        _ = cm % sm1
    with assert_raises(contains='not supported'):
        _ = cm % sm2
    # scalar op sparse
    with assert_raises(contains='not supported'):
        _ = z % sm1
    with assert_raises(contains='not supported'):
        _ = z % sm2
    with assert_raises(contains='not supported'):
        _ = s % sm1
    with assert_raises(contains='not supported'):
        _ = s % sm2
    # sparse op= sparse
    with assert_raises(contains='not supported'):
        sm1 %= sm2
    # sparse op= matrix
    with assert_raises(contains='not supported'):
        sm1 %= cm
    # sparse op= scalar
    si = sm1
    si %= s
    _assert_matrix_almost_equal(si.to_dense(), sm1.to_dense() % s, 'imod')

def test_matmul():
    m1 = CSRCMatrix[type](4, 6, 
        List[Int, True](0, 0, 1, 1, 2, 2, 2, 3),
        List[Int, True](0, 1, 1, 3, 2, 3, 4, 5),
        List[ComplexScalar[type], True](10, 20, 30, 40, 50, 60, (70, 3), 80),
    )
    vec = CMatrix[type](6, 1)
    for i in range(6):
        vec[i] = (int(hash(i)) % 20 - 10, int(hash(i)) % 20 - 10)
    _assert_matrix_almost_equal(m1._dense_vec_matmul(vec), m1.to_dense() @ vec, 'sparse_vec_mul')
    _assert_matrix_almost_equal(m1._dense_vec_matmul(CMatrix[type](6, 1)), CMatrix[type](4, 1), 'sparse_vec_mul')
    cm = CMatrix[type](6, 8)
    for r in range(6):
        for c in range(8):
            cm[r, c] = (int(hash(r)) % 20 - 10, int(hash(c + r)) % 20 - 10)
    _assert_matrix_almost_equal(m1._dense_mat_matmul(cm), m1.to_dense() @ cm, 'sparse_dense_mul')
    _assert_matrix_almost_equal(m1._dense_mat_matmul(cm.zeros_like()), CMatrix[type](4, 8), 'sparse_dense_mul')
    m2 = CSRCMatrix[type](6, 8)
    for r in range(6):
        for c in range(8):
            m2[r, c] = (int(hash(r)) % 20 - 10, int(hash(c + r)) % 20 - 10)
    _assert_matrix_almost_equal(m1._sparse_matmul_gustavson(m2).to_dense(), m1.to_dense() @ m2.to_dense(), 'sparse_sparse_mul_gustavson')
    _assert_matrix_almost_equal(m1._sparse_matmul_gustavson(CSRCMatrix[type](6, 8)).to_dense(), CMatrix[type](4, 8), 'sparse_sparse_mul_gustavson')
    _assert_matrix_almost_equal(m1._sparse_matmul_hash(m2).to_dense(), m1.to_dense() @ m2.to_dense(), 'sparse_sparse_mul_hash')
    _assert_matrix_almost_equal(m1._sparse_matmul_hash(CSRCMatrix[type](6, 8)).to_dense(), CMatrix[type](4, 8), 'sparse_sparse_mul_hash')
    _assert_matrix_almost_equal(m1._sparse_matmul_tranposed(m2.transpose()).to_dense(), m1.to_dense() @ m2.to_dense(), 'sparse_matmul_tranposed')
    _assert_matrix_almost_equal(m1._sparse_matmul_tranposed(CSRCMatrix[type](6, 8).transpose()).to_dense(), CMatrix[type](4, 8), 'sparse_matmul_tranposed')
    m1 = CSRCMatrix[type](4, 6, 
        List[Int, True](0, 0, 0, 0, 2, 2, 2, 2),
        List[Int, True](0, 1, 1, 3, 2, 3, 4, 5),
        List[ComplexScalar[type], True](10, 20, 30, 40, 50, 60, (70, 3), 80),
    )
    _assert_matrix_almost_equal(m1._sparse_matmul_tranposed(m2.transpose()).to_dense(), m1.to_dense() @ m2.to_dense(), 'sparse_matmul_tranposed')
    _assert_matrix_almost_equal(m1._sparse_matmul_tranposed(CSRCMatrix[type](6, 8).transpose()).to_dense(), CMatrix[type](4, 8), 'sparse_matmul_tranposed')

def test_shape():
    m1 = CSRCMatrix[type](4, 6,
        (1, 1, ComplexScalar[type](30, 5)),
        (0, 0, ComplexScalar[type](10)),
        (0, 1, ComplexScalar[type](20)),
        (1, 3, ComplexScalar[type](40)),
        (3, 5, ComplexScalar[type](80)),
        (2, 2, ComplexScalar[type](50)),
        (1, 3, ComplexScalar[type](60)),
        (2, 4, ComplexScalar[type](70)),
    )
    _assert_matrix_equal(m1.transpose().to_dense(), m1.to_dense().transpose(), 'transpose')
    _assert_matrix_equal(CSRCMatrix[type](6, 6).transpose().to_dense(), CSRCMatrix[type](6, 6).to_dense(), 'transpose')
    _assert_matrix_equal(CSRCMatrix[type](0, 0).transpose().to_dense(), CSRCMatrix[type](0, 0).to_dense(), 'transpose')
    _assert_matrix_equal(m1.reshape(4, 6).to_dense(), m1.to_dense().reshape(4, 6), 'reshape')
    _assert_matrix_equal(m1.reshape(6, 4).to_dense(), m1.to_dense().reshape(6, 4), 'reshape')
    _assert_matrix_equal(m1.reshape(3, 8).to_dense(), m1.to_dense().reshape(3, 8), 'reshape')
    _assert_matrix_equal(m1.reshape(8, 3).to_dense(), m1.to_dense().reshape(8, 3), 'reshape')
    _assert_matrix_equal(m1.reshape(2, 12).to_dense(), m1.to_dense().reshape(2, 12), 'reshape')
    _assert_matrix_equal(m1.reshape(12, 2).to_dense(), m1.to_dense().reshape(12, 2), 'reshape')
    _assert_matrix_equal(m1.reshape(1, 24).to_dense(), m1.to_dense().reshape(1, 24), 'reshape')
    _assert_matrix_equal(m1.reshape(24, 1).to_dense(), m1.to_dense().reshape(24, 1), 'reshape')
    _assert_matrix_equal(CSRCMatrix[type](200, 54).reshape(20, 540).to_dense(), CMatrix[type](20, 540), 'reshape')
    with assert_raises(contains='reshape'):
        _ = m1.reshape(4, 7)

def test_compare():
    m1 = CSRCMatrix[type](4, 6,
        (1, 1, ComplexScalar[type](30, 5)),
        (0, 0, ComplexScalar[type](10)),
        (0, 1, ComplexScalar[type](20)),
        (1, 3, ComplexScalar[type](40)),
        (3, 5, ComplexScalar[type](80)),
        (2, 2, ComplexScalar[type](50)),
        (1, 3, ComplexScalar[type](60)),
        (2, 4, ComplexScalar[type](70)),
    )
    assert_true(m1.matrix_equals(m1), 'matrix_equals')
    assert_true(m1.is_close(m1), 'is_close')
    m2 = CSRCMatrix[type](5, 6,
        (1, 1, ComplexScalar[type](30, 5)),
        (0, 0, ComplexScalar[type](10)),
        (0, 1, ComplexScalar[type](20)),
        (1, 3, ComplexScalar[type](40)),
        (3, 5, ComplexScalar[type](80)),
        (2, 2, ComplexScalar[type](50)),
        (1, 3, ComplexScalar[type](60)),
        (2, 4, ComplexScalar[type](70)),
    )
    assert_false(m1.matrix_equals(m2), 'matrix_equals')
    assert_false(m1.is_close(m2), 'is_close')
    m2 = CSRCMatrix[type](4, 5,
        (1, 1, ComplexScalar[type](30, 5)),
        (0, 0, ComplexScalar[type](10)),
        (0, 1, ComplexScalar[type](20)),
        (1, 3, ComplexScalar[type](40)),
        (3, 3, ComplexScalar[type](80)),
        (2, 2, ComplexScalar[type](50)),
        (1, 3, ComplexScalar[type](60)),
        (2, 4, ComplexScalar[type](70)),
    )
    assert_false(m1.matrix_equals(m2), 'matrix_equals')
    assert_false(m1.is_close(m2), 'is_close')
    m2 = CSRCMatrix[type](4, 6,
        (1, 1, ComplexScalar[type](30, 5.00000001)),
        (0, 0, ComplexScalar[type](10)),
        (0, 1, ComplexScalar[type](20)),
        (1, 3, ComplexScalar[type](40)),
        (3, 5, ComplexScalar[type](80)),
        (2, 2, ComplexScalar[type](50)),
        (1, 3, ComplexScalar[type](60)),
        (2, 4, ComplexScalar[type](70)),
    )
    assert_false(m1.matrix_equals(m2), 'matrix_equals')
    assert_false(m1.is_close[1e-15](m2), 'is_close')
    assert_true(m1.is_close[1e-5](m2), 'is_close')
    m2 = CSRCMatrix[type](4, 6,
        (1, 1, ComplexScalar[type](30, 5)),
        (0, 0, ComplexScalar[type](10)),
        (0, 1, ComplexScalar[type](20)),
        (1, 3, ComplexScalar[type](40)),
        (3, 5, ComplexScalar[type](80)),
        (2, 2, ComplexScalar[type](50)),
        (1, 3, ComplexScalar[type](60)),
        (2, 4, ComplexScalar[type](70)),
        (3, 4, ComplexScalar[type](0.0000001))
    )
    assert_false(m1.matrix_equals(m2), 'matrix_equals')
    assert_false(m1.is_close[1e-15](m2), 'is_close')
    assert_true(m1.is_close[1e-5](m2), 'is_close')
    assert_false(m1.is_close(0))
    assert_true((CSRCMatrix[type](1000, 32456).is_close(0)), 'is_close')
    m2 = CSRCMatrix[type](12, 56)
    for r in range(12):
        for c in range(56):
            m2[r, c] = (6, 9)
    assert_true(m2.is_close((6, 9)))
    assert_true(m2.is_close[1.5]((7, 8)))
    assert_false(m2.is_close[1.4]((7, 8)))
    assert_equal((m1 == m1).sum().re, m1.rows * m1.cols, 'eq')
    assert_equal((m1 == m1).sum().im, 0, 'eq')
    m2 = CSRCMatrix[type](4, 6,
        (1, 1, ComplexScalar[type](30, 5)),
        (0, 0, ComplexScalar[type](10)),
        (0, 1, ComplexScalar[type](20)),
        (1, 3, ComplexScalar[type](40)),
        (3, 4, ComplexScalar[type](80)),
        (2, 2, ComplexScalar[type](50)),
        (1, 3, ComplexScalar[type](60)),
        (2, 4, ComplexScalar[type](70)),
    )
    for r in range(4):
        for c in range(6):
            if (r == 3 and c == 4) or (r == 3 and c == 5):
                assert_equal((m1 == m2)[r, c].re, 0, 'eq')
                assert_equal((m1 == m2)[r, c].im, 0, 'eq')
            else:
                assert_equal((m1 == m2)[r, c].re, 1, 'eq')
                assert_equal((m1 == m2)[r, c].im, 0, 'eq')
    with assert_raises():
        _ = m1 == CSRCMatrix[type](4, 5)
    cm1 = CMatrix[type](8, 7)
    cm1.fill_range()
    cm2 = CMatrix[type](8, 7)
    cm2.fill_range(0, -1)
    for r in range(8):
        for c in range(7):
            if r == 0 and c == 0:
                assert_equal((CSRCMatrix(cm1) == CSRCMatrix(cm2))[r, c].re, 1, 'eq')
                assert_equal((CSRCMatrix(cm1) == CSRCMatrix(cm2))[r, c].im, 0, 'eq')
            else:
                assert_equal((CSRCMatrix(cm1) == CSRCMatrix(cm2))[r, c].re, 0, 'eq')
                assert_equal((CSRCMatrix(cm1) == CSRCMatrix(cm2))[r, c].im, 0, 'eq')
    _assert_matrix_equal((m1 == m1).to_dense(), m1.to_dense().ones_like(), 'eq')
    _assert_matrix_equal((m1 == CSRCMatrix[type](m1.rows, m1.cols)).to_dense(), m1.to_dense() == m1.to_dense().zeros_like(), 'eq')
    _assert_matrix_equal((m1 == (30, 5)).to_dense(), m1.to_dense() == (30, 5), 'eq')
    _assert_matrix_equal((m1 == 0).to_dense(), m1.to_dense() == 0, 'eq')
    _assert_matrix_equal((CSRCMatrix[type](5, 6) == 0).to_dense(), CMatrix[type](5, 6).ones_like(), 'eq')
    _assert_matrix_equal((CSRCMatrix[type](5, 6) == 2).to_dense(), CMatrix[type](5, 6).zeros_like(), 'eq')
    m2 = CSRCMatrix[type](4, 6,
        (1, 1, ComplexScalar[type](35, 5)),
        (0, 0, ComplexScalar[type](10)),
        (0, 1, ComplexScalar[type](22)),
        (1, 3, ComplexScalar[type](40)),
        (1, 4, ComplexScalar[type](80)),
        (3, 2, ComplexScalar[type](60)),
        (1, 3, ComplexScalar[type](60)),
        (2, 2, ComplexScalar[type](70)),
    )
    _assert_matrix_equal((m1 != m1).to_dense(), m1.to_dense().zeros_like(), 'ne')
    _assert_matrix_equal((m1 != CSRCMatrix[type](m1.rows, m1.cols)).to_dense(), m1.to_dense() != m1.to_dense().zeros_like(), 'ne')
    _assert_matrix_equal((m1 != m2).to_dense(), m1.to_dense() != m2.to_dense(), 'ne')
    _assert_matrix_equal((m1 != (30, 5)).to_dense(), m1.to_dense() != (30, 5), 'ne')
    _assert_matrix_equal((m1 != 0).to_dense(), m1.to_dense() != 0, 'ne')
    _assert_matrix_equal((CSRCMatrix[type](5, 6) != 0).to_dense(), CMatrix[type](5, 6).zeros_like(), 'ne')
    _assert_matrix_equal((CSRCMatrix[type](5, 6) != 2).to_dense(), CMatrix[type](5, 6).ones_like(), 'ne')

    _assert_matrix_equal((m1 > m1).to_dense(), m1.to_dense().zeros_like(), 'gt')
    _assert_matrix_equal((m1 > CSRCMatrix[type](m1.rows, m1.cols)).to_dense(), m1.to_dense() > m1.to_dense().zeros_like(), 'gt')
    _assert_matrix_equal((m1 > m2).to_dense(), m1.to_dense() > m2.to_dense(), 'gt')
    _assert_matrix_equal((m1 > (30, 5)).to_dense(), m1.to_dense() > (30, 5), 'gt')
    _assert_matrix_equal((m1 > 0).to_dense(), m1.to_dense() > 0, 'gt')
    _assert_matrix_equal((CSRCMatrix[type](5, 6) > 0).to_dense(), CMatrix[type](5, 6).zeros_like(), 'gt')
    _assert_matrix_equal((CSRCMatrix[type](5, 6) > 2).to_dense(), CMatrix[type](5, 6).zeros_like(), 'gt')

    _assert_matrix_equal((m1 >= m1).to_dense(), m1.to_dense().ones_like(), 'ge')
    _assert_matrix_equal((m1 >= CSRCMatrix[type](m1.rows, m1.cols)).to_dense(), m1.to_dense() >= m1.to_dense().zeros_like(), 'ge')
    _assert_matrix_equal((m1 >= m2).to_dense(), m1.to_dense() >= m2.to_dense(), 'ge')
    _assert_matrix_equal((m1 >= (30, 5)).to_dense(), m1.to_dense() >= (30, 5), 'ge')
    _assert_matrix_equal((m1 >= 0).to_dense(), m1.to_dense() >= 0, 'ge')
    _assert_matrix_equal((CSRCMatrix[type](5, 6) >= 0).to_dense(), CMatrix[type](5, 6).ones_like(), 'ge')
    _assert_matrix_equal((CSRCMatrix[type](5, 6) >= 2).to_dense(), CMatrix[type](5, 6).zeros_like(), 'ge')

    _assert_matrix_equal((m1 < m1).to_dense(), m1.to_dense().zeros_like(), 'lt')
    _assert_matrix_equal((m1 < CSRCMatrix[type](m1.rows, m1.cols)).to_dense(), m1.to_dense() < m1.to_dense().zeros_like(), 'lt')
    _assert_matrix_equal((m1 < m2).to_dense(), m1.to_dense() < m2.to_dense(), 'lt')
    _assert_matrix_equal((m1 < (30, 5)).to_dense(), m1.to_dense() < (30, 5), 'lt')
    _assert_matrix_equal((m1 < 0).to_dense(), m1.to_dense() < 0, 'lt')
    _assert_matrix_equal((CSRCMatrix[type](5, 6) < 0).to_dense(), CMatrix[type](5, 6).zeros_like(), 'lt')
    _assert_matrix_equal((CSRCMatrix[type](5, 6) < 2).to_dense(), CMatrix[type](5, 6).ones_like(), 'lt')

    _assert_matrix_equal((m1 <= m1).to_dense(), m1.to_dense().ones_like(), 'le')
    _assert_matrix_equal((m1 <= CSRCMatrix[type](m1.rows, m1.cols)).to_dense(), m1.to_dense() <= m1.to_dense().zeros_like(), 'le')
    _assert_matrix_equal((m1 <= m2).to_dense(), m1.to_dense() <= m2.to_dense(), 'le')
    _assert_matrix_equal((m1 <= (30, 5)).to_dense(), m1.to_dense() <= (30, 5), 'le')
    _assert_matrix_equal((m1 <= 0).to_dense(), m1.to_dense() <= 0, 'le')
    _assert_matrix_equal((CSRCMatrix[type](5, 6) <= 0).to_dense(), CMatrix[type](5, 6).ones_like(), 'le')
    _assert_matrix_equal((CSRCMatrix[type](5, 6) <= 2).to_dense(), CMatrix[type](5, 6).ones_like(), 'le')

def test_fill():
    sm = CSRCMatrix[type](100, 100)
    cm = CMatrix[type](40, 30)
    cm.fill_range()
    cm[20, 20] = 0
    cm[10, 5] = 0
    cm[10, 6] = 0
    cm[11, 6] = 0
    i = sm.inset(cm, 0, 0)
    assert_equal(i.v.size, cm.rows * cm.cols - 5, 'inset_into_zero')
    for r in range(cm.rows):
        for c in range(cm.cols):
            assert_equal(i[r, c].re, cm[r, c].re, 'inset_into_zero')
            assert_equal(i[r, c].im, cm[r, c].im, 'inset_into_zero')
    i = sm.inset(cm, 8, 7)
    assert_equal(i.v.size, cm.rows * cm.cols - 5, 'inset_into_zero')
    for r in range(cm.rows):
        for c in range(cm.cols):
            assert_equal(i[r + 8, c + 7].re, cm[r, c].re, 'inset_into_zero')
            assert_equal(i[r + 8, c + 7].im, cm[r, c].im, 'inset_into_zero')
    i = sm.inset(cm, 60, 70)
    assert_equal(i.v.size, cm.rows * cm.cols - 5, 'inset_into_zero')
    for r in range(cm.rows):
        for c in range(cm.cols):
            assert_equal(i[r + 60, c + 70].re, cm[r, c].re, 'inset_into_zero')
            assert_equal(i[r + 60, c + 70].im, cm[r, c].im, 'inset_into_zero')
    with assert_raises():
        _ = sm.inset(cm, 61, 70)
    with assert_raises():
        _ = sm.inset(cm, 60, 71)
    with assert_raises():
        _ = sm.inset(cm, -1, 70)
    with assert_raises():
        _ = sm.inset(cm, 0, -1)
    i = sm.inset(cm, 0, 0, 2, 2)
    assert_equal(i.v.size, cm.rows * cm.cols - 5, 'inset_into_zero')
    for r in range(cm.rows):
        for c in range(cm.cols):
            assert_equal(i[2 * r, 2 * c].re, cm[r, c].re, 'inset_into_zero')
            assert_equal(i[2 * r, 2 * c].im, cm[r, c].im, 'inset_into_zero')
    i = sm.inset(cm, 8, 7, 2, 3)
    assert_equal(i.v.size, cm.rows * cm.cols - 5, 'inset_into_zero')
    for r in range(cm.rows):
        for c in range(cm.cols):
            assert_equal(i[2 * r + 8, 3 * c + 7].re, cm[r, c].re, 'inset_into_zero')
            assert_equal(i[2 * r + 8, 3 * c + 7].im, cm[r, c].im, 'inset_into_zero')
    i = sm.inset(cm, 21, 12, 2, 3)
    assert_equal(i.v.size, cm.rows * cm.cols - 5, 'inset_into_zero')
    for r in range(cm.rows):
        for c in range(cm.cols):
            assert_equal(i[2 * r + 21, 3 * c + 12].re, cm[r, c].re, 'inset_into_zero')
            assert_equal(i[2 * r + 21, 3 * c + 12].im, cm[r, c].im, 'inset_into_zero')
    with assert_raises():
        _ = sm.inset(cm, 22, 12, 2, 3)
    with assert_raises():
        _ = sm.inset(cm, 21, 13, 2, 3)
    sm = CSRCMatrix[type](100, 100)
    for r in range(100):
        for c in range(100):
            sm[hash(r * c + 5) % 100, hash(c + 7) % 100] = ComplexScalar[type](r * r + 1, c // r + c - 5)
    cm = CMatrix[type](20, 10)
    cm.fill_range()
    cm -= cm.ones_like() * 3 + cm.i_like() * 8
    cm[18, 8] = 0
    cm[10, 5] = 0
    cm[10, 6] = 0
    cm[11, 6] = 0
    i = sm.inset(cm, 0, 0)
    for r in range(100):
        for c in range(100):
            if r < 20 and c < 10:
                assert_equal(i[r, c].re, cm[r, c].re, 'inset')
                assert_equal(i[r, c].im, cm[r, c].im, 'inset')
            else:
                assert_equal(i[r, c].re, sm[r, c].re, 'inset')
                assert_equal(i[r, c].im, sm[r, c].im, 'inset')
    i = sm.inset(cm, 17, 51)
    for r in range(100):
        for c in range(100):
            if r >= 17 and r < 37 and c >= 51 and c < 61:
                assert_equal(i[r, c].re, cm[r - 17, c - 51].re, 'inset')
                assert_equal(i[r, c].im, cm[r - 17, c - 51].im, 'inset')
            else:
                assert_equal(i[r, c].re, sm[r, c].re, 'inset')
                assert_equal(i[r, c].im, sm[r, c].im, 'inset')
    i = sm.inset(cm, 80, 90)
    for r in range(100):
        for c in range(100):
            if r >= 80 and r < 100 and c >= 90 and c < 100:
                assert_equal(i[r, c].re, cm[r - 80, c - 90].re, 'inset')
                assert_equal(i[r, c].im, cm[r - 80, c - 90].im, 'inset')
            else:
                assert_equal(i[r, c].re, sm[r, c].re, 'inset')
                assert_equal(i[r, c].im, sm[r, c].im, 'inset')
    with assert_raises():
        _ = sm.inset(cm, 90, 90)
    with assert_raises():
        _ = sm.inset(cm, 80, 95)
    i = sm.inset(cm, 2, 3, 4, 5)
    for r in range(100):
        for c in range(100):
            if r >= 2 and c >= 3 and (r - 2) % 4 == 0 and (c - 3) % 5 == 0 and r <= (2 + (20 - 1) * 4) and c <= (3 + (10 - 1) * 5):
                assert_equal(i[r, c].re, cm[(r - 2) // 4, (c - 3) // 5].re, 'inset')
                assert_equal(i[r, c].im, cm[(r - 2) // 4, (c - 3) // 5].im, 'inset')
            else:
                assert_equal(i[r, c].re, sm[r, c].re, 'inset')
                assert_equal(i[r, c].im, sm[r, c].im, 'inset')
    i = sm.inset(cm, 2, 3, 2, 3).inset(cm, 55, 71, 1, 2)
    for r in range(100):
        for c in range(100):
            if r >= 2 and c >= 3 and (r - 2) % 2 == 0 and (c - 3) % 3 == 0 and r <= (2 + (20 - 1) * 2) and c <= (3 + (10 - 1) * 3):
                assert_equal(i[r, c].re, cm[(r - 2) // 2, (c - 3) // 3].re, 'inset')
                assert_equal(i[r, c].im, cm[(r - 2) // 2, (c - 3) // 3].im, 'inset')
            elif r >= 55 and c >= 71 and (r - 55) % 1 == 0 and (c - 71) % 2 == 0 and r <= (55 + (20 - 1) * 1) and c <= (71 + (10 - 1) * 2):
                assert_equal(i[r, c].re, cm[(r - 55) // 1, (c - 71) // 2].re, 'inset')
                assert_equal(i[r, c].im, cm[(r - 55) // 1, (c - 71) // 2].im, 'inset')
            else:
                assert_equal(i[r, c].re, sm[r, c].re, 'inset')
                assert_equal(i[r, c].im, sm[r, c].im, 'inset')
    cm = CMatrix[type](22, 1).range_like(-13, 2)
    i = sm.inset(cm, 0, 0)
    for r in range(100):
        for c in range(100):
            if r < 22 and c == 0:
                assert_equal(i[r, c].re, cm[r, c].re, 'inset_column')
                assert_equal(i[r, c].im, cm[r, c].im, 'inset_column')
            else:
                assert_equal(i[r, c].re, sm[r, c].re, 'inset_column')
                assert_equal(i[r, c].im, sm[r, c].im, 'inset_column')
    i = sm.inset(cm, 78, 99)
    for r in range(100):
        for c in range(100):
            if r >= 78 and c == 99:
                assert_equal(i[r, c].re, cm[r - 78, c - 99].re, 'inset_column')
                assert_equal(i[r, c].im, cm[r - 78, c - 99].im, 'inset_column')
            else:
                assert_equal(i[r, c].re, sm[r, c].re, 'inset_column')
                assert_equal(i[r, c].im, sm[r, c].im, 'inset_column')
    i = sm.inset(cm, 4, 7, 3, 895)
    for r in range(100):
        for c in range(100):
            if r >= 4 and r <= (4 + (22 - 1) * 3) and c == 7 and (r - 4) % 3 == 0:
                assert_equal(i[r, c].re, cm[(r - 4) // 3, c - 7].re, 'inset_column')
                assert_equal(i[r, c].im, cm[(r - 4) // 3, c - 7].im, 'inset_column')
            else:
                assert_equal(i[r, c].re, sm[r, c].re, 'inset_column')
                assert_equal(i[r, c].im, sm[r, c].im, 'inset_column')
    with assert_raises():
        _ = sm.inset(cm, 80, 0)
    with assert_raises():
        _ = sm.inset(cm, 40, 0, 6)
    with assert_raises():
        _ = sm.inset(cm, 50, 100)
    c = ComplexScalar[type]((-12.5, 6.6))
    smc = sm
    smc.fill_diag(c)
    cmc = sm.to_dense()
    cmc.fill_diag(c)
    _assert_matrix_equal(smc.to_dense(), cmc, 'fill_diag')
    smc = sm
    smc.fill_diag(c, 14)
    cmc = sm.to_dense()
    cmc.fill_diag(c, 14)
    _assert_matrix_equal(smc.to_dense(), cmc, 'fill_diag')
    smc = sm
    smc.fill_diag(c, -14)
    cmc = sm.to_dense()
    cmc.fill_diag(c, -14)
    _assert_matrix_equal(smc.to_dense(), cmc, 'fill_diag')

def test_static_constructors():
    for r in range(1, 50):
        for c in range(1, 50):
            _assert_matrix_equal(CSRCMatrix[type].zeros(r, c).to_dense(), CMatrix[type](r, c), 'zeros')
            _assert_matrix_equal(CSRCMatrix[type].eye(r, c).to_dense(), CMatrix[type].eye(r, c), 'zeros')
