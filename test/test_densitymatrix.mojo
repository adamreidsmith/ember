from collections import Dict
from math import sqrt

from testing import assert_equal, assert_raises, assert_almost_equal, assert_true, assert_false
from ._testing import _assert_matrix_equal, _assert_matrix_almost_equal

from ember import DensityMatrix, CMatrix, is_positive_semidefinite
from ember import ComplexScalar, Statevector

alias type = DType.float64

def run_densitymatrix_tests():
    print('Running densitymatrix tests')
    test_init()
    test_purity()
    test_partial_trace()
    test_to_sv_data()
    print('All tests passed')

def test_init():
    dm1 = DensityMatrix[type]()
    assert_equal(dm1.n_qubits, 0, 'init')
    assert_equal(dm1.matrix.rows, 0, 'init')
    assert_equal(dm1.matrix.cols, 0, 'init')

    dm_mat = CMatrix[type](4, 4,
        (0.20925406, 1.47808531e-18),  (0.02749093, -1.56398637e-02), (-0.11303181, -8.06957426e-02), (-0.01738071, -1.27452033e-01),
        (0.02749093, 1.56398637e-02),  (0.25263993, -6.81872229e-19), (-0.03743478, -1.41255212e-01), (-0.09342171, 7.14286538e-02),
        (-0.11303181, 8.06957426e-02), (-0.03743478, 1.41255212e-01), (0.27195739, 1.49465348e-18), (0.02339803, -5.69018290e-02),
        (-0.01738071, 1.27452033e-01), (-0.09342171, -7.14286538e-02), (0.02339803, 5.69018290e-02),  (0.26614863, 6.78572149e-18),
    )
    dm2 = DensityMatrix[type, 1e-6](dm_mat)
    assert_equal(dm2.n_qubits, 2, 'init')
    assert_equal(dm2.matrix.rows, 4, 'init')
    assert_equal(dm2.matrix.cols, 4, 'init')

    dm_mat = CMatrix[type](4, 4,
        (0.40925406, 1.47808531e-18),  (0.02749093, -1.56398637e-02), (-0.11303181, -8.06957426e-02), (-0.01738071, -1.27452033e-01),
        (0.02749093, 1.56398637e-02),  (0.25263993, -6.81872229e-19), (-0.03743478, -1.41255212e-01), (-0.09342171, 7.14286538e-02),
        (-0.11303181, 8.06957426e-02), (-0.03743478, 1.41255212e-01), (0.27195739, 1.49465348e-18), (0.02339803, -5.69018290e-02),
        (-0.01738071, 1.27452033e-01), (-0.09342171, -7.14286538e-02), (0.02339803, 5.69018290e-02),  (0.26614863, 6.78572149e-18),
    )
    with assert_raises(contains='trace 1'):
        dm2 = DensityMatrix[type, 1e-6](dm_mat)
    
    dm_mat = CMatrix[type](5, 4,
        (0.20925406, 1.47808531e-18),  (0.02749093, -1.56398637e-02), (-0.11303181, -8.06957426e-02), (-0.01738071, -1.27452033e-01),
        (0.02749093, 1.56398637e-02),  (0.25263993, -6.81872229e-19), (-0.03743478, -1.41255212e-01), (-0.09342171, 7.14286538e-02),
        (-0.11303181, 8.06957426e-02), (-0.03743478, 1.41255212e-01), (0.27195739, 1.49465348e-18), (0.02339803, -5.69018290e-02),
        (-0.01738071, 1.27452033e-01), (-0.09342171, -7.14286538e-02), (0.02339803, 5.69018290e-02),  (0.26614863, 6.78572149e-18),
        (-0.01738071, 1.27452033e-01), (-0.09342171, -7.14286538e-02), (0.02339803, 5.69018290e-02),  (0.26614863, 6.78572149e-18),
    )
    with assert_raises(contains='power of 2'):
        dm2 = DensityMatrix[type, 1e-6](dm_mat)
    
    dm_mat = CMatrix[type](4, 5,
        (0.20925406, 1.47808531e-18),  (0.02749093, -1.56398637e-02), (-0.11303181, -8.06957426e-02), (-0.01738071, -1.27452033e-01), (-0.01738071, -1.27452033e-01),
        (0.02749093, 1.56398637e-02),  (0.25263993, -6.81872229e-19), (-0.03743478, -1.41255212e-01), (-0.09342171, 7.14286538e-02), (-0.01738071, -1.27452033e-01),
        (-0.11303181, 8.06957426e-02), (-0.03743478, 1.41255212e-01), (0.27195739, 1.49465348e-18), (0.02339803, -5.69018290e-02), (-0.01738071, -1.27452033e-01),
        (-0.01738071, 1.27452033e-01), (-0.09342171, -7.14286538e-02), (0.02339803, 5.69018290e-02),  (0.26614863, 6.78572149e-18), (-0.01738071, -1.27452033e-01),
    )
    with assert_raises(contains='not Hermitian'):
        dm2 = DensityMatrix[type, 1e-6](dm_mat)
    
    dm_mat = CMatrix[type](4, 4,
        (0.20925406, 1.47808531e-18),  (0.02749093, -1.56398637e-02), (-0.11303181, -8.06957426e-02), (-0.01738071, 1.27452033e-01),
        (0.02749093, 1.56398637e-02),  (0.25263993, -6.81872229e-19), (-0.03743478, -1.41255212e-01), (-0.09342171, 7.14286538e-02),
        (-0.11303181, 8.06957426e-02), (-0.03743478, 1.41255212e-01), (0.27195739, 1.49465348e-18), (0.02339803, -5.69018290e-02),
        (-0.01738071, 1.27452033e-01), (-0.09342171, -7.14286538e-02), (0.02339803, 5.69018290e-02),  (0.26614863, 6.78572149e-18),
    )
    with assert_raises(contains='not Hermitian'):
        dm2 = DensityMatrix[type, 1e-6](dm_mat)
    
    dm_mat = CMatrix[type](2, 2,
        1,  (0.0, 1.0),
        (0.0, -1.0),  0,
    )
    with assert_raises(contains='not positive semi-definite'):
        dm2 = DensityMatrix[type, 1e-6](dm_mat)
    
    sv = Dict[Int, ComplexScalar[type]]()
    sv[0] = (0.5, 0.5)
    sv[2] = (0.5, 0.5)

    with assert_raises(contains='Invalid statevector size'):
        dm = DensityMatrix[type](sv, -1)
    with assert_raises(contains='power of 2'):
        dm = DensityMatrix[type](sv, 5)
    with assert_raises(contains='power of 2'):
        dm = DensityMatrix[type](sv, 3)
    dm = DensityMatrix[type](sv, 0)
    assert_equal(dm.n_qubits, 0, 'init')
    assert_equal(dm.matrix.rows, 0, 'init')
    assert_equal(dm.matrix.cols, 0, 'init')

    sv[1] = (0.5, 0.5)
    with assert_raises(contains='not normalized'):
        dm = DensityMatrix[type](sv, 4)
    _ = sv.pop(1)

    dm = DensityMatrix(sv, 4)
    assert_equal(dm.n_qubits, 2, 'init')
    assert_equal(dm.matrix.rows, 4, 'init')
    assert_equal(dm.matrix.cols, 4, 'init')

def test_purity():
    svl = List[ComplexScalar[DType.float64], True](
        (1.0, 3.1235),
        (1.21352, -5.123),
        (0.4532, 1.12),
        (-4.123, 0.0),
        (3.1234, 1.34),
        (1.5234, -1.2341),
        (-0.3452, 5.97),
        (2.51234, 2.123),
        (-1.12354, -3.245),
        (2.25, -0.00123),
        (2.235144445, 2.2),
        (6.45, 4.1234),
        (-4.2134, -2.1),
        (1.354, 5.4),
        (0.000123, -0.234),
        (0.1, 0.1)
    )
    norm = 0.0
    for item in svl:
        norm += item[].squared_norm()
    norm = sqrt(norm)
    for i in range(len(svl)):
        svl[i] /= norm
    sv = Statevector(svl)
    dm = DensityMatrix(sv._data, sv.size)
    assert_almost_equal(dm.purity(), 1, 'purity')
    assert_true(dm.is_pure(), 'is_pure')

    dm_mat = CMatrix[type].eye(8, 8) / 8
    dm = DensityMatrix(dm_mat)
    assert_false(dm.is_pure(), 'is_pure')
    assert_almost_equal(dm.purity(), 1/8, 'purity')

def test_partial_trace():
    sv = Statevector[type].zero(4)
    dm = DensityMatrix[type](sv._data, sv.size)
    with assert_raises(contains='Invalid qubit'):
        _ = dm.partial_trace(-1)
    with assert_raises(contains='Invalid qubit'):
        _ = dm.partial_trace(4)
    with assert_raises(contains='repeated qubit'):
        _ = dm.partial_trace(0, 0)
    _assert_matrix_almost_equal(dm.partial_trace(0).matrix, DensityMatrix[type](Statevector[type].zero(3)._data, 8).matrix, 'partial_trace')
    _assert_matrix_almost_equal(dm.partial_trace(1).matrix, DensityMatrix[type](Statevector[type].zero(3)._data, 8).matrix, 'partial_trace')
    _assert_matrix_almost_equal(dm.partial_trace(2).matrix, DensityMatrix[type](Statevector[type].zero(3)._data, 8).matrix, 'partial_trace')
    _assert_matrix_almost_equal(dm.partial_trace(3).matrix, DensityMatrix[type](Statevector[type].zero(3)._data, 8).matrix, 'partial_trace')
    _assert_matrix_almost_equal(dm.partial_trace(0, 1).matrix, DensityMatrix[type](Statevector[type].zero(2)._data, 4).matrix, 'partial_trace')
    _assert_matrix_almost_equal(dm.partial_trace(0, 2).matrix, DensityMatrix[type](Statevector[type].zero(2)._data, 4).matrix, 'partial_trace')
    _assert_matrix_almost_equal(dm.partial_trace(0, 3).matrix, DensityMatrix[type](Statevector[type].zero(2)._data, 4).matrix, 'partial_trace')
    _assert_matrix_almost_equal(dm.partial_trace(2, 0).matrix, DensityMatrix[type](Statevector[type].zero(2)._data, 4).matrix, 'partial_trace')
    _assert_matrix_almost_equal(dm.partial_trace(2, 1).matrix, DensityMatrix[type](Statevector[type].zero(2)._data, 4).matrix, 'partial_trace')
    _assert_matrix_almost_equal(dm.partial_trace(2, 3).matrix, DensityMatrix[type](Statevector[type].zero(2)._data, 4).matrix, 'partial_trace')
    _assert_matrix_almost_equal(dm.partial_trace(0, 1, 2).matrix, DensityMatrix[type](Statevector[type].zero(1)._data, 2).matrix, 'partial_trace')
    _assert_matrix_almost_equal(dm.partial_trace(2, 1, 3).matrix, DensityMatrix[type](Statevector[type].zero(1)._data, 2).matrix, 'partial_trace')
    _assert_matrix_almost_equal(dm.partial_trace(0, 3, 1).matrix, DensityMatrix[type](Statevector[type].zero(1)._data, 2).matrix, 'partial_trace')

    svl = List[ComplexScalar[DType.float64], True](
        (1.0, 3.1235),
        (1.21352, -5.123),
        (0.4532, 1.12),
        (-4.123, 0.0),
        (3.1234, 1.34),
        (1.5234, -1.2341),
        (-0.3452, 5.97),
        (2.51234, 2.123),
        (-1.12354, -3.245),
        (2.25, -0.00123),
        (2.235144445, 2.2),
        (6.45, 4.1234),
        (-4.2134, -2.1),
        (1.354, 5.4),
        (0.000123, -0.234),
        (0.1, 0.1)
    )
    norm = 0.0
    for item in svl:
        norm += item[].squared_norm()
    norm = sqrt(norm)
    for i in range(len(svl)):
        svl[i] /= norm
    sv = Statevector(svl)
    dm = DensityMatrix(sv._data, sv.size)

    pdm = dm.partial_trace(0)
    tdm = CMatrix[type](8, 8,
        ( 0.14945227, +1.50660358e-18), (-0.00408581, +8.31970387e-02),
        ( 0.06013151, +8.19305681e-03), ( 0.04068912, -8.73827956e-02),
        (-0.03310594, -4.57968278e-02), (-0.01627656, -1.29220608e-01),
        (-0.14292566, -9.53648666e-02), (-0.00435734, -1.55095555e-03),
        (-0.00408581, -8.31970387e-02), ( 0.07170366, +1.93270356e-19),
        (-0.01306999, -8.53531060e-03), (-0.0148715 , +2.19897663e-02),
        (-0.0521313 , +8.04860147e-04), (-0.08979557, +7.18908179e-02),
        (-0.03823922, +7.18512427e-02), (-0.00261941, +2.01406071e-03),
        ( 0.06013151, -8.19305681e-03), (-0.01306999, +8.53531060e-03),
        ( 0.05980178, -2.85182619e-20), ( 0.03157673, -9.88369024e-02),
        (-0.0172021 , +2.27438920e-02), ( 0.05697198, -7.03790446e-02),
        (-0.07992587, -3.48989148e-02), (-0.00110415, +1.76857534e-03),
        ( 0.04068912, +8.73827956e-02), (-0.0148715 , -2.19897663e-02),
        ( 0.03157673, +9.88369024e-02), ( 0.18093627, -8.75142585e-19),
        (-0.05179844, -1.18394398e-02), ( 0.1449734 , +6.77348682e-02),
        ( 0.01469651, -1.42060072e-01), (-0.00362614, -4.62164098e-04),
        (-0.03310594, +4.57968278e-02), (-0.0521313 , -8.04860147e-04),
        (-0.0172021 , -2.27438920e-02), (-0.05179844, +1.18394398e-02),
        ( 0.06547273, -6.13806303e-19), ( 0.01886761, -5.46426279e-02),
        ( 0.05666817, -3.25760494e-03), ( 0.00382261, -1.89730794e-03),
        (-0.01627656, +1.29220608e-01), (-0.08979557, -7.18908179e-02),
        ( 0.05697198, +7.03790446e-02), ( 0.1449734 , -6.77348682e-02),
        ( 0.01886761, +5.46426279e-02), ( 0.26585824, +6.78533057e-18),
        ( 0.06588928, -1.31383841e-01), ( 0.00210856, +1.12896914e-03),
        (-0.14292566, +9.53648666e-02), (-0.03823922, -7.18512427e-02),
        (-0.07992587, +3.48989148e-02), ( 0.01469651, +1.42060072e-01),
        ( 0.05666817, +3.25760494e-03), ( 0.06588928, +1.31383841e-01),
        ( 0.20648466, +2.10845978e-18), ( 0.00453042, -2.25920113e-03),
        (-0.00435734, +1.55095555e-03), (-0.00261941, -2.01406071e-03),
        (-0.00110415, -1.76857534e-03), (-0.00362614, +4.62164098e-04),
        ( 0.00382261, +1.89730794e-03), ( 0.00210856, -1.12896914e-03),
        ( 0.00453042, +2.25920113e-03), ( 0.00029039, +3.90919901e-22)
    )
    _assert_matrix_almost_equal(pdm.matrix, tdm, 'partial_trace', atol=1e-7)

    pdm = dm.partial_trace(3, 1)
    tdm = CMatrix[type](4, 4,
        ( 0.13146808, +7.53717237e-19), ( 0.01673238, +7.63942409e-03),
        ( 0.09661816, +6.66582237e-02), (-0.06765134, +3.69570828e-02),
        ( 0.01673238, -7.63942409e-03), ( 0.42101881, +7.11768096e-18),
        (-0.04697487, +5.13802910e-02), ( 0.00741858, -3.86040364e-02),
        ( 0.09661816, -6.66582237e-02), (-0.04697487, -5.13802910e-02),
        ( 0.27008452, +2.02385889e-18), (-0.00837804, +1.61257303e-01),
        (-0.06765134, -3.69570828e-02), ( 0.00741858, +3.86040364e-02),
        (-0.00837804, -1.61257303e-01), ( 0.17742858, -8.18669037e-19),
    )
    _assert_matrix_almost_equal(pdm.matrix, tdm, 'partial_trace', atol=1e-7)

    pdm = dm.partial_trace(2, 3, 0)
    tdm = CMatrix[type](2, 2,
        (0.48121145, +2.97273879e-18), (0.05088896, -7.25416927e-02),
        (0.05088896, +7.25416927e-02), (0.51878855, +6.10384926e-18),
    )
    _assert_matrix_almost_equal(pdm.matrix, tdm, 'partial_trace', atol=1e-7)

    _assert_matrix_almost_equal(dm.partial_trace(0).partial_trace(1).matrix, dm.partial_trace(0, 2).matrix, 'partial_trace')
    _assert_matrix_almost_equal(dm.partial_trace(0).partial_trace(0).partial_trace(1).matrix, dm.partial_trace(0, 1, 3).matrix, 'partial_trace')

def test_to_sv_data():
    dm_mat = CMatrix[type].eye(8, 8) / 8
    dm = DensityMatrix(dm_mat)
    with assert_raises(contains='not pure'):
        _ = dm._to_sv_data()

    sv = Statevector[type].zero(5)
    dm = DensityMatrix(sv._data, sv.size)
    sv_data = dm._to_sv_data()
    assert_equal(sv_data[1], sv.size, '_to_sv_data')
    assert_equal(len(sv_data[0]), len(sv._data), '_to_sv_data')
    for idx in sv._data:
        assert_true(idx[] in sv_data[0])
        assert_almost_equal(sv_data[0][idx[]].re, sv._data[idx[]].re, '_to_sv_data')
        assert_almost_equal(sv_data[0][idx[]].im, sv._data[idx[]].im, '_to_sv_data')

    svl = List[ComplexScalar[DType.float64], True](
        (1.0, 3.1235),
        (1.21352, -5.123),
        (0.4532, 1.12),
        (-4.123, 0.0),
        (3.1234, 1.34),
        (1.5234, -1.2341),
        (-0.3452, 5.97),
        (2.51234, 2.123),
        (-1.12354, -3.245),
        (2.25, -0.00123),
        (2.235144445, 2.2),
        (6.45, 4.1234),
        (-4.2134, -2.1),
        (1.354, 5.4),
        (0.000123, -0.234),
        (0.1, 0.1)
    )
    norm = 0.0
    for item in svl:
        norm += item[].squared_norm()
    norm = sqrt(norm)
    for i in range(len(svl)):
        svl[i] /= norm
    sv = Statevector(svl)
    dm = DensityMatrix(sv._data, sv.size)
    sv_data = dm._to_sv_data()
    # Statevectors are only defined up to a phase, so we must adjust the phase of the new
    # statevector so it matcehs the old before we can compare them
    sv_matrix = CMatrix[type](16, 1)
    for i in range(16):
        sv_matrix[i] = sv_data[0].get(i, 0)
    sv_matrix *= ComplexScalar[type](0.062325751924004975, 0.19467448613462954) / sv_matrix[0]
    assert_equal(sv_matrix.size, sv.size, '_to_sv_data')
    for idx in range(16):
        assert_almost_equal(sv_matrix[idx].re, sv._data[idx].re, '_to_sv_data')
        assert_almost_equal(sv_matrix[idx].im, sv._data[idx].im, '_to_sv_data')
