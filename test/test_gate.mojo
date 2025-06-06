from math import sqrt, pi, sin, cos

from testing import assert_equal, assert_raises
from ._testing import _assert_matrix_equal

from ember import Gate, Measure, X, Y, Z, H, S, T, SX, I, RX, RY, RZ, CX, CCX, ECR, SWAP, U, PHASE, R, RXX, RYY, RZZ, RZX, XXMinusYY, XXPlusYY
from ember import CMatrix, ComplexScalar


alias type = DType.float64

def run_gate_tests():
    print('Running gate tests')
    test_init()
    test_measure()
    test_control()
    test_other()
    test_unparameterized_single_qubit_gates()
    test_unparameterized_multi_qubit_gates()
    test_parameterized_single_qubit_gates()
    test_parameterized_multi_qubit_gates()
    print('All tests passed')

def test_init():
    with assert_raises(contains='Duplicate'):
        gate = Gate[type]('test_gate', CMatrix[type].eye(4, 4), List[Int, True](0, 0))
    with assert_raises(contains='Invalid gate matrix'):
        gate = Gate[type]('test_gate', CMatrix[type].eye(4, 5), List[Int, True](0, 1))
    with assert_raises(contains='Invalid gate matrix'):
        gate = Gate[type]('test_gate', CMatrix[type].eye(5, 5), List[Int, True](0, 1))
    with assert_raises(contains='not unitary'):
        gate = Gate[type]('test_gate', CMatrix[type].arange(8, 8), List[Int, True](0, 1, 2))
    with assert_raises(contains='number of qubits'):
        gate = Gate[type]('test_gate', CMatrix[type].eye(4, 4), List[Int, True](0, 1, 2))
    with assert_raises(contains='number of qubits'):
        gate = Gate[type]('test_gate', CMatrix[type].eye(4, 4), List[Int, True](0))
    gate = Gate[type]('test_gate', CMatrix[type].eye(8, 8), List[Int, True](0, 1, 2))
    assert_equal(len(gate.controls), 0, 'init')
    assert_equal(len(gate._measure_targs), 0, 'init')
    assert_equal(len(gate.classical_controls), 0, 'init')
    assert_equal(len(gate.qubits), 3, 'init')
    assert_equal(len(gate.params), 0, 'init')
    gate = Gate[type]('test_gate', CMatrix[type].eye(8, 8), List[Int, True](0, 1, 2), List[Scalar[type], True](0.2, 0.876, 3))
    assert_equal(len(gate.params), 3, 'init')
    assert_equal(gate._is_measure, False, 'init')

def test_measure():
    gate = Gate[type]._measure(List[Int, True](1, 2, 3), List[Int, True](9, 5, 1))
    assert_equal(len(gate.qubits), 3, '_measure')
    assert_equal(len(gate.controls), 0, '_measure')
    assert_equal(len(gate._measure_targs), 3, '_measure')
    assert_equal(len(gate.params), 0, '_measure')
    assert_equal(len(gate.classical_controls), 0, '_measure')
    assert_equal(gate.name, 'Measure', '_measure')
    assert_equal(gate._is_measure, True, '_measure')
    assert_equal(gate.matrix.size, 0, '_measure')
    # Measure gate
    gate = Measure[type](0, 0)
    assert_equal(gate._is_measure, True, 'Measure')
    with assert_raises(contains='no qubit specifiers'):
        gate = Measure[type](List[Int, True](), List[Int, True]())
    with assert_raises(contains='qubits must be unique'):
        gate = Measure[type](List[Int, True](0, 0), List[Int, True](0, 1))
    with assert_raises(contains='classical bits must be unique'):
        gate = Measure[type](List[Int, True](0, 1), List[Int, True](1, 1))
    with assert_raises(contains='same length'):
        gate = Measure[type](List[Int, True](0, 1, 2), List[Int, True](0, 1))
    gate = Measure[type](List[Int, True](0, 1, 2), List[Int, True](0, 3, 1))
    assert_equal(gate._is_measure, True, 'Measure')

def test_control():
    gate = Gate[type]('test_gate', CMatrix[type].eye(8, 8), List[Int, True](0, 1, 2))
    with assert_raises(contains='on which it is applied'):
        gate.control(0, 3, 4)
    gate.control(3, 4, 5, 5)
    assert_equal(gate.n_qubits, 6, 'control')
    assert_equal(len(gate.controls), 3, 'control')
    gate = Gate[type]('test_gate', CMatrix[type].eye(8, 8), List[Int, True](0, 1, 2))
    gate.control(List[Int, True](4, 5), List[Int, True](1, 3, 4, 3))
    assert_equal(gate.n_qubits, 5, 'control')
    assert_equal(len(gate.classical_controls), 3, 'control')
    gate = Gate[type]._measure(List[Int, True](1, 2, 3), List[Int, True](9, 5, 1))
    with assert_raises(contains='measurement'):
        gate.control(0)

def test_other():
    g1 = Gate[type]('test_gate', CMatrix[type].eye(8, 8), List[Int, True](0, 1, 2))
    assert_equal(g1.n_qubits, len(g1), 'len')
    g2 = Gate[type]('test_gate', CMatrix[type].eye(8, 8), List[Int, True](0, 1, 2))
    assert_equal(g1 == g2, True, '__eq__')
    g2 = Gate[type]('test_gate', CMatrix[type].eye(8, 8), List[Int, True](0, 1, 3))
    assert_equal(g1 == g2, False, '__eq__')
    g2 = Gate[type]('other_gate', CMatrix[type].eye(8, 8), List[Int, True](0, 1, 2))
    assert_equal(g1 == g2, False, '__eq__')
    g2 = Gate[type]('test_gate', CMatrix[type].eye(4, 4), List[Int, True](0, 2))
    assert_equal(g1 == g2, False, '__eq__')
    g2 = Gate[type]('test_gate', CMatrix[type].eye(8, 8), List[Int, True](0, 1, 2))
    g2.control(3)
    assert_equal(g1 == g2, False, '__eq__')
    g2 = Gate[type]('test_gate', CMatrix[type].eye(8, 8), List[Int, True](0, 1, 2))
    g2.control(clbits=List[Int, True](3))
    assert_equal(g1 == g2, False, '__eq__')
    g1 = Gate[type]('test_gate', CMatrix[type].eye(2, 2), List[Int, True](0))
    m = CMatrix[type].eye(2, 2)
    m[1, 1] = -1
    g2 = Gate[type]('test_gate', m, List[Int, True](0))
    assert_equal(g1 == g2, False, '__eq__')
    g2 = Gate[type]('test_gate', CMatrix[type].eye(2, 2) + 1e-15, List[Int, True](0))
    assert_equal(g1 == g2, True, '__eq__')

def test_unparameterized_single_qubit_gates():
    # X
    x = CMatrix[type](2, 2,
        0, 1, 
        1, 0,
    )
    gate = X[type](1)
    _assert_matrix_equal(gate.matrix, x, 'X')
    assert_equal(gate.name, 'X', 'X')
    assert_equal(gate.n_qubits, 1, 'X')
    assert_equal(len(gate.qubits), 1, 'X')
    assert_equal(gate.qubits[0], 1, 'X')
    assert_equal(len(gate.params), 0, 'X')
    assert_equal(len(gate.controls), 0, 'X')
    assert_equal(len(gate._measure_targs), 0, 'X')
    assert_equal(len(gate.classical_controls), 0, 'X')
    assert_equal(gate._is_measure, False, 'X')
    # Y
    y = CMatrix[type](2, 2,
        0, ComplexScalar[type](0, -1), 
        ComplexScalar[type](0, 1), 0,
    )
    gate = Y[type](1)
    _assert_matrix_equal(gate.matrix, y, 'Y')
    assert_equal(gate.name, 'Y', 'Y')
    assert_equal(gate.n_qubits, 1, 'Y')
    assert_equal(len(gate.qubits), 1, 'Y')
    assert_equal(gate.qubits[0], 1, 'Y')
    assert_equal(len(gate.params), 0, 'Y')
    assert_equal(len(gate.controls), 0, 'Y')
    assert_equal(len(gate._measure_targs), 0, 'Y')
    assert_equal(len(gate.classical_controls), 0, 'Y')
    assert_equal(gate._is_measure, False, 'Y')
    # Z
    z = CMatrix[type](2, 2,
        1, 0, 
        0, -1,
    )
    gate = Z[type](1)
    _assert_matrix_equal(gate.matrix, z, 'Z')
    assert_equal(gate.name, 'Z', 'Z')
    assert_equal(gate.n_qubits, 1, 'Z')
    assert_equal(len(gate.qubits), 1, 'Z')
    assert_equal(gate.qubits[0], 1, 'Z')
    assert_equal(len(gate.params), 0, 'Z')
    assert_equal(len(gate.controls), 0, 'Z')
    assert_equal(len(gate._measure_targs), 0, 'Z')
    assert_equal(len(gate.classical_controls), 0, 'Z')
    assert_equal(gate._is_measure, False, 'Z')
    # H
    inv_sqrt2 = 1 / sqrt(2.0).cast[type]()
    h = CMatrix[type](2, 2,
        inv_sqrt2, inv_sqrt2, 
        inv_sqrt2, -inv_sqrt2,
    )
    gate = H[type](1)
    _assert_matrix_equal(gate.matrix, h, 'H')
    assert_equal(gate.name, 'H', 'H')
    assert_equal(gate.n_qubits, 1, 'H')
    assert_equal(len(gate.qubits), 1, 'H')
    assert_equal(gate.qubits[0], 1, 'H')
    assert_equal(len(gate.params), 0, 'H')
    assert_equal(len(gate.controls), 0, 'H')
    assert_equal(len(gate._measure_targs), 0, 'H')
    assert_equal(len(gate.classical_controls), 0, 'H')
    assert_equal(gate._is_measure, False, 'H')
    # S
    s = CMatrix[type](2, 2,
        1, 0, 
        0, ComplexScalar[type](0, 1),
    )
    gate = S[type](1)
    _assert_matrix_equal(gate.matrix, s, 'S')
    assert_equal(gate.name, 'S', 'S')
    assert_equal(gate.n_qubits, 1, 'S')
    assert_equal(len(gate.qubits), 1, 'S')
    assert_equal(gate.qubits[0], 1, 'S')
    assert_equal(len(gate.params), 0, 'S')
    assert_equal(len(gate.controls), 0, 'S')
    assert_equal(len(gate._measure_targs), 0, 'S')
    assert_equal(len(gate.classical_controls), 0, 'S')
    assert_equal(gate._is_measure, False, 'S')
    # T
    t = CMatrix[type](2, 2,
        1, 0, 
        0, ComplexScalar[type](0, pi / 4).exp(),
    )
    gate = T[type](1)
    _assert_matrix_equal(gate.matrix, t, 'T')
    assert_equal(gate.name, 'T', 'T')
    assert_equal(gate.n_qubits, 1, 'T')
    assert_equal(len(gate.qubits), 1, 'T')
    assert_equal(gate.qubits[0], 1, 'T')
    assert_equal(len(gate.params), 0, 'T')
    assert_equal(len(gate.controls), 0, 'T')
    assert_equal(len(gate._measure_targs), 0, 'T')
    assert_equal(len(gate.classical_controls), 0, 'T')
    assert_equal(gate._is_measure, False, 'T')
    # I
    i = CMatrix[type](2, 2,
        1, 0, 
        0, 1,
    )
    gate = I[type](1)
    _assert_matrix_equal(gate.matrix, i, 'I')
    assert_equal(gate.name, 'I', 'I')
    assert_equal(gate.n_qubits, 1, 'I')
    assert_equal(len(gate.qubits), 1, 'I')
    assert_equal(gate.qubits[0], 1, 'I')
    assert_equal(len(gate.params), 0, 'I')
    assert_equal(len(gate.controls), 0, 'I')
    assert_equal(len(gate._measure_targs), 0, 'I')
    assert_equal(len(gate.classical_controls), 0, 'I')
    assert_equal(gate._is_measure, False, 'I')
    # SX
    sx = CMatrix[type](2, 2,
        ComplexScalar[type](1, 1), ComplexScalar[type](1, -1), 
        ComplexScalar[type](1, -1), ComplexScalar[type](1, 1),
    ) / 2
    gate = SX[type](1)
    _assert_matrix_equal(gate.matrix, sx, 'SX')
    assert_equal(gate.name, 'SX', 'SX')
    assert_equal(gate.n_qubits, 1, 'SX')
    assert_equal(len(gate.qubits), 1, 'SX')
    assert_equal(gate.qubits[0], 1, 'SX')
    assert_equal(len(gate.params), 0, 'SX')
    assert_equal(len(gate.controls), 0, 'SX')
    assert_equal(len(gate._measure_targs), 0, 'SX')
    assert_equal(len(gate.classical_controls), 0, 'SX')
    assert_equal(gate._is_measure, False, 'SX')
    

def test_unparameterized_multi_qubit_gates():
    # I
    i = CMatrix[type].eye(16, 16)
    gate = I[type](List[Int, True](0, 3, 6, 77))
    _assert_matrix_equal(gate.matrix, i, 'I')
    assert_equal(gate.name, 'I', 'I')
    assert_equal(gate.n_qubits, 4, 'I')
    assert_equal(len(gate.qubits), 4, 'I')
    assert_equal(len(gate.params), 0, 'I')
    assert_equal(len(gate.controls), 0, 'I')
    assert_equal(len(gate._measure_targs), 0, 'I')
    assert_equal(len(gate.classical_controls), 0, 'I')
    assert_equal(gate._is_measure, False, 'I')
    # CX
    from ember import swap_rows_inplace
    x = CMatrix[type].eye(2, 2)
    swap_rows_inplace(x, 0, 1)
    gate = CX[type](0, 3)
    _assert_matrix_equal(gate.matrix, x, 'X')
    assert_equal(gate.name, 'X', 'X')
    assert_equal(gate.n_qubits, 2, 'X')
    assert_equal(len(gate.qubits), 1, 'X')
    assert_equal(gate.qubits[0], 3, 'X')
    assert_equal(len(gate.params), 0, 'X')
    assert_equal(len(gate.controls), 1, 'X')
    assert_equal(gate.controls[0], 0, 'X')
    assert_equal(len(gate._measure_targs), 0, 'X')
    assert_equal(len(gate.classical_controls), 0, 'X')
    assert_equal(gate._is_measure, False, 'X')
    # CCX
    gate = CCX[type](0, 1, 3)
    _assert_matrix_equal(gate.matrix, x, 'X')
    assert_equal(gate.name, 'X', 'X')
    assert_equal(gate.n_qubits, 3, 'X')
    assert_equal(len(gate.qubits), 1, 'X')
    assert_equal(gate.qubits[0], 3, 'X')
    assert_equal(len(gate.params), 0, 'X')
    assert_equal(len(gate.controls), 2, 'X')
    assert_equal(gate.controls[0], 0, 'X')
    assert_equal(gate.controls[1], 1, 'X')
    assert_equal(len(gate._measure_targs), 0, 'X')
    assert_equal(len(gate.classical_controls), 0, 'X')
    assert_equal(gate._is_measure, False, 'X')
    # ECR
    gate = ECR[type](0, 1)
    ecr = CMatrix[type](4, 4,
    0, 1, 0, ComplexScalar[type].i(),
    1, 0, -ComplexScalar[type].i(), 0,
    0, ComplexScalar[type].i(), 0, 1,
    -ComplexScalar[type].i(), 0, 1, 0,
    ) / sqrt(2.0)
    _assert_matrix_equal(gate.matrix, ecr, 'ECR')
    assert_equal(gate.name, 'ECR', 'ECR')
    assert_equal(gate.n_qubits, 2, 'ECR')
    assert_equal(len(gate.qubits), 2, 'ECR')
    assert_equal(gate.qubits[0], 0, 'ECR')
    assert_equal(gate.qubits[1], 1, 'ECR')
    assert_equal(len(gate.params), 0, 'ECR')
    assert_equal(len(gate.controls), 0, 'ECR')
    assert_equal(len(gate._measure_targs), 0, 'ECR')
    assert_equal(len(gate.classical_controls), 0, 'ECR')
    assert_equal(gate._is_measure, False, 'ECR')

def test_parameterized_single_qubit_gates():
    # RX
    t = 247.55
    a = ComplexScalar[type](cos(t / 2), 0)
    b = ComplexScalar[type](0, -sin(t / 2))
    rx = CMatrix[type](2, 2,
        a, b,
        b, a,
    )
    gate = RX[type](12, t)
    _assert_matrix_equal(gate.matrix, rx, 'RX')
    assert_equal(gate.name, 'RX', 'RX')
    assert_equal(gate.n_qubits, 1, 'RX')
    assert_equal(len(gate.qubits), 1, 'RX')
    assert_equal(gate.qubits[0], 12, 'RX')
    assert_equal(len(gate.params), 1, 'RX')
    assert_equal(gate.params[0], t, 'RX')
    assert_equal(len(gate.controls), 0, 'RX')
    assert_equal(len(gate._measure_targs), 0, 'RX')
    assert_equal(len(gate.classical_controls), 0, 'RX')
    assert_equal(gate._is_measure, False, 'RX')
    # RY
    a = ComplexScalar[type](cos(t / 2), 0)
    b = ComplexScalar[type](sin(t / 2), 0)
    ry = CMatrix[type](2, 2,
        a, -b,
        b, a,
    )
    gate = RY[type](12, t)
    _assert_matrix_equal(gate.matrix, ry, 'RY')
    assert_equal(gate.name, 'RY', 'RY')
    assert_equal(gate.n_qubits, 1, 'RY')
    assert_equal(len(gate.qubits), 1, 'RY')
    assert_equal(gate.qubits[0], 12, 'RY')
    assert_equal(len(gate.params), 1, 'RY')
    assert_equal(gate.params[0], t, 'RY')
    assert_equal(len(gate.controls), 0, 'RY')
    assert_equal(len(gate._measure_targs), 0, 'RY')
    assert_equal(len(gate.classical_controls), 0, 'RY')
    assert_equal(gate._is_measure, False, 'RY')
    # RZ
    c = cos(t / 2)
    s = sin(t / 2)
    rz = CMatrix[type](2, 2,
        ComplexScalar[type](c, -s), 0, 
        0, ComplexScalar[type](c, s),
    )
    gate = RZ[type](12, t)
    _assert_matrix_equal(gate.matrix, rz, 'RZ')
    assert_equal(gate.name, 'RZ', 'RZ')
    assert_equal(gate.n_qubits, 1, 'RZ')
    assert_equal(len(gate.qubits), 1, 'RZ')
    assert_equal(gate.qubits[0], 12, 'RZ')
    assert_equal(len(gate.params), 1, 'RZ')
    assert_equal(gate.params[0], t, 'RZ')
    assert_equal(len(gate.controls), 0, 'RZ')
    assert_equal(len(gate._measure_targs), 0, 'RZ')
    assert_equal(len(gate.classical_controls), 0, 'RZ')
    assert_equal(gate._is_measure, False, 'RZ')
    # U
    t = 52.22
    p = -442.2
    l = 0.345
    ct = cos(t / 2)
    st = sin(t / 2)
    u = CMatrix[type](2, 2,
        ComplexScalar[type](ct, 0), 
        ComplexScalar[type](-cos(l) * st, -sin(l) * st), 
        ComplexScalar[type](cos(p) * st, sin(p) * st), 
        ComplexScalar[type](cos(p + l) * ct, sin(p + l) * ct),
    )
    gate = U[type](12, t, p, l)
    _assert_matrix_equal(gate.matrix, u, 'U')
    assert_equal(gate.name, 'U', 'U')
    assert_equal(gate.n_qubits, 1, 'U')
    assert_equal(len(gate.qubits), 1, 'U')
    assert_equal(gate.qubits[0], 12, 'U')
    assert_equal(len(gate.params), 3, 'U')
    assert_equal(gate.params[0], t, 'U')
    assert_equal(gate.params[1], p, 'U')
    assert_equal(gate.params[2], l, 'U')
    assert_equal(len(gate.controls), 0, 'U')
    assert_equal(len(gate._measure_targs), 0, 'U')
    assert_equal(len(gate.classical_controls), 0, 'U')
    assert_equal(gate._is_measure, False, 'U')
    # PHASE
    phase = CMatrix[type](2, 2,
        1, 0,
        0, ComplexScalar[type](0, t).exp(),
    )
    gate = PHASE[type](12, t)
    _assert_matrix_equal(gate.matrix, phase, 'PHASE')
    assert_equal(gate.name, 'PHASE', 'PHASE')
    assert_equal(gate.n_qubits, 1, 'PHASE')
    assert_equal(len(gate.qubits), 1, 'PHASE')
    assert_equal(gate.qubits[0], 12, 'PHASE')
    assert_equal(len(gate.params), 1, 'PHASE')
    assert_equal(gate.params[0], t, 'PHASE')
    assert_equal(len(gate.controls), 0, 'PHASE')
    assert_equal(len(gate._measure_targs), 0, 'PHASE')
    assert_equal(len(gate.classical_controls), 0, 'PHASE')
    assert_equal(gate._is_measure, False, 'PHASE')
    # R
    a = cos(t / 2)
    b = -ComplexScalar[type].i() * sin(t / 2)
    r = CMatrix[type](2, 2,
        a, b * ComplexScalar(0, -p).exp(),
        b * ComplexScalar(0, p).exp(), a,
    )
    gate = R[type](12, t, p)
    _assert_matrix_equal(gate.matrix, r, 'R')
    assert_equal(gate.name, 'R', 'R')
    assert_equal(gate.n_qubits, 1, 'R')
    assert_equal(len(gate.qubits), 1, 'R')
    assert_equal(gate.qubits[0], 12, 'R')
    assert_equal(len(gate.params), 2, 'R')
    assert_equal(gate.params[0], t, 'R')
    assert_equal(gate.params[1], p, 'R')
    assert_equal(len(gate.controls), 0, 'R')
    assert_equal(len(gate._measure_targs), 0, 'R')
    assert_equal(len(gate.classical_controls), 0, 'R')
    assert_equal(gate._is_measure, False, 'R')

def test_parameterized_multi_qubit_gates():
    # RXX
    t = -224.13
    a = cos(t / 2)
    b = -ComplexScalar[type].i() * sin(t / 2)
    rxx = CMatrix[type](4, 4,
        a, 0, 0, b,
        0, a, b, 0,
        0, b, a, 0,
        b, 0, 0, a,
    )
    gate = RXX[type](12, 19, t)
    _assert_matrix_equal(gate.matrix, rxx, 'RXX')
    assert_equal(gate.name, 'RXX', 'RXX')
    assert_equal(gate.n_qubits, 2, 'RXX')
    assert_equal(len(gate.qubits), 2, 'RXX')
    assert_equal(gate.qubits[0], 12, 'RXX')
    assert_equal(gate.qubits[1], 19, 'RXX')
    assert_equal(len(gate.params), 1, 'RXX')
    assert_equal(gate.params[0], t, 'RXX')
    assert_equal(len(gate.controls), 0, 'RXX')
    assert_equal(len(gate._measure_targs), 0, 'RXX')
    assert_equal(len(gate.classical_controls), 0, 'RXX')
    assert_equal(gate._is_measure, False, 'RXX')
    # RYY
    t = -224.13
    a = cos(t / 2)
    b = ComplexScalar[type].i() * sin(t / 2)
    ryy = CMatrix[type](4, 4,
        a, 0, 0, b,
        0, a, -b, 0,
        0, -b, a, 0,
        b, 0, 0, a,
    )
    gate = RYY[type](12, 19, t)
    _assert_matrix_equal(gate.matrix, ryy, 'RYY')
    assert_equal(gate.name, 'RYY', 'RYY')
    assert_equal(gate.n_qubits, 2, 'RYY')
    assert_equal(len(gate.qubits), 2, 'RYY')
    assert_equal(gate.qubits[0], 12, 'RYY')
    assert_equal(gate.qubits[1], 19, 'RYY')
    assert_equal(len(gate.params), 1, 'RYY')
    assert_equal(gate.params[0], t, 'RYY')
    assert_equal(len(gate.controls), 0, 'RYY')
    assert_equal(len(gate._measure_targs), 0, 'RYY')
    assert_equal(len(gate.classical_controls), 0, 'RYY')
    assert_equal(gate._is_measure, False, 'RYY')
    # RZZ
    t = -224.13
    e = ComplexScalar[type](0, t / 2).exp()
    b = ComplexScalar[type](0, -t / 2).exp()
    rzz = CMatrix[type](4, 4,
        b, 0, 0, 0,
        0, e, 0, 0,
        0, 0, e, 0,
        0, 0, 0, b,
    )
    gate = RZZ[type](12, 19, t)
    _assert_matrix_equal(gate.matrix, rzz, 'RZZ')
    assert_equal(gate.name, 'RZZ', 'RZZ')
    assert_equal(gate.n_qubits, 2, 'RZZ')
    assert_equal(len(gate.qubits), 2, 'RZZ')
    assert_equal(gate.qubits[0], 12, 'RZZ')
    assert_equal(gate.qubits[1], 19, 'RZZ')
    assert_equal(len(gate.params), 1, 'RZZ')
    assert_equal(gate.params[0], t, 'RZZ')
    assert_equal(len(gate.controls), 0, 'RZZ')
    assert_equal(len(gate._measure_targs), 0, 'RZZ')
    assert_equal(len(gate.classical_controls), 0, 'RZZ')
    assert_equal(gate._is_measure, False, 'RZZ')
    # RZX
    t = -224.13
    a = cos(t / 2)
    b = ComplexScalar[type].i() * sin(t / 2)
    rzx = CMatrix[type](4, 4,
        a, 0, -b, 0,
        0, a, 0, b,
        -b, 0, a, 0,
        0, b, 0, a,
    )
    gate = RZX[type](12, 19, t)
    _assert_matrix_equal(gate.matrix, rzx, 'RZX')
    assert_equal(gate.name, 'RZX', 'RZX')
    assert_equal(gate.n_qubits, 2, 'RZX')
    assert_equal(len(gate.qubits), 2, 'RZX')
    assert_equal(gate.qubits[0], 12, 'RZX')
    assert_equal(gate.qubits[1], 19, 'RZX')
    assert_equal(len(gate.params), 1, 'RZX')
    assert_equal(gate.params[0], t, 'RZX')
    assert_equal(len(gate.controls), 0, 'RZX')
    assert_equal(len(gate._measure_targs), 0, 'RZX')
    assert_equal(len(gate.classical_controls), 0, 'RZX')
    assert_equal(gate._is_measure, False, 'RZX')