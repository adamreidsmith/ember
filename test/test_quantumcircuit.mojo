from collections import Dict, Set
from math import sqrt

from testing import assert_equal, assert_raises, assert_almost_equal
from ._testing import _assert_matrix_equal

from ember import Gate, Measure, X, Y, T, CCX, RZ, H
from ember import QuantumCircuit, CMatrix, CSRCMatrix


alias type = DType.float64

def run_quantumcircuit_tests():
    print('Running quantumcircuit tests')
    test_init()
    test_set_cbits()
    test_apply()
    test_set_initial_state()
    test_join()
    print('All tests passed')

def test_init():
    with assert_raises(contains='at least one qubit'):
        qc = QuantumCircuit(0, 0)
    qc = QuantumCircuit(12, 15)
    assert_equal(qc.n_qubits, 12, 'init')
    assert_equal(qc.n_clbits, 15, 'init')
    assert_equal(len(qc.clbits), 15, 'init')
    for c in qc.clbits:
        assert_equal(c[], 0, 'init')
    assert_equal(len(qc._data), 0, 'init')

def test_set_cbits():
    qc = QuantumCircuit(2, 3)
    cbits_dict = Dict[Int, Int]()
    cbits_dict[0] = 0
    cbits_dict[1] = 2
    cbits_dict[2] = 1
    with assert_raises(contains='set to 0 or 1'):
        qc.set_clbits(cbits_dict)
    cbits_dict[1] = 1
    cbits_dict[3] = 0
    with assert_raises(contains='not in the circuit'):
        qc.set_clbits(cbits_dict)
    _ = cbits_dict.pop(3)
    qc.set_clbits(cbits_dict)
    assert_equal(qc.clbits[0], cbits_dict[0], 'set_clbits')
    assert_equal(qc.clbits[1], cbits_dict[1], 'set_clbits')
    assert_equal(qc.clbits[2], cbits_dict[2], 'set_clbits')
    qc = QuantumCircuit(2, 3)
    cbits_dict = Dict[Int, Int]()
    cbits_dict[1] = 1
    qc.set_clbits(cbits_dict)
    assert_equal(qc.clbits[0], 0, 'set_clbits')
    assert_equal(qc.clbits[1], 1, 'set_clbits')
    assert_equal(qc.clbits[2], 0, 'set_clbits')
    qc = QuantumCircuit(2, 3)
    cbits = List[Int, True](0, 1, 2)
    with assert_raises(contains='set to 0 or 1'):
        qc.set_clbits(cbits)
    _ = cbits.pop()
    with assert_raises(contains='same length'):
        qc.set_clbits(cbits)
    cbits.append(0)
    qc.set_clbits(cbits)
    assert_equal(qc.clbits[0], cbits[0], 'set_clbits')
    assert_equal(qc.clbits[1], cbits[1], 'set_clbits')
    assert_equal(qc.clbits[2], cbits[2], 'set_clbits')

def test_apply():
    qc = QuantumCircuit(3, 3)
    gate = X(-1)
    with assert_raises(contains='invalid qubit specifier'):
        qc.apply(gate)
    gate = Y(4)
    with assert_raises(contains='invalid qubit specifier'):
        qc.apply(gate)
    m = Measure(0, -1)
    with assert_raises(contains='invalid classical bit specifier'):
        qc.apply(m)
    m = Measure(0, 4)
    with assert_raises(contains='invalid classical bit specifier'):
        qc.apply(m)
    gate = X(0)
    _ = gate.control(-1)
    with assert_raises(contains='invalid qubit specifier'):
        qc.apply(gate)
    gate = X(0)
    _ = gate.control(5)
    with assert_raises(contains='invalid qubit specifier'):
        qc.apply(gate)
    gate = X(0)
    _ = gate.control(clbits=List[Int, True](-1))
    with assert_raises(contains='invalid classical bit specifier'):
        qc.apply(gate)
    gate = X(0)
    _ = gate.control(clbits=List[Int, True](5))
    with assert_raises(contains='invalid classical bit specifier'):
        qc.apply(gate)
    gates = List[Gate](X(0), T(2), CCX(0, 2, 1))
    for g in gates:
        qc.apply(g[])
    for i in range(len(gates)):
        assert_equal(gates[i], qc._data[i], 'apply')

def test_set_initial_state():
    qc = QuantumCircuit(3)
    statevector = CSRCMatrix[type](CMatrix[type].arange(3, 3))
    with assert_raises(contains='Expected 1D'):
        qc.set_initial_state(statevector)
    statevector = CSRCMatrix[type](CMatrix[type].arange(1, 9))
    with assert_raises(contains='Statevector must have'):
        qc.set_initial_state(statevector)
    statevector = CSRCMatrix[type](CMatrix[type].arange(1, 8))
    with assert_raises(contains='not normalized'):
        qc.set_initial_state(statevector)
    qc.set_initial_state(statevector, normalize=True)
    var t: Scalar[type] = 0
    for i in range(statevector.size):
        t += statevector[0, i].squared_norm()
    svn = statevector / sqrt(t)
    for i in range(svn.size):
        assert_almost_equal(svn[0, i].re, qc._initial_state[i].re, 'set_initial_state', atol=1e-20, rtol=1e-20)
        assert_almost_equal(svn[0, i].im, qc._initial_state[i].im, 'set_initial_state', atol=1e-20, rtol=1e-20)
    qc.set_initial_state(svn.to_list(), normalize=True)
    for i in range(svn.size):
        assert_almost_equal(svn[0, i].re, qc._initial_state[i].re, 'set_initial_state', atol=1e-20, rtol=1e-20)
        assert_almost_equal(svn[0, i].im, qc._initial_state[i].im, 'set_initial_state', atol=1e-20, rtol=1e-20)
    qc.set_initial_state(svn.to_dense(), normalize=True)
    for i in range(svn.size):
        assert_almost_equal(svn[0, i].re, qc._initial_state[i].re, 'set_initial_state', atol=1e-20, rtol=1e-20)
        assert_almost_equal(svn[0, i].im, qc._initial_state[i].im, 'set_initial_state', atol=1e-20, rtol=1e-20)

def test_join():
    qc1 = QuantumCircuit(4, 4)
    rz = RZ(1, 0.665)
    rz.control(3)
    rz.control(clbits=List[Int, True](2, 1, 3))
    qc1.apply(X(0), CCX(0, 1, 3), Measure(List[Int, True](2, 3, 0), List[Int, True](2, 1, 3)), rz)
    qc2 = QuantumCircuit(5, 3)
    with assert_raises(contains='Cannot join'):
        qc1.join(qc2, List[Int, True]())
    qc2 = QuantumCircuit(3, 5)
    with assert_raises(contains='Cannot join'):
        qc1.join(qc2, List[Int, True]())
    qc2 = QuantumCircuit(3, 3)
    with assert_raises(contains='Number of qubits'):
        qc1.join(qc2, List[Int, True](0, 1))
    with assert_raises(contains='Number of qubits'):
        qc1.join(qc2, List[Int, True](0, 1, 2, 3))
    with assert_raises(contains='Number of classical'):
        qc1.join(qc2, List[Int, True](0, 1, 3), List[Int, True](0))
    with assert_raises(contains='Number of classical'):
        qc1.join(qc2, List[Int, True](0, 1, 3), List[Int, True](0, 1, 2, 3))
    with assert_raises(contains='not unique'):
        qc1.join(qc2, List[Int, True](0, 1, 1), List[Int, True](0, 1, 2))
    with assert_raises(contains='not unique'):
        qc1.join(qc2, List[Int, True](0, 1, 2), List[Int, True](0, 0, 2))
    with assert_raises(contains='Invalid qubit'):
        qc1.join(qc2, List[Int, True](0, 1, 5), List[Int, True](0, 1, 2))
    with assert_raises(contains='Invalid classical'):
        qc1.join(qc2, List[Int, True](0, 1, 2), List[Int, True](0, 1, 4))
    qc2 = QuantumCircuit(3, 3)
    cx = X(0)
    cx.control(1, 2)
    cx.control(clbits=List[Int, True](2, 1))
    qc2.apply(Measure(2, 2), H(0), cx, Measure(List[Int, True](0, 1, 2), List[Int, True](0, 1, 2)))
    qc = qc1  # Copy qc1
    qc.join(qc2, List[Int, True](3, 2, 1), List[Int, True](1, 2, 3))
    for i in range(len(qc1._data)):
        assert_equal(qc._data[i], qc1._data[i], 'join')
    assert_equal(qc._data[4], Measure(1, 3), 'join')
    assert_equal(qc._data[5], H(3), 'join')
    x = X(3)
    x.control(2, 1)
    x.control(clbits=List[Int, True](3, 2)) 
    assert_equal(qc._data[6], x, 'join')
    assert_equal(qc._data[7], Measure(List[Int, True](3, 2, 1), List[Int, True](1, 2, 3)), 'join')
