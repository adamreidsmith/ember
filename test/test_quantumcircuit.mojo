from collections import Dict

from testing import assert_equal, assert_raises
from ._testing import _assert_matrix_equal

from src import Gate, Measure, X, Y, T, CCX
from src import QuantumCircuit


alias type = DType.float64

def run_quantumcircuit_tests():
    print('Running quantumcircuit tests')
    test_init()
    test_set_cbits()
    test_apply()
    print('All tests passed')

def test_init():
    with assert_raises(contains='at least one qubit'):
        qc = QuantumCircuit[type](0, 0)
    qc = QuantumCircuit[type](12, 15)
    assert_equal(qc.n_qubits, 12, 'init')
    assert_equal(qc.n_clbits, 15, 'init')
    assert_equal(len(qc.clbits), 15, 'init')
    for c in qc.clbits:
        assert_equal(c[], 0, 'init')
    assert_equal(len(qc._data), 0, 'init')

def test_set_cbits():
    qc = QuantumCircuit[type](2, 3)
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
    qc = QuantumCircuit[type](2, 3)
    cbits_dict = Dict[Int, Int]()
    cbits_dict[1] = 1
    qc.set_clbits(cbits_dict)
    assert_equal(qc.clbits[0], 0, 'set_clbits')
    assert_equal(qc.clbits[1], 1, 'set_clbits')
    assert_equal(qc.clbits[2], 0, 'set_clbits')
    qc = QuantumCircuit[type](2, 3)
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
    qc = QuantumCircuit[type](3, 3)
    gate = X[type](-1)
    with assert_raises(contains='invalid qubit specifier'):
        qc.apply(gate)
    gate = Y[type](4)
    with assert_raises(contains='invalid qubit specifier'):
        qc.apply(gate)
    m = Measure[type](0, -1)
    with assert_raises(contains='invalid classical bit specifier'):
        qc.apply(m)
    m = Measure[type](0, 4)
    with assert_raises(contains='invalid classical bit specifier'):
        qc.apply(m)
    gate = X[type](0)
    _ = gate.control(-1)
    with assert_raises(contains='invalid qubit specifier'):
        qc.apply(gate)
    gate = X[type](0)
    _ = gate.control(5)
    with assert_raises(contains='invalid qubit specifier'):
        qc.apply(gate)
    gate = X[type](0)
    _ = gate.control(clbits=List[Int, True](-1))
    with assert_raises(contains='invalid classical bit specifier'):
        qc.apply(gate)
    gate = X[type](0)
    _ = gate.control(clbits=List[Int, True](5))
    with assert_raises(contains='invalid classical bit specifier'):
        qc.apply(gate)
    gates = List[Gate[type]](X[type](0), T[type](2), CCX[type](0, 2, 1))
    for g in gates:
        qc.apply(g[])
    for i in range(len(gates)):
        assert_equal(gates[i], qc._data[i], 'apply')