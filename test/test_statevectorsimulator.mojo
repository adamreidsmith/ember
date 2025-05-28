from collections import Dict
from math import sqrt, cos, sin, pi

from testing import assert_equal, assert_raises, assert_almost_equal
from ._testing import _assert_matrix_equal, _assert_matrix_almost_equal

from ember import Gate, Measure, I, X, Y, Z, H, S, T, CX, CCX, RX, RY, RZ, U
from ember import QuantumCircuit, CMatrix, ComplexScalar
from ember import StatevectorSimulator, Statevector

alias type = DType.float64

def run_statevectorsimulator_tests():
    print('Running statevectorsimulator tests')
    test_init()
    test_run()
    test_measure()
    print('All tests passed')

def test_init():
    sim = StatevectorSimulator()
    assert_equal(sim._sv.size, 0, 'init')
    assert_equal(len(sim._cb), 0, 'init')
    assert_equal(sim._qc.n_qubits, 0, 'init')
    assert_equal(sim._qc.n_clbits, 0, 'init')
    sim = StatevectorSimulator(15)
    assert_equal(sim._sv.size, 0, 'init')
    assert_equal(len(sim._cb), 0, 'init')
    assert_equal(sim._qc.n_qubits, 0, 'init')
    assert_equal(sim._qc.n_clbits, 0, 'init')

def test_run():
    # Uniform superposition state
    for n_qubits in range(1, 11):
        qc = QuantumCircuit(n_qubits)
        for i in range(n_qubits):
            qc.apply(H(i))
        sim = StatevectorSimulator()
        sim.run[False](qc)
        sv = sim.get_statevector()
        amp = sqrt(1.0 / (2.0 ** n_qubits))
        for i in range(len(sv)):
            assert_almost_equal(amp, sv[i].re, 'run')
            assert_almost_equal(0, sv[i].im, 'run')
        sim.run[True](qc)
        sv = sim.get_statevector()
        amp = sqrt(1.0 / (2.0 ** n_qubits))
        for i in range(len(sv)):
            assert_almost_equal(amp, sv[i].re, 'run')
            assert_almost_equal(0, sv[i].im, 'run')

    # Random circuit
    qc = QuantumCircuit(4)
    qc.apply(H(0), H(3))
    qc.apply(RY(1, 0.44), RY(2, 2.25))
    qc.apply(CCX(1, 2, 0), CX(3, 2))
    qc.apply(T(2), S(1), H(0), Z(3))
    rz = RZ(2, 1.5)
    rz.control(0, 1, 3)
    rx = RX(1, 5.55)
    rx.control(0, 3)
    qc.apply(rz, rx)
    true_sv = List[ComplexScalar[DType.float64]](
        (2.97539264e-1, +1.38777878e-17),
        (0.0, +0.0),
        (0.0, +6.65355582e-2),
        (0.0, +0.0),
        (4.40260322e-1, +4.40260322e-1),
        (0.0, +0.0),
        (-9.84507585e-2, +9.84507585e-2),
        (0.0, +0.0),
        (-6.22622118e-1, +1.38777878e-17),
        (0.0, +0.0),
        (3.46944695e-18, -1.39230398e-1),
        (0.0, +0.0),
        (-2.10392031e-1, -2.10392031e-1),
        (0.0, +0.0),
        (4.70477444e-2, -4.70477444e-2),
        (0.0, +0.0)
    )
    sim = StatevectorSimulator()
    sim.run[False](qc)
    sv = sim.get_statevector()
    assert_equal(len(true_sv), len(sv), 'run')
    for i in range(len(sv)):
        assert_almost_equal(true_sv[i].re, sv[i].re, 'run')
        assert_almost_equal(true_sv[i].im, sv[i].im, 'run')
    sim.run[True](qc)
    sv = sim.get_statevector()
    assert_equal(len(true_sv), len(sv), 'run')
    for i in range(len(sv)):
        assert_almost_equal(true_sv[i].re, sv[i].re, 'run')
        assert_almost_equal(true_sv[i].im, sv[i].im, 'run')
    
    # Random circuit
    qc = QuantumCircuit(7, 3)
    qc.set_clbits(List[Int, True](1, 0, 1))
    for i in range(7):
        qc.apply(H(i))
    u = U(0, 0.2, -3.1, 1.7)
    u.control(qubits=List[Int, True](1, 2, 3, 4, 5, 6), clbits=List[Int, True](0))
    qc.apply(u)
    rzz_matrix = CMatrix(4, 4, fill_zeros=True)
    rzz_matrix[0, 0] = (1, -1)
    rzz_matrix[1, 1] = (1, 1)
    rzz_matrix[2, 2] = (1, 1)
    rzz_matrix[3, 3] = (1, -1)
    rzz_matrix /= sqrt(2.0)
    rzz = Gate('RZZ', rzz_matrix, List[Int, True](1, 3), List[Scalar[DType.float64], True](pi / 2))
    rzz.control(2, 4)
    qc.apply(rzz)
    theta = 1.95
    st = sin(theta / 2)
    ct = cos(theta / 2)
    ryy_matrix = CMatrix(4, 4,
        ct, 0, 0, (0.0, st),
        0, ct, (0.0, -st), 0,
        0, (0.0, -st), ct, 0,
        (0.0, st), 0, 0, ct,
    )
    ryy = Gate('RYY', ryy_matrix, List[Int, True](3, 5), List[Scalar[DType.float64], True](theta))
    ryy.control(2, 4, 6)
    ryy.control(clbits=List[Int, True](2))
    qc.apply(ryy)
    ryy.control(clbits=List[Int, True](1))
    qc.apply(ryy)
    qc.apply(T(0), T(5), S(3), Z(2), X(3), H(4), Y(1), S(6))
    theta = 0.7
    st = sin(theta / 2)
    ct = cos(theta / 2)
    ryy_matrix = CMatrix(4, 4,
        ct, 0, 0, (0.0, st),
        0, ct, (0.0, -st), 0,
        0, (0.0, -st), ct, 0,
        (0.0, st), 0, 0, ct,
    )
    ryy = Gate('RYY', ryy_matrix, List[Int, True](1, 0), List[Scalar[DType.float64], True](theta))
    qc.apply(ryy)
    sim = StatevectorSimulator()
    true_sv = List[ComplexScalar[DType.float64]](
        ( 1.47729760e-1, -3.03081706e-2),
        (8.30296019e-2, +1.25891828e-1),
        (-8.71134185e-2, -3.03081706e-2),
        (-8.30296019e-2, -4.01673760e-2),
        (-1.36810794e-1, +5.66688863e-2),
        (-8.50715102e-2, -7.80999992e-2),
        (8.50715102e-2, +7.80999992e-2),
        (2.63607157e-2, +6.36403972e-2),
        (-3.03081706e-2, -1.47729760e-1),
        (1.25891828e-1, -8.30296019e-2),
        (-3.03081706e-2, +8.71134185e-2),
        (-4.01673760e-2, +8.30296019e-2),
        (-4.92960268e-3, +1.15379681e-1),
        (-1.36810794e-1, +5.66688863e-2),
        (-2.63607157e-2, -6.36403972e-2),
        (4.92960268e-3, -1.15379681e-1),
        (1.22706187e-16, -6.24696129e-17),
        (5.66689765e-17, +5.92836552e-17),
        (-5.19970725e-19, -4.71150545e-17),
        (-2.89039309e-17, +7.85759765e-18),
        (-1.09189659e-2, -2.63607157e-2),
        (2.04190829e-3, -4.77918286e-2),
        (2.04190829e-3, -4.77918286e-2),
        (5.66688863e-2, -2.34730213e-2),
        (-2.42420977e-17, -8.30769703e-17),
        (5.17255062e-17, -9.20961270e-17),
        (-3.10061067e-17, +1.21466613e-16),
        (-6.57020466e-17, +1.09338646e-16),
        (3.52377733e-2, +3.23500789e-2),
        (1.09189659e-2, +2.63607157e-2),
        (5.66688863e-2, -2.34730213e-2),
        (3.52377733e-2, +3.23500789e-2),
        (1.25891828e-1, +8.30296019e-2),
        (-3.03081706e-2, +1.47729760e-1),
        (-4.01673760e-2, -8.30296019e-2),
        (-3.03081706e-2, -8.71134185e-2),
        (-1.36810794e-1, -5.66688863e-2),
        (-4.92960268e-3, -1.15379681e-1),
        (4.92960268e-3, +1.15379681e-1),
        (-2.63607157e-2, +6.36403972e-2),
        (8.30296019e-2, -1.25891828e-1),
        (1.47729760e-1, +3.03081706e-2),
        (-8.30296019e-2, +4.01673760e-2),
        (-8.71134185e-2, +3.03081706e-2),
        (-8.50715102e-2, +7.80999992e-2),
        (-1.36810794e-1, -5.66688863e-2),
        (2.63607157e-2, -6.36403972e-2),
        (8.50715102e-2, -7.80999992e-2),
        (1.18966575e-17, -6.51820759e-17),
        (5.53088255e-17, +2.43132858e-17),
        (-5.89842066e-18, -2.01892973e-17),
        (2.37933151e-17, -3.25910380e-17),
        (1.09189659e-2, -2.63607157e-2),
        (3.52377733e-2, -3.23500789e-2),
        (3.52377733e-2, -3.23500789e-2),
        (5.66688863e-2, +2.34730213e-2),
        (9.67153230e-18, -6.88009813e-17),
        (1.22186216e-16, -1.30792353e-17),
        (-7.86416074e-17, +5.14260576e-17),
        (-2.65928011e-17, +1.39768006e-17),
        (2.04190829e-3, +4.77918286e-2),
        (-1.09189659e-2, +2.63607157e-2),
        (5.66688863e-2, +2.34730213e-2),
        (2.04190829e-3, +4.77918286e-2),
        (3.03081706e-2, +1.47729760e-1),
        (-1.25891828e-1, +8.30296019e-2),
        (3.03081706e-2, -8.71134185e-2),
        (4.01673760e-2, -8.30296019e-2),
        (-7.28127445e-2, -1.25811422e-1),
        (5.89068375e-2, -1.27095535e-1),
        (-3.38540148e-2, +3.24916102e-2),
        (1.66117788e-2, +2.04677474e-2),
        (1.47729760e-1, -3.03081706e-2),
        (8.30296019e-2, +1.25891828e-1),
        (-8.71134185e-2, -3.03081706e-2),
        (-8.30296019e-2, -4.01673760e-2),
        (-1.35195739e-1, -4.50851505e-2),
        (-8.53261735e-2, -9.52245430e-2),
        (1.92346128e-2, +4.36859618e-2),
        (4.57703399e-2, -2.30382335e-3),
        (6.75794685e-17, +7.40796150e-17),
        (-7.71785496e-17, +5.25687771e-17),
        (-2.25965613e-17, -9.00153324e-17),
        (-1.96544390e-17, -4.91988822e-17),
        (4.25045739e-2, -2.19183374e-2),
        (6.69849904e-2, +4.40659333e-2),
        (3.54584422e-3, +5.46218083e-2),
        (-5.67791548e-2, +6.25618545e-2),
        (9.83929008e-17, -2.59482552e-17),
        (3.17498172e-17, +2.25537412e-17),
        (-2.13141674e-17, -1.42011919e-17),
        (-4.96219986e-17, +1.37560183e-17),
        (-1.25340203e-2, +7.53933211e-2),
        (2.29657161e-3, -3.06672848e-2),
        (6.78788057e-2, -1.33777912e-2),
        (3.72592620e-2, +4.24711993e-2),
        (-8.30296019e-2, +1.25891828e-1),
        (-1.47729760e-1, -3.03081706e-2),
        (8.30296019e-2, -4.01673760e-2),
        (8.71134185e-2, -3.03081706e-2),
        (5.93000966e-2, -6.86037504e-2),
        (1.81792998e-2, +5.47250532e-2),
        (-1.19477012e-1, +4.04105114e-2),
        (-1.06832871e-1, -4.36952085e-2),
        (1.25891828e-1, +8.30296019e-2),
        (-3.03081706e-2, +1.47729760e-1),
        (-4.01673760e-2, -8.30296019e-2),
        (-3.03081706e-2, -8.71134185e-2),
        (-8.39929675e-2, -4.81933418e-3),
        (4.08898621e-3, -9.25648094e-2),
        (1.05664422e-1, +4.55538774e-2),
        (6.71006278e-2, +1.01215368e-1),
        (-2.85914805e-17, +1.83949019e-16),
        (-1.39981294e-16, -5.60567966e-17),
        (5.34049858e-17, -5.97038098e-17),
        (1.39681845e-16, -8.35616151e-17),
        (2.37295054e-2, -5.72880775e-2),
        (1.29550460e-1, -2.44168825e-2),
        (3.64474096e-2, -2.43135464e-4),
        (1.97194527e-2, +7.40033791e-2),
        (1.30364152e-17, +3.72152459e-17),
        (1.53089725e-17, +1.42260809e-16),
        (1.49955922e-17, -8.68422565e-17),
        (-3.97006645e-17, +4.75866301e-18),
        (-4.18988603e-2, -7.82102677e-2),
        (2.62191844e-2, -5.51649503e-2),
        (-6.54970463e-2, +3.74757245e-2),
        (-3.67924572e-2, -1.41019499e-2),
    )
    sim.run[False](qc)
    sv = sim.get_statevector()
    assert_equal(len(sv), len(true_sv), 'run')
    for i in range(len(sv)):
        assert_almost_equal(true_sv[i].re, sv[i].re, 'run')
        assert_almost_equal(true_sv[i].im, sv[i].im, 'run')
    sim.run[True](qc)
    sv = sim.get_statevector()
    assert_equal(len(sv), len(true_sv), 'run')
    for i in range(len(sv)):
        assert_almost_equal(true_sv[i].re, sv[i].re, 'run')
        assert_almost_equal(true_sv[i].im, sv[i].im, 'run')

def test_measure():
    trials = 50
    for _ in range(trials):
        qc = QuantumCircuit(5, 5)
        qc.apply(Measure(List[Int, True](0, 1, 2, 3, 4), List[Int, True](0, 1, 2, 3, 4)))
        sim = StatevectorSimulator()
        sim.run(qc)
        for bit in sim._cb:
            assert_equal(bit[], 0, 'Measure')
    for _ in range(trials):
        qc = QuantumCircuit(5, 5)
        for q in range(5):
            qc.apply(X(q))
        qc.apply(Measure(List[Int, True](0, 1, 2, 3, 4), List[Int, True](0, 1, 2, 3, 4)))
        sim = StatevectorSimulator()
        sim.run(qc)
        for bit in sim._cb:
            assert_equal(bit[], 1, 'Measure')
    for i in range(trials):
        qc = QuantumCircuit(5, 5)
        # Set classical bits randomly to ensure they get overwritten by measurements
        clbits = List[Int, True]()
        for q in range(5):
            clbits.append((i >> q) & 1)
        qc.set_clbits(clbits)
        qc.apply(Measure(List[Int, True](0, 1, 2, 3, 4), List[Int, True](0, 1, 2, 3, 4)))
        sim = StatevectorSimulator()
        sim.run(qc)
        for bit in sim._cb:
            assert_equal(bit[], 0, 'Measure')
    
    # Random statevector
    SV = List[ComplexScalar[DType.float64], True](
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
    for item in SV:
        norm += item[].squared_norm()
    norm = sqrt(norm)
    for i in range(len(SV)):
        SV[i] /= norm
    qc = QuantumCircuit(4, 1)
    qc.set_initial_state(SV)
    true_probs = List[Tuple[Float64, Float64]]()
    true_probs.append((0.401465, 0.598535)) #0
    true_probs.append((0.481679, 0.518321)) #1
    true_probs.append((0.552895, 0.447105)) #2
    true_probs.append((0.461693, 0.538307)) #3
    true_svs = List[Tuple[Statevector[tol=1e-4], Statevector[tol=1e-4]]]()
    true_svs.append(( #0
        Statevector[tol=1e-4](CMatrix(16, 1,
            (9.83549677e-2, 0.30721174),
            (0.0, 0.0),
            (4.45744714e-2, 0.11015756),
            (0.0, 0.0),
            (3.07201906e-1, 0.13179566),
            (0.0, 0.0),
            (-3.39521349e-2, 0.58717916),
            (0.0, 0.0),
            (-1.10505740e-1, -0.31916187),
            (0.0, 0.0),
            (2.19837560e-1, 0.21638093),
            (0.0, 0.0),
            (-4.14408821e-1, -0.20654543),
            (0.0, 0.0),
            (1.20976610e-5, -0.02301506),
            (0.0, 0.0)
        )),
        Statevector[tol=1e-4](CMatrix(16, 1,
            (0.0, 0.0), 
            (0.09776907, -4.12742218e-1),
            (0.0, 0.0), 
            (-0.33217571, 0.0),
            (0.0, 0.0), 
            (0.12273502, -9.94271270e-2),
            (0.0, 0.0), 
            (0.20241046, 1.71042696e-1),
            (0.0, 0.0), 
            (0.18127464, -9.90968043e-5),
            (0.0, 0.0), 
            (0.51965397, 3.32207937e-1),
            (0.0, 0.0), 
            (0.10908705, 4.35059141e-1),
            (0.0, 0.0), 
            (0.00805665, 8.05665076e-3)
        )),
    ))
    true_svs.append(( #1
        Statevector[tol=1e-4](CMatrix(16, 1,
            (0.08984617, 2.80634501e-1),
            (0.10903012, -4.60281911e-1),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.28062552, 1.20393863e-1),
            (0.13687165, -1.10879154e-1),
            (0.0, 0.0),
            (0.0, 0.0),
            (-0.10094576, -2.91550811e-1),
            (0.20215387, -1.10510785e-4),
            (0.0, 0.0),
            (0.0, 0.0),
            (-0.37855784, -1.88676950e-1),
            (0.12165171, 4.85169300e-1),
            (0.0, 0.0),
            (0.0, 0.0)
        )),
        Statevector[tol=1e-4](CMatrix(16, 1,
            (0.0, 0.0),
            (0.0, 0.0),
            (3.92159044e-2, 0.09691486),
            (-3.56767815e-1, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (-2.98705432e-2, 0.5165908) ,
            (2.17395599e-1, 0.18370557),
            (0.0, 0.0),
            (0.0, 0.0),
            (1.93409556e-1, 0.19036847),
            (5.58125736e-1, 0.35680243),
            (0.0, 0.0),
            (0.0, 0.0),
            (1.06433280e-5, -0.02024828),
            (8.65311218e-3, 0.00865311)
        )),
    ))
    true_svs.append(( #2
        Statevector[tol=1e-4](CMatrix(16, 1,
            (0.08385067, 2.61907567e-1),
            (0.10175446, -4.29566980e-1),
            (0.03800112, 9.39127500e-2),
            (-0.34571631, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (-0.09420958, -2.72095423e-1),
            (0.18866401, -1.03136324e-4),
            (0.18741836, 1.84471473e-1),
            (0.54083682, 3.45749851e-1),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0)
        )),
        Statevector[tol=1e-4](CMatrix(16, 1,
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (2.90999508e-1, 0.12484451),
            (1.41931437e-1, -0.11497807),
            (-3.21614363e-2, 0.55621024), 
            (2.34068549e-1, 0.1977947),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (-3.92552132e-1, -0.19565184),
            (1.26148855e-1, 0.50310474),
            (1.14596080e-5, -0.02180121),
            (9.31675444e-3, 0.00931675),
        )),
    ))
    true_svs.append(( #3
        Statevector[tol=1e-4](CMatrix(16, 1,
            (0.09170571, 0.28644278),
            (0.11128671, -0.46980835),
            (0.04156103, 0.10271039),
            (-0.37810264, 0.0),
            (0.28643361, 0.12288565),
            (0.13970448, -0.11317402),
            (-0.03165681, 0.54748308),
            (0.23039592, 0.19469122),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        )),
        Statevector[tol=1e-4](CMatrix(16, 1,
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (-9.54601564e-2, -2.75707325e-1),
            (1.91168407e-1, -1.04505396e-4),
            (1.89906223e-1, 1.86920220e-1),
            (5.48016100e-1, 3.50339471e-1),
            (-3.57986207e-1, -1.78423846e-1),
            (1.15040899e-1, 4.58804177e-1),
            (1.04505396e-5, -1.98815143e-2),
            (8.49637364e-3, 8.49637364e-3),
        )),
    ))
    for q in range(4):
        qcn = qc
        qcn.apply(Measure(List[Int, True](q), List[Int, True](0)))
        sim = StatevectorSimulator()
        trials = 100_000
        counts = Dict[Int, Int]()
        counts[0] = 0
        counts[1] = 0
        sv0 = Statevector()
        sv1 = Statevector()
        for _ in range(trials):
            sim.run[False](qcn)
            clbit = sim.get_clbits()[0]
            counts[clbit] += 1
            if sv0.size == 0 and clbit == 0:
                sv0 = sim.get_statevector()
            if sv1.size == 0 and clbit == 1:
                sv1 = sim.get_statevector()
        prob0 = Float64(counts[0]) / Float64(trials)
        prob1 = Float64(counts[1]) / Float64(trials)
        assert_almost_equal(prob0, true_probs[q][0], 'Measure q=' + String(q), atol=1e-2)
        assert_almost_equal(prob1, true_probs[q][1], 'Measure q=' + String(q), atol=1e-2)
        for i in range(len(sv0)):
            assert_almost_equal(sv0[i].re, true_svs[q][0][i].re, 'Measure q=' + String(q), atol=1e-7)
            assert_almost_equal(sv0[i].im, true_svs[q][0][i].im, 'Measure q=' + String(q), atol=1e-7)
        for i in range(len(sv1)):
            assert_almost_equal(sv1[i].re, true_svs[q][1][i].re, 'Measure q=' + String(q), atol=1e-7)
            assert_almost_equal(sv1[i].im, true_svs[q][1][i].im, 'Measure q=' + String(q), atol=1e-7)

    qc = QuantumCircuit(4, 2)
    qc.set_initial_state(SV)
    qc.apply(Measure(List[Int, True](1, 3), List[Int, True](0, 1)))
    sim = StatevectorSimulator()
    trials = 100_000
    cts = Dict[String, Int]()
    cts['00'] = 0
    cts['01'] = 0
    cts['10'] = 0
    cts['11'] = 0
    sv00 = Statevector()
    sv01 = Statevector()
    sv10 = Statevector()
    sv11 = Statevector()
    true_pbs = List[Float64]()
    true_pbs.append(0.208773) #00
    true_pbs.append(0.252718) #01
    true_pbs.append(0.271921) #10
    true_pbs.append(0.266588) #11
    sv00_t = Statevector[tol=1e-4](CMatrix(16, 1,
        (0.13624814, 0.42557107),
        (0.16533984, -0.69799922),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.42555744, 0.18257251),
        (0.20756042, -0.16814383),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
    ))
    sv01_t = Statevector[tol=1e-4](CMatrix(16, 1,
        (0.0, 0.0),
        (0.0, 0.0),
        (0.05619613, 0.13887835), 
        (-0.51124593, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (-0.04280429, 0.74027122),
        (0.31152646, 0.26324888),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
    ))
    sv10_t = Statevector[tol=1e-4](CMatrix(16, 1,
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (-0.13427824, -3.87821423e-1),
        (0.26890546, -1.47001649e-4),
        (0.0, 0.0),
        (0.0, 0.0),
        (-0.50355833, -2.50978425e-1),
        (0.16182133, 6.45373093e-1),
        (0.0, 0.0),
        (0.0, 0.0),
    ))
    sv11_t = Statevector[tol=1e-4](CMatrix(16, 1,
        (0.0, 0.0), 
        (0.0, 0.0),
        (0.0, 0.0), 
        (0.0, 0.0),
        (0.0, 0.0), 
        (0.0, 0.0),
        (0.0, 0.0), 
        (0.0, 0.0),
        (0.0, 0.0), 
        (0.0, 0.0),
        (2.70029372e-1, 0.26578355),
        (7.79229036e-1, 0.49815085),
        (0.0, 0.0), 
        (0.0, 0.0),
        (1.48597165e-5, -0.0282697),
        (1.20810703e-2, 0.01208107),
    ))
    for _ in range(trials):
        sim.run(qc)
        clbits = sim.get_clbits()
        state_str = String(clbits[1]) + String(clbits[0])  # Reverse ordering of qubits
        cts[state_str] += 1
        if sv00.size == 0 and state_str == '00':
            sv00 = sim.get_statevector()
        if sv01.size == 0 and state_str == '01':
            sv01 = sim.get_statevector()
        if sv10.size == 0 and state_str == '10':
            sv10 = sim.get_statevector()
        if sv11.size == 0 and state_str == '11':
            sv11 = sim.get_statevector()
    prob00 = Float64(cts['00']) / Float64(trials)
    prob01 = Float64(cts['01']) / Float64(trials)
    prob10 = Float64(cts['10']) / Float64(trials)
    prob11 = Float64(cts['11']) / Float64(trials)
    assert_almost_equal(prob00, true_pbs[0], 'Measure c=00', atol=1e-2)
    assert_almost_equal(prob01, true_pbs[1], 'Measure c=01', atol=1e-2)
    assert_almost_equal(prob10, true_pbs[2], 'Measure c=10', atol=1e-2)
    assert_almost_equal(prob11, true_pbs[3], 'Measure c=11', atol=1e-2)
    for i in range(len(sv00)):
        assert_almost_equal(sv00[i].re, sv00_t[i].re, 'Measure c=00', atol=1e-7)
        assert_almost_equal(sv00[i].im, sv00_t[i].im, 'Measure c=00', atol=1e-7)
    for i in range(len(sv01)):
        assert_almost_equal(sv01[i].re, sv01_t[i].re, 'Measure c=00', atol=1e-7)
        assert_almost_equal(sv01[i].im, sv01_t[i].im, 'Measure c=00', atol=1e-7)
    for i in range(len(sv10)):
        assert_almost_equal(sv10[i].re, sv10_t[i].re, 'Measure c=00', atol=1e-7)
        assert_almost_equal(sv10[i].im, sv10_t[i].im, 'Measure c=00', atol=1e-7)
    for i in range(len(sv11)):
        assert_almost_equal(sv11[i].re, sv11_t[i].re, 'Measure c=00', atol=1e-7)
        assert_almost_equal(sv11[i].im, sv11_t[i].im, 'Measure c=00', atol=1e-7)