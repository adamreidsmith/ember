'''
Quantum Teleportation

Suppose there are two researchers, Alice and Bob, and Alice wants to send a quantum state
to Bob. Quantum telepportation allows Alice to do just this by only sending classical
information, at the expense of a shared entangled state.  The protocol is as follows:

    1. Alice initializes a qubit to the state |ψ⟩ she wishes to teleport to Bob.

    2. An entangled state (Bell state) is generated, with one qubit sent to Alice and
        the other sent to Bob.
    
    3. A Bell measurment of Alice's Bell state qubit and the qubit to be teleported is
        performed by Alice. This yields one of 4 measurement outcomes which can be sent
        to Bob using two classical bits.
    
    4. Alice sends the two classical bits to Bob. As a result of Alice's measurement, Bob's
        qubit is in one of four possible states. Based on the classical bits received from
        Alice, Bob performs an operation on his qubit, which places it in a state identical
        to the state |ψ⟩ chosen for teleportation.

Below we implement the quantum teleportation protocol to teleport the state
(1 / sqrt(2) + 0.3i)|0⟩ + (0.4 - 0.5i)|1⟩
'''

from math import sqrt

from ember import QuantumCircuit, StatevectorSimulator, Statevector, DensityMatrix, Gate
from ember import ComplexScalar, CMatrix
from ember import X, Z, H, CX, Measure


fn state_preparation(state: Tuple[ComplexScalar, ComplexScalar]) raises -> QuantumCircuit:
    '''Prepare the single qubit state we wish to teleport.'''
    var statevector = Statevector(state[0], state[1])
    var qc = QuantumCircuit(1)
    qc.set_initial_state(statevector)
    return qc^

fn entangle_generation() raises -> QuantumCircuit:
    '''Generate an entangled Bell state.'''
    var qc = QuantumCircuit(2)
    qc.apply(H(0), CX(0, 1))
    return qc^

fn change_of_basis() raises -> QuantumCircuit:
    '''Change the basis so we can perform a Bell measurement.'''
    var qc = QuantumCircuit(2)
    qc.apply(CX(0, 1), H(0))
    return qc^

fn measure_and_update() raises -> QuantumCircuit:
    '''Measure Alice's qubits and perform operations on Bob's qubit based on the measurment.'''
    var qc = QuantumCircuit(3, 2)
    qc.apply(Measure(0, 0), Measure(1, 1))
    var x: Gate = X(2)
    x.control(clbits=List[Int, True](1))
    var z: Gate = Z(2)
    z.control(clbits=List[Int, True](0))
    qc.apply(x, z)
    return qc^

fn teleport(state: Tuple[ComplexScalar, ComplexScalar]) raises -> QuantumCircuit:
    '''Put it all together to create the quantum teleportation circuit.'''
    var teleportation_circuit = QuantumCircuit(3, 2)
    teleportation_circuit.join(state_preparation(state), 0)
    teleportation_circuit.join(entangle_generation(), 1, 2)
    teleportation_circuit.join(change_of_basis(), 0, 1)
    teleportation_circuit.join(
        measure_and_update(),
        List[Int, True](0, 1, 2),
        List[Int, True](0, 1),
    )
    return teleportation_circuit^

fn quantum_teleportation() raises:
    # The state we will teleport is (1 / sqrt(2) + 0.3i)|0⟩ + (0.4 - 0.5i)|1⟩
    var state = (ComplexScalar(1.0 / sqrt(2.0), 0.3), ComplexScalar(0.4, -0.5))

    print('State to teleport:', Statevector(state[0], state[1]), sep='\n')

    # Call our teleportation function to generate the quantum teleporation circuit
    var teleportation_circuit: QuantumCircuit = teleport(state)

    # Initialize a statevector simultor and run the teleportation circuit
    var sim = StatevectorSimulator()
    sim.run(teleportation_circuit)

    # Get the resulting statevector and trace out Alice's qubits to isolate Bob's qubits
    var statevector: Statevector = sim.get_statevector()
    statevector = statevector.partial_trace(0, 1)

    print('Bob\'s teleported state:', statevector, sep='\n')

    # The result may not look exactly like the original state, so what's going on?
    # An important note is that statevectors are only defined up to a global phase,
    # so the statevector for Bob's telported qubit may differ from the original state
    # by a phase factor. If this is the case, each amplitude should differ from the
    # original by the same factor, and that factor should have norm 1.
    var amp_0_factor: ComplexScalar = state[0] / statevector[0]
    var amp_1_factor: ComplexScalar = state[1] / statevector[1]
    
    # Check if the phase factors agree and have norm 1
    if not amp_0_factor.is_close(amp_1_factor):
        raise Error('Statevector amplitude factors differ!')
    if not abs(amp_0_factor).is_close(1):
        raise Error('Statevector is not normalized!')
    
    # Adjust Bob's teleported state by the global phase to confirm the teleportation was successful
    var phase_corrected_statevector = statevector.apply_global_phase_factor(amp_0_factor)
    print('Bob\'s global phase-corrected teleported state:', phase_corrected_statevector, sep='\n')
