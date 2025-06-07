'''
Quantum Fourier Transform

The quantum Fourier transform (QFT) is the quantum analogue of the discrete Fourier transform (DFT),
and is a crucial component in many quantum algorithms, including Shor's factoring algorithm
and quantum phase estimation.

Given a quantum state |x⟩ represented on n qubits, the QFT maps this state to

    (1 / sqrt(2^n)) ∑_y exp(2πi * x * y / 2^n) |y⟩

which is the representation of the state |x⟩ in the Fourier basis.

The QFT can be implemented efficiently using a series of Hadamard gates and controlled phase
shifts. The QFT circuit can be approximated by ignoring small phase rotations and omitting
final swaps, which can reduce the gate depth for practical purposes.

Below we implement the QFT and apply it to the computational basis state |1⟩ on a 3-qubit register.
'''

from math import pi

from ember import QuantumCircuit, H, PHASE, SWAP, X, Gate
from ember import Statevector, StatevectorSimulator


fn build_qft_circuit(
    n_qubits: Int, approximation_degree: Int = 0, do_swaps: Bool = True
) raises -> QuantumCircuit:
    '''Build the quantum Fourier transform circuit.

    Args:
        n_qubits: The number of qubits in the QFT circuit.
        approximation_degree: The degree of approximation. A value of 0 means no approximation.
        do_swaps: Whether to include the final swaps in the QFT.

    Returns:
        The quantum circuit representing the QFT.
    '''
    # Initialize a quantum circuit with the specified number of qubits.
    var qc = QuantumCircuit(n_qubits)

    # The QFT algorithm consists of a series of Hadamard and controlled-phase gates.
    # We iterate through each qubit, starting from the most significant one.
    for j in range(n_qubits - 1, -1, -1):
        # Apply a Hadamard gate to the current qubit. This puts the qubit into a
        # superposition, which is the first step of creating the Fourier basis states.
        qc.apply(H(j))

        # Apply a series of controlled-phase rotations. The number of these rotations
        # can be reduced by the approximation_degree to create an approximate QFT,
        # which can be useful for reducing gate count and circuit depth.
        var num_entanglements: Int = max(0, j - max(0, approximation_degree - (n_qubits - j - 1)))
        for k in range(j - 1, j - num_entanglements - 1, -1):
            # The angle of rotation depends on the distance between the control and target qubits.
            var theta: Float64 = pi * 2.0 ** (k - j)
            var cp: Gate = PHASE(k, theta)
            cp.control(j)
            qc.apply(cp)

    # After the rotations, the order of the qubits is reversed relative to the
    # desired output. If do_swaps is true, we apply a series of SWAP gates
    # to reverse the qubits to the correct order.
    if do_swaps:
        for i in range(n_qubits // 2):
            qc.apply(SWAP(i, n_qubits - i - 1))
    
    return qc

fn quantum_fourier_transform() raises:
    '''This example demonstrates the use of the Quantum Fourier Transform (QFT) on a
    simple quantum state. The steps are as follows:
        1. An initial state is prepared. For this example, we'll use a 4-qubit state
        representing the number 5 in the computational basis, i.e., |0101⟩.

        2. The QFT circuit is constructed and applied to this initial state.

        3. To verify the transformation, we also construct the inverse QFT (QFT†).
        Applying the inverse QFT to the transformed state should return it to the
        initial state.

        4. A simulator is used to run the circuits and observe the resulting statevectors
        at each step.
    '''
    # Initialize a quantum circuit to implement state preparation and the quantum 
    # Fourier transform.
    var n_qubits: Int = 4
    var full_circuit = QuantumCircuit(n_qubits)

    # Let's prepare the state |5⟩ = |0101⟩ using Pauli X gates.
    # Note: Ember uses a little-endian convention, so |0101⟩ is prepared by applying
    # Pauli X gates to the qubits at index 0 and 2 to flip each of them from the |0⟩
    # to the |1⟩ state.
    full_circuit.apply(X(0), X(2))

    # Build the QFT circuit for 4 qubits.
    var qft_circuit: QuantumCircuit = build_qft_circuit(n_qubits)

    # Combine the preparation and QFT circuits.
    full_circuit.join(qft_circuit, 0, 1, 2, 3)

    # Initialize a statevector simulator and run the QFT circuit.
    var sim = StatevectorSimulator()
    sim.run(full_circuit)
    var qft_state: Statevector = sim.get_statevector()
    
    print('State after Quantum Fourier Transform:', qft_state, sep='\n')

    # To verify the QFT, we can apply its inverse. The inverse of the QFT circuit
    # should transform the state back to the original |0101⟩.
    var iqft_circuit: QuantumCircuit = qft_circuit.inverse()
    
    # Apply the inverse QFT to the transformed state.
    iqft_circuit.set_initial_state(qft_state)
    sim.run(iqft_circuit)
    var final_state: Statevector = sim.get_statevector()

    print('State after applying the inverse QFT:', final_state, sep='\n')

    # If the QFT and its inverse were successful, the element at index 5 of the
    # final statevector should have norm one (representing the |0101⟩ state)
    if not final_state.is_normalized():
        raise Error('Verification failed: The final state is not normalized.')
    if not abs(final_state[5]).is_close(1):
        raise Error('Verification failed: The final state does not match the initial state.')


    print('Successfully verified the QFT and its inverse.')