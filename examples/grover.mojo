'''
Grover's Search Algorithm

Grover's algorithm is a quantum algorithm for unstructured search, providing a
quadratic speedup over the best possible classical algorithm. Given a "black box"
function, or an "oracle," that identifies a specific item in an unsorted list,
Grover's algorithm can find that item much faster than a classical computer could.

The protocol is as follows:

    1. All qubits are initialized into a uniform superposition using Hadamard gates.
       This represents a state where every possible item in the search space is
       equally likely.

    2. The "Oracle" is applied. This quantum circuit marks the desired item by
       inverting its phase. For example, if we are searching for the state |101⟩,
       the oracle will transform the state into -|101⟩, leaving all other
       states unchanged.

    3. The "Diffuser" (or amplitude amplification) operator is applied. This
       operator reflects the state of the system about the average amplitude, which
       has the effect of significantly increasing the amplitude of the marked item
       while decreasing the amplitudes of all other items.
    
    4. Steps 2 and 3 are repeated a specific number of times to maximize the
       amplitude (and thus the measurement probability) of the marked item. The
       optimal number of iterations is approximately (π/4) * sqrt(N), where N
       is the total number of items in the search space.

    5. The qubits are measured. Due to the amplification, the measurement outcome
       will be the marked state with a very high probability.

Below we implement Grover's algorithm to find the state |101⟩ in a 3-qubit system.
'''

from math import pi, sqrt, floor

from ember import QuantumCircuit, StatevectorSimulator, Statevector, H, X, Z, Gate


fn build_oracle(n_qubits: Int, marked_state: Int) raises -> QuantumCircuit:
    '''Builds the oracle circuit that marks a specific computational basis state.
    
    The oracle flips the phase of the `marked_state` and leaves all other states
    unchanged. This is achieved by applying X-gates to the qubits that are |0⟩ in the
    marked state, performing a multi-controlled Z gate, and then undoing the X-gates.

    Args:
        n_qubits: The total number of qubits in the circuit.
        marked_state: The integer representation of the state to be marked.
    
    Returns:
        The quantum circuit for the oracle.
    '''
    var qc = QuantumCircuit(n_qubits)
    
    # Flip the qubits that are 0 in the binary representation of marked_state
    for i in range(n_qubits):
        if (marked_state >> i) & 1 == 0:
            qc.apply(X(i))

    # Apply a multi-controlled Z gate. This is a Z gate on the target qubit (n_qubits - 1)
    # controlled by all other qubits (0 to n_qubits - 2).
    var mcz: Gate = Z(n_qubits - 1)
    var controls = List[Int, True]()
    for i in range(n_qubits - 1):
        controls.append(i)
    mcz.control(controls)
    qc.apply(mcz)

    # Un-flip the qubits that were flipped before
    for i in range(n_qubits):
        if (marked_state >> i) & 1 == 0:
            qc.apply(X(i))
            
    return qc


fn build_diffuser(n_qubits: Int) raises -> QuantumCircuit:
    '''Builds the Grover diffuser (amplitude amplification) circuit.

    The diffuser circuit reflects the state vector about the mean amplitude.
    This operation amplifies the amplitude of the marked state.

    Args:
        n_qubits: The total number of qubits.

    Returns:
        The quantum circuit for the diffuser.
    '''
    var qc = QuantumCircuit(n_qubits)

    # Apply Hadamard gates to all qubits
    for i in range(n_qubits):
        qc.apply(H(i))

    # Apply X gates to all qubits
    for i in range(n_qubits):
        qc.apply(X(i))

    # Apply a multi-controlled Z gate
    var mcz: Gate = Z(n_qubits - 1)
    var controls = List[Int, True]()
    for i in range(n_qubits - 1):
        controls.append(i)
    mcz.control(controls)
    qc.apply(mcz)

    # Apply X gates to all qubits
    for i in range(n_qubits):
        qc.apply(X(i))

    # Apply Hadamard gates to all qubits
    for i in range(n_qubits):
        qc.apply(H(i))
        
    return qc


fn grovers_algorithm() raises:
    var n_qubits: Int = 3
    # The state we will search for is |101⟩, which is 5 in integer representation
    var marked_state: Int = 5

    print('Grover\'s Algorithm Example')
    print('==========================')
    print('Number of qubits:', n_qubits)
    print('Marked state:', marked_state, '(|101⟩)')

    # The total number of states in the search space is 2^n_qubits
    var n_states: Int = 2 ** n_qubits

    # The optimal number of iterations to run the algorithm.
    var optimal_iterations = Int(floor(pi / 4.0 * sqrt(Float64(n_states))))
    print('Optimal number of iterations:', optimal_iterations)
    
    # --- Step 1: Create the uniform superposition ---
    var grover_circuit = QuantumCircuit(n_qubits)
    for i in range(n_qubits):
        grover_circuit.apply(H(i))

    # --- Step 2 & 3: Repeatedly apply the Oracle and Diffuser ---
    var oracle: QuantumCircuit = build_oracle(n_qubits, marked_state)
    var diffuser: QuantumCircuit = build_diffuser(n_qubits)
    for _ in range(optimal_iterations):
        grover_circuit.join(oracle)
        grover_circuit.join(diffuser)

    # --- Step 4: Simulate the circuit ---
    var sim = StatevectorSimulator()
    sim.run(grover_circuit)

    # --- Step 5: Get results and verify ---
    var final_statevector: Statevector = sim.get_statevector()
    print('\nFinal Statevector:', final_statevector, sep='\n')
    
    var probabilities: List[Float64, True] = final_statevector.get_probabilities()
    print('\nFinal Probabilities:')
    for i in range(n_states):
        print('    ' + bin(i)[2:].rjust(n_qubits, '0') + ': ' + String(probabilities[i]))
    
    # Find the outcome with the highest probability
    var max_prob: Float64 = -1.0
    var measured_state: Int = -1
    for i in range(len(probabilities)):
        if probabilities[i] > max_prob:
            max_prob = probabilities[i]
            measured_state = i

    print('\n--- Verification ---')
    print(
        'Most likely measured state:',
        measured_state,
        '(|' + bin(measured_state)[2:].rjust(n_qubits, '0') + '⟩)',
    )
    print('Probability: ' + String(round(max_prob * 100, 2)) + '%')
    
    if measured_state == marked_state:
        print('\nSuccess! The algorithm correctly identified the marked state.')
    else:
        raise Error('Failure! The algorithm did not find the marked state.')
