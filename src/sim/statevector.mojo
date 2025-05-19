'''
A quantum circuit statevector simulator.

Algorithms taken/adapted from [1].

[1] Jones, Tyson, Bálint Koczor, and Simon C. Benjamin. “Distributed Simulation of
    Statevectors and Density Matrices.” arXiv.org, November 2, 2023. 
    https://arxiv.org/abs/2311.01512. 
'''

import random
from algorithm import parallelize, vectorize
from math import sqrt

from ..quantum import QuantumCircuit, Gate
from ..cplx import CMatrix, CSRCMatrix, CSRBuilder
from ..config import DEFAULT_TOL


## Bit-twiddling functions of unsigned integers ###################################################

@always_inline
fn get_bit(n: UInt, t: UInt) -> UInt:
    '''Retrieves the value of the bit at a specific position t within an unsigned integer n.
    
    Args:
        n: The integer to get a bit from.
        t: The bit position to get.
    
    Returns:
        The specified bit.
    '''
    return (n >> t) & 1

@always_inline
fn flip_bit(n: UInt, t: UInt) -> UInt:
    '''Inverts the bit at position t in the unsigned integer n.
    
    Args:
        n: The integer to flip a bit in.
        t: The bit position to flip.
    
    Returns:
        The integer with the specified bit flipped.
    '''
    return n ^ (1 << t)

@always_inline
fn flip_bits(owned n: UInt, t: List[UInt, True]) -> UInt:
    '''Flips all the bits in the unsigned integer n at the positions specified in t.
    
    Args:
        n: The integer to flip bits in.
        t: The bit positions to flip.
    
    Returns:
        The integer with the specified bits flipped.
    '''
    for q in t:
        n = flip_bit(n, q[])
    return n

@always_inline
fn insert_bit(n: UInt, t: UInt, b: UInt) -> UInt:
    '''Sets the bit at position t of integer n to the value b.
    
    Args:
        n: The integer to set a bit in.
        t: The bit position to set.
        b: The value to set the bit to.
    
    Returns:
        The integer with the specified bit set to the specified value.
    '''
    # This function actually *sets* the bit at position t of n to b rather than inserting it
    var l: UInt = (n >> t) << (t + 1)
    var m: UInt = b << t
    var r: UInt = n & ((1 << t) - 1)
    return l | m | r

@always_inline
fn insert_bits(owned n: UInt, t: List[UInt, True], b: UInt) -> UInt:
    '''Sets the bits of n at all positions specified in t to the value b.
    
    Args:
        n: The integer to set bits in.
        t: The bit positions to set.
        b: The value to set the bits to.
    
    Returns:
        The integer with the specified bits set to the specified value.
    '''
    # This function actually *sets* the bits at positions t of n to b rather than inserting them
    for q in t:
        n = insert_bit(n, q[], b)
    return n

@always_inline
fn set_bits(owned n: UInt, t: List[UInt, True], v: UInt) -> UInt:
    '''Sets specific bits within the integer n based on the values from another integer v.
    
    Args:
        n: The integer to set bits in.
        t: The bit positions to set.
        v: The integer to get bit values from.
    
    Returns:
        The integer with the specified bits set to the bit values in v.
    '''
    var b: UInt
    for q in range(len(t)):
        b = get_bit(v, q) << t[q]
        n = (n & (~b)) | b
    return n

###################################################################################################

@value
struct StatevectorSimulator[type: DType, tol: Scalar[type] = DEFAULT_TOL]:
    '''A quantum circuit statevector simulator.
    
    Parameters:
        type: A type for the circuit and statevector data.
        tol: A tolerance for unitarity and closeness checks.
    '''

    var _sv: List[ComplexScalar[Self.type], True]
    '''A list to store the statevector elements.'''
    var _cb: List[Int, True]
    '''The values of the classical bits in the circuit.'''
    var _qc: QuantumCircuit[Self.type, Self.tol]
    '''The quantum circuit that is being simulated.'''

    fn __init__(out self):
        '''Initialize a statevector simulator.'''
        self.set_seed()
        # Placeholders that are overwritten when `run` is called
        self._sv = List[ComplexScalar[Self.type], True]()
        self._cb = List[Int, True]()
        # Initialize directly to avoid the raising __init__ method
        self._qc = QuantumCircuit[Self.type, Self.tol](
            n_qubits=0,
            n_clbits=0,
            _cbits=List[Int, True](),
            _data=List[Gate[Self.type, Self.tol]](),
        )
    
    fn __init__(out self, seed: Int):
        '''Initialize a statevector simulator.
        
        Args:
            seed: An integer to seed the random number generator.
        '''
        # Placeholders that are overwritten when `run` is called
        self._sv = List[ComplexScalar[Self.type], True]()
        self._cb = List[Int, True]()
        # Initialize directly to avoid the raising __init__ method
        self._qc = QuantumCircuit[Self.type, Self.tol](
            n_qubits=0,
            n_clbits=0,
            _cbits=List[Int, True](),
            _data=List[Gate[Self.type, Self.tol]](),
        )
        self.set_seed(seed)
    
    @staticmethod
    @always_inline
    fn set_seed():
        '''Set the seed of the random number generator with the current time.'''
        random.seed()
    
    @staticmethod
    @always_inline
    fn set_seed(seed: Int):
        '''Set the seed of the random number generator.

        Args:
            seed: An integer seed.
        '''
        random.seed(seed)

    fn run(mut self, qc: QuantumCircuit[Self.type, Self.tol], *, parallel: Bool = True) -> Self:
        '''Run the statevector simulator on the quantum circuit.
        
        Args:
            qc: The quantum circuit to simulate.
            parallel: Whether or not to parallelize statevector simulation.
        
        Returns:
            Self.
        '''
        self._qc = qc
        self._cb = self._qc._cbits

        # Initialize the statevector to the |0> state
        self._sv = List[ComplexScalar[Self.type], True](length=2 ** self._qc.n_qubits, fill=0)
        self._sv[0] = 1

        for gate_ref in qc._data:
            var gate: Gate[Self.type, Self.tol] = gate_ref[]

            if gate._is_measure:
                self._measure(gate)
                continue

            # Check classical controls
            for c in gate.classical_controls:
                if not self._cb[c[]]:
                    break
            else:
                if gate.n_qubits == 1:
                    # Single-qubit uncontrolled gate
                    self._apply_one_qubit_gate(gate, parallel)
                elif len(gate.qubits) == 1:
                    # Single-qubit controlled gate
                    self._apply_multi_control_one_qubit_gate(gate, parallel)
                elif len(gate.controls) == 0:
                    # Multi-qubit uncontrolled gate
                    self._apply_multi_qubit_gate(gate, parallel)
                else:
                    # Multi-qubit controlled gate
                    self._apply_multi_control_multi_qubit_gate(gate, parallel)
        
        return self

    fn _apply_one_qubit_gate(mut self, gate: Gate[Self.type, Self.tol], parallel: Bool):
        '''Apply a single qubit gate to the statevector.

        Implements Alg. 2 of [1].

        Args:
            gate: The gate to apply. Must be a single-qubit uncontrolled gate.
            parallel: Whether or not to parallelize statevector simulation.
        '''
        var t: UInt = gate.qubits[0]

        if parallel:
            @parameter
            fn op_sv_elem(n: Int):
                var ib: UInt = insert_bit(UInt(n), t, 0)
                var ig: UInt = flip_bit(ib, t)
                var b: ComplexScalar[Self.type] = self._sv[ib]
                var g: ComplexScalar[Self.type] = self._sv[ig]
                self._sv[ib] = gate.matrix.load_idx[1](0) * b + gate.matrix.load_idx[1](1) * g
                self._sv[ig] = gate.matrix.load_idx[1](2) * b + gate.matrix.load_idx[1](3) * g
            parallelize[op_sv_elem](2 ** (self._qc.n_qubits - 1))
        else:
            var ib: UInt
            var ig: UInt
            var b: ComplexScalar[Self.type]
            var g: ComplexScalar[Self.type]
            for n in range(2 ** (self._qc.n_qubits - 1)):
                ib = insert_bit(UInt(n), t, 0)
                ig = flip_bit(ib, t)
                b = self._sv[ib]
                g = self._sv[ig]
                self._sv[ib] = gate.matrix.load_idx[1](0) * b + gate.matrix.load_idx[1](1) * g
                self._sv[ig] = gate.matrix.load_idx[1](2) * b + gate.matrix.load_idx[1](3) * g
    
    fn _apply_multi_control_one_qubit_gate(
        mut self, 
        gate: Gate[Self.type, Self.tol], 
        parallel: Bool
    ):
        '''Apply a single-qubit gate controlled on one or more qubits to the statevector.

        Implements Alg. 3 of [1].

        Args:
            gate: The gate to apply. Must be a single-qubit controlled gate.
            parallel: Whether or not to parallelize statevector simulation.
        '''
        var s: UInt = len(gate.controls)

        # Get a list of the control qubit indices
        var q = List[UInt, True](capacity=s + 1)
        for qubit in gate.controls:
            q.append(qubit[])
        
        # Target qubit index
        var t: UInt = gate.qubits[0]
        
        # Union the controlled and target qubits and sort the list
        q.append(t)
        sort[T=UInt](q)

        if parallel:
            @parameter
            fn op_sv_elem(j: Int):
                var ig: UInt = insert_bits(UInt(j), q, 1)
                var ib: UInt = flip_bit(ig, t)
                var b: ComplexScalar[Self.type] = self._sv[ib]
                var g: ComplexScalar[Self.type] = self._sv[ig]
                self._sv[ib] = gate.matrix.load_idx[1](0) * b + gate.matrix.load_idx[1](1) * g
                self._sv[ig] = gate.matrix.load_idx[1](2) * b + gate.matrix.load_idx[1](3) * g
            parallelize[op_sv_elem](2 ** (self._qc.n_qubits - s))
        else:
            var ib: UInt
            var ig: UInt
            var b: ComplexScalar[Self.type]
            var g: ComplexScalar[Self.type]
            for j in range(2 ** (self._qc.n_qubits - s)):
                ig = insert_bits(j, q, 1)
                ib = flip_bit(ig, t)
                b = self._sv[ib]
                g = self._sv[ig]
                self._sv[ib] = gate.matrix.load_idx[1](0) * b + gate.matrix.load_idx[1](1) * g
                self._sv[ig] = gate.matrix.load_idx[1](2) * b + gate.matrix.load_idx[1](3) * g

    fn _apply_multi_qubit_gate(mut self, gate: Gate[Self.type, Self.tol], parallel: Bool):
        '''Apply a multi-qubit gate with no controls to the statevector.

        Implements Alg. 4 of [1].

        Args:
            gate: The gate to apply. Must be a multi-qubit uncontrolled gate.
            parallel: Whether or not to parallelize statevector simulation.
        '''
        var n: UInt = len(gate.qubits)

        # Form a list of qubit indices and a sorted version
        var t = List[UInt, True](capacity=n)
        for qubit in gate.qubits:
            t.append(qubit[])
        var q: List[UInt, True] = t.copy()
        sort[T=UInt](q)

        if parallel:
            @parameter
            fn op_sv_elem(k: Int):
                var i: UInt
                var v = List[ComplexScalar[Self.type], True](length=2 ** n, fill=0)
                var z: UInt = insert_bits(UInt(k), q, 0)

                for j in range(2 ** n):
                    i = set_bits(z, t, UInt(j))
                    v[j] = self._sv[i]
                
                for j in range(2 ** n):
                    i = set_bits(z, t, UInt(j))
                    self._sv[i] = 0

                    for l in range(2 ** n):
                        self._sv[i] += gate.matrix.load_crd[1](j, l) * v[l]
            parallelize[op_sv_elem](2 ** (self._qc.n_qubits - n))
        else:
            var v = List[ComplexScalar[Self.type], True](length=2 ** n, fill=0)
            var z: UInt
            var i: UInt
            for k in range(2 ** (self._qc.n_qubits - n)):
                z = insert_bits(UInt(k), q, 0)

                for j in range(2 ** n):
                    i = set_bits(z, t, UInt(j))
                    v[j] = self._sv[i]
                
                for j in range(2 ** n):
                    i = set_bits(z, t, UInt(j))
                    self._sv[i] = 0

                    for l in range(2 ** n):
                        self._sv[i] += gate.matrix.load_crd[1](j, l) * v[l]
    
    fn _apply_multi_control_multi_qubit_gate(
        mut self, 
        gate: Gate[Self.type, Self.tol], 
        parallel: Bool
    ):
        '''Apply a multi-controlled multi-qubit gate to the statevector.

        Adaptation of Alg. 4 of [1].

        Args:
            gate: The gate to apply. Must be a multi-qubit controlled gate.
            parallel: Whether or not to parallelize statevector simulation.
        '''
        alias zero = ComplexScalar[Self.type](0)
        alias one = ComplexScalar[Self.type](1)

        var n: UInt = len(gate.qubits) + len(gate.controls)
        var N: Int = 2 ** n

        # Form a list of qubit indices and a sorted version
        var t = List[UInt, True](capacity=n)
        for qubit in gate.qubits:
            t.append(qubit[])
        for qubit in gate.controls:
            t.append(qubit[])
        var q: List[UInt, True] = t.copy()
        sort[T=UInt](q)
        
        var rc_start: Int = N - 2 ** len(gate.qubits)
        @parameter
        fn get_full_unitary_elem(r: Int, c: Int) -> ComplexScalar[Self.type]:
            if r >= rc_start and c >= rc_start:
                return gate.matrix.load_crd[1](r - rc_start, c - rc_start)
            return one if r == c else zero
        
        if parallel:
            @parameter
            fn op_sv_elem(k: Int):
                var i: UInt
                var v = List[ComplexScalar[Self.type], True](length=N, fill=0)
                var z: UInt = insert_bits(UInt(k), q, 0)

                for j in range(N):
                    i = set_bits(z, t, UInt(j))
                    v[j] = self._sv[i]
                
                for j in range(rc_start, N):
                    i = set_bits(z, t, UInt(j))
                    self._sv[i] = 0
                    for l in range(rc_start, N):
                        self._sv[i] += gate.matrix.load_crd[1](j - rc_start, l - rc_start) * v[l]
            parallelize[op_sv_elem](2 ** (self._qc.n_qubits - n))
        else:
            var v = List[ComplexScalar[Self.type], True](length=N, fill=0)
            var z: UInt
            var i: UInt
            for k in range(2 ** (self._qc.n_qubits - n)):
                z = insert_bits(UInt(k), q, 0)

                for j in range(N):
                    i = set_bits(z, t, UInt(j))
                    v[j] = self._sv[i]
                
                for j in range(rc_start, N):
                    i = set_bits(z, t, UInt(j))
                    self._sv[i] = 0
                    for l in range(rc_start, N):
                        self._sv[i] += gate.matrix.load_crd[1](j - rc_start, l - rc_start) * v[l]
    
    fn _measure(mut self, gate: Gate[Self.type, Self.tol]):
        '''Measure a subset of qubits and update the statevector.
        
        Args:
            gate: The measurement gate to simulate.
        '''
        var n: Int = self._qc.n_qubits
        var N: Int = 2 ** n
        var m: Int = len(gate.qubits)
        var M: Int = 2 ** m
        var probabilities = List[Scalar[Self.type], True](length=M, fill=0)

        # Compute probabilities for each possible outcome
        var outcome: Int
        var bit: Int
        for i in range(N):
            # Extract the bits at the measured positions
            outcome = 0
            for j in range(m):
                bit = (i >> gate.qubits[j]) & 1
                outcome |= (bit << j)
            probabilities[outcome] += self._sv[i].squared_norm()
        
        # Normalize probabilities to ensure they sum to one
        var total_prob: Scalar[Self.type] = 0
        for p in probabilities:
            total_prob += p[]
        if abs(total_prob - 1) >= Self.tol:
            for i in range(M):
                probabilities[i] /= total_prob
        
        # Select outcome based on probabilities
        var r: Scalar[Self.type] = random.random_float64().cast[Self.type]()
        var cumulative_prob: Scalar[Self.type] = 0
        var selected_outcome: Int = 0
        for i in range(M):
            cumulative_prob += probabilities[i]
            if r <= cumulative_prob:
                selected_outcome = i
                break
        
        norm_factor = sqrt(probabilities[selected_outcome])
        new_statevector = List[ComplexScalar[Self.type], True](length=N, fill=0)

        # For each basis state, check if it's consistent with the selected measurement outcome
        var consistent: Bool
        var expected_bit: Int
        for i in range(N):
            consistent = True
            for j in range(m):
                bit = (i >> gate.qubits[j]) & 1
                expected_bit = (selected_outcome >> j) & 1
                if bit != expected_bit:
                    consistent = False
                    break
            
            # Update statevector
            if consistent:
                new_statevector[i] = self._sv[i] / norm_factor

        # The position i bit (where position 0 is the rightmost [LSB] bit) in 
        # `selected_outcome` corresponds to the ith qubit in gate.applied_to
        for i in range(m):
            bit = (selected_outcome >> i) & 1
            # self._cb is a list of classical bits in the circuit
            self._cb[gate._measure_targs[i]] = bit

        self._sv = new_statevector^

    fn get_statevector(self) -> CMatrix[Self.type]:
        '''Get the statevector as a CMatrix.
        
        Returns:
            The statevector.
        '''
        # Normalize probabilities to ensure they sum to one
        var statevector = CMatrix[Self.type](rows=len(self._sv), cols=1, fill_zeros=False)
        var total_prob: Scalar[type] = 0
        for i in range(len(self._sv)):
            total_prob += self._sv[i].squared_norm()
            statevector.store_idx[1](i, self._sv[i])

        # Ensure the statevector is normalized
        if abs(total_prob - 1) >= Self.tol:
            statevector /= sqrt(total_prob)

        return statevector^
    
    # fn set_clbits(mut self, bit_values: Dict[Int, Int]) raises:
    #     '''Set the classical bits to the values specified in bit_values.

    #     Args:
    #         bit_values: A mapping of classical bit specifiers to bit values.
    #     '''
    #     # Check that the bits are in the circuit and that all values are binary
    #     for kv in bit_values.items():
    #         var bit: Int = kv[].key
    #         var value: Int = kv[].value
    #         if bit < 0 or bit >= self._qc.n_clbits:
    #             raise Error('Cannot set bit ' + String(bit) + ' as it is not in the circuit.')
    #         if value not in Tuple[Int, Int](0, 1):
    #             raise Error(
    #                 'Bits can only be set to 0 or 1. Received ' + String(value) 
    #                 + ' for bit ' + String(bit) + '.'
    #             )
        
    #     for kv in bit_values.items():
    #         self._cb[kv[].key] = kv[].value
        