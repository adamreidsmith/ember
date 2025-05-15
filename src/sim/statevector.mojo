'''
[1] https://arxiv.org/abs/2311.01512
'''

import random
from collections import Dict, Set
from algorithm import parallelize, vectorize
from memory import memcpy
from sys import simdwidthof
from math import sqrt

from ..quantum import QuantumCircuit, Gate
from ..cplx import CMatrix, CSRCMatrix, CSRBuilder
from ..config import DEFAULT_TOL

## Bit-twiddling functions of unsigned integers ###################################################

@always_inline
fn get_bit(n: UInt, t: UInt) -> UInt:
    '''Retrieves the value of the bit at a specific position t within an unsigned integer n.'''
    return (n >> t) & 1

@always_inline
fn flip_bit(n: UInt, t: UInt) -> UInt:
    '''Inverts the bit at position t in the unsigned integer n.'''
    return n ^ (1 << t)

fn flip_bits(owned n: UInt, t: List[UInt, True]) -> UInt:
    '''Flips all the bits in the unsigned integer n at the positions specified in t.'''
    for q in t:
        n = flip_bit(n, q[])
    return n

fn insert_bit(n: UInt, t: UInt, b: UInt) -> UInt:
    '''Sets the bit at position t of integer n to the value b.'''
    # This function actually *sets* the bit at position t of n to b rather than inserting it
    var l: UInt = (n >> t) << (t + 1)
    var m: UInt = b << t
    var r: UInt = n & ((1 << t) - 1)
    return l | m | r

fn insert_bits(owned n: UInt, t: List[UInt, True], b: UInt) -> UInt:
    '''Sets the bits of n at all positions specified in t to the value b.'''
    # This function actually *sets* the bits at positions t of n to b rather than inserting them
    for q in t:
        n = insert_bit(n, q[], b)
    return n

fn set_bits(owned n: UInt, t: List[UInt, True], v: UInt) -> UInt:
    '''Sets specific bits within the integer n based on the values from another integer v.'''
    var b: UInt
    for q in range(len(t)):
        b = get_bit(v, q) << t[q]
        n = (n & (~b)) | b
    return n

###################################################################################################

@value
struct StatevectorSimulator[type: DType, tol: Scalar[type] = DEFAULT_TOL]:
    var sv: List[ComplexScalar[Self.type], True]
    var _cb: List[Bool, True]  # Classical bit values
    var _qc: QuantumCircuit[Self.type, Self.tol]

    fn __init__(out self):
        # Placeholders that are overwritten when `run` is called
        random.seed()  # Seed with the current time
        self.sv = List[ComplexScalar[Self.type], True]()
        self._cb = List[Bool, True]()
        # Initialize directly to avoid raising __init__ method
        self._qc = QuantumCircuit[Self.type, Self.tol](
            n_qubits=0,
            n_clbits=0,
            _data=List[Gate[Self.type, Self.tol]]()
        )
    
    fn __init__(out self, seed: Int):
        # Placeholders that are overwritten when `run` is called
        random.seed(seed)
        self.sv = List[ComplexScalar[Self.type], True]()
        self._cb = List[Bool, True]()
        # Initialize directly to avoid raising __init__ method
        self._qc = QuantumCircuit[Self.type, Self.tol](
            n_qubits=0,
            n_clbits=0,
            _data=List[Gate[Self.type, Self.tol]]()
        )

    fn run(mut self, qc: QuantumCircuit[Self.type, Self.tol], *, parallel: Bool = True) raises:
        '''Run the statevector simulator on the quantum circuit.'''
        self._qc = qc
        self._cb = List[Bool, True](length=self._qc.n_clbits, fill=0)

        # Initialize the statevector to the |0> state
        self.sv = List[ComplexScalar[Self.type], True](length=2 ** self._qc.n_qubits, fill=0)
        self.sv[0] = 1

        for gate_ref in qc._data:
            var gate: Gate[Self.type, Self.tol] = gate_ref[]
            if gate._is_measure:
                self._measure(gate)
            elif gate.n_qubits == 1:
                # Single-qubit uncontrolled gate
                self._apply_one_qubit_gate(gate, parallel)
            elif len(gate.applied_to) == 1:
                # Single-qubit controlled gate
                self._apply_multi_control_one_qubit_gate(gate, parallel)
            elif len(gate.controlled_on) == 0:
                # Multi-qubit uncontrolled gate
                self._apply_multi_qubit_gate(gate, parallel)
            else:
                # Multi-qubit controlled gate
                self._apply_multi_control_multi_qubit_gate(gate, parallel)

    fn _apply_one_qubit_gate(mut self, gate: Gate[Self.type, Self.tol], parallel: Bool):
        '''Apply a single qubit gate to the statevector.
        Alg. 2 of [1].
        '''
        var t: UInt = gate.applied_to[0]

        if parallel:
            @parameter
            fn op_sv_elem(n: Int):
                var ib: UInt = insert_bit(UInt(n), t, 0)
                var ig: UInt = flip_bit(ib, t)
                var b: ComplexScalar[Self.type] = self.sv[ib]
                var g: ComplexScalar[Self.type] = self.sv[ig]
                self.sv[ib] = gate.matrix.load_idx[1](0) * b + gate.matrix.load_idx[1](1) * g
                self.sv[ig] = gate.matrix.load_idx[1](2) * b + gate.matrix.load_idx[1](3) * g
            parallelize[op_sv_elem](2 ** (self._qc.n_qubits - 1))
        else:
            var ib: UInt
            var ig: UInt
            var b: ComplexScalar[Self.type]
            var g: ComplexScalar[Self.type]
            for n in range(2 ** (self._qc.n_qubits - 1)):
                ib = insert_bit(UInt(n), t, 0)
                ig = flip_bit(ib, t)
                b = self.sv[ib]
                g = self.sv[ig]
                self.sv[ib] = gate.matrix.load_idx[1](0) * b + gate.matrix.load_idx[1](1) * g
                self.sv[ig] = gate.matrix.load_idx[1](2) * b + gate.matrix.load_idx[1](3) * g
    
    fn _apply_multi_control_one_qubit_gate(
        mut self, 
        gate: Gate[Self.type, Self.tol], 
        parallel: Bool
    ):
        '''Apply a single-qubit gate controlled on one or more qubits to the statevector.
        Alg. 3 of [1].
        '''
        var s: UInt = len(gate.controlled_on)

        # Get a list of the control qubit indices
        var q = List[UInt, True](capacity=s + 1)
        for qubit in gate.controlled_on:
            q.append(qubit[])
        
        # Target qubit index
        var t: UInt = gate.applied_to[0]
        
        # Union the controlled and target qubits and sort the list
        q.append(t)
        sort[T=UInt](q)

        if parallel:
            @parameter
            fn op_sv_elem(j: Int):
                var ig: UInt = insert_bits(UInt(j), q, 1)
                var ib: UInt = flip_bit(ig, t)
                var b: ComplexScalar[Self.type] = self.sv[ib]
                var g: ComplexScalar[Self.type] = self.sv[ig]
                self.sv[ib] = gate.matrix.load_idx[1](0) * b + gate.matrix.load_idx[1](1) * g
                self.sv[ig] = gate.matrix.load_idx[1](2) * b + gate.matrix.load_idx[1](3) * g
            parallelize[op_sv_elem](2 ** (self._qc.n_qubits - s))
        else:
            var ib: UInt
            var ig: UInt
            var b: ComplexScalar[Self.type]
            var g: ComplexScalar[Self.type]
            for j in range(2 ** (self._qc.n_qubits - s)):
                ig = insert_bits(j, q, 1)
                ib = flip_bit(ig, t)
                b = self.sv[ib]
                g = self.sv[ig]
                self.sv[ib] = gate.matrix.load_idx[1](0) * b + gate.matrix.load_idx[1](1) * g
                self.sv[ig] = gate.matrix.load_idx[1](2) * b + gate.matrix.load_idx[1](3) * g

    fn _apply_multi_qubit_gate(mut self, gate: Gate[Self.type, Self.tol], parallel: Bool):
        '''Apply a multi-qubit gate with no controls to the statevector.
        Alg. 4 of [1].
        '''
        var n: UInt = len(gate.applied_to)

        # Form a list of qubit indices and a sorted version
        var t = List[UInt, True](capacity=n)
        for qubit in gate.applied_to:
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
                    v[j] = self.sv[i]
                
                for j in range(2 ** n):
                    i = set_bits(z, t, UInt(j))
                    self.sv[i] = 0

                    for l in range(2 ** n):
                        self.sv[i] += gate.matrix.load_crd[1](j, l) * v[l]
            parallelize[op_sv_elem](2 ** (self._qc.n_qubits - n))
        else:
            var v = List[ComplexScalar[Self.type], True](length=2 ** n, fill=0)
            var z: UInt
            var i: UInt
            for k in range(2 ** (self._qc.n_qubits - n)):
                z = insert_bits(UInt(k), q, 0)

                for j in range(2 ** n):
                    i = set_bits(z, t, UInt(j))
                    v[j] = self.sv[i]
                
                for j in range(2 ** n):
                    i = set_bits(z, t, UInt(j))
                    self.sv[i] = 0

                    for l in range(2 ** n):
                        self.sv[i] += gate.matrix.load_crd[1](j, l) * v[l]
    
    fn _apply_multi_control_multi_qubit_gate(
        mut self, 
        gate: Gate[Self.type, Self.tol], 
        parallel: Bool
    ):
        '''Apply a multi-controlled multi-qubit gate to the statevector.
        Adaptation of Alg. 4 of [1].
        '''
        alias zero = ComplexScalar[Self.type](0)
        alias one = ComplexScalar[Self.type](1)

        var n: UInt = len(gate.applied_to) + len(gate.controlled_on)
        var N: Int = 2 ** n

        # Form a list of qubit indices and a sorted version
        var t = List[UInt, True](capacity=n)
        for qubit in gate.applied_to:
            t.append(qubit[])
        for qubit in gate.controlled_on:
            t.append(qubit[])
        var q: List[UInt, True] = t.copy()
        sort[T=UInt](q)
        
        var rc_start: Int = N - 2 ** len(gate.applied_to)
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
                    v[j] = self.sv[i]
                
                for j in range(rc_start, N):
                    i = set_bits(z, t, UInt(j))
                    self.sv[i] = 0
                    for l in range(rc_start, N):
                        self.sv[i] += gate.matrix.load_crd[1](j - rc_start, l - rc_start) * v[l]
            parallelize[op_sv_elem](2 ** (self._qc.n_qubits - n))
        else:
            var v = List[ComplexScalar[Self.type], True](length=N, fill=0)
            var z: UInt
            var i: UInt
            for k in range(2 ** (self._qc.n_qubits - n)):
                z = insert_bits(UInt(k), q, 0)

                for j in range(N):
                    i = set_bits(z, t, UInt(j))
                    v[j] = self.sv[i]
                
                for j in range(rc_start, N):
                    i = set_bits(z, t, UInt(j))
                    self.sv[i] = 0
                    for l in range(rc_start, N):
                        self.sv[i] += gate.matrix.load_crd[1](j - rc_start, l - rc_start) * v[l]
    
    fn _measure(mut self, gate: Gate[Self.type, Self.tol]) -> String:
        '''Measure a subset of qubits and update the statevector.'''
        var n: Int = self._qc.n_qubits
        var N: Int = 2 ** n
        var m: Int = len(gate.applied_to)
        var M: Int = 2 ** m
        var probabilities = List[Scalar[Self.type], True](length=N, fill=0)

        var measured_qubits_idx = List[Int, True]()
        for q in gate.applied_to:
            measured_qubits_idx.append(q[])

        # Compute probabilities for each possible outcome
        var outcome: Int
        var bit: Int
        for i in range(N):
            # Extract the bits at the measured positions
            outcome = 0
            for j in range(m):
                bit = (i >> measured_qubits_idx[j]) & 1
                outcome |= (bit << j)
            probabilities[outcome] += self.sv[i].squared_norm()
        
        # Normalize probabilities to ensure they sum to one
        var total_prob: Scalar[Self.type] = 0
        for p in probabilities:
            total_prob += p[]
        if abs(total_prob - 1) >= Self.tol:
            for i in range(N):
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
                bit = (i >> measured_qubits_idx[j]) & 1
                expected_bit = (selected_outcome >> j) & 1
                if bit != expected_bit:
                    consistent = False
                    break
            
            # Update statevector
            if consistent:
                new_statevector[i] = self.sv[i] / norm_factor

        # TODO: Check if this measure function gives correct results
        
        self.sv = new_statevector^
        print(bin(selected_outcome))
        return bin(selected_outcome)


    # fn _expand_controls(self, gate: Gate[Self.type, Self.tol]) -> Gate[Self.type, Self.tol]:
    #     '''Expand the gate's unitary to include control qubits, and return the corresponding Gate.
    #     '''

    #     var n: Int = 2 ** len(gate.applied_to)  # Dimension of uncontrolled unitary
    #     var N: Int = 2 ** (len(gate.controlled_on) + len(gate.applied_to))  # Dimension of controlled unitary 

    #     var full_unitary = CMatrix[Self.type](rows=N, cols=N, fill_zeros=True)

    #     # Copy the uncontrolled unitary into the bottom right corner of the controlled unitary
    #     var idx_base: Int = N * (N + 1 - n) - n
    #     for r in range(gate.matrix.rows):
    #         var idx: Int = idx_base + N * r
    #         memcpy(full_unitary.re + idx, gate.matrix.re + r * n, n)
    #         memcpy(full_unitary.im + idx, gate.matrix.im + r * n, n)
        
    #     # Fill the diagonal above the inset unitary with ones
    #     var stride: Int = N + 1
    #     @parameter
    #     fn fill_diag[simd_width: Int](p: Int):
    #         (full_unitary.re + p * stride).strided_store[width=simd_width](1, stride)
    #     vectorize[fill_diag, simdwidthof[Self.type]()](N - n)

    #     return Gate[Self.type, Self.tol](
    #         n_qubits=len(gate.controlled_on) + len(gate.applied_to),
    #         name=gate.name + '_expanded',
    #         matrix=full_unitary^, 
    #         applied_to=gate.applied_to + gate.controlled_on,
    #         params=gate.params,
    #         controlled_on = List[Qubit, True]()
    #     )
    