from math import ceil, log10
from collections import Set

from ..cplx import CMatrix
from .gates import *
from ..config import DEFAULT_TOL

# TODO: Add classical bit control

@value
struct QuantumCircuit[type: DType, tol: Scalar[type] = DEFAULT_TOL](Stringable, Writable, Movable):
    '''A quantum circuit.

    Parameters:
        type: A type for the circuit data.
        tol: A tolerance for unitarity and closeness checks.
    '''

    var n_qubits: Int
    '''The number of qubits in the quantum circuit.'''
    var n_clbits: Int
    '''The number of classical bits in the quantum circuit.'''
    var _data: List[Gate[Self.type, Self.tol]]
    '''The gates applied to the qubits in the quantum circuit.'''
    
    fn __init__(out self, n_qubits: Int, n_clbits: Int = 0) raises:
        '''Initialize a QuantumCircuit.
        
        Args:
            n_qubits: The number of qubits in the quantum circuit.
            n_clbits: The number of classical bits in the quantum circuit.
        '''
        if n_qubits < 1:
            raise Error('Quantum circuit must have at least one qubit')
        self.n_qubits = n_qubits
        self.n_clbits = n_clbits
        self._data = List[Gate[Self.type, Self.tol]]()

    fn apply(mut self, gate: Gate[Self.type, Self.tol]) raises:
        '''Apply a gate to the quantum circuit.
        
        Args:
            gate: The gate to apply.
        '''
        for q in gate.applied_to:
            if q[] < 0 or q[] >= self.n_qubits:
                raise Error(
                    'Gate ' + gate.name + ' contains invalid qubit specifier: ' + String(q[])
                    + '. Quantum circuit has ' + String(self.n_qubits) + ' qubits indexed 0 to ' 
                    + String(self.n_qubits - 1) + '.'
                )
        if gate._is_measure:
            for q in gate._measure_targs:
                if q[] < 0 or q[] >= self.n_clbits:
                    raise Error(
                        'Gate ' + gate.name + ' contains invalid classical bit specifier: '
                        + String(q[]) + '. Quantum circuit has ' + String(self.n_clbits)
                        + ' classical bits labeled 0 to ' + String(self.n_clbits - 1) + '.'
                    )
        else:
            for q in gate.controlled_on:
                if q[] < 0 or q[] >= self.n_qubits:
                    raise Error(
                        'Gate ' + gate.name + ' contains invalid control qubit specifier: '
                        + String(q[]) + '. Quantum circuit has ' + String(self.n_qubits)
                        + ' qubits labeled 0 to ' + String(self.n_qubits - 1) + '.'
                    )
        self._data.append(gate)

    @no_inline
    fn __str__(self) -> String:
        '''Convert the quantum circuit to a string.

        Returns:
            A string representation of the quantum circuit.
        '''
        alias max_width: Int = 120

        fn max(l: List[Int, True]) -> Int:
            var mx: Int = l[0]
            for e in l[1:]:
                if e[] > mx:
                    mx = e[]
            return mx

        var lines = List[List[String]](List[String]('|0⟩ -')) * self.n_qubits
        for gate in self._data:
            var gate_str = String(gate[])
            for i in range(self.n_qubits):
                var mod_gate_str: String
                if gate[]._is_measure:
                    var measure_specifier_len: Int = len(String(max(gate[]._measure_targs))) + 2
                    var measure_to: Int
                    try:
                        measure_to = gate[]._measure_targs[gate[].applied_to.index(i)]
                    except:
                        measure_to = 0
                    var n_dashes: Int = measure_specifier_len - len(String(measure_to)) - 1
                    mod_gate_str = gate_str + '-' * n_dashes + '>' + String(measure_to)
                else:
                    mod_gate_str = gate_str
                if i in gate[].applied_to:
                    lines[i].append(mod_gate_str)
                elif i in gate[].controlled_on and not gate[]._is_measure:
                    var n_dashes: Int = len(mod_gate_str) - 1
                    var left_dashes: Int = n_dashes // 2
                    var right_dashes: Int = n_dashes - left_dashes
                    lines[i].append('-' * left_dashes + '*' + '-' * right_dashes)
                else:
                    lines[i].append('-' * len(mod_gate_str))
            for line in lines:
                line[][-1] += '-'

        var lines_split = List[List[String]]()
        for line in lines:
            var line_str = List[String]('')
            for gate_str in line[]:
                var gs: String = gate_str[]
                if len(gate_str[]) > max_width:
                    gs = gate_str[][:max_width]
                if len(line_str[-1]) + len(gs) + 1 > max_width:
                    line_str.append('')
                line_str[-1] += gs
            lines_split.append(line_str)
        
        var ret: String = ''
        for i in range(len(lines_split[0])):
            for line in lines_split:
                ret += line[][i] + '\n'
            ret += '\n'
        ret = ret[:-2]
        return ret

    @no_inline
    fn write_to[W: Writer](self, mut writer: W):
        '''Write the quantum circuit to a writer.

        Args:
            writer: The writer to write to.
        '''
        writer.write(String(self))
