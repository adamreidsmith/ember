from math import ceil, log10
from collections import Set

from ..cplx import CMatrix
from .gates import *
from .bit import Clbit, Qubit


@value
struct QuantumCircuit[type: DType](Stringable, Writable):
    var n_qubits: Int
    '''The number of qubits in the quantum circuit.'''
    var qubits: List[Qubit, True]
    '''The Qubits in the quantum circuit.'''
    var clbits: List[Clbit, True]
    '''The classical bits in the quantum circuit.'''
    var _data: List[Gate[Self.type]]
    '''The gates applied to the qubits in the quantum circuit.'''
    
    @always_inline
    fn __init__(out self, n_qubits: Int, n_clbits: Int = 0) raises:
        '''Initialize a QuantumCircuit with `n_qubits` qubits and `n_clbits` classical bits.'''
        if n_qubits < 1:
            raise Error('Quantum circuit must have at least one qubit')
        self.n_qubits = n_qubits
        self.qubits = List[Qubit, True](Qubit(0))
        self.clbits = List[Clbit, True]()
        self._data = List[Gate[Self.type]]()
        for q in range(1, n_qubits):
            self.qubits.append(Qubit(q))
        for c in range(n_clbits):
            self.clbits.append(Clbit(c))
    
    @always_inline
    fn __init__(out self, owned qubits: List[Qubit, True], owned clbits: List[Clbit, True] = List[Clbit, True]()) raises:
        '''Initialize a QuantumCircuit from lists of qubits and classical bits.'''
        if len(qubits) < 1:
            raise Error('Quantum circuit must have at least one qubit')
        self.n_qubits = len(qubits)
        self.qubits = qubits^
        self.clbits = clbits^
        self._data = List[Gate[Self.type]]()

    @always_inline
    fn apply(mut self, gate: Gate[Self.type]) raises:
        '''Apply `gate` to the quantum circuit.'''
        for q in gate.applied_to:
            if q[] not in self.qubits:
                raise Error(
                    'Gate ' + gate.name + ' is applied to qubit ' + String(q[]) 
                    + ', but qubit ' + String(q[]) + ' is not in the circuit'
                )
        for q in gate.controlled_on:
            if q[] not in self.qubits:
                raise Error(
                    'Gate ' + gate.name + ' is controlled on qubit ' + String(q[]) 
                    + ', but qubit ' + String(q[]) + ' is not in the circuit'
                )
        self._data.append(gate)

    @no_inline
    fn __str__(self) -> String:
        alias max_width: Int = 120
        
        var lines = List[List[String]](List[String]('|0⟩ -')) * self.n_qubits
        for gate in self._data:
            var gate_str = String(gate[])
            for i in range(self.n_qubits):
                if i in gate[].applied_to:
                    lines[i].append(gate_str)
                elif i in gate[].controlled_on:
                    var n_dashes: Int = len(gate_str) - 1
                    var left_dashes: Int = n_dashes // 2
                    var right_dashes: Int = n_dashes - left_dashes
                    lines[i].append('-' * left_dashes + '*' + '-' * right_dashes)
                else:
                    lines[i].append('-' * len(gate_str))
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
        writer.write(String(self))
