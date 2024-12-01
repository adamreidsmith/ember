from math import ceil, log10
from collections import Set

from ..cplx import CMatrix
from .gates import *
from .bit import Clbit, Qubit


@value
struct QuantumCircuit[type: DType](Stringable, Formattable):
    var qubits: List[Qubit, True]
    var clbits: List[Clbit, True]
    var _data: List[Gate[Self.type]]
    
    @always_inline
    fn __init__(inout self, n_qubits: Int, n_clbits: Int = 0) raises:
        if n_qubits < 1:
            raise Error('Quantum circuit must have at least one qubit')
        self.qubits = List[Qubit, True](Qubit(0))
        for q in range(1, n_qubits):
            self.qubits.append(Qubit(q))
        self.clbits = List[Clbit, True]()
        for c in range(n_clbits):
            self.clbits.append(Clbit(c))
        self._data = List[Gate[Self.type]]()
    
    @always_inline
    fn __init__(inout self, owned qubits: List[Qubit, True], owned clbits: List[Clbit, True] = List[Clbit, True]()) raises:
        if len(qubits) < 1:
            raise Error('Quantum circuit must have at least one qubit')
        self.qubits = qubits^
        self.clbits = clbits^
        self._data = List[Gate[Self.type]]()

    @always_inline
    fn apply(inout self, gate: Gate[Self.type]) raises:
        for q in gate.applied_to:
            if q[] not in self.qubits:
                raise Error(
                    'Gate ' + gate.name + ' is applied to qubit ' + str(q[]) 
                    + ', but qubit ' + str(q[]) + ' is not in the circuit'
                )
        for q in gate.controlled_on:
            if q[] not in self.qubits:
                raise Error(
                    'Gate ' + gate.name + ' is controlled on qubit ' + str(q[]) 
                    + ', but qubit ' + str(q[]) + ' is not in the circuit'
                )
        self._data.append(gate)

    @no_inline
    fn __str__(self) -> String:
        return ''
    
    @no_inline
    fn format_to(self, inout writer: Formatter):
        writer.write(self.__str__())
