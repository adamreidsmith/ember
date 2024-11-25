from cplx import CMatrix
from .gates import *


@value
struct AppliedGate[type: DType]():
    var gate: String
    var qubits: List[Int, True]
    var matrix: CMatrix[Self.type]

    @always_inline
    fn __init__(inout self, gate: String, qubits: List[Int, True], matrix: CMatrix[Self.type]):
        self.gate = gate
        self.qubits = qubits
        self.matrix = matrix

@value
struct QuantumCircuit[type: DType]:
    var n_qubits: Int
    var n_clbits: Int
    var _data: List[AppliedGate[Self.type]]
    
    fn __init__(inout self, n_qubits: Int, n_clbits: Int = 0):
        self.n_qubits = n_qubits
        self.n_clbits = n_clbits
        self._data = List[AppliedGate[Self.type]]()

    fn apply(inout self, gate: Gate[Self.type], *qubits: Int) raises:
        if len(qubits) != gate.n_qubits:
            raise Error(
                'Gate ' + gate.name + ' is a ' + str(gate.n_qubits) 
                + '-qubit gate, but ' + str(len(qubits)) + ' qubits were provided'
            )
        var qubit_list = List[Int, True]()
        for q in qubits:
            qubit_list.append(q)
        self._data.append(AppliedGate[Self.type](gate.name, qubit_list, gate.matrix)) 
