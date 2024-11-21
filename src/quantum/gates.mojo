from cplx import CMatrix


@value
struct Gate[type: DType, qubits: Int]:
    alias dim = 2 ** qubits

    var matrix: CMatrix[type, Self.dim, Self.dim]

    fn __init__(inout self):
        self.matrix = CMatrix[type, Self.dim, Self.dim]()
