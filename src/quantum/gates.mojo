from math import sqrt, cos, sin, pi
from cplx import CMatrix, ComplexScalar


alias _gate_function = fn[type: DType](params: VariadicList[Scalar[type]]) -> CMatrix[type]

trait GateLike(Formattable, Sized, StringableCollectionElement): ...

@value
struct Gate[
    type: DType,
    name: String,
    n_qubits: Int,
    n_params: Int,
    gate_func: _gate_function
](GateLike):
    var params: VariadicList[Scalar[Self.type]]
    var matrix: CMatrix[Self.type]
    
    @always_inline
    fn __init__(inout self, *params: Scalar[Self.type]) raises:
        if len(params) != Self.n_params:
            raise Error(
                'Gate ' + Self.name + ' requires ' + str(Self.n_params) 
                + ' parameters, but ' + str(len(params)) + ' were provided'
            )
        self.params = params
        self.matrix = gate_func[Self.type](params)
    
    @no_inline
    fn __str__(self) -> String:
        return Self.name
    
    @no_inline
    fn format_to(self, inout writer: Formatter):
        writer.write(self.__str__())
    
    @always_inline
    fn __len__(self) -> Int:
        return Self.n_qubits

# TODO: Change `params: VariadicList[Scalar[type]]` argument to `*params: Scalar[type]` once argument unpacking is supported

# Unparameterized single-qubit gate functions ###

fn x[type: DType](params: VariadicList[Scalar[type]]) -> CMatrix[type]:
    return CMatrix[type](2, 2,
        0, 1, 
        1, 0,
    )

fn y[type: DType](params: VariadicList[Scalar[type]]) -> CMatrix[type]:
    return CMatrix[type](2, 2,
        0, ComplexScalar[type](0, -1), 
        ComplexScalar[type](0, 1), 0,
    )

fn z[type: DType](params: VariadicList[Scalar[type]]) -> CMatrix[type]:
    return CMatrix[type](2, 2,
        1, 0, 
        0, -1,
    )

fn h[type: DType](params: VariadicList[Scalar[type]]) -> CMatrix[type]:
    return CMatrix[type](2, 2,
        1, 1, 
        1, -1,
    ) / sqrt(2)

fn s[type: DType](params: VariadicList[Scalar[type]]) -> CMatrix[type]:
    return CMatrix[type](2, 2,
        1, 0, 
        0, ComplexScalar[type](0, 1),
    )

fn t[type: DType](params: VariadicList[Scalar[type]]) -> CMatrix[type]:
    return CMatrix[type](2, 2,
        1, 0, 
        0, ComplexScalar[type](0, pi / 4).exp(),
    )

# Unparameterized multi-qubit gate functions ####

fn cx[type: DType](params: VariadicList[Scalar[type]]) -> CMatrix[type]:
    r'''\ket{0}\bra{0} \otimes I + \ket{0}\bra{0} \otimes X.'''
    return CMatrix[type](4, 4,
        1, 0, 0, 0, 
        0, 1, 0, 0, 
        0, 0, 0, 1, 
        0, 0, 1, 0,
    )

fn ccx[type: DType](params: VariadicList[Scalar[type]]) -> CMatrix[type]:
    r'''\ket{0}\bra{0} \otimes I^2 + \ket{0}\bra{0} \otimes CX.'''
    return CMatrix[type](8, 8,
        1, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 0, 0, 0, 1, 0,
    )

# Parameterized single-qubit gate functions #####

fn rx[type: DType](params: VariadicList[Scalar[type]]) -> CMatrix[type]:
    var a = ComplexScalar[type](cos(params[0] / 2), 0)
    var b = ComplexScalar[type](0, -sin(params[0] / 2))
    return CMatrix[type](2, 2,
        a, b, 
        b, a,
    )

fn ry[type: DType](params: VariadicList[Scalar[type]]) -> CMatrix[type]:
    var a = ComplexScalar[type](cos(params[0] / 2), 0)
    var b = ComplexScalar[type](sin(params[0] / 2), 0)
    return CMatrix[type](2, 2,
        a, -b, 
        b, a,
    )

fn rz[type: DType](params: VariadicList[Scalar[type]]) -> CMatrix[type]:
    var c: Scalar[type] = cos(params[0] / 2)
    var s: Scalar[type] = sin(params[0] / 2)
    return CMatrix[type](2, 2,
        ComplexScalar[type](c, -s), 0, 
        0, ComplexScalar[type](c, s),
    )

fn u[type: DType](params: VariadicList[Scalar[type]]) -> CMatrix[type]:
    var ct: Scalar[type] = cos(params[0] / 2)
    var st: Scalar[type] = sin(params[0] / 2)
    return CMatrix[type](2, 2,
        ComplexScalar[type](ct, 0), 
        ComplexScalar[type](-cos(params[2]) * st, -sin(params[2]) * st), 
        ComplexScalar[type](cos(params[1]) * st, sin(params[1]) * st), 
        ComplexScalar[type](cos(params[1] + params[2]) * ct, sin(params[1] + params[2]) * ct),
    )

# Parameterized multi-qubit gate functions ######
...


# Unparameterized single-qubit gates ############

alias X = Gate[name='X', n_qubits=1, n_params=0, gate_func=x]
alias Y = Gate[name='Y', n_qubits=1, n_params=0, gate_func=y]
alias Z = Gate[name='Z', n_qubits=1, n_params=0, gate_func=z]
alias H = Gate[name='H', n_qubits=1, n_params=0, gate_func=h]
alias S = Gate[name='S', n_qubits=1, n_params=0, gate_func=s]
alias T = Gate[name='T', n_qubits=1, n_params=0, gate_func=t]

# Unparameterized multi-qubit gates #############

alias CX = Gate[name='CX', n_qubits=2, n_params=0, gate_func=cx]
alias CCX = Gate[name='CCX', n_qubits=3, n_params=0, gate_func=ccx]

# Parameterized single-qubit gates ##############

alias RX = Gate[name='RX', n_qubits=1, n_params=1, gate_func=rx]
alias RY = Gate[name='RY', n_qubits=1, n_params=1, gate_func=ry]
alias RZ = Gate[name='RZ', n_qubits=1, n_params=1, gate_func=rz]
alias U = Gate[name='U', n_qubits=1, n_params=3, gate_func=u]

# Parameterized multi-qubit gates ###############
...
