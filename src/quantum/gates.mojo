from math import sqrt, cos, sin, pi
from collections import Set
from utils import Variant

from ..cplx import CMatrix, ComplexScalar
from ..config import DEFAULT_TOL


@value
struct Gate[type: DType, tol: Scalar[type] = DEFAULT_TOL](Writable, Sized, Stringable, CollectionElement):
    '''A quantum gate.'''

    var name: String
    '''An identifier for the gate.'''
    var n_qubits: Int
    '''The number of qubits the gate acts on.'''
    var matrix: CMatrix[Self.type]
    '''The unitary matrix that implements the gate.'''
    var applied_to: List[Int, True]
    '''The qubits the gate is applied to.'''
    var controlled_on: List[Int, True]
    '''The qubits the gate is controlled on.'''
    var params: List[Scalar[Self.type], True]
    '''The parameter values applied to the gate.'''

    # Special fields for the measurement gate
    var _is_measure: Bool
    var _measure_targs: List[Int, True]

    fn __init__(
        out self,
        owned name: String,
        owned matrix: CMatrix[Self.type], 
        owned qubits: List[Int, True] = List[Int, True](),
        owned params: List[Scalar[Self.type], True] = List[Scalar[Self.type], True](),
    ) raises:
        '''Initialize a Gate with a name, unitary matrix, set of qubits, and parameters.'''
        if len(Set(qubits)) != len(qubits):
            raise Error('Duplicate qubit specified for gate ' + name)
        # Ensure the matrix is square and its dimension is a power of 2
        if not matrix.is_square() or (matrix.rows & (matrix.rows - 1)) != 0:
            raise Error('Invalid gate matrix for gate ' + name)
        if not matrix.is_unitary[Self.tol]():
            raise Error('Gate matrix is not unitary for gate ' + name)
        var n_qubits: Int = len(bin(matrix.rows)) - 3
        if len(qubits) != n_qubits:
            raise Error('Invalid number of qubits for gate ' + name)

        self.n_qubits = n_qubits
        self.name = name^ 
        self.matrix = matrix^
        self.applied_to = qubits^
        self.params = params^
        self.controlled_on = List[Int, True]()

        # Placeholders. Only applicable to the measurement gate
        self._is_measure = False
        self._measure_targs = List[Int, True]()

    @staticmethod
    fn _measure(owned qubits: List[Int, True], clbits: List[Int, True]) -> Self:
        '''Creates a measurement gate. `qubits` and `clbits` must be lists of the same length.'''
        return Self(
            name='Measure',
            n_qubits=len(qubits),
            matrix=CMatrix[Self.type](rows=0, cols=0, fill_zeros=False),
            applied_to=qubits^,
            controlled_on=List[Int, True](),
            params=List[Scalar[Self.type], True](),
            _is_measure=True,
            _measure_targs=clbits,
        )
    
    @no_inline
    fn __str__(self) -> String:
        if len(self.params) == 0:
            return self.name
        var str_rep: String = self.name + '('
        for p in self.params:
            str_rep += String(p[]) + ', '
        str_rep = str_rep[:-2]
        return str_rep + ')'
    
    # TODO: Debug seg fault with RX 
    # @no_inline
    # fn __repr__(self) -> String:
    #     var str_rep: String = '<Gate ' + self.name + ': '
    #     if not self.applied_to and not self.controlled_on and not self.params:
    #         return str_rep[:-2] + '>'
    #     var qubits_str: String = ''
    #     if self.applied_to:
    #         qubits_str = 'qubtis=(' + str(self.applied_to[0])
    #         for qubit in self.applied_to[1:]:
    #             qubits_str += ', ' + str(qubit[])
    #         qubits_str += ')'
    #     var controls_str: String = ''
    #     if self.applied_to:
    #         controls_str = 'controls=(' + str(self.controlled_on[0])
    #         for qubit in self.controlled_on[1:]:
    #             controls_str += ', ' + str(qubit[])
    #         controls_str += ')'
    #     var params_str: String = ''
    #     if self.params:
    #         params_str = 'params=(' + str(self.params[0])
    #         for param in self.params[1:]:
    #             params_str += ', ' + str(param[])
    #         params_str += ')'
    #     if qubits_str:
    #         str_rep += qubits_str if str_rep[-2:] == ': ' else ', ' + qubits_str
    #     if controls_str:
    #         str_rep += controls_str if str_rep[-2:] == ': ' else ', ' + controls_str
    #     if params_str:
    #         str_rep += params_str if str_rep[-2:] == ': ' else ', ' + params_str
    #     return str_rep + '>'
        
    @no_inline
    fn write_to[W: Writer](self, mut writer: W):
        writer.write(String(self))
    
    @always_inline
    fn __len__(self) -> Int:
        return self.n_qubits

    fn control(mut self, *qubits: Int) raises -> Self:
        '''Control the gate on a qubit or set of qubits.'''
        if self._is_measure:
            raise Error('Cannot control a measurement gate.')
        for q in qubits:
            if q in self.applied_to:
                raise Error(
                    'Cannot control gate ' + self.name + ' on qubit ' 
                    + String(q) + ' on which it is applied'
                )
            if q not in self.controlled_on:
                self.controlled_on.append(q)
                self.n_qubits += 1
        return self
    
    fn __eq__(self, other: Self) -> Bool:
        '''Gates are considered equal if they have the same name, they are applied to and controlled
        on the same qubits, and their matrices are (close to) equal.'''
        if (
            self.name != other.name 
            or self.n_qubits != other.n_qubits
            or self.applied_to != other.applied_to
            or Set(self.controlled_on) != Set(other.controlled_on)
            # or not self.matrix.matrix_equals(other.matrix)
            or not self.matrix.is_close[Self.tol](other.matrix)
        ):
            return False
        return True

# Measurement

fn Measure[type: DType, tol: Scalar[type] = DEFAULT_TOL](
    qubits: List[Int, True], clbits: List[Int, True]
) raises -> Gate[type, tol]:
    '''Create a measurement on a set of qubits.'''
    if len(qubits) == 0:
        raise Error('Measure received no qubit specifiers')
    if len(qubits) != len(Set(qubits)):
        raise Error('List of qubits must be unique')
    if len(clbits) != len(Set(clbits)):
        raise Error('List of clbits must be unique')
    if len(qubits) != len(clbits):
        raise Error('qubits and clbits must be lists of the same length')
    return Gate[type, tol]._measure(qubits, clbits)

fn Measure[type: DType, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int, clbit: Int
) raises -> Gate[type, tol]:
    '''Create a measurement on a qubit.'''
    return Gate[type, tol]._measure(List[Int, True](qubit), List[Int, True](clbit))

# Unparameterized single-qubit gates ############

fn x[type: DType]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        0, 1, 
        1, 0,
    )
fn X[type: DType, tol: Scalar[type] = DEFAULT_TOL](q: Int) raises -> Gate[type, tol]:
    '''Create an X gate applied to qubit `q`.'''
    return Gate[type, tol]('X', x[type](), List[Int, True](q))

fn y[type: DType]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        0, ComplexScalar[type](0, -1), 
        ComplexScalar[type](0, 1), 0,
    )
fn Y[type: DType, tol: Scalar[type] = DEFAULT_TOL](q: Int) raises -> Gate[type, tol]:
    '''Create a Y gate applied to qubit `q`.'''
    return Gate[type, tol]('Y', y[type](), List[Int, True](q))

fn z[type: DType]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        1, 0, 
        0, -1,
    )
fn Z[type: DType, tol: Scalar[type] = DEFAULT_TOL](q: Int) raises -> Gate[type, tol]:
    '''Create a Z gate applied to qubit `q`.'''
    return Gate[type, tol]('Z', z[type](), List[Int, True](q))

fn h[type: DType]() raises -> CMatrix[type]:
    var inv_sqrt2: Scalar[type] = 1 / sqrt(2.0).cast[type]()
    return CMatrix[type](2, 2,
        inv_sqrt2, inv_sqrt2, 
        inv_sqrt2, -inv_sqrt2,
    )
fn H[type: DType, tol: Scalar[type] = DEFAULT_TOL](q: Int) raises -> Gate[type, tol]:
    '''Create a Hadamard gate applied to qubit `q`.'''
    return Gate[type, tol]('H', h[type](), List[Int, True](q))     

fn s[type: DType]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        1, 0, 
        0, ComplexScalar[type](0, 1),
    )
fn S[type: DType, tol: Scalar[type] = DEFAULT_TOL](q: Int) raises -> Gate[type, tol]:
    '''Create an S gate applied to qubit `q`.'''
    return Gate[type, tol]('S', s[type](), List[Int, True](q))

fn t[type: DType]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        1, 0, 
        0, ComplexScalar[type](0, pi / 4).exp(),
    )
fn T[type: DType, tol: Scalar[type] = DEFAULT_TOL](q: Int) raises -> Gate[type, tol]:
    '''Create a T gate applied to qubit `q`.'''
    return Gate[type, tol]('T', t[type](), List[Int, True](q))

# Unparameterized multi-qubit gates #############

fn CX[type: DType, tol: Scalar[type] = DEFAULT_TOL](control: Int, target: Int) raises -> Gate[type, tol]:
    '''Create a controlled-X gate applied to qubit `target` and controlled on qubit `control`.'''
    var cx = Gate[type, tol]('CX', x[type](), List[Int, True](target))
    return cx.control(control)

fn CCX[type: DType, tol: Scalar[type] = DEFAULT_TOL](control1: Int, control2: Int, target: Int) raises -> Gate[type, tol]:
    '''Create a controlled-controlled-X gate applied to qubit `target` and controlled on qubits
     `control1` and `control2`.'''
    var ccx = Gate[type, tol]('CCX', x[type](), List[Int, True](target))
    return ccx.control(control1, control2)

# Parameterized single-qubit gates ##############

fn rx[type: DType](t: Scalar[type]) raises -> CMatrix[type]:
    var a = ComplexScalar[type](cos(t / 2), 0)
    var b = ComplexScalar[type](0, -sin(t / 2))
    return CMatrix[type](2, 2,
        a, b,
        b, a,
    )
fn RX[type: DType, tol: Scalar[type] = DEFAULT_TOL](q: Int, theta: Scalar[type]) raises -> Gate[type, tol]:
    '''Create an RX rotation gate with angle `theta` applied to qubit `q`.'''
    return Gate[type, tol]('RX', rx[type](theta), List[Int, True](q), List[Scalar[type], True](theta))

fn ry[type: DType](t: Scalar[type]) raises -> CMatrix[type]:
    var a = ComplexScalar[type](cos(t / 2), 0)
    var b = ComplexScalar[type](sin(t / 2), 0)
    return CMatrix[type](2, 2,
        a, -b,
        b, a,
    )
fn RY[type: DType, tol: Scalar[type] = DEFAULT_TOL](q: Int, theta: Scalar[type]) raises -> Gate[type, tol]:
    '''Create an RY rotation gate with angle `theta` applied to qubit `q`.'''
    return Gate[type, tol]('RY', ry[type](theta), List[Int, True](q), List[Scalar[type], True](theta))

fn rz[type: DType](t: Scalar[type]) raises -> CMatrix[type]:
    var c: Scalar[type] = cos(t / 2)
    var s: Scalar[type] = sin(t / 2)
    return CMatrix[type](2, 2,
        ComplexScalar[type](c, -s), 0, 
        0, ComplexScalar[type](c, s),
    )
fn RZ[type: DType, tol: Scalar[type] = DEFAULT_TOL](q: Int, theta: Scalar[type]) raises -> Gate[type, tol]:
    '''Create an RZ rotation gate with angle `theta` applied to qubit `q`.'''
    return Gate[type, tol]('RZ', rz[type](theta), List[Int, True](q), List[Scalar[type], True](theta))

fn u[type: DType](t: Scalar[type], p: Scalar[type], l: Scalar[type]) raises -> CMatrix[type]:
    var ct: Scalar[type] = cos(t / 2)
    var st: Scalar[type] = sin(t / 2)
    return CMatrix[type](2, 2,
        ComplexScalar[type](ct, 0), 
        ComplexScalar[type](-cos(l) * st, -sin(l) * st), 
        ComplexScalar[type](cos(p) * st, sin(p) * st), 
        ComplexScalar[type](cos(p + l) * ct, sin(p + l) * ct),
    )
fn U[type: DType, tol: Scalar[type] = DEFAULT_TOL](q: Int, theta: Scalar[type], phi: Scalar[type], lbda: Scalar[type]) raises -> Gate[type, tol]:
    '''Create a U gate applied to qubit `q`.'''
    return Gate[type, tol]('U', u[type](theta, phi, lbda), List[Int, True](q), List[Scalar[type], True](theta, phi, lbda))

# Parameterized multi-qubit gates ###############

# TODO: Add more gates
