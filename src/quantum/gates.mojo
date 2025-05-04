from math import sqrt, cos, sin, pi
from collections import Set
from utils import Variant

from ..cplx import CMatrix, ComplexScalar
from .bit import Qubit
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
    var applied_to: List[Qubit, True]
    '''The qubits the gate is applied to.'''
    var controlled_on: List[Qubit, True]
    '''The qubits the gate is controlled on.'''
    var params: List[Scalar[Self.type], True]
    '''The parameter values applied to the gate.'''

    @always_inline
    fn __init__(
        out self,
        owned name: String,
        owned matrix: CMatrix[Self.type], 
        owned qubits: List[Qubit, True] = List[Qubit, True](),
        owned params: List[Scalar[Self.type], True] = List[Scalar[Self.type], True](),
    ) raises:
        '''Initialize a Gate with a name, unitary matrix, set of qubits, and parameters.'''
        if len(Set[Qubit](qubits)) != len(qubits):
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
        self.controlled_on = List[Qubit, True]()
    
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

    @always_inline
    fn control(mut self, *qubits: Qubit) raises -> Self:
        '''Control the gate on a qubit or set of qubits.'''
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
    
    @always_inline
    fn __eq__(self, other: Self) -> Bool:
        '''Gates are considered equal if they have the same name, they are applied to and controlled
        on the same qubits, and their matrices are (close to) equal.'''
        if (
            self.name != other.name 
            or self.n_qubits != other.n_qubits
            or self.applied_to != other.applied_to
            or Set[Qubit](self.controlled_on) != Set[Qubit](other.controlled_on)
            # or not self.matrix.matrix_equals(other.matrix)
            or not self.matrix.is_close[Self.tol](other.matrix)
        ):
            return False
        return True

# Unparameterized single-qubit gates ############

fn x[type: DType]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        0, 1, 
        1, 0,
    )
fn X[type: DType](q: Qubit) raises -> Gate[type]:
    '''Create an X gate applied to qubit `q`.'''
    return Gate[type]('X', x[type](), List[Qubit, True](q))

fn y[type: DType]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        0, ComplexScalar[type](0, -1), 
        ComplexScalar[type](0, 1), 0,
    )
fn Y[type: DType](q: Qubit) raises -> Gate[type]:
    '''Create a Y gate applied to qubit `q`.'''
    return Gate[type]('Y', y[type](), List[Qubit, True](q))

fn z[type: DType]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        1, 0, 
        0, -1,
    )
fn Z[type: DType](q: Qubit) raises -> Gate[type]:
    '''Create a Z gate applied to qubit `q`.'''
    return Gate[type]('Z', z[type](), List[Qubit, True](q))

fn h[type: DType]() raises -> CMatrix[type]:
    var inv_sqrt2: Scalar[type] = 1 / sqrt(2.0).cast[type]()
    return CMatrix[type](2, 2,
        inv_sqrt2, inv_sqrt2, 
        inv_sqrt2, -inv_sqrt2,
    )
fn H[type: DType](q: Qubit) raises -> Gate[type]:
    '''Create a Hadamard gate applied to qubit `q`.'''
    return Gate[type]('H', h[type](), List[Qubit, True](q))     

fn s[type: DType]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        1, 0, 
        0, ComplexScalar[type](0, 1),
    )
fn S[type: DType](q: Qubit) raises -> Gate[type]:
    '''Create an S gate applied to qubit `q`.'''
    return Gate[type]('S', s[type](), List[Qubit, True](q))

fn t[type: DType]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        1, 0, 
        0, ComplexScalar[type](0, pi / 4).exp(),
    )
fn T[type: DType](q: Qubit) raises -> Gate[type]:
    '''Create a T gate applied to qubit `q`.'''
    return Gate[type]('T', t[type](), List[Qubit, True](q))

# Unparameterized multi-qubit gates #############

fn CX[type: DType](control: Qubit, target: Qubit) raises -> Gate[type]:
    '''Create a controlled-X gate applied to qubit `target` and controlled on qubit `control`.'''
    var cx = Gate[type]('CX', x[type](), List[Qubit, True](target))
    return cx.control(control)

fn CCX[type: DType](control1: Qubit, control2: Qubit, target: Qubit) raises -> Gate[type]:
    '''Create a controlled-controlled-X gate applied to qubit `target` and controlled on qubits
     `control1` and `control2`.'''
    var ccx = Gate[type]('CCX', x[type](), List[Qubit, True](target))
    return ccx.control(control1, control2)

# Parameterized single-qubit gates ##############

fn rx[type: DType](t: Scalar[type]) raises -> CMatrix[type]:
    var a = ComplexScalar[type](cos(t / 2), 0)
    var b = ComplexScalar[type](0, -sin(t / 2))
    return CMatrix[type](2, 2,
        a, b,
        b, a,
    )
fn RX[type: DType](q: Qubit, theta: Scalar[type]) raises -> Gate[type]:
    '''Create an RX rotation gate with angle `theta` applied to qubit `q`.'''
    return Gate[type]('RX', rx[type](theta), List[Qubit, True](q), List[Scalar[type], True](theta))

fn ry[type: DType](t: Scalar[type]) raises -> CMatrix[type]:
    var a = ComplexScalar[type](cos(t / 2), 0)
    var b = ComplexScalar[type](sin(t / 2), 0)
    return CMatrix[type](2, 2,
        a, -b,
        b, a,
    )
fn RY[type: DType](q: Qubit, theta: Scalar[type]) raises -> Gate[type]:
    '''Create an RY rotation gate with angle `theta` applied to qubit `q`.'''
    return Gate[type]('RY', ry[type](theta), List[Qubit, True](q), List[Scalar[type], True](theta))

fn rz[type: DType](t: Scalar[type]) raises -> CMatrix[type]:
    var c: Scalar[type] = cos(t / 2)
    var s: Scalar[type] = sin(t / 2)
    return CMatrix[type](2, 2,
        ComplexScalar[type](c, -s), 0, 
        0, ComplexScalar[type](c, s),
    )
fn RZ[type: DType](q: Qubit, theta: Scalar[type]) raises -> Gate[type]:
    '''Create an RZ rotation gate with angle `theta` applied to qubit `q`.'''
    return Gate[type]('RZ', rz[type](theta), List[Qubit, True](q), List[Scalar[type], True](theta))

fn u[type: DType](t: Scalar[type], p: Scalar[type], l: Scalar[type]) raises -> CMatrix[type]:
    var ct: Scalar[type] = cos(t / 2)
    var st: Scalar[type] = sin(t / 2)
    return CMatrix[type](2, 2,
        ComplexScalar[type](ct, 0), 
        ComplexScalar[type](-cos(l) * st, -sin(l) * st), 
        ComplexScalar[type](cos(p) * st, sin(p) * st), 
        ComplexScalar[type](cos(p + l) * ct, sin(p + l) * ct),
    )
fn U[type: DType](q: Qubit, theta: Scalar[type], phi: Scalar[type], lbda: Scalar[type]) raises -> Gate[type]:
    '''Create a U gate applied to qubit `q`.'''
    return Gate[type]('U', u[type](theta, phi, lbda), List[Qubit, True](q), List[Scalar[type], True](theta, phi, lbda))

# Parameterized multi-qubit gates ###############

# TODO: Add more gates
