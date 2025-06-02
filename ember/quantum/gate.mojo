from math import sqrt, cos, sin, pi
from collections import Set
from utils import Variant

from ..cplx import CMatrix, ComplexScalar
from ..config import DEFAULT_TOL, DEFAULT_TYPE


@value
struct Gate[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    Writable, Sized, Stringable, CollectionElement, EqualityComparable
):
    '''A quantum gate.
    
    Parameters:
        type: A type for the gate matrix and parameter data.
        tol: A tolerance for unitarity and closeness checks.
    '''

    var name: String
    '''An identifier for the gate.'''
    var n_qubits: Int
    '''The number of qubits the gate acts on.'''
    # TODO: Once Mojo supports dynamic traits, we can define a `CMatrixLike` trait to represent
    # a CMatrix or a CSRCMatrix, allowing gate matrices to be optionally sparse
    var matrix: CMatrix[Self.type]
    '''The unitary matrix that implements the gate.'''
    var qubits: List[Int, True]
    '''The qubits the gate is applied to.'''
    var controls: List[Int, True]
    '''The qubits the gate is controlled on.'''
    var classical_controls: List[Int, True]
    '''The classical bits the gate is controlled on.'''
    var params: List[Scalar[Self.type], True]
    '''The parameter values applied to the gate.'''

    # Special fields for the measurement gate
    # TODO: We need measurement to be a Gate so it can be added to the circuit and treated
    # like the other operators. One Mojo supports adding fields to traits, we can define
    # a GateLike trait and have measurement be handled separately.
    var _is_measure: Bool
    '''Whether or not the gate is the measurement gate.'''
    var _measure_targs: List[Int, True]
    '''The target qubits for the measurement gate.'''

    fn __init__(
        out self,
        owned name: String,
        owned matrix: CMatrix[Self.type], 
        owned qubits: List[Int, True],
        owned params: List[Scalar[Self.type], True] = List[Scalar[Self.type], True](),
    ) raises:
        '''Initialize a gate with a name, unitary matrix, set of qubits, and parameters.
        
        Args:
            name: A name for the gate.
            matrix: The unitary matrix implementing the gate.
            qubits: The qubits the gate acts on.
            params: Any parameters applied to the gate.
        '''
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
        self.qubits = qubits^
        self.params = params^
        self.controls = List[Int, True]()
        self.classical_controls = List[Int, True]()

        # Placeholders only applicable to the measurement gate, which must be 
        # instantiated with the _measure method
        self._is_measure = False
        self._measure_targs = List[Int, True]()

    @staticmethod
    @always_inline
    fn _measure(owned qubits: List[Int, True], clbits: List[Int, True]) -> Self:
        '''Creates a measurement gate. `qubits` and `clbits` must be lists of the same length.
        
        Args:
            qubits: The qubits to be measured.
            clbits: The classical bits storing the corresponding measurement results.

        Returns:
            A measurement gate.
        '''
        return Self(
            name='Measure',
            n_qubits=len(qubits),
            matrix=CMatrix[Self.type](rows=0, cols=0, fill_zeros=False),
            qubits=qubits^,
            controls=List[Int, True](),
            classical_controls=List[Int, True](),
            params=List[Scalar[Self.type], True](),
            _is_measure=True,
            _measure_targs=clbits,
        )
    
    @no_inline
    fn __str__(self) -> String:
        '''Convert the gate to a string.
        
        Returns:
            A string representation of the gate.
        '''
        if len(self.params) == 0:
            return self.name
        var str_rep: String = self.name + '('
        for p in self.params:
            str_rep += String(p[]) + ', '
        str_rep = str_rep[:-2]
        return str_rep + ')'
    
    @no_inline
    fn __repr__(self) -> String:
        '''Convert the gate to a string.
        
        Returns:
            A string representation of the gate.
        '''
        var str_rep: String = '<Gate ' + self.name + ': '
        if not self.qubits and not self.controls and not self.params:
            return str_rep[:-2] + '>'
        if self._is_measure:
            var meas_str: String = 'qubits=(q' + String(self.qubits[0])
            for qubit in self.qubits[1:]:
                meas_str += ', q' + String(qubit[])
            meas_str += '), clbits=(c' + String(self._measure_targs[0])
            for clbit in self._measure_targs[1:]:
                meas_str += ', c' + String(clbit[])
            return str_rep + meas_str + ')>'
        var qubits_str: String = ''
        if self.qubits:
            qubits_str = 'qubits=(q' + String(self.qubits[0])
            for qubit in self.qubits[1:]:
                qubits_str += ', q' + String(qubit[])
            qubits_str += ')'
        var controls_str: String = ''
        if self.controls:
            controls_str = 'controls=(q' + String(self.controls[0])
            for qubit in self.controls[1:]:
                controls_str += ', q' + String(qubit[])
            controls_str += ')'
        var classical_controls_str: String = ''
        if self.classical_controls:
            classical_controls_str = 'classical controls=(c' + String(self.classical_controls[0])
            for clbit in self.classical_controls[1:]:
                classical_controls_str += ', c' + String(clbit[])
            classical_controls_str += ')'
        var params_str: String = ''
        if self.params:
            params_str = 'params=(' + String(self.params[0])
            for param in self.params[1:]:
                params_str += ', ' + String(param[])
            params_str += ')'
        if qubits_str:
            str_rep += qubits_str if str_rep[-2:] == ': ' else ', ' + qubits_str
        if controls_str:
            str_rep += controls_str if str_rep[-2:] == ': ' else ', ' + controls_str
        if classical_controls_str:
            str_rep += classical_controls_str if str_rep[-2:] == ': ' else ', ' + classical_controls_str
        if params_str:
            str_rep += params_str if str_rep[-2:] == ': ' else ', ' + params_str
        return str_rep + '>'
        
    @no_inline
    fn write_to[W: Writer](self, mut writer: W):
        '''Write the gate to a writer.

        Args:
            writer: The writer to write to.
        '''
        writer.write(String(self))
    
    @always_inline
    fn __len__(self) -> Int:
        '''The size of the circuit.

        Returns:
            The number of qubits in the circuit.
        '''
        return self.n_qubits

    fn control(mut self, *qubits: Int) raises:
        '''Control the gate on a qubit or set of qubits.
        
        Args:
            qubits: The qubits to set as controls.
        '''
        if self._is_measure:
            raise Error('Cannot control a measurement gate.')
        for q in qubits:
            if q in self.qubits:
                raise Error(
                    'Cannot control gate ' + self.name + ' on qubit ' 
                    + String(q) + ' on which it is applied'
                )
            if q not in self.controls:
                self.controls.append(q)
                self.n_qubits += 1
    
    fn control(
        mut self, 
        qubits: List[Int, True] = List[Int, True](), 
        clbits: List[Int, True] = List[Int, True](),
    ) raises:
        '''Control the gate on specified qubits or classical bits.
        
        Args:
            qubits: The qubits to set as controls.
            clbits: The classical bits to set as controls.
        '''
        if self._is_measure:
            raise Error('Cannot control a measurement gate.')
        for q in qubits:
            if q[] in self.qubits:
                raise Error(
                    'Cannot control gate ' + self.name + ' on qubit ' 
                    + String(q) + ' on which it is applied'
                )
            if q[] not in self.controls:
                self.controls.append(q[])
                self.n_qubits += 1
        for c in clbits:
            if c[] not in self.classical_controls:
                self.classical_controls.append(c[])
    
    @always_inline
    fn __eq__(self, other: Self) -> Bool:
        '''Check two gates for equality. Gates are considered equal if all of their fields are the
        same.
        
        The complex and floating point values in the gate matrix and list of parameters are 
        considered equal if they are within tol of eachother.
        
        Args:
            other: The gate to check for equality.
        
        Returns:
            True if the gates are equal, False otherwise.
        '''
        if (
            self.name != other.name 
            or self.n_qubits != other.n_qubits
            or self.qubits != other.qubits
            or Set(self.controls) != Set(other.controls)
            or Set(self.classical_controls) != Set(other.classical_controls)
            # or not self.matrix.matrix_equals(other.matrix)
            or not self.matrix.is_close[Self.tol](other.matrix)
            or self._is_measure != other._is_measure
            or self._measure_targs != other._measure_targs
        ):
            return False
        if len(self.params) != len(other.params):
            return False
        for i in range(len(self.params)):
            if abs(self.params[i] - other.params[i]) >= Self.tol:
                return False
        return True
    
    @always_inline
    fn __ne__(self, other: Self) -> Bool:
        '''Check two gates for inequality. Gates are considered inequal if any of their fields
        differ.

        The complex and floating point values in the gate matrix and list of parameters are 
        considered inequal if they differ by at least tol.
        
        Args:
            other: The gate to check for inequality.
        
        Returns:
            True if the gates are inequal, False otherwise.
        '''
        return not self.__eq__(other)


# Measurement


fn Measure[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubits: List[Int, True], clbits: List[Int, True]
) raises -> Gate[type, tol]:
    '''Create a measurement on a set of qubits.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubits: The qubits to measure.
        clbits: The classical bits to store the measurements.
    
    Returns:
        A gate representing a measurement.
    '''
    if len(qubits) == 0:
        raise Error('Measure received no qubit specifiers')
    if len(qubits) != len(Set(qubits)):
        raise Error('List of qubits must be unique')
    if len(clbits) != len(Set(clbits)):
        raise Error('List of classical bits must be unique')
    if len(qubits) != len(clbits):
        raise Error('qubits and clbits must be lists of the same length')
    return Gate[type, tol]._measure(qubits, clbits)


fn Measure[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int, clbit: Int
) raises -> Gate[type, tol]:
    '''Create a measurement on a qubit.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit to measure.
        clbit: The classical bit to store the measurements.
    
    Returns:
        A gate representing a measurement.
    '''
    return Gate[type, tol]._measure(List[Int, True](qubit), List[Int, True](clbit))


# Unparameterized single-qubit gates ############


@always_inline
fn x[type: DType = DEFAULT_TYPE]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        0, 1, 
        1, 0,
    )
fn X[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int
) raises -> Gate[type, tol]:
    '''Create an Pauli X gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit on which the gate applies.
    
    Returns:
        An X gate.
    '''
    return Gate[type, tol]('X', x[type](), List[Int, True](qubit))


@always_inline
fn y[type: DType = DEFAULT_TYPE]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        0, ComplexScalar[type](0, -1), 
        ComplexScalar[type](0, 1), 0,
    )
fn Y[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int
) raises -> Gate[type, tol]:
    '''Create a Pauli Y gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit on which the gate applies.
    
    Returns:
        A Y gate.
    '''
    return Gate[type, tol]('Y', y[type](), List[Int, True](qubit))


@always_inline
fn z[type: DType = DEFAULT_TYPE]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        1, 0, 
        0, -1,
    )
fn Z[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int
) raises -> Gate[type, tol]:
    '''Create a Pauli Z gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit on which the gate applies.
    
    Returns:
        A Z gate.
    '''
    return Gate[type, tol]('Z', z[type](), List[Int, True](qubit))


@always_inline
fn h[type: DType = DEFAULT_TYPE]() raises -> CMatrix[type]:
    var inv_sqrt2: Scalar[type] = 1 / sqrt(2.0).cast[type]()
    return CMatrix[type](2, 2,
        inv_sqrt2, inv_sqrt2, 
        inv_sqrt2, -inv_sqrt2,
    )
fn H[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int
) raises -> Gate[type, tol]:
    '''Create a Hadamard gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit on which the gate applies.
    
    Returns:
        A Hadamard gate.
    '''
    return Gate[type, tol]('H', h[type](), List[Int, True](qubit))     


@always_inline
fn s[type: DType = DEFAULT_TYPE]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        1, 0, 
        0, ComplexScalar[type](0, 1),
    )
fn S[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int
) raises -> Gate[type, tol]:
    '''Create an S gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit on which the gate applies.
    
    Returns:
        An S gate.
    '''
    return Gate[type, tol]('S', s[type](), List[Int, True](qubit))


@always_inline
fn t[type: DType = DEFAULT_TYPE]() raises -> CMatrix[type]:
    return CMatrix[type](2, 2,
        1, 0, 
        0, ComplexScalar[type](0, pi / 4).exp(),
    )
fn T[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int
) raises -> Gate[type, tol]:
    '''Create a T gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit on which the gate applies.
    
    Returns:
        A T gate.
    '''
    return Gate[type, tol]('T', t[type](), List[Int, True](qubit))


@always_inline
fn sx[type: DType = DEFAULT_TYPE]() raises -> CMatrix[type]:
    var opi = ComplexScalar[type](0.5, 0.5)
    var omi: ComplexScalar[type] = opi.conj()
    return CMatrix[type](2, 2, 
        opi, omi,
        omi, opi,
    )
fn SX[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int
) raises -> Gate[type, tol]:
    '''Create an SX (square root fo X) gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit on which the gate applies.
    
    Returns:
        An SX gate.
    '''
    return Gate[type, tol]('SX', sx[type](), List[Int, True](qubit))


# Unparameterized multi-qubit gates #############


fn I[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int
) raises -> Gate[type, tol]:
    '''Create an identity gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit on which the gate applies.
    
    Returns:
        An identity gate.
    '''
    return Gate[type, tol]('I', CMatrix[type].eye(2, 2), List[Int, True](qubit))


fn I[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubits: List[Int, True]
) raises -> Gate[type, tol]:
    '''Create an identity gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubits: The qubit on which the gate applies.
    
    Returns:
        An identity gate.
    '''
    var dim: Int = 2 ** len(qubits)
    return Gate[type, tol]('I', CMatrix[type].eye(dim, dim), qubits)


fn CX[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    control: Int, target: Int
) raises -> Gate[type, tol]:
    '''Create a controlled X gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        control: The qubit on which the gate is controlled.
        target: The qubit on which the gate is applied.
    
    Returns:
        A controlled X gate.
    '''
    var cx = Gate[type, tol]('X', x[type](), List[Int, True](target))
    cx.control(control)
    return cx


fn CCX[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    control1: Int, control2: Int, target: Int
) raises -> Gate[type, tol]:
    '''Create a controlled controlled X (Toffoli) gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        control1: The first control qubit.
        control2: The second control qubit.
        target: The qubit on which the gate is applied.
    
    Returns:
        A controlled controlled X (Toffoli) gate.
    '''
    var ccx = Gate[type, tol]('X', x[type](), List[Int, True](target))
    ccx.control(control1, control2)
    return ccx


@always_inline
fn ecr[type: DType = DEFAULT_TYPE]() raises -> CMatrix[type]:
    alias i = ComplexScalar[type].i()
    return CMatrix[type](4, 4,
        0,  1,  0,  i,
        1,  0, -i,  0,
        0,  i,  0,  1,
        -i, 0,  1,  0,
    ) / sqrt(2.0).cast[type]()
fn ECR[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit0: Int, qubit1: Int
) raises -> Gate[type, tol]:
    '''Create an ECR (echoed cross-resonance) gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit0: The first qubit on which the gate applies.
        qubit1: The second qubit on which the gate applies.
    
    Returns:
        An ECR gate.
    '''
    return Gate[type, tol]('ECR', ecr[type](), List[Int, True](qubit0, qubit1))


@always_inline
fn swap[type: DType = DEFAULT_TYPE]() raises -> CMatrix[type]:
    return CMatrix[type](4, 4,
        1, 0, 0, 0,
        0, 0, 1, 0,
        0, 1, 0, 0,
        0, 0, 0, 1,
    )
fn SWAP[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit0: Int, qubit1: Int
) raises -> Gate[type, tol]:
    '''Create a swap gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit0: The first qubit on which the gate applies.
        qubit1: The second qubit on which the gate applies.
    
    Returns:
        A swap gate.
    '''
    return Gate[type, tol]('SWAP', swap[type](), List[Int, True](qubit0, qubit1))


# Parameterized single-qubit gates ##############


@always_inline
fn rx[type: DType = DEFAULT_TYPE](t: Scalar[type]) raises -> CMatrix[type]:
    var a = ComplexScalar[type](cos(t / 2.0), 0)
    var b = ComplexScalar[type](0, -sin(t / 2.0))
    return CMatrix[type](2, 2,
        a, b,
        b, a,
    )
fn RX[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int, theta: Scalar[type]
) raises -> Gate[type, tol]:
    '''Create an RX rotation gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit on which the gate applies.
        theta: A rotation angle for the gate.
    
    Returns:
        An RX rotation gate.
    '''
    return Gate[type, tol](
        'RX', rx[type](theta), List[Int, True](qubit), List[Scalar[type], True](theta)
    )


@always_inline
fn ry[type: DType = DEFAULT_TYPE](t: Scalar[type]) raises -> CMatrix[type]:
    var a = ComplexScalar[type](cos(t / 2.0), 0)
    var b = ComplexScalar[type](sin(t / 2.0), 0)
    return CMatrix[type](2, 2,
        a, -b,
        b, a,
    )
fn RY[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int, theta: Scalar[type]
) raises -> Gate[type, tol]:
    '''Create an RY rotation gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit on which the gate applies.
        theta: A rotation angle for the gate.
    
    Returns:
        An RY rotation gate.
    '''
    return Gate[type, tol](
        'RY', ry[type](theta), List[Int, True](qubit), List[Scalar[type], True](theta)
    )


@always_inline
fn rz[type: DType = DEFAULT_TYPE](t: Scalar[type]) raises -> CMatrix[type]:
    var c: Scalar[type] = cos(t / 2.0)
    var s: Scalar[type] = sin(t / 2.0)
    return CMatrix[type](2, 2,
        ComplexScalar[type](c, -s), 0, 
        0, ComplexScalar[type](c, s),
    )
fn RZ[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int, theta: Scalar[type]
) raises -> Gate[type, tol]:
    '''Create an RZ rotation gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit on which the gate applies.
        theta: A rotation angle for the gate.
    
    Returns:
        An RZ rotation gate.
    '''
    return Gate[type, tol](
        'RZ', rz[type](theta), List[Int, True](qubit), List[Scalar[type], True](theta)
    )


@always_inline
fn u[type: DType = DEFAULT_TYPE](
    t: Scalar[type], p: Scalar[type], l: Scalar[type]
) raises -> CMatrix[type]:
    var ct: Scalar[type] = cos(t / 2.0)
    var st: Scalar[type] = sin(t / 2.0)
    return CMatrix[type](2, 2,
        ComplexScalar[type](ct, 0), 
        ComplexScalar[type](-cos(l) * st, -sin(l) * st), 
        ComplexScalar[type](cos(p) * st, sin(p) * st), 
        ComplexScalar[type](cos(p + l) * ct, sin(p + l) * ct),
    )
fn U[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int, theta: Scalar[type], phi: Scalar[type], lbda: Scalar[type]
) raises -> Gate[type, tol]:
    '''Create a U gate applied.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit on which the gate applies.
        theta: The first Euler angle.
        phi: The second Euler angle.
        lbda: The third Euler angle.

    Returns:
        A U gate.
    '''
    return Gate[type, tol](
        'U', 
        u[type](theta, phi, lbda), 
        List[Int, True](qubit), 
        List[Scalar[type], True](theta, phi, lbda),
    )


@always_inline
fn phase[type: DType = DEFAULT_TYPE](t: Scalar[type]) raises -> CMatrix[type]:
    var exp_it = ComplexScalar[type](cos(t), sin(t))
    return CMatrix[type](2, 2,
        1, 0,
        0, exp_it,
    )
fn PHASE[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int, theta: Scalar[type]
) raises -> Gate[type, tol]:
    '''Create a phase gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit on which the gate applies.
        theta: The parameter value for the gate.
    
    Returns:
        A phase gate.
    '''
    return Gate[type, tol](
        'PHASE', phase[type](theta), List[Int, True](qubit), List[Scalar[type], True](theta)
    )


@always_inline
fn r[type: DType = DEFAULT_TYPE](t: Scalar[type], p: Scalar[type]) raises -> CMatrix[type]:
    var ct: Scalar[type] = cos(t / 2.0)
    var st: Scalar[type] = sin(t / 2.0)
    return CMatrix[type](2, 2,
        ct, -ComplexScalar[type].i() * (-ComplexScalar[type](0, p)).exp() * st,
        -ComplexScalar[type].i() * (ComplexScalar[type](0, p)).exp() * st, ct,
    )
fn R[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit: Int, theta: Scalar[type], phi: Scalar[type]
) raises -> Gate[type, tol]:
    '''Create an R gate (a rotation by angle θ around the cos(φ)x + sin(φ)y axis).
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit: The qubit on which the gate applies.
        theta: A rotation angle for the gate.
        phi: The angle defining the rotation axis.
    
    Returns:
        An R gate.
    '''
    return Gate[type, tol](
        'R', r[type](theta, phi), List[Int, True](qubit), List[Scalar[type], True](theta, phi)
    )


# Parameterized multi-qubit gates ###############


@always_inline
fn rxx[type: DType = DEFAULT_TYPE](t: Scalar[type]) raises -> CMatrix[type]:
    var ct: Scalar[type] = cos(t / 2.0)
    var mist = ComplexScalar[type](0, -sin(t / 2.0)) 
    return CMatrix[type](4, 4,
        ct, 0, 0, mist,
        0, ct, mist, 0,
        0, mist, ct, 0,
        mist, 0, 0, ct,
    )
fn RXX[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit0: Int, qubit1: Int, theta: Scalar[type]
) raises -> Gate[type, tol]:
        '''Create an RXX rotation gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit0: The first qubit on which the gate applies.
        qubit1: The second qubit on which the gate applies.
        theta: A rotation angle for the gate.
    
    Returns:
        An RXX rotation gate.
    '''
    return Gate[type, tol](
        'RXX', rxx[type](theta), List[Int, True](qubit0, qubit1), List[Scalar[type], True](theta)
    )


@always_inline
fn ryy[type: DType = DEFAULT_TYPE](t: Scalar[type]) raises -> CMatrix[type]:
    var ct: Scalar[type] = cos(t / 2.0)
    var ist = ComplexScalar[type](0, sin(t / 2.0)) 
    return CMatrix[type](4, 4,
        ct, 0, 0, ist,
        0, ct, -ist, 0,
        0, -ist, ct, 0,
        ist, 0, 0, ct,
    )
fn RYY[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit0: Int, qubit1: Int, theta: Scalar[type]
) raises -> Gate[type, tol]:
        '''Create an RYY rotation gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit0: The first qubit on which the gate applies.
        qubit1: The second qubit on which the gate applies.
        theta: A rotation angle for the gate.
    
    Returns:
        An RYY rotation gate.
    '''
    return Gate[type, tol](
        'RYY', ryy[type](theta), List[Int, True](qubit0, qubit1), List[Scalar[type], True](theta)
    )


@always_inline
fn rzz[type: DType = DEFAULT_TYPE](t: Scalar[type]) raises -> CMatrix[type]:
    var expit = ComplexScalar[type](cos(t / 2.0), sin(t / 2.0))
    return CMatrix[type](4, 4,
        expit.conj(), 0, 0, 0,
        0, expit, 0, 0,
        0, 0, expit, 0,
        0, 0, 0, expit.conj(),
    )
fn RZZ[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit0: Int, qubit1: Int, theta: Scalar[type]
) raises -> Gate[type, tol]:
        '''Create an RZZ rotation gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit0: The first qubit on which the gate applies.
        qubit1: The second qubit on which the gate applies.
        theta: A rotation angle for the gate.
    
    Returns:
        An RZZ rotation gate.
    '''
    return Gate[type, tol](
        'RZZ', rzz[type](theta), List[Int, True](qubit0, qubit1), List[Scalar[type], True](theta)
    )


@always_inline
fn rzx[type: DType = DEFAULT_TYPE](t: Scalar[type]) raises -> CMatrix[type]:
    var ct: Scalar[type] = cos(t / 2.0)
    var ist = ComplexScalar[type](0, sin(t / 2.0)) 
    return CMatrix[type](4, 4,
        ct, 0, -ist, 0,
        0, ct, 0, ist,
        -ist, 0, ct, 0,
        0, ist, 0, ct,
    )
fn RZX[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit0: Int, qubit1: Int, theta: Scalar[type]
) raises -> Gate[type, tol]:
        '''Create an RZX rotation gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit0: The first qubit on which the gate applies.
        qubit1: The second qubit on which the gate applies.
        theta: A rotation angle for the gate.
    
    Returns:
        An RZX rotation gate.
    '''
    return Gate[type, tol](
        'RZX', rzx[type](theta), List[Int, True](qubit0, qubit1), List[Scalar[type], True](theta)
    )

@always_inline
fn xxminusyy[type: DType = DEFAULT_TYPE](t: Scalar[type], b: Scalar[type]) raises -> CMatrix[type]:
    var ct: Scalar[type] = cos(t / 2.0)
    var mist = ComplexScalar[type](0, -sin(t / 2.0))
    var expib = ComplexScalar[type](cos(b), sin(b))
    return CMatrix[type](4, 4,
        ct, 0, 0, mist * expib.conj(),
        0, 1, 0, 0,
        0, 0, 1, 0,
        mist * expib, 0, 0, ct,
    )
fn XXMinusYY[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit0: Int, qubit1: Int, theta: Scalar[type], beta: Scalar[type] = 0.0,
) raises -> Gate[type, tol]:
    '''Create an XX-YY gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit0: The first qubit on which the gate applies.
        qubit1: The second qubit on which the gate applies.
        theta: A rotation angle for the gate.
        beta: A phase angle for the gate.
    
    Returns:
        An XX-YY gate.
    '''
    return Gate[type, tol](
        'XXMinusYY',
        xxminusyy[type](theta, beta),
        List[Int, True](qubit0, qubit1),
        List[Scalar[type], True](theta, beta),
    )


@always_inline
fn xxplusyy[type: DType = DEFAULT_TYPE](t: Scalar[type], b: Scalar[type]) raises -> CMatrix[type]:
    var ct: Scalar[type] = cos(t / 2.0)
    var mist = ComplexScalar[type](0, -sin(t / 2.0))
    var expib = ComplexScalar[type](cos(b), sin(b))
    return CMatrix[type](4, 4,
        1, 0, 0, 0,
        0, ct, mist * expib.conj(), 0,
        0, mist * expib, ct, 0,
        0, 0, 0, 1,
    )
fn XXPlusYY[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    qubit0: Int, qubit1: Int, theta: Scalar[type], beta: Scalar[type] = 0.0,
) raises -> Gate[type, tol]:
    '''Create an XX+YY gate.
    
    Parameters:
        type: A type for the gate data.
        tol: A tolerance for unitarity and closeness checks.
    
    Args:
        qubit0: The first qubit on which the gate applies.
        qubit1: The second qubit on which the gate applies.
        theta: A rotation angle for the gate.
        beta: A phase angle for the gate.
    
    Returns:
        An XX+YY gate.
    '''
    return Gate[type, tol](
        'XXPlusYY',
        xxplusyy[type](theta, beta),
        List[Int, True](qubit0, qubit1),
        List[Scalar[type], True](theta, beta),
    )
