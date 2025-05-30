from math import sqrt
from collections import Set, Dict

from ..cplx import CMatrix, CSRCMatrix, kron
from .gate import Gate
from .statevector import Statevector
from .densitymatrix import DensityMatrix
from ..config import DEFAULT_TOL, DEFAULT_TYPE


@value
struct QuantumCircuit[type: DType = DEFAULT_TYPE, tol: Scalar[type] = DEFAULT_TOL](
    Stringable, Writable, Movable
):
    '''A quantum circuit.

    Parameters:
        type: A type for the circuit data.
        tol: A tolerance for unitarity and closeness checks.
    '''

    var n_qubits: Int
    '''The number of qubits in the quantum circuit.'''
    var n_clbits: Int
    '''The number of classical bits in the quantum circuit.'''
    var clbits: List[Int, True]
    '''The initial values of the classical bits in the circuit.'''
    var _data: List[Gate[Self.type, Self.tol]]
    '''The gates applied to the qubits in the quantum circuit.'''
    var _initial_state: Statevector[Self.type, Self.tol]
    '''The initial state of the quantum circuit. Defaults to |0⟩ on all qubits.'''

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
        self.clbits = List[Int, True](length=n_clbits, fill=0)
        self._data = List[Gate[Self.type, Self.tol]]()
        self._initial_state = Statevector[Self.type, Self.tol]()
    
    fn apply(mut self, owned *gates: Gate[Self.type, Self.tol]) raises:
        '''Apply a gate to the quantum circuit.
        
        Args:
            gates: The gates to apply.
        '''
        for gate in gates:
            self.apply(gate[])

    fn apply(mut self, owned gate: Gate[Self.type, Self.tol]) raises:
        '''Apply a gate to the quantum circuit.
        
        Args:
            gate: The gate to apply.
        '''
        fn raise_invalid_bit(bit: Int, quantum: Bool) raises:
            var n_bits: Int = self.n_qubits if quantum else self.n_clbits
            var bit_type: String = 'qubit' if quantum else 'classical bit'
            raise Error(
                'Gate ' + gate.name + ' contains invalid ' + bit_type + ' specifier: ' 
                + String(bit) + '. Quantum circuit has ' + String(n_bits) + ' ' + bit_type 
                + 's indexed 0 to ' + String(n_bits - 1) + '.'
            )
        # Make sure all qubits the gate is applied to are in the circuit
        for q in gate.qubits:
            if q[] < 0 or q[] >= self.n_qubits:
                raise_invalid_bit(q[], True)
        if gate._is_measure:
            # Make sure all qubits being measured are in the circuit
            for c in gate._measure_targs:
                if c[] < 0 or c[] >= self.n_clbits:
                    raise_invalid_bit(c[], False)
        else:
            # Make sure all qubits the gate is controlled on are in the circuit
            for q in gate.controls:
                if q[] < 0 or q[] >= self.n_qubits:
                    raise_invalid_bit(q[], True)
            # Make sure all classical bits the gate is controlled on are in the circuit
            for c in gate.classical_controls:
                if c[] < 0 or c[] >= self.n_clbits:
                    raise_invalid_bit(c[], False)
        self._data.append(gate)
    
    fn set_clbits(mut self, bit_values: Dict[Int, Int]) raises:
        '''Set the classical bits to the values specified in bit_values.

        Args:
            bit_values: A mapping of classical bit specifiers to bit values.
        '''
        # Check that the bits are in the circuit and that all values are binary
        for kv in bit_values.items():
            var bit: Int = kv[].key
            var value: Int = kv[].value
            if bit < 0 or bit >= self.n_clbits:
                raise Error('Cannot set bit ' + String(bit) + ' as it is not in the circuit.')
            if value not in Tuple[Int, Int](0, 1):
                raise Error(
                    'Bits can only be set to 0 or 1. Received ' + String(value) 
                    + ' for bit ' + String(bit) + '.'
                )

        for kv in bit_values.items():
            self.clbits[kv[].key] = kv[].value
    
    fn set_clbits(mut self, owned bit_values: List[Int, True]) raises:
        '''Set the classical bits to the values specified in bit_values.

        Args:
            bit_values: A list of bit values corresponding to the classical bits. Must have length
                self.n_clbits.
        '''
        if len(bit_values) != self.n_clbits:
            raise Error('bit_values must have the same length as the number of classical bits')
        for bit in bit_values:
            if bit[] not in Tuple[Int, Int](0, 1):
                raise Error('Bits can only be set to 0 or 1. Received ' + String(bit) + '.')
        self.clbits = bit_values^
    
    @always_inline
    fn set_initial_state(
        mut self, owned statevector: CSRCMatrix[Self.type], normalize: Bool = False
    ) raises:
        '''Set the initial circuit state to the provided statevector.

        Args:
            statevector: The statevector to use as the initial state. Should be a sparse matrix 
                with shape (1, 2^n_qubits) or (2^n_qubits, 1).
            normalize: If True, the statevector will be automatically normalized. If False and the
                statevector is unnormalized, an error will be raised.
        '''
        self._initial_state = Statevector[Self.type, Self.tol](
            statevector=statevector^, normalize=normalize, enforce_n_elements=2 ** self.n_qubits
        )
    
    @always_inline
    fn set_initial_state(
        mut self, owned statevector: CMatrix[Self.type], normalize: Bool=False
    ) raises:
        '''Set the initial circuit state to the provided statevector.

        Args:
            statevector: The statevector to use as the initial state. Should be a matrix 
                with shape (1, 2^n_qubits) or (2^n_qubits, 1).
            normalize: If True, the statevector will be automatically normalized. If False and the
                statevector is unnormalized, an error will be raised.
        '''
        self._initial_state = Statevector[Self.type, Self.tol](
            statevector=statevector^, normalize=normalize, enforce_n_elements=2 ** self.n_qubits
        )
    
    @always_inline
    fn set_initial_state(
        mut self, owned statevector: List[ComplexScalar[Self.type], True], normalize: Bool = False
    ) raises:
        '''Set the initial circuit state to the provided statevector.

        Args:
            statevector: The statevector to use as the initial state. Should be a list with 
                2^n_qubits values.
            normalize: If True, the statevector will be automatically normalized. If False and the
                statevector is unnormalized, an error will be raised.
        '''
        self._initial_state = Statevector[Self.type, Self.tol](
            statevector=statevector^, normalize=normalize, enforce_n_elements=2 ** self.n_qubits
        )
    
    @always_inline
    fn set_initial_state(
        mut self, owned statevector: Statevector[Self.type, Self.tol], normalize: Bool = False
    ) raises:
        '''Set the initial circuit state to the provided statevector.

        Args:
            statevector: The statevector to use as the initial state. Should be a list with 
                2^n_qubits values.
            normalize: If True, the statevector will be automatically normalized. If False and the
                statevector is unnormalized, an error will be raised.
        '''
        self._initial_state = statevector^
    
    @always_inline
    fn set_initial_state(
        mut self,
        owned statevector: Dict[Int, ComplexScalar[Self.type]],
        normalize: Bool = False
    ) raises:
        '''Set the initial circuit state to the provided statevector.
        
        Args:
            statevector: A dictionary of statevector index element pairs.
            normalize: If True, the statevector will be automatically normalized. If False and the
                statevector is unnormalized, an error will be raised.
        '''
        self._initial_state = Statevector[Self.type, Self.tol](
            statevector=statevector^, normalize=normalize, n_elements=2 ** self.n_qubits
        )
    
    fn join(
        mut self,
        owned other: QuantumCircuit[Self.type, Self.tol],
        *qubits: Int,
    ) raises:
        '''Apply the instructions from one circuit onto the qubits of self.

        Args:
            other: The circuit to join with self.
            qubits: The qubits in self to join onto.
        '''
        var self_qubits = List[Int, True](capacity=len(qubits))
        for q in qubits:
            self_qubits.append(q)
        self.join(other, self_qubits)

    fn join(
        mut self,
        owned other: QuantumCircuit[Self.type, Self.tol],
        qubits: List[Int, True],
        clbits: List[Int, True] = List[Int, True]()
    ) raises:
        '''Apply the instructions from one circuit onto the qubits/clbits of self.

        Args:
            other: The circuit to join with self.
            qubits: The qubits in self to join onto. Must be a list of length other.n_qubits.
            clbits: The classical bits in self to join onto. Must be a list of length 
                other.n_clbits.
        '''
        if self.n_qubits < other.n_qubits:
            raise Error(
                'Cannot join a circuit with ' + String(other.n_qubits) 
                + ' qubits onto a circuit with ' + String(self.n_qubits) + ' qubits.'
            )
        if self.n_clbits < other.n_clbits:
            raise Error(
                'Cannot join a circuit with ' + String(other.n_clbits) 
                + ' classical bits onto a circuit with ' + String(self.n_clbits) 
                + ' classical bits.'
            )
        if len(qubits) != len(Set(qubits)):
            raise Error('List of qubits is not unique.')
        if len(clbits) != len(Set(clbits)):
            raise Error('List of classical bits is not unique.')
        if len(qubits) != other.n_qubits:
            raise Error(
                'Number of qubits specified does not match the number of qubits in the circuit.'
            )
        if len(clbits) != other.n_clbits:
            raise Error(
                'Number of classical bits specified does not match the number of classical bits'
                ' in the circuit.'
            )
        for q in qubits:
            if q[] < 0 or q[] >= self.n_qubits:
                raise Error('Invalid qubit specifier: ' + String(q[]) + '.')
        for c in clbits:
            if c[] < 0 or c[] >= self.n_clbits:
                raise Error('Invalid classical bit specifier: ' + String(c[]) + '.')
        if other._initial_state.size != 0:
            if self._initial_state.size != 0:
                raise Error(
                    'Cannot join a circuit with qubit initialization onto another circuit '
                    'with qubit initialization.'
                )
            # If any qubits `other` is applied to in self have gates, raise an error
            for q in qubits:
                for gate in self._data:
                    if q[] in gate[].qubits or q[] in gate[].controls:
                        raise Error('Cannot apply mid-circuit qubit initialization at this time.')
            
            var new_initial_state = Statevector[Self.type, Self.tol]._empty(self.n_qubits)
            for other_sv_idx in range(2 ** other.n_qubits):
                var other_sv_idx_arr: List[Int, True] = self._int_to_bin_array(other_sv_idx, other.n_qubits)
                
                var self_sv_idx_arr = List[Int, True](length=self.n_qubits, fill=0)
                for qubit_idx in range(other.n_qubits):
                    var new_idx: Int = qubits[qubit_idx]
                    self_sv_idx_arr[new_idx] = other_sv_idx_arr[qubit_idx]
                
                var self_sv_idx: Int = self._bin_array_to_int(self_sv_idx_arr)
                new_initial_state[self_sv_idx] = other._initial_state[other_sv_idx]
            self._initial_state = new_initial_state^
        
        # Map classical bit values
        for i in range(other.n_clbits):
            self.clbits[clbits[i]] = other.clbits[i]

        fn get_new_bits(bits: List[Int, True], quantum: Bool) raises -> List[Int, True]:
            var new_bits = List[Int, True]()
            for b in bits:
                new_bits.append((qubits if quantum else clbits)[b[]])
            return new_bits
        
        for gate_ref in other._data:
            var gate: Gate[Self.type, Self.tol] = gate_ref[]
            var new_qubits: List[Int, True] = get_new_bits(gate.qubits, True)
            var new_controls: List[Int, True] = get_new_bits(gate.controls, True)
            var new_cl_controls: List[Int, True] = get_new_bits(gate.classical_controls, False)
            var new_measure_targs: List[Int, True] = get_new_bits(gate._measure_targs, False)
            var new_gate = Gate[Self.type, Self.tol](
                name=gate.name,
                n_qubits=gate.n_qubits,
                matrix=gate.matrix,
                qubits=new_qubits^,
                controls=new_controls^,
                classical_controls=new_cl_controls^,
                params=gate.params,
                _is_measure=gate._is_measure,
                _measure_targs=new_measure_targs,
            )
            self.apply(new_gate^)
    
    @staticmethod
    fn _int_to_bin_array(i: Int, length: Int) -> List[Int, True]:
        '''Convert an integer to a binary array of specified length.'''
        var bin_arr = List[Int, True](capacity=length)
        for j in range(length):
            bin_arr.append((i >> j) & 1)
        return bin_arr
    
    @staticmethod
    fn _bin_array_to_int(arr: List[Int, True]) -> Int:
        '''Convert binary array to integer.'''
        var i: Int = 0
        for j in range(len(arr)):
            i += arr[j] * (2 ** j)
        return i

    @no_inline
    fn __str__(self) -> String:
        '''Convert the quantum circuit to a string.

        Returns:
            A string representation of the quantum circuit.
        '''
        alias max_width: Int = 120

        @parameter
        fn max(l: List[Int, True]) -> Int:
            var mx: Int = l[0]
            for e in l[1:]:
                if e[] > mx:
                    mx = e[]
            return mx
        
        var max_bit_num_len: Int = len(String(max(List[Int, True](self.n_qubits, self.n_clbits))))
        # List of qubit and clbit lines
        # Each qubit line is a list of gate strings
        var lines = List[List[String]]()
        for i in range(self.n_qubits):
            var i_str: String = String(i)
            var spaces: Int = max_bit_num_len - len(i_str) + 1
            lines.append(List[String]('q' + i_str + ':' + ' ' * spaces + '|0⟩ -'))
        for i in range(self.n_clbits):
            var i_str: String = String(i)
            var spaces: Int = max_bit_num_len - len(i_str) + 2
            lines.append(
                List[String]('c' + i_str + ':' + ' ' * spaces + String(self.clbits[i]) + '  -')
            )
        
        if len(self._initial_state) > 0:
            for i in range(self.n_qubits):
                lines[i].append('Initialize-')
            for i in range(self.n_qubits, self.n_qubits + self.n_clbits):
                lines[i].append('-----------')

        for gate in self._data:
            var gate_str = String(gate[])
            for i in range(self.n_qubits + self.n_clbits):
                var is_clbit: Bool = False
                if i >= self.n_qubits:
                    i -= self.n_qubits
                    is_clbit = True
                var mod_gate_str: String
                if gate[]._is_measure:
                    var measure_specifier_len: Int = len(String(max(gate[]._measure_targs))) + 2
                    var measure_to: Int
                    try:
                        measure_to = gate[]._measure_targs[gate[].qubits.index(i)]
                    except:
                        measure_to = 0
                    var n_dashes: Int = measure_specifier_len - len(String(measure_to)) - 1
                    mod_gate_str = gate_str + '-' * n_dashes + '>' + String(measure_to)
                else:
                    mod_gate_str = gate_str
                if not is_clbit and i in gate[].qubits:
                    lines[i].append(mod_gate_str)
                elif not gate[]._is_measure and (
                    (not is_clbit and i in gate[].controls) 
                    or (is_clbit and i in gate[].classical_controls)
                ):
                    var n_dashes: Int = len(mod_gate_str) - 1
                    var left_dashes: Int = n_dashes // 2
                    var right_dashes: Int = n_dashes - left_dashes
                    lines[i + self.n_qubits if is_clbit else i].append(
                        '-' * left_dashes + '*' + '-' * right_dashes
                    )
                else:
                    lines[i + self.n_qubits if is_clbit else i].append('-' * len(mod_gate_str))
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
