
<div align="center">

# Ember

</div>

## Overview

Ember is a quantum computing framework that provides tools for building, simulating, and analyzing quantum circuits. To facilitate efficient statevector computation, Ember additionally contains a comprehensive complex number SIMD vector implementation, efficient complex dense and sparse matrices, and numerical linear algebra processing tools.

This project serves as one of my passion projects for me to learn about the Mojo language.

---

## Features

- **🤖 Feature-Rich Quantum Circuits:** Supporting arbitrary quantum gates, classical bit control, measurement, inversion, and merging.
- **🧪 Comprehensive Quantum Gate Library:** Utilize preconfigured quantum gates or configure custom gates for complex operations.
- **🧮 Statevector Handling:** Efficient quantum state representations and operations through sparse statevectors.
- **🧠 Density Matrix Support:** Accurate representation of mixed states with purity validation. 
- **🚀 Statevector Simulation:** Efficient, parallelizable statevector simulation for large-circuit simulation.
- **🔢 Complex Number SIMD Vectors:** A complex number type supporting vectorized SIMD operations.
- **📚 Complex Matrix Support:** Sparse and dense complex matrices supporting efficient numerical linear algebra.

---

## Project Structure

```sh
└── ember/
    ├── LICENSE
    ├── README.md
    ├── ember
    │   ├── __init__.mojo
    │   ├── config.mojo
    │   ├── cplx
    │   │   ├── __init__.mojo
    │   │   ├── complexsimd.mojo
    │   │   ├── cmatrix.mojo
    │   │   ├── csrcmatrix.mojo
    │   │   ├── mmath.mojo
    │   │   └── qr.mojo
    │   ├── quantum
    │   │   ├── __init__.mojo
    │   │   ├── gate.mojo
    │   │   ├── quantumcircuit.mojo
    │   │   ├── statevector.mojo
    │   │   └── densitymatrix.mojo
    │   └── sim
    │       ├── __init__.mojo
    │       └── statevectorsimulator.mojo
    ├── examples
    │   ├── fourier.mojo
    │   ├── grover.mojo
    │   └── teleportation.mojo
    └── test
```

---

## Getting Started

### Prerequisites

This project requires the following dependencies:

- **Programming Language:** Mojo 25.3.0

### Installation

Install [Mojo](https://docs.modular.com/mojo/manual/get-started/) and clone the repository:

```sh
❯ git clone https://github.com/adamreidsmith/ember
```

---

## Usage

The following example implements [quantum teleportation](https://en.wikipedia.org/wiki/Quantum_teleportation) in Ember:

```mojo
from ember import QuantumCircuit, Statevector, StatevectorSimulator
from ember import X, Z, H, CX, Measure

fn main() raises:
	# Create a quantum circuit with 3 qubits and 2 classical bits
	var qc = QuantumCircuit(3, 2)

	# Apply a Hadamard gate to qubit 0 to prepare the state to teleport (we will teleport the |+⟩ state)
	qc.apply(H(0))

	# Generate the entangled Bell state on qubits 1 and 2
	qc.apply(H(1), CX(1, 2))

	# Change the basis so we can perform a Bell measurement
	qc.apply(CX(0, 1), H(0))

	# Measure qubit 0 to classicao bit 0 and qubit 1 to classical bit 1
	qc.apply(Measure(0, 0), Measure(0, 1))

	# Create gates controlled on the classical bits containing the measurement outcomes
	var x: Gate = X(2)
	x.control(clbits=List[Int, True](1))
	var z: Gate = Z(2)
	z.control(clbits=List[Int, True](0))
	qc.apply(x, z)

	# Initialize a simulator and run the quantum circuit
	sim = StatevectorSimulator()
	sim.run(qc)
	var statevector = sim.get_statevector()

	# Trace out qubits 0 and 1 to obtain the state of the teleported qubit
	statevector = statevector.partial_trace(0, 1)

	print('Teleported state:')
	print(statevector)
```

See [examples](./examples) for more detailed examples of [quantum teleportation](./examples/teleportation.mojo), [Grover's algorithm](./examples/grover.mojo), and the [quantum Fourier transform](./examples/fourier.mojo).

---

## License

Ember is protected under the [MIT](https://choosealicense.com/licenses/mit/) License. For more details, refer to the [LICENSE](./LICENSE) file.

---

## Sources

- Jones, Tyson, Bálint Koczor, and Simon C. Benjamin. "Distributed Simulation of Statevectors and Density Matrices." arXiv.org, November 2, 2023. https://arxiv.org/abs/2311.01512.
