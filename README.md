<div id="top">

<!-- HEADER STYLE: CLASSIC -->
<div align="center">

<img src="readmeai/assets/logos/purple.svg" width="30%" style="position: relative; top: 0; right: 0;" alt="Project Logo"/>

# EMBER.GIT

<em></em>

<!-- BADGES -->
<img src="https://img.shields.io/github/license/adamreidsmith/ember.git?style=default&logo=opensourceinitiative&logoColor=white&color=0080ff" alt="license">
<img src="https://img.shields.io/github/last-commit/adamreidsmith/ember.git?style=default&logo=git&logoColor=white&color=0080ff" alt="last-commit">
<img src="https://img.shields.io/github/languages/top/adamreidsmith/ember.git?style=default&color=0080ff" alt="repo-top-language">
<img src="https://img.shields.io/github/languages/count/adamreidsmith/ember.git?style=default&color=0080ff" alt="repo-language-count">

<!-- default option, no dependency badges. -->


<!-- default option, no dependency badges. -->

</div>
<br>

---

## Table of Contents

- [Table of Contents](#table-of-contents)
- [Overview](#overview)
- [Features](#features)
- [Project Structure](#project-structure)
    - [Project Index](#project-index)
- [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Installation](#installation)
    - [Usage](#usage)
    - [Testing](#testing)
- [Roadmap](#roadmap)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)

---

## Overview



---

## Features

<code>❯ REPLACE-ME</code>

---

## Project Structure

```sh
└── ember.git/
    ├── LICENSE
    ├── README.md
    ├── ember
    │   ├── __init__.mojo
    │   ├── config.mojo
    │   ├── cplx
    │   ├── quantum
    │   └── sim
    ├── examples
    │   ├── fourier.mojo
    │   ├── grover.mojo
    │   └── teleportation.mojo
    └── test
        ├── __init__.mojo
        ├── _testing.mojo
        ├── test_cmatrix.mojo
        ├── test_complexsimd.mojo
        ├── test_csrcmatrix.mojo
        ├── test_densitymatrix.mojo
        ├── test_gate.mojo
        ├── test_mmath.mojo
        ├── test_qr.mojo
        ├── test_quantumcircuit.mojo
        └── test_statevectorsimulator.mojo
```

### Project Index

<details open>
	<summary><b><code>EMBER.GIT/</code></b></summary>
	<!-- __root__ Submodule -->
	<details>
		<summary><b>__root__</b></summary>
		<blockquote>
			<div class='directory-path' style='padding: 8px 0; color: #666;'>
				<code><b>⦿ __root__</b></code>
			<table style='width: 100%; border-collapse: collapse;'>
			<thead>
				<tr style='background-color: #f8f9fa;'>
					<th style='width: 30%; text-align: left; padding: 8px;'>File Name</th>
					<th style='text-align: left; padding: 8px;'>Summary</th>
				</tr>
			</thead>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/LICENSE'>LICENSE</a></b></td>
					<td style='padding: 8px;'>- The LICENSE file specifies the projects open-source licensing terms under the MIT License<br>- It grants users broad permissions to use, modify, and distribute the software, disclaiming any warranty or liability<br>- This ensures legal clarity and facilitates community contribution and wider adoption of the project.</td>
				</tr>
			</table>
		</blockquote>
	</details>
	<!-- test Submodule -->
	<details>
		<summary><b>test</b></summary>
		<blockquote>
			<div class='directory-path' style='padding: 8px 0; color: #666;'>
				<code><b>⦿ test</b></code>
			<table style='width: 100%; border-collapse: collapse;'>
			<thead>
				<tr style='background-color: #f8f9fa;'>
					<th style='width: 30%; text-align: left; padding: 8px;'>File Name</th>
					<th style='text-align: left; padding: 8px;'>Summary</th>
				</tr>
			</thead>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/test/test_qr.mojo'>test_qr.mojo</a></b></td>
					<td style='padding: 8px;'>- Tests for QR decomposition algorithms are implemented<br>- The code verifies the complex Schur decomposition, eigenvalue calculation, and right eigenvector computation functions<br>- Various test matrices, including identity, random, and custom matrices, are used to assess the accuracy and robustness of these numerical linear algebra routines within the Ember library<br>- Test results confirm the correctness of the implemented algorithms.</td>
				</tr>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/test/_testing.mojo'>_testing.mojo</a></b></td>
					<td style='padding: 8px;'>- The <code>_testing.mojo</code> file provides custom assertion functions for testing the <code>CMatrix</code> and <code>ComplexScalar</code> types within the Ember project<br>- It offers functions to compare matrices and scalars for exact equality and approximate equality, handling potential NaN values<br>- These functions enhance the testing framework by providing tailored comparison methods for the projects core data structures.</td>
				</tr>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/test/test_csrcmatrix.mojo'>test_csrcmatrix.mojo</a></b></td>
					<td style='padding: 8px;'>- The file <code>test/test_csrcmatrix.mojo</code> contains unit tests for the <code>CSRCMatrix</code> class within the Ember project<br>- It verifies the correctness of <code>CSRCMatrix</code> functionality, including initialization, property access, arithmetic operations, matrix multiplication, comparisons, and static constructors<br>- These tests ensure the <code>CSRCMatrix</code> class behaves as expected within the larger Ember library, likely contributing to a robust and reliable linear algebra component.</td>
				</tr>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/test/test_cmatrix.mojo'>test_cmatrix.mojo</a></b></td>
					<td style='padding: 8px;'>- The file <code>test/test_cmatrix.mojo</code> contains unit tests for the <code>CMatrix</code> class within the Ember project<br>- It verifies the correctness of <code>CMatrix</code> functionality, including initialization, arithmetic operations, matrix multiplication, comparisons, and various other methods, ensuring the <code>CMatrix</code> class behaves as expected within the broader Ember library<br>- The tests cover a comprehensive range of scenarios to validate the robustness and accuracy of the <code>CMatrix</code> implementation.</td>
				</tr>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/test/test_quantumcircuit.mojo'>test_quantumcircuit.mojo</a></b></td>
					<td style='padding: 8px;'>- Tests comprise the <code>test_quantumcircuit.mojo</code> file, verifying the <code>QuantumCircuit</code> class functionality within the Ember quantum computing library<br>- It validates initialization, classical bit manipulation, gate application, statevector setting, circuit joining, and error handling for invalid operations<br>- The tests ensure correct behavior and data integrity across various <code>QuantumCircuit</code> methods.</td>
				</tr>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/test/test_gate.mojo'>test_gate.mojo</a></b></td>
					<td style='padding: 8px;'>- The file <code>test/test_gate.mojo</code> contains unit tests for the <code>ember</code> librarys quantum gate functionality<br>- It verifies the correct operation of various gates (single-qubit, multi-qubit, parameterized and unparameterized) within the larger <code>ember</code> project, ensuring the accuracy and reliability of the gate implementations<br>- The tests cover initialization, measurement, control operations, and other gate-specific behaviors.</td>
				</tr>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/test/test_statevectorsimulator.mojo'>test_statevectorsimulator.mojo</a></b></td>
					<td style='padding: 8px;'>- The file <code>test/test_statevectorsimulator.mojo</code> contains unit tests for the <code>StatevectorSimulator</code> class within the Ember quantum computing library<br>- It verifies the correctness of the simulators functionality by comparing its output (statevector calculations) against expected results using assertions<br>- This ensures the accuracy and reliability of the <code>StatevectorSimulator</code> component within the broader Ember project, which likely involves a larger quantum circuit simulation framework.</td>
				</tr>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/test/test_complexsimd.mojo'>test_complexsimd.mojo</a></b></td>
					<td style='padding: 8px;'>- Unit tests validate the <code>ComplexSIMD</code> and <code>ComplexScalar</code> classes<br>- The tests encompass initialization, static constructors, arithmetic operations, comparisons, element-wise access, and various mathematical functions like exponentiation, logarithms, and trigonometric functions<br>- Successful execution confirms the correctness and functionality of these classes within the broader Ember project.</td>
				</tr>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/test/test_mmath.mojo'>test_mmath.mojo</a></b></td>
					<td style='padding: 8px;'>- The file <code>test/test_mmath.mojo</code> contains unit tests for the <code>ember</code> librarys matrix math functions<br>- It verifies the correctness of functions like matrix multiplication (kron, sparse_kron), matrix manipulations (swap_rows, swap_cols), linear algebra operations (solve, expm), matrix norms (one_norm), and other specialized functions (mmax, mmin, augmented_ref, matrix_power, hstack, vstack, is_positive_semidefinite)<br>- These tests are crucial for ensuring the reliability and accuracy of the core matrix math capabilities within the larger <code>ember</code> project.</td>
				</tr>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/test/__init__.mojo'>__init__.mojo</a></b></td>
					<td style='padding: 8px;'>- It orchestrates the execution of unit tests for various components within a quantum computing simulation library<br>- These components include matrix operations, quantum gates, circuits, simulators (statevector and density matrix), and related mathematical functions<br>- The file acts as a central entry point for running the complete test suite, ensuring comprehensive validation of the librarys functionality.</td>
				</tr>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/test/test_densitymatrix.mojo'>test_densitymatrix.mojo</a></b></td>
					<td style='padding: 8px;'>- The file <code>test/test_densitymatrix.mojo</code> contains unit tests for the <code>DensityMatrix</code> class within the <code>ember</code> project<br>- It verifies the correctness of the <code>DensityMatrix</code> classs initialization and purity checks, among other functionalities, contributing to the overall quality assurance of the <code>ember</code> library<br>- The tests utilize assertion functions to compare expected and actual results, ensuring the reliability of the density matrix calculations within the larger <code>ember</code> codebase.</td>
				</tr>
			</table>
		</blockquote>
	</details>
	<!-- examples Submodule -->
	<details>
		<summary><b>examples</b></summary>
		<blockquote>
			<div class='directory-path' style='padding: 8px 0; color: #666;'>
				<code><b>⦿ examples</b></code>
			<table style='width: 100%; border-collapse: collapse;'>
			<thead>
				<tr style='background-color: #f8f9fa;'>
					<th style='width: 30%; text-align: left; padding: 8px;'>File Name</th>
					<th style='text-align: left; padding: 8px;'>Summary</th>
				</tr>
			</thead>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/examples/teleportation.mojo'>teleportation.mojo</a></b></td>
					<td style='padding: 8px;'>- The <code>teleportation.mojo</code> example demonstrates quantum teleportation using the Ember quantum computing framework<br>- It simulates Alice sending a quantum state to Bob via classical communication and a pre-shared entangled state<br>- The code constructs a quantum circuit implementing the teleportation protocol, simulates its execution, and verifies the successful transfer of the quantum state, accounting for global phase differences.</td>
				</tr>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/examples/grover.mojo'>grover.mojo</a></b></td>
					<td style='padding: 8px;'>- Grovers algorithm implementation demonstrates quantum search<br>- It constructs and simulates a quantum circuit to locate a specific state within a larger search space, leveraging oracle and diffuser components for amplitude amplification<br>- The algorithm achieves a quadratic speedup compared to classical search methods by iteratively increasing the probability of measuring the target state<br>- The example finds the state |101⟩ in a 3-qubit system.</td>
				</tr>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/examples/fourier.mojo'>fourier.mojo</a></b></td>
					<td style='padding: 8px;'>- The <code>examples/fourier.mojo</code> file provides an implementation of the Quantum Fourier Transform (QFT) algorithm<br>- It constructs a QFT circuit, allowing for approximation to reduce gate count, and applies it to a sample quantum state<br>- The code then verifies the transformations correctness by applying the inverse QFT, demonstrating the algorithms functionality within the Ember quantum computing framework<br>- The example showcases state preparation, QFT application, and statevector simulation for verification.</td>
				</tr>
			</table>
		</blockquote>
	</details>
	<!-- ember Submodule -->
	<details>
		<summary><b>ember</b></summary>
		<blockquote>
			<div class='directory-path' style='padding: 8px 0; color: #666;'>
				<code><b>⦿ ember</b></code>
			<table style='width: 100%; border-collapse: collapse;'>
			<thead>
				<tr style='background-color: #f8f9fa;'>
					<th style='width: 30%; text-align: left; padding: 8px;'>File Name</th>
					<th style='text-align: left; padding: 8px;'>Summary</th>
				</tr>
			</thead>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/config.mojo'>config.mojo</a></b></td>
					<td style='padding: 8px;'>- Embers <code>config.mojo</code> defines project-wide constants<br>- It establishes default data types, tolerances for numerical comparisons, and zero thresholds for sparse data structures<br>- Crucially, it also sets machine epsilon and numeric limits for various floating-point precisions (Float64, Float32, Float16), ensuring consistent and predictable numerical behavior across the entire Ember application.</td>
				</tr>
				<tr style='border-bottom: 1px solid #eee;'>
					<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/__init__.mojo'>__init__.mojo</a></b></td>
					<td style='padding: 8px;'>Code>❯ REPLACE-ME</code></td>
				</tr>
			</table>
			<!-- sim Submodule -->
			<details>
				<summary><b>sim</b></summary>
				<blockquote>
					<div class='directory-path' style='padding: 8px 0; color: #666;'>
						<code><b>⦿ ember.sim</b></code>
					<table style='width: 100%; border-collapse: collapse;'>
					<thead>
						<tr style='background-color: #f8f9fa;'>
							<th style='width: 30%; text-align: left; padding: 8px;'>File Name</th>
							<th style='text-align: left; padding: 8px;'>Summary</th>
						</tr>
					</thead>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/sim/statevectorsimulator.mojo'>statevectorsimulator.mojo</a></b></td>
							<td style='padding: 8px;'>Code>❯ REPLACE-ME</code></td>
						</tr>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/sim/__init__.mojo'>__init__.mojo</a></b></td>
							<td style='padding: 8px;'>- StatevectorSimulator provides a quantum simulation capability within the Ember project<br>- Its a core component of the <code>ember/sim</code> module, offering a crucial simulation engine for the broader application<br>- The modules role is to facilitate quantum computation simulations, enabling the rest of the Ember codebase to perform and analyze quantum algorithms.</td>
						</tr>
					</table>
				</blockquote>
			</details>
			<!-- cplx Submodule -->
			<details>
				<summary><b>cplx</b></summary>
				<blockquote>
					<div class='directory-path' style='padding: 8px 0; color: #666;'>
						<code><b>⦿ ember.cplx</b></code>
					<table style='width: 100%; border-collapse: collapse;'>
					<thead>
						<tr style='background-color: #f8f9fa;'>
							<th style='width: 30%; text-align: left; padding: 8px;'>File Name</th>
							<th style='text-align: left; padding: 8px;'>Summary</th>
						</tr>
					</thead>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/cplx/qr.mojo'>qr.mojo</a></b></td>
							<td style='padding: 8px;'>Code>❯ REPLACE-ME</code></td>
						</tr>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/cplx/mmath.mojo'>mmath.mojo</a></b></td>
							<td style='padding: 8px;'>Code>❯ REPLACE-ME</code></td>
						</tr>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/cplx/complexsimd.mojo'>complexsimd.mojo</a></b></td>
							<td style='padding: 8px;'>Code>❯ REPLACE-ME</code></td>
						</tr>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/cplx/csrcmatrix.mojo'>csrcmatrix.mojo</a></b></td>
							<td style='padding: 8px;'>Code>❯ REPLACE-ME</code></td>
						</tr>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/cplx/__init__.mojo'>__init__.mojo</a></b></td>
							<td style='padding: 8px;'>Code>❯ REPLACE-ME</code></td>
						</tr>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/cplx/cmatrix.mojo'>cmatrix.mojo</a></b></td>
							<td style='padding: 8px;'>Code>❯ REPLACE-ME</code></td>
						</tr>
					</table>
				</blockquote>
			</details>
			<!-- quantum Submodule -->
			<details>
				<summary><b>quantum</b></summary>
				<blockquote>
					<div class='directory-path' style='padding: 8px 0; color: #666;'>
						<code><b>⦿ ember.quantum</b></code>
					<table style='width: 100%; border-collapse: collapse;'>
					<thead>
						<tr style='background-color: #f8f9fa;'>
							<th style='width: 30%; text-align: left; padding: 8px;'>File Name</th>
							<th style='text-align: left; padding: 8px;'>Summary</th>
						</tr>
					</thead>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/quantum/gate.mojo'>gate.mojo</a></b></td>
							<td style='padding: 8px;'>Code>❯ REPLACE-ME</code></td>
						</tr>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/quantum/statevector.mojo'>statevector.mojo</a></b></td>
							<td style='padding: 8px;'>Code>❯ REPLACE-ME</code></td>
						</tr>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/quantum/densitymatrix.mojo'>densitymatrix.mojo</a></b></td>
							<td style='padding: 8px;'>Code>❯ REPLACE-ME</code></td>
						</tr>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/quantum/__init__.mojo'>__init__.mojo</a></b></td>
							<td style='padding: 8px;'>Code>❯ REPLACE-ME</code></td>
						</tr>
						<tr style='border-bottom: 1px solid #eee;'>
							<td style='padding: 8px;'><b><a href='https://github.com/adamreidsmith/ember.git/blob/master/ember/quantum/quantumcircuit.mojo'>quantumcircuit.mojo</a></b></td>
							<td style='padding: 8px;'>Code>❯ REPLACE-ME</code></td>
						</tr>
					</table>
				</blockquote>
			</details>
		</blockquote>
	</details>
</details>

---

## Getting Started

### Prerequisites

This project requires the following dependencies:

- **Programming Language:** unknown

### Installation

Build ember.git from the source and intsall dependencies:

1. **Clone the repository:**

    ```sh
    ❯ git clone https://github.com/adamreidsmith/ember.git
    ```

2. **Navigate to the project directory:**

    ```sh
    ❯ cd ember.git
    ```

3. **Install the dependencies:**

echo 'INSERT-INSTALL-COMMAND-HERE'

### Usage

Run the project with:

echo 'INSERT-RUN-COMMAND-HERE'

### Testing

Ember.git uses the {__test_framework__} test framework. Run the test suite with:

echo 'INSERT-TEST-COMMAND-HERE'

---

## Roadmap

- [X] **`Task 1`**: <strike>Implement feature one.</strike>
- [ ] **`Task 2`**: Implement feature two.
- [ ] **`Task 3`**: Implement feature three.

---

## Contributing

- **💬 [Join the Discussions](https://github.com/adamreidsmith/ember.git/discussions)**: Share your insights, provide feedback, or ask questions.
- **🐛 [Report Issues](https://github.com/adamreidsmith/ember.git/issues)**: Submit bugs found or log feature requests for the `ember.git` project.
- **💡 [Submit Pull Requests](https://github.com/adamreidsmith/ember.git/blob/main/CONTRIBUTING.md)**: Review open PRs, and submit your own PRs.

<details closed>
<summary>Contributing Guidelines</summary>

1. **Fork the Repository**: Start by forking the project repository to your github account.
2. **Clone Locally**: Clone the forked repository to your local machine using a git client.
   ```sh
   git clone https://github.com/adamreidsmith/ember.git
   ```
3. **Create a New Branch**: Always work on a new branch, giving it a descriptive name.
   ```sh
   git checkout -b new-feature-x
   ```
4. **Make Your Changes**: Develop and test your changes locally.
5. **Commit Your Changes**: Commit with a clear message describing your updates.
   ```sh
   git commit -m 'Implemented new feature x.'
   ```
6. **Push to github**: Push the changes to your forked repository.
   ```sh
   git push origin new-feature-x
   ```
7. **Submit a Pull Request**: Create a PR against the original project repository. Clearly describe the changes and their motivations.
8. **Review**: Once your PR is reviewed and approved, it will be merged into the main branch. Congratulations on your contribution!
</details>

<details closed>
<summary>Contributor Graph</summary>
<br>
<p align="left">
   <a href="https://github.com{/adamreidsmith/ember.git/}graphs/contributors">
      <img src="https://contrib.rocks/image?repo=adamreidsmith/ember.git">
   </a>
</p>
</details>

---

## License

Ember.git is protected under the [LICENSE](https://choosealicense.com/licenses) License. For more details, refer to the [LICENSE](https://choosealicense.com/licenses/) file.

---

## Acknowledgments

- Credit `contributors`, `inspiration`, `references`, etc.

<div align="right">

[![][back-to-top]](#top)

</div>


[back-to-top]: https://img.shields.io/badge/-BACK_TO_TOP-151515?style=flat-square


---
