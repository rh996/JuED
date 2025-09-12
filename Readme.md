# JuED: A Julia-based Exact Diagonalization Package

JuED is a powerful and flexible Julia package for performing exact diagonalization (ED) of quantum Hamiltonians, with a focus on condensed matter physics models. It provides a framework for constructing Hilbert spaces with various symmetries, building sparse Hamiltonian matrices, and using iterative solvers to find low-lying eigenvalues and eigenvectors.

## Features

*   **Versatile Model Support:** JuED is designed to handle a variety of quantum models, including spinless and spinful fermions in one and two dimensions.
*   **Symmetry Exploitation:** The package can work with different Hilbert spaces to exploit symmetries, such as momentum and spin, significantly reducing the computational cost of diagonalization.
*   **High-Performance:** JuED is written in Julia and leverages its high-performance capabilities. It uses multi-threading to parallelize key computations, such as Hamiltonian construction and the calculation of reduced density matrices.
*   **Sparse Matrices:** Hamiltonians are constructed as sparse matrices, which is essential for handling the large dimensions of the Hilbert spaces encountered in many-body systems.
*   **Iterative Solvers:** The package uses the `KrylovKit.jl` library, a powerful iterative solver, to efficiently find the lowest-energy eigenstates of the Hamiltonian.
*   **Reduced Density Matrices:** JuED can compute one-, two-, and three-body reduced density matrices (RDMs), which are crucial for calculating various physical observables.
*   **Caching:** A caching mechanism is implemented for RDM calculations to speed up repeated computations.

## Supported Models and Hilbert Spaces

JuED supports a range of models and Hilbert spaces:

*   **Models:**
    *   1D and 2D Hubbard-like models
    *   Spinless and spinful fermion systems
    *   Two-band models (e.g., for Transition Metal Dichalcogenides - TMDs)

*   **Hilbert Spaces:**
    *   General fermion Hilbert spaces
    *   Momentum Hilbert spaces (1D and 2D)
    *   Spin Hilbert spaces
    *   Combined spin-momentum Hilbert spaces (1D and 2D)
    *   Two-band momentum Hilbert spaces

## Usage

The main entry point for using JuED is the `EDMain.jl` module. Here's a basic overview of the workflow:

1.  **Define Model Parameters:** Create a model parameter object (e.g., `ModelParams2DSpinlessList`) with information about the system, such as the number of particles, lattice size, and one- and two-body interaction terms.

2.  **Diagonalize the Hamiltonian:** Use the `DiagonalizeOneMomentum` or `DiagonalizeAllMomentum` functions to find the eigenvalues and eigenvectors of the Hamiltonian. You can specify the number of eigenvalues to compute.

3.  **Calculate Observables:** Once you have the ground state eigenvector, you can use the `RDM1`, `RDM2`, and `RDM3` functions to calculate the reduced density matrices.

### Example

```julia
using JuED.EDMod

# 1. Define model parameters
nparticle = 4
Nkx = 2
Nky = 2
OneBody = ... # Define your one-body integrals
TwoBody = ... # Define your two-body integrals
model = InputModel(nparticle, Nkx, Nky, OneBody, TwoBody)

# 2. Diagonalize the Hamiltonian for a specific momentum sector
momentum = 0
neigenv = 1
eigenvalues, eigenvectors = DiagonalizeOneMomentum(model, momentum, neigenv)

# 3. Calculate the 1-RDM
rdm1 = RDM1(model, eigenvectors[1], momentum)
```

## Installation

To use JuED, you need to have Julia installed. You can then clone this repository and install the required dependencies.

```bash
# (Assuming you have Julia and Git installed)
git clone https://github.com/your-username/JuED.git
cd JuED
julia -e 'using Pkg; Pkg.instantiate()'
```

## Running Tests

The tests for JuED are located in the `tests/` directory. Each test file focuses on a specific part of the codebase, such as a particular model, Hilbert space, or feature. To run a test, you can execute the corresponding test file using Julia:

```bash
julia tests/test_hubbard_momentum2d.jl
```

For a comprehensive check of the package, you can run all the test files. A main test runner script could be created to automate this process.


## Contributing

Contributions to JuED are welcome! If you have ideas for new features, bug fixes, or improvements, please open an issue or submit a pull request.
