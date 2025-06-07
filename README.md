# Fisher-Kolmogorov Equation for Neurodegenerative Diseases

This project implements a parallel solver for a nonlinear anisotropic reaction-diffusion model describing the spread of misfolded proteins in brain tissue. The model is based on the generalized **Fisher-Kolmogorov equation**, adapted to simulate neurodegenerative processes with biologically inspired parameters and geometries.

## 📘 Model Overview

The model captures the autocatalytic transformation of healthy proteins into misfolded ones and their spatio-temporal spread through diffusion and reaction terms. It includes:

- Production and clearance of both healthy and misfolded proteins.
- Autocatalytic reaction converting healthy into misfolded forms.
- Spatial diffusion with anisotropic transport along preferred directions.

The normalized governing equation reads:

$$
\frac{\partial c}{\partial t} = \nabla \cdot (D \nabla c) + \alpha c(1 - c)
$$

with homogeneous Neumann boundary conditions.

## 🧮 Numerical Method

The model is solved using:

- **Finite Element Method (FEM)** in space
- **Backward Euler scheme** in time
- **Newton's method** for nonlinearity
- **Galerkin projection** for weak form discretization

Linear systems are solved using:

- **Conjugate Gradient (1D)** with SSOR preconditioning
- **GMRES (3D)** with ILU preconditioning (for nonsymmetric matrices)

## ⚙️ Computational Features

- Written in C++ with support for MPI parallelism
- FEM code supports 1D and 3D simulations
- Parallel assembly and solver using **Trilinos**

## 🧪 Benchmarks

Performance tests were run on:

- A local workstation (Intel i7, 4 cores, 8 GB RAM)
- **MareNostrum 5** supercomputer (Xeon 8480+, 112 cores, 256 GB RAM)

## 🚀 Usage

The project is compiled using **CMake**. After building, two executables are available:

### 🔹 1D Simulation

Run the executable without any input arguments:

    ./simulation1d

- Parameters for the 1D simulation are **hardcoded** in the source code (meaning that compilation must be done again every time the parameters are changed).
- No input files are required at runtime.

### 🔹 3D Simulation

Run the executable optionally providing a CSV parameter file path:

    ./simulation3d [/path/to/custom_param_file.csv]

- The 3D simulation reads parameters from a CSV file.
- If no command-line argument is given, the default file `../data/params.csv` is used.
### Input Data

- The CSV parameter file should contain one `<key>,<value>` pair per line. A default example is available at `data/params.csv`.
- The simulation mesh can be downloaded [here](https://polimi365-my.sharepoint.com/:u:/g/personal/10461512_polimi_it/EY9ZPoq279JArvbXLPR1pNcB-wjU5tPZLClfO_4O9EYbtg?e=C1aIRH). It is provided in `.stl` format and must be preprocessed with `gmsh` to generate a `.msh` file.
- To preprocess the mesh, use the `convert_mesh.geo` script in the `/scripts` directory.

The converter script contains two strings specifying the input `.stl` file path and the output `.msh` file path. The default values are `./brain-h3.0.stl` and `./data/volume.msh`.

If the input file with the default name is present, you can generate the mesh by running:

    gmsh scripts/convert_mesh.geo -