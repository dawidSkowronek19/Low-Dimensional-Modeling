This repository contains source code and scripts developed during the 
Low-Dimensional Modeling (Modelowanie Niskowymiarowe) course. 

About the Project
The codes included here are numerical implementations of physical models and simulations created 
as part of the laboratory assignments.


##Important Note on Running the Code

Please note that some codes require you to manually create an output directory before execution 
(details on which specific folders are needed are provided below). 
Since these programs were developed as university assignments, they were not designed with strict 
user-friendliness or robust error-handling in mind. Therefore, the code does not automatically check for or 
create missing output folders, and failing to create them may result in a segmentation fault.

Additionally, each project folder contains a PDF file, which serves as a comprehensive report detailing the completed 
assignment, methodology, and results.

************************************
project01

* **Matrix Diagonalization:** Uses the Finite Difference Method (FDM) to map the Hamiltonian into a tridiagonal 
matrix, computing energy levels and wavefunctions via `SciPy` eigensolvers.
* **Shooting Method:** Implements the Numerov algorithm as an alternative numerical approach to find energy eigenvalues
 and wavefunctions.
* **Grid Optimization:** Automatically determines the minimum number of grid points (N) required to achieve a target 
accuracy when compared to analytical solutions.
* **Domain Analysis:** Investigates how the width of the simulated system (L) impacts the computed energy levels.
* **Data Export & Visualization:** Automatically normalizes wavefunctions, saves probability density data to text 
files, and generates plots for visual analysis.

************************************
project02

This repository contains a numerical solver for the 1D Schrödinger equation designed to 
find the ground state of a quantum system using the Imaginary Time Evolution method.


* **Supported Potentials:** Models both the infinite square well and the finite square well 
* **Imaginary Time Evolution:** Iteratively propagates the wavefunction in imaginary time. 
The state is normalized at each step, and its energy is evaluated using a discrete Hamiltonian operator.
* **Tolerance-Based Convergence:** The solver automatically stops when the energy difference between 
consecutive iterations drops below a user-defined tolerance threshold.
* **Domain Size Analysis:** Allows investigating how the size of the computational box affects the numerical 
ground state of a finite well with a fixed physical width.

************************************
project03

in order to compile: g++ -fopenmp -O3 *.cpp  -o ./test

This repository contains an C++ object-oriented numerical solver for the 2D Poisson equation using the 
Finite Difference Method (FDM).

* **Object-Oriented Design:** The architecture cleanly separates the physical domain definition (`Geometry` class) 
 from the iterative numerical calculations (`Solver` class).
* **Iterative Algorithms:** Computes spatial potential distributions using **Gauss-Seidel** and 
**Successive Over-Relaxation (SOR)** solvers, implemented via a checkerboard update scheme.
* **Physical Systems:** Includes built-in test environments for three distinct nanostructures:
  * Uniform disk
  * Quantum wire
  * Quantum dot controlled by external gating electrodes
* **Boundary Conditions:** Fully handles grid-edge behavior by supporting 
  **Dirichlet**, **Neumann**, and **Born**boundary conditions.
* **Data Export:** Automatically saves the converged numerical results of the spatial potential into text 
  files for further analysis and visualization.

************************************
project04



This repository contains a numerical simulator for tracking the time evolution of a 1D wave packet by solving 
the Time-Dependent Schrödinger Equation.

Require ./outdir

Key Features

* **Gaussian Initialization:** The initial quantum state is defined as a Gaussian wave packet with customizable 
parameters for spatial width and initial momentum.
* **Physical Scenarios:** Simulates various quantum phenomena, including free particle propagation, 
scattering/tunneling through a potential barrier, and boundary handling using a Complex Absorbing Potential.
* **Crank-Nicolson Scheme:** Employs the implicit, unconditionally stable Crank-Nicolson method for robust 
spatial and temporal discretization.
* **Thomas Algorithm:** Efficiently solves the resulting tridiagonal system of linear equations at every time step.

************************************
project05


This repository contains a numerical solver designed to calculate the transmission coefficient of a quantum 
particle passing through various 1D potential barriers by solving the Time-Independent Schrödinger Equation.

requires ./outdir

Key Features

* **Barrier Geometries:** Simulates quantum tunneling and scattering across rectangular, Gaussian, and double barrier 
potentials.
* **Numerical Core:** Discretizes the stationary Schrödinger equation using the 
Finite Difference Method, reducing the problem to a tridiagonal system of linear equations efficiently solved 
via the Thomas algorithm.
* **Open Boundary Conditions:** Implements specialized boundary conditions to accurately model incoming and outgoing 
plane waves without artificial reflections.
* **Energy Sweeps:** Iteratively calculates the transmission coefficient T as a function of the incident energy E,
 exporting the results to text files for further analysis.
* **Resonant Tunneling:** Identifies and analyzes specific energy levels where perfect transmission T(E)=1 occurs, 
which is particularly relevant for double-barrier systems.

************************************
project06

This repository contains a numerical solver for calculating the quantum conductance of a 2D 
nanostructure using the Quantum Transmitting Boundary Method.

requires ./outdir

Key Features

* **Nanostructure Geometry:** Defines a 2D potential profile V(x, y) designed to simulate systems such as a 
Quantum Point Contact (QPC).
* **Transverse Modes:** Diagonalizes the transverse Hamiltonian for each cross-section along the transport axis 
to determine the active propagating quantum modes.
* **Numerical Core:** Solves the discretized Schrödinger equation using the Thomas algorithm for rapid 
inversion of tridiagonal matrices.
* **Open Boundary Conditions:** Implements specialized boundary conditions to accurately model the injection and 
transmission of electrons through the system leads.
* **Landauer Formalism:** Computes the total quantum conductance by summing the transmission coefficients of 
all active transverse modes at a given Fermi energy using the Landauer formula.
* **Data Export:** Automatically saves the conductance vs. energy characteristics and the 2D potential map to text 
files for further visualization and analysis.

