# The DSE-BSE solver code

This code uses the expressions generated in FunExp.nb to solve the Dyson-Schwinger and Bethe-Salpeter equations. 
The template provided solves the Quark Propagator DSE in the Rainbow-Ladder truncation using the Maris-Tandy model and provides also a routine to solve the homogoneous and inhomogeneous BSE in this truncation. 

## Compiling

The code was created using C and C++ functions and executed on Ubuntu 18.04. The make file already contains all the necessary flags required, but a user employing a different OS may need to install libraries in order for it to run. Moreover, the libraries were the most optimal for the hardware at hand, feel free to use recent versions.

Some of the most important features include:

- QMDSE.cc contains the main file, where the DSE or BSE can be chosen.
- All the headers are available in the headers.hpp file.
- The code is compiled using clang++-6.0 (faster and more optimal than g++).
- The flag -O3 is used for more efficient running, it is recommended using lower optimisation flags before final run.
- Parametrisation is implemented using openmp (might require library installation).
- Some libraries (Eigen for instance) are not installed and are called from the MathTools folder.

Check the path of the libraries in the headers file to make sure that every library and file is included correctly.


## DSE

The user is allowed to chose between real-valued or complex-valued Euclidean momentum. 

#### Real DSEs

The real DSEs part generates momentum grid for integration and performs:

- Calculation of renormalisation constants via solving the DSEs (implemented in momintreal.hpp) at the ren. scale.
- Calculation of every dressing at every momentum squared grid point using the Zs computed in the previous step.
- Iteration of previous steps until desired converenge.

The results are saved in the folder Data/RealDSE.

#### Complex DSEs

Using the results of the Zs computed in the real DSEs (as they are equivalent) the routine performs the following:

- Generation of complex momentum parabola according to the initial conditions.
- Implementation of complex interpolation using Cauchy's formula.
- Recalculation of the complex-valued dressings via redefinition of complex parameters.
- Projection to the real axis for comparison with Real DSEs results.

The generalisation to other DSE's is straightforward: one needs only to generate new expressions in FunExp.nb and keep them in the expression folders as txt files. The renormalisation conditions need to be implemented manually, while with the new expressions only extra parameters and more for() functions are necessary.

The results are saved in the folder Data/CompDSE.

## BSE

The implementation of the BSE is found directly in the QMDSE.cc file. Using Eigen, the code generates a very large matrix that is filled with the values of the evaluated expressions obtained using the FunExp.nb file. Hence, this routines allows the user to:

- Create the BSE-related matrix.
- Solve the inhomogeneous BSE with matrix manipulation for each of the four possible vectors.
- Solve the homogeneous BSE obtaining the correct Eigenvalues and Eigenvectors.

The results are saved in the folder Data/BSE-Solutions.
