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
- Some libraries are not installed and called from the MathTools folder.


## DSE

The generalisation to other DSE's is straightforward: one needs only to generate new expressions in FunExp.nb and keep them in the expressions folders as txt files
