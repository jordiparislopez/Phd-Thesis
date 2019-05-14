/*
    Headers file
    Includes all relevant packages and external files.
*/


// Include omp parallelisation package
#include <omp.h>


// Include C and C++ packages
#include <float.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <complex>
#include <ctime>
#include <string>


// Include files from the MathTools folder
#include "../../MathTools/Interpolation/spline.hpp"
#include "../../MathTools/Interpolation/cauchy.hpp"
#include "../../MathTools/Interpolation/cauchy1.hpp"
#include "../../MathTools/Removedir/removedir.cc"
#include "../../MathTools/Complexvector/complexvector.hpp"


// Include used Eigen packages and 'using' terms
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>

using Eigen::VectorXcd;
using Eigen::ComplexEigenSolver;
using Eigen::MatrixXcd;


// Global definitions of interpolating functions
tk::spline splineA;
tk::spline splineB;
cauchyint::cauchy complex_parabola_A;
cauchyint::cauchy complex_parabola_B;


// Redefinition of short functions and constants
typedef std::vector  < double > state_type;
typedef std::vector  < double > vdouble;
typedef std::vector  < std::complex < double > > vcdouble;
typedef std::complex < double > comp;
const comp I(0.0,1.0);

# define PI 3.141592653589793238462643
# define E  2.718281828459045235360287
# define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()


// Std using
using namespace std;


// Include all the functions in the proper order
#include "globalvariables.hpp"
#include "simpfunctionsDSE.hpp"
#include "simpfunctionscomplexDSE.hpp"
#include "simpfunctionsBSE.hpp"
#include "beginend.hpp"
#include "momintreal.hpp"
#include "momintcomp.hpp"
#include "propcomp.hpp"
#include "propreal.hpp"
#include "DSE-rc.hpp"
#include "DSE.hpp"
#include "BSE-Solution.hpp"
#include "BSE.hpp"
