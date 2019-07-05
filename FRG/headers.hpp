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

// Include boost for ODE solving
#include <boost/numeric/odeint.hpp>

// Include files from the MathTools folder
#include "../../MathTools/Removedir/removedir.cc"
#include "../../MathTools/Eigen/Eigen/Dense"
#include "../../MathTools/Eigen/Eigen/Eigenvalues"
#include "../../MathTools/Eigen/Eigen/Core"
#include "../../MathTools/Interpolation/spline.h"
#include "../../MathTools/Cubature/cubature.h"
#include "../../MathTools/Cubature/hcubature.c"


// Global definitions of interpolating functions
tk::spline splineH;
tk::spline splineFP;
tk::spline splineA;
tk::spline splineFS;

// Redefinition of short functions and constants
using namespace std;
using namespace boost::numeric::odeint;
typedef std::vector< double  > state_type;
typedef std::vector< double  > Vector;

# define PI    3.141592653589793238462643383279502884197169399375105820974944
# define E     2.718281828459045235360287471352662497757247093699959574966968
# define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

// Include all the functions in the proper order
#include "globalvariables.hpp"
#include "regulators.hpp"
#include "simpfunctions.hpp"
#include "momintegration.hpp"
#include "minimization.hpp"
#include "integratefunctions.hpp"
#include "beginend.hpp"
#include "output.hpp"
#include "rg_integration.hpp"
