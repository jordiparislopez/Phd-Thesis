#include <iostream>
#include <stdio.h>      /* printf */
#include <math.h>       /* Pow */
#include "cubature.h"
#include "hcubature.c"

using namespace std;
/*
int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    double sigma = *((double *) fdata); // we can pass σ via fdata argument
    double sum = 0;
    unsigned i;
    for (i = 0; i < 1000; ++i)
    {fval[i] = i*exp(x[0]*x[0]*i*sigma);}
    //cout << x[i] << endl;
    // compute the output value: note that fdim should == 1 from below


    return 0; // success
}

int main()
{
    int i;
    double xmin[1] = {-1}, xmax[1] = {1}, sigma = 0.5, val[1000], err;
hcubature(1000, f, &sigma, 1, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, val, &err);
cout.precision(25);
for (i = 0; i < 1000; ++i)
{cout << val[i] << endl;}
    return 0;
}

*/


int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    fval[0] = x[0]*x[1];
    return 0; // success
}

int main()
{
    double xmin[2] = {0,-1}, xmax[2] = {200,1}, sigma = 0.5, val, err;
hcubature(1, f, NULL, 2, xmin, xmax, 0, 1e-12, 1e-4, ERROR_INDIVIDUAL, &val, &err);
printf("Computed integral = %0.10g +/- %g\n", val, err);
    return 0;
}
