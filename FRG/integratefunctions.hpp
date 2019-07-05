/*
   File containing functions executed in every step during the ODE integration
*/

struct ode
{
  // Transference of variables from rg_integration file
  double * P2;
  ode(double* PPP2) : P2( PPP2 ) {}
	void operator()( const state_type  x0 , state_type  &dxdt , const double k )
	{

    // Definition of integration limit and interpolation variables
    double LUV = pow(k*2.5,2);
    state_type PP2(NL2), HQ2(NL2), AQ2(NL2), FP2(NL2),FS2(NL2);


    // Definition of hcubature parameters
    double y[8];
    double xmin[2]={lir,-1}, xmax[2]={LUV,1}, val[par], err[par];
    double fdata[10 + NL2];


    // Assignation of hcubature external parameters
    fdata[0] = k;
    for(unsigned int i = 0; i < 9; i++){
      y[i] = x0[i];
      fdata[1 + i] = y[i];
    }
    for(unsigned int i = 0; i < NL2; i++){
      fdata[10 + i] = P2[i];
    }


    // Assignation of momentum and dressing grid values
    for(unsigned int i = 0; i < NL2; i++){

        PP2[i] = P2[i];
        HQ2[i] = x0[9 + i];
        FP2[i] = x0[9 + NL2 + i];
        FS2[i] = x0[9 + 2*NL2 + i];
        AQ2[i] = x0[9 + 3*NL2 + i];

       }

    // Computing interpolation values to compute at evaluation time
    splineH.set_points(PP2,HQ2);
    splineFP.set_points(PP2,FP2);
    splineFS.set_points(PP2,FS2);
    splineA.set_points(PP2,AQ2);


    // Momentum integration of 2 variables
    hcubature(par, flowequations, &fdata, 2, xmin, xmax, 0, abs_err, rel_err, ERROR_INDIVIDUAL, val, err);


    // Assigning derivative values from the hcubature result
    for(unsigned int i = 0; i < 9; i++){dxdt[i] = val[i];}
    for(unsigned int i = 9 ; i < 9 + NL2; i++){dxdt[i] = val[i] - x0[1]*val[4*NL2 + i]/HQ2[i - 9];}
    for(unsigned int i = 9 + NL2; i < 9 + 4*NL2; i++){dxdt[i] = val[i];}
    for(unsigned int i = 9 + 4*NL2; i < par; i++){dxdt[i] = val[i]/HQ2[i - 4*NL2 - 9];}


	}
};
