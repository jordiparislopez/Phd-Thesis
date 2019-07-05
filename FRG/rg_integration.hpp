/*
    rg_integration.hpp file.
    File to define parameters to execute the ODE solver
*/

void rg_integration()
{

    // Definition arrays containing the values of the global variables
    double inval0[12];
    double rgvalues[7];

    // Functions to read global variables
    initialconditions(inval0);
    initialrgvalues(rgvalues);

    // Definition of the necessary parameters
    NL2 = inval0[0];
    par = 9 + 5*NL2;
    double inval[3] = {pow(inval0[1],2),inval0[2],inval0[3]};

    // Assignation of the global variables from the array
    Rho0 = pow(inval0[4],2)/2.0;
    c0 = inval0[5];
    luv = 4*inval0[6]*inval0[6];
    lir = inval0[7];
    Nf = inval0[8];
    Nc = inval0[9];
    abs_err = inval0[10];
    rel_err = inval0[11];

    // Definition and generation of the momentum grid
    double P2[NL2];
    mom_gen(P2);

    // Assignation of the relevant ODE values from the RG-global variables
    double initial_scale = rgvalues[0];
    double final_scale = rgvalues[1];
    double initial_step = rgvalues[2];
    double abs_err = rgvalues[3];
    double rel_err = rgvalues[4];
    double a_x = rgvalues[5];
    double a_dxdt = rgvalues[6];

    // Boost definition of the RKCK54
    typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
    typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type; controlled_stepper_type controlled_stepper(default_error_checker< double , range_algebra , default_operations >(abs_err,rel_err,a_x,a_dxdt ));

    // Choosing initial values for every relevant quantity.
    state_type x0(par);
    x0[0] = 1.e+8;                  // Potential 0
    x0[1] = inval[0];               // Potential 1
    x0[2] = inval[1];               // Potential 2
    x0[3] = 0.0;                    // Potential 3
    x0[4] = 0.0;                    // Potential 4
    x0[5] = 0.0;                    // Potential 5
    x0[6] = 0.0;                    // Potential 6
    x0[7] = 0.0;                    // Potential 7
    x0[8] = 0.0;

    for(unsigned int i = 0; i < NL2; i++){
        x0[9 + i] = inval[2];       // Yukawa
        x0[9 + i + NL2] = 1;        // Zpion
        x0[9 + i + 2*NL2] = 1;      // Zsigma
        x0[9 + i + 3*NL2] = 1;      // Zq
        x0[9 + i + 4*NL2] = 0;      // A
    }


    cout << "\n" << "Solving the system for " << NL2 << " momentum points." << "\n"<< "\n";

    // Calling ODE solver
    integrate_adaptive( controlled_stepper, ode(P2), x0 , initial_scale , final_scale , initial_step, write_cout(P2) );

    // Printing final values to external vile once the system converges
    file_output(P2);

}
