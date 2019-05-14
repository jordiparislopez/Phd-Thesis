/*
  File for calling Real or Complex momentum DSE routines
*/


// Real DSE routine
void dse_real(){

    // Momentum grid
    vdouble P2(np);

    //Generation of external momentum
    mom_gen_0(np,P2);

    // Generation of dressings
    dress_gen(np,P2);

    // Routine of real propagator
    real_prop(P2);
}



// Complex DSE routine
void dse_complex(){

    // Real and complex grid
    vdouble  P2(np);
    vcdouble Pc2(ncp);
    vcdouble jac(ncp);

    // Remove data from complex momentum folders
    comp_file_remover();

    // Generates complex momentum
    comp_mom_gen_0(Pc2,jac);

    // Generates complex dressings
    comp_dress_gen(Pc2,jac);

    // Calls complex DSE routine
    comp_prop(Pc2, jac);

    // Generates real grid
    mom_gen_0(np,P2);

    // Projects results to real axis
    comp_real_projection(P2);
}
