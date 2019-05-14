/*
  File to run the BSE routine using the complex propagator data
*/

void bse_solving()
{
    // Read BSE initial conditions for the matrix
    global_variable_BSE_reading();

    // Read complex propagator dressings
    dressing_read();

    // Remove any previous BSE result
    bse_remover();


    // Generate real momentum grid
    vdouble P2(NP2);
    vcdouble a(NP2);
    mom_gen_0(NP2,P2);



    int j;
    for(unsigned int i = 0; i < a.size(); i++){

      // Shift to work with negative real values
      a[i] = 0.25 - i*0.001;
        if(real(a[i]) <= 0){
        j = 1;
        }
        else{j=-1;}
      a[i] =  complex <double> (j,0)*a[i]*a[i];

      // Start BSE routine
      bse_solution(a[i]);
    }

}
