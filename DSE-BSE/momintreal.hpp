/*
   Execution of momentum integral file for the calculation of Z's and dressings

   momintreal() performs the integral given P2
   Zs1() calls momintreal() at renormalisation scale
   realpropagator1() calls momintreal for proper P2 with updated dressings
*/


void momintreal(const vdouble P2, vdouble & intA, vdouble & intB){

  int np2 = P2.size();
  vdouble sqrtp2(np2);
  vdouble k2(nk);
  vdouble xk(nk);
  vdouble wk(nk);
  vdouble jac(nk);
  double zpk(nz);
  double wz(nz);

  double mtk;
  double q2, sqrtq2;
  double aq, bq, aqterm, bqterm, z;

  // Generate quadrature points
  mom_gen_1(nk, k2, jac, wk);
  cheb2xw(nz,zpk,wz);


  vdouble res0(np2);
  vdouble res1(np2);

  for(unsigned int i = 0; i < np2; i++){

    // Initial values for every p^2 dependent quantity.
    sqrtp2[i] = sqrt(P2[i]);
    res0[i] = 0;
    res1[i] = 0;

		for(unsigned int j = 0; j < nk; j++){

      // Call Maris-Tandy coupling function
      mtk = 4.0*PI*alphamt(k2[j])/k2[j];

        for(unsigned int k = 0; k < nz; k++){

          // q^2 dependent quantities
    			q2 = k2[j] + P2[i] + 2.0*zpk[k]*sqrt(k2[j]*P2[i]);
    			sqrtq2 = sqrt(q2);
    			z = 1.0/sqrtq2*(sqrt(k2[j])*zpk[k] + sqrtp2[i]);
	        aq = splineA(q2);
    			bq = splineB(q2);
          aqterm = aq/(pow(aq,2)*q2 + pow(bq,2));
          bqterm = bq/(pow(aq,2)*q2 + pow(bq,2));

          // Add every quadrature value to the total result
          res0[i] = res0[i] + jac[j]*wk[j]*wz[k]*k2[j]*
          #include "Expressions/Exp-DSE/expA.txt"
          ;

          res1[i] = res1[i] + jac[j]*wk[j]*wz[k]*k2[j]*
          #include "Expressions/Exp-DSE/expB.txt"
          ;
          }
  	    }
     intA[i] = res0[i];
     intB[i] = res1[i];
    }
}



// Function to determine new renormalisation constants
int Zs1(double &Z2, double &Zm) {

    int j1;
    vdouble intA(1);  // Definition 1-element valued vectors
    vdouble intB(1);
    vdouble MU(1);
    MU[0] = mu;

    // Call previous Z's values for convergence conditions
    double Z2pre = Z2;
    double Zmpre = Zm;

    // Solve the integral at renormalisation scale
    momintreal(MU, intA, intB);


    // Determine values of Z's
    Z2 = ( - 1.0 + sqrt( 1.0 + 4.0*intA[0]))/(2.0*(intA[0]));
    if(mq != 0 ){
        Zm = (mq - Z2*Z2*(intB[0]))/(mq*Z2);
    }else{
        Zm = 0;
    }

    // convergence condition through relative error
    if(abs(Z2-Z2pre)/ Z2 > err || abs(Zm-Zmpre)/ Zm > err){
        j1 = 0;
    }else{
        j1 = 1;
    }

  // use j1 as integer boolean
  return j1;
}



// Function to update dressings
int realpropagator1(vdouble P2) {

    vdouble intA(np);
    vdouble intB(np);
    vdouble AP(np);
    vdouble BP(np);
    vdouble APpre(np);
    vdouble BPpre(np);
    int j1 = 1;

    // Solve DSE for every p^2
    momintreal(P2, intA, intB);


    for(unsigned int i = 0; i < np; i++)
    {
        // Keep old values and update new ones for the dressings
        APpre[i] = splineA(P2[i]);
        BPpre[i] = splineB(P2[i]);
        AP[i] = Z2 + Z2*Z2*intA[i];
        BP[i] = Z2*Zm*mq + Z2*Z2*intB[i];

        // Apply convergence condition
        if(abs(AP[i]-APpre[i])/ AP[i] > err || abs(BP[i]-BPpre[i])/ BP[i] > err){
            j1 = j1 * 0;
        }else{
            j1 = j1 * 1;
        }
    }

    // Erase previous written values
    file_remover();

    // Write these and set new grid for interpolation
    dress_setting(np,P2,AP,BP);

    return j1;
}
