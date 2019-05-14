/*
  File to perform the complex plane integration
*/

void momintcomp(const vcdouble P2, const vcdouble JAC ,vcdouble & APpre, vcdouble &BPpre){

    int NP2 = P2.size();
    int np2 = NP2/2;
    vcdouble AP(NP2);
    vcdouble BP(NP2);

    vcdouble sqrtp2(np2);
    vdouble k2(nk);
    vdouble sqrtk2(nk);
    vdouble xk(nk);
    vdouble wk(nk);
    vdouble jac(nk);
    vdouble mtk(nk);
    vdouble zpk(nz);
    vdouble wz(nz);
    vcdouble q2(nk*nz*np2);
    vcdouble sqrtq2(nk*nz*np2);
    vcdouble aq(nk*nz*np2);
    vcdouble bq(nk*nz*np2);
    vcdouble aqterm(nk*nz*np2);
    vcdouble bqterm(nk*nz*np2);
    vcdouble denterm(nk*nz*np2);
    vcdouble z(nk*nz*np2);

    mom_gen_1(nk, k2, jac, wk);
    cheb2xw(nz,zpk,wz);

    vcdouble res0(np2);
    vcdouble res1(np2);

    // j1 is an integer boolean to determine convergence condition
    int j1 = 0;
    while( j1 == 0){
      j1 = 1;

      // We paralelise using openmp
      omp_set_num_threads(8);
      #pragma omp parallel for private(mtk, q2, sqrtq2, z, aq, bq, denterm, aqterm, bqterm)
      for(unsigned int i = 0; i < np2; i++){

        sqrtp2[i] = sqrt(P2[i]);
        res0[i] = 0;
        res1[i] = 0;

        	for(unsigned int j = 0; j < nk; j++){

            mtk = 4.0*PI*alphamt(k2[j])/k2[j];
            sqrtk2[j] = sqrt(k2[j]);

            for(unsigned int k = 0; k < nz; k++){

        			q2 = k2[j] + P2[i] + 2.0*zpk[k]*sqrtk2[j]*sqrtp2[i];
        			sqrtq2 = sqrt(q2);			// q²=(k-p)²
        			z = 1.0/sqrtq2*(sqrtk2[j]*zpk[k] + sqrtp2[i]);
    	        aq = complex_parabola_A(q2);
        			bq = complex_parabola_B(q2);
              denterm = (pow(aq,2)*q2 + pow(bq,2));

              // Condition to avoid NaN from 0 denominator
              if(denterm == 0.0){  aqterm = 0.0;     bqterm = 0.0;  }
              else{
                  aqterm = aq/denterm;
                  bqterm = bq/denterm;
              }

    	        res0[i] = res0[i] + jac[j]*wk[j]*wz[k]*k2[j]*
    	        #include "Expressions/Exp-DSE/expA.txt"
    	        ;

  		        res1[i] = res1[i] + jac[j]*wk[j]*wz[k]*k2[j]*
    	        #include "Expressions/Exp-DSE/expB.txt"
    	        ;

                   }
           }

           // Update of the dressings
           AP[i] = Z2 + Z2*Z2*res0[i];
           BP[i] = Z2*Zm*mq + Z2*Z2*res1[i];

           // Check convergence condition
           if(abs(1.0 - APpre[i]/AP[i] ) > err || abs(1.0 - BPpre[i]/BP[i]) > err){j1 = j1 * 0;}
           else{j1 = j1 * 1;}

            APpre[i] = AP[i];
            BPpre[i] = BP[i];
            AP[NP2 - i - 1] = conj(AP[i]);
            BP[NP2 - i - 1] = conj(BP[i]);

        }

        // Remove previous files and update grid
        comp_file_remover();
        compdress_setting(P2,JAC,AP,BP);
    }

    // Write obtained values
    compdress_print(NP2,P2,AP,BP);
}
