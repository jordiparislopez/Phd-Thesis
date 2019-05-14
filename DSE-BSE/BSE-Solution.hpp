/*
  BSE file

  This file uses the Eigen package to work with the matrix and
  its eigenvalues
*/

void bse_solution(complex <double> P2){

  MatrixXcd m(4*nQ*nZ,4*nq*nz);   // Generates relevant BSE matrices using eigen
  MatrixXcd n(4*nQ*nZ,4*nq*nz);
  MatrixXcd M(4*nQ*nZ,4*nq*nz);
  VectorXcd V(4*nq*nz);

  // Definition of relevant arrays in form of real and complex vectors
  vdouble q2(nq) , jac(nq) , wl(nq), sqrtq2(nq);
  vdouble Q2(nQ) , jacQ(nQ) , wL(nQ);
  vdouble zqP(nz) , wz(nz) , zQP(nZ) , wZ(nZ), xy(ny) , wy(ny);
  vcdouble p1(nq*nz) , p2(nq*nz), DT1(nq*nz) , DT2(nq*nz);
  vcdouble AP1(nq*nz) , AP2(nq*nz), BP1(nq*nz) , BP2(nq*nz);

  // Definition of complext vectors indicating matrix terms
  vcdouble bse11(nQ*nZ*nq*nz), bse12(nQ*nZ*nq*nz), bse13(nQ*nZ*nq*nz), bse14(nQ*nZ*nq*nz);
  vcdouble bse21(nQ*nZ*nq*nz), bse22(nQ*nZ*nq*nz), bse23(nQ*nZ*nq*nz), bse24(nQ*nZ*nq*nz);
  vcdouble bse31(nQ*nZ*nq*nz), bse32(nQ*nZ*nq*nz), bse33(nQ*nZ*nq*nz), bse34(nQ*nZ*nq*nz);
  vcdouble bse41(nQ*nZ*nq*nz), bse42(nQ*nZ*nq*nz), bse43(nQ*nZ*nq*nz), bse44(nQ*nZ*nq*nz);

  // Generation of weights and abscissas for the quadrature
  mom_gen_1_bse(nq, q2, jac,  wl);
  mom_gen_1_bse(nQ, Q2, jacQ, wL);
  cheb2xw(nz, zqP, wz);
  cheb2xw(nZ, zQP, wZ);
  legxw(ny, xy, wy);

  double zqQ, k2, mtk;
  complex <double> sqrtP2 = sqrt(P2);
  unsigned int iQ;
  unsigned int jQ;
  unsigned int k;
  unsigned int i;
  unsigned int j;

  // Parallelised loop
  omp_set_num_threads(8);
  #pragma omp parallel for private(zqQ, k2, mtk, iQ, jQ, k, i, j)
  // Loops for columns
  for(i = 0; i < nq; i++){

    sqrtq2[i] = sqrt(q2[i]);

    for(j = 0; j < nz; j++){

      p1[i + j*nq] = q2[i] + P2/4.0 + sqrtq2[i]*sqrtP2*zqP[j];
      p2[i + j*nq] = q2[i] + P2/4.0 - sqrtq2[i]*sqrtP2*zqP[j];
      AP1[i + j*nq] = complex_parabola_A( p1[i + j*nq] );
      AP2[i + j*nq] = conj(AP1[i + j*nq]);
      BP1[i + j*nq] = complex_parabola_B( p1[i + j*nq] );
      BP2[i + j*nq] = conj(BP1[i + j*nq]);
      DT1[i + j*nq] = 1.0/(p1[i + j*nq]*AP1[i + j*nq]*AP1[i + j*nq] + BP1[i + j*nq]*BP1[i + j*nq]);
      DT2[i + j*nq] = 1.0/(p2[i + j*nq]*AP2[i + j*nq]*AP2[i + j*nq] + BP2[i + j*nq]*BP2[i + j*nq]);

      // Loops for rows
      for(iQ = 0; iQ < nQ; iQ++){
        for(jQ = 0; jQ < nZ; jQ++){

          bse11[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse12[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse13[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse14[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse21[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse22[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse23[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse24[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse31[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse32[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse33[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse34[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse41[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse42[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse43[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;
          bse44[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = 0.0;

          for(k = 0; k < ny; k++){

            zqQ = zqP[j]*zQP[jQ] + xy[k]*sqrt(1.0 - zqP[j]*zqP[j])*sqrt(1.0 - zQP[jQ]*zQP[jQ]);
            k2 = q2[i] + Q2[iQ] - 2.0*sqrtq2[i]*sqrt(Q2[iQ])*zqQ;
            mtk = alphamt(k2);

            // Matrix elements computation

            bse11[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse11[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse11.txt"
            ;

            bse12[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse12[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse12.txt"
            ;

            bse13[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse13[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse13.txt"
            ;

            bse14[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse14[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse14.txt"
            ;

            bse21[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse21[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse21.txt"
            ;

            bse22[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse22[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse22.txt"
            ;

            bse23[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse23[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse23.txt"
            ;

            bse24[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse24[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse24.txt"
            ;

            bse31[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse31[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse31.txt"
            ;

            bse32[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse32[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse32.txt"
            ;

            bse33[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse33[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse33.txt"
            ;

            bse34[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse34[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse34.txt"
            ;

            bse41[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse41[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse41.txt"
            ;

            bse42[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse42[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse42.txt"
            ;

            bse43[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse43[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse43.txt"
            ;

            bse44[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] = bse44[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i] +
            #include "Exp-BSE/bse44.txt"
            ;
        }

        // Matrix elements assignation
        m(jQ + iQ*nZ + 0*(nQ*nZ) , j + i*nz + 0*(nq*nz)) = bse11[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 0*(nQ*nZ) , j + i*nz + 1*(nq*nz)) = bse12[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 0*(nQ*nZ) , j + i*nz + 2*(nq*nz)) = bse13[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 0*(nQ*nZ) , j + i*nz + 3*(nq*nz)) = bse14[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 1*(nQ*nZ) , j + i*nz + 0*(nq*nz)) = bse21[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 1*(nQ*nZ) , j + i*nz + 1*(nq*nz)) = bse22[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 1*(nQ*nZ) , j + i*nz + 2*(nq*nz)) = bse23[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 1*(nQ*nZ) , j + i*nz + 3*(nq*nz)) = bse24[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 2*(nQ*nZ) , j + i*nz + 0*(nq*nz)) = bse31[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 2*(nQ*nZ) , j + i*nz + 1*(nq*nz)) = bse32[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 2*(nQ*nZ) , j + i*nz + 2*(nq*nz)) = bse33[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 2*(nQ*nZ) , j + i*nz + 3*(nq*nz)) = bse34[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 3*(nQ*nZ) , j + i*nz + 0*(nq*nz)) = bse41[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 3*(nQ*nZ) , j + i*nz + 1*(nq*nz)) = bse42[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 3*(nQ*nZ) , j + i*nz + 2*(nq*nz)) = bse43[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        m(jQ + iQ*nZ + 3*(nQ*nZ) , j + i*nz + 3*(nq*nz)) = bse44[jQ + iQ*nZ + nQ*nZ*j + nQ*nZ*nz*i];
        }
      }
    }
  }

  // Definition vector for Inhomogeneous BSE
  for( j = 0; j < 4*nq*nz; j++){
    V(j) = 0;
    if(j < nq*nz){V(j) = Z2;}
  }

  // Inhomogeneous BSE matrix
  n = MatrixXcd::Identity(4*nQ*nZ, 4*nq*nz);
  M = n - m;
  VectorXcd Inh =  M.completeOrthogonalDecomposition().solve(V);
    //cout << P2 << "\t"<< Inh(nz/2) << endl;


  // Using Eigen to obtain eigenvalues and eigenvectors
  VectorXcd EigVec(4*nQ*nZ);
  ComplexEigenSolver<MatrixXcd> ces;
  ces.compute(m);
  double EigVal = abs(ces.eigenvalues()[4*nQ*nZ-1]);
  EigVec = ces.eigenvectors().col(4*nQ*nZ-1);
  // double EigVal = 0;



  // Writing data in external files
  for(unsigned int i = 0; i < nQ; i++){
      for(unsigned int j = 0; j < nZ; j++){

      ofstream BSE1("Data/BSE-Solutions/Homogeneous/F1.txt", ios_base::app | ios_base::out);
      BSE1 << real(P2) << "\t" <<  Q2[i] << "\t" << zQP[j] << "\t"  << real(EigVec(j + i*nZ + 0*nZ*nQ))
           << "\t"  << imag(EigVec[j + i*nZ + 0*nZ*nQ]) << "\t"  << EigVal << endl;
      ofstream BSE2("Data/BSE-Solutions/Homogeneous/F2.txt", ios_base::app | ios_base::out);
      BSE2 << real(P2) << "\t" <<  Q2[i] << "\t" << zQP[j] << "\t"  << real(EigVec[j + i*nZ + 1*nZ*nQ])
           << "\t"  << imag(EigVec[j + i*nZ + 1*nZ*nQ])<< "\t"  << EigVal << endl;
      ofstream BSE3("Data/BSE-Solutions/Homogeneous/F3.txt", ios_base::app | ios_base::out);
      BSE3 << real(P2) << "\t" <<  Q2[i] << "\t" << zQP[j] << "\t"  << real(EigVec[j + i*nZ + 2*nZ*nQ])
           << "\t" << imag(EigVec[j + i*nZ + 2*nZ*nQ])<< "\t"  << EigVal << endl;
      ofstream BSE4("Data/BSE-Solutions/Homogeneous/F4.txt", ios_base::app | ios_base::out);
      BSE4 << real(P2) << "\t" <<  Q2[i] << "\t" << zQP[j] << "\t"  << real(EigVec[j + i*nZ + 3*nZ*nQ])
           << "\t"  << imag(EigVec[j + i*nZ + 3*nZ*nQ]) << "\t"  << EigVal << endl;


      ofstream BSE1i("Data/BSE-Solutions/Inhomogeneous/F1.txt", ios_base::app | ios_base::out);
      BSE1i << real(P2) << "\t" <<  Q2[i] << "\t" << zQP[j] << "\t" << real(Inh(j + i*nZ + 0*nZ*nQ))
            << "\t"  << imag(Inh[j + i*nZ + 0*nZ*nQ])<< "\t"  << EigVal << endl;
      ofstream BSE2i("Data/BSE-Solutions/Inhomogeneous/F2.txt", ios_base::app | ios_base::out);
      BSE2i << real(P2) << "\t" <<  Q2[i] << "\t" << zQP[j] << "\t"  << real(Inh[j + i*nZ + 1*nZ*nQ])
            << "\t"  << imag(Inh[j + i*nZ + 1*nZ*nQ]) << "\t"  << EigVal << endl;
      ofstream BSE3i("Data/BSE-Solutions/Inhomogeneous/F3.txt", ios_base::app | ios_base::out);
      BSE3i << real(P2) << "\t" <<  Q2[i] << "\t" << zQP[j] << "\t"  << real(Inh[j + i*nZ + 2*nZ*nQ])
            << "\t"  << imag(Inh[j + i*nZ + 2*nZ*nQ]) << "\t"  << EigVal << endl;
      ofstream BSE4i("Data/BSE-Solutions/Inhomogeneous/F4.txt", ios_base::app | ios_base::out);
      BSE4i << real(P2) << "\t" <<  Q2[i] << "\t" << zQP[j] << "\t"  << real(Inh[j + i*nZ + 3*nZ*nQ])
            << "\t"  << imag(Inh[j + i*nZ + 3*nZ*nQ]) << "\t"  << EigVal << endl;

        }
    }
}
