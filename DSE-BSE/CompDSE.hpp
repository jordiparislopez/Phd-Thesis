/*
  Necessary functions for complex DSE
*/

// Executes complex DSE routine
void comp_prop(vcdouble P2,vcdouble JAC) {

    vcdouble APpre(ncp);
    vcdouble BPpre(ncp);

    // Read real data to use as initial contidion for the complex plane
    std::ifstream Asr("Data/RealDSE/A.txt");
    std::ifstream Bsr("Data/RealDSE/B.txt");
    complex <double> a, b;

    unsigned int i = 0; while (Asr >> a >> APpre[i]){i++;}
                 i = 0; while (Bsr >> b >> BPpre[i]){i++;}

    // Call the complex DSE routine using real axis Z's
    momintcomp(P2, JAC, APpre, BPpre);
}

// Projects obtained data to real axis to compare with real results
void comp_real_projection(const vdouble P2){

    // Project real axis with the values obtained using Cauchy's formula
    for(unsigned int i = 0; i < np; i++)
    {
        ofstream aP("ComplexDSE/Projection/Ap.txt", ios_base::app | ios_base::out);
            aP << P2[i] << "\t" <<  complex_parabola_A(P2[i])<< "\t" <<  complex_parabola_A(P2[i]).imag() << endl;
        ofstream bP("ComplexDSE/Projection/Bp.txt", ios_base::app | ios_base::out);
            bP << P2[i] << "\t" <<  complex_parabola_B(P2[i]).real()<< "\t" <<  complex_parabola_B(P2[i]).imag() << endl;
        ofstream mP("ComplexDSE/Projection/Mp.txt", ios_base::app | ios_base::out);
            mP << P2[i] << "\t" <<  real(complex_parabola_B(P2[i])/complex_parabola_A(P2[i]))<< "\t" << imag(complex_parabola_B(P2[i])/complex_parabola_A(P2[i])) << endl;

    }
cout << "Complex propagator solved with M2 = " << M2 << " GeV."<< endl;
}
