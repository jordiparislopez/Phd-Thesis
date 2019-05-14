void dressing_read(){

    vcdouble P2(ncp), jac(ncp), AP(ncp), BP(ncp);
    comp_mom_gen_0(P2,jac);

    std::ifstream Asr("Complex_ABs/As.txt");
    std::ifstream Bsr("Complex_ABs/Bs.txt");

    complex <double> a, b;

    unsigned int i = 0; while (Asr >> P2[i] >> AP[i]){i++;}
                 i = 0; while (Bsr >> P2[i] >> BP[i]){i++;}

     complex_parabola_A.set_points(AP,P2,jac);
     complex_parabola_B.set_points(BP,P2,jac);

}



template<typename T, size_t N>
void initialconditions_BSE(T (&arr)[N])
{
    ifstream file("initialconditions-BSE.txt");
    if(file.is_open())
    {
        for(unsigned int i = 0; i < N; ++i)
        {
            file >> arr[i];
        }
    }
}

void global_variable_BSE_reading(){

    double initcond[6];
    initialconditions_BSE(initcond);
    nq = initcond[0];
    nz = initcond[1];
    ny = initcond[2];
    nQ = initcond[3];
    nZ = initcond[4];
    NP2 = initcond[5];

}

void mom_gen_1_bse(int ndum, vdouble & P2, vdouble & jac, vdouble & wl2){
    vdouble xl2(ndum);
    legxw(ndum,xl2,wl2);
    for(unsigned int i = 0; i < ndum; i++){
        P2[i] = (1.0 + xl2[i])/(1.0 - xl2[i]);
        jac[i] = 2.0*pow(1.0 - xl2[i],-2);
    }
}

void bse_remover(){
removedirectoryrecursively("BSE-Solutions/Homogeneous");
removedirectoryrecursively("BSE-Solutions/Inhomogeneous");
}
