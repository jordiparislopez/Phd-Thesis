void bse_solving()
{
    global_variable_BSE_reading();
    dressing_read_real();
    bse_remover();
    vdouble P2(NP2);
    mom_gen_0(NP2,P2);
    vdouble a(NP2);


    // vcdouble a(15);
    // a[0] = -0.15*0.15;
    // a[1] = -0.145*0.145;
    // a[2] = -0.143*0.143;
    // a[3] = -0.1415*0.1415;
    // a[4] = -0.14*0.14;
    // a[5] = -0.139*0.139;
    // a[6] = -0.138*0.138;
    // a[7] = -0.137*0.137;
    // a[8] = -0.136*0.136;
    // a[9] = -0.135*0.135;
    // a[10] = -0.134*0.134;
    // a[11] = -0.133*0.133;
    // a[12] = -0.132*0.132;
    // a[13] = -0.131*0.131;
    // a[14] = -0.130*0.130;

    for(unsigned int i = 0; i < a.size(); i++){
        a[i] = P2[i];
        bse_solution(a[i]);
    }

}
