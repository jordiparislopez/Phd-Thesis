/*
    File containing the explicit form of the regulators


    The definitions in this file are related to the general
    definition of the regulators used in the thesis
*/

// Regulator rk
double Rk( double p2,const  double k1){
    double res;
    double x = p2/(k1*k1);

    res = x/(exp(x*x)-1.0);
    return res;
}

// Derivative of the regulator over the scale
double dRk( double p2,const  double k1){
    double res;
    double x = p2/(k1*k1);

    res=2*x*(1 + exp(x*x)*(-1 + 2*x*x))*pow(-1 + exp(x*x),-2)/k1;
    return res;
}
