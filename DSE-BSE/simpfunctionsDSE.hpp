/*
  simpfunctionsDSE.hpp file

  This file contains custom made functions that were useful during the
  coding time. Some of them were defined to match with Mathematica C++
  exporting type or others to practice numerical routines.

*/


// Gauss-Legendre quadrature points and weights generator (Numerical Recipes)
void legxw(int n, vdouble &xleg, vdouble &wleg)
{
    vdouble x(n); // Abscissas
    vdouble w(n); // Weights
    double error = 3.0e-15; // Desired numerical error
    double p1, p2, p3, pp;
    double z, z1;
    int i,j,m = (n + 1)/2;

    for(i = 1;i <= m; i++){

      z1 = 2.0; // Initial values for z1 and z in function of number of points
      z = cos(PI*(i - 0.25)/(n + 0.5));

        while(abs((z-z1))>error){ // Iteration loop until z and z1 are close

          p1 = 1.0;
          p2 = 0.0;

          for(j = 1;j <= n; j++){ // Reassignation of values for p_i
            p3 = p2;
            p2 = p1;
            p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j;
          }

          pp = n*(z*p1 - p2)/(z*z - 1.0);
          z1 = z;
          z = z1 - p1*1.0/pp;
        }

        x[i-1] = -z;  // Abscissas and weights assigned following z and pp.
        x[n-i] = z;
        w[i-1] = 2.0/((1.0-z*z)*pp*pp);
        w[n-i] = w[i-1];
     }
     for(i = 0; i < n; i++){

        if(abs(x[i]) < error){
          xleg[i] = 0;
        }
        else{
          xleg[i] = x[i];
        };

        wleg[i] = w[i];
	}
}



// Gauss-Chebyshev quadrature generator
void cheb2xw(int n, vdouble &xcheb2, vdouble & wcheb2)
{
    for(unsigned int i = 1;  i <= n ; i++){
        xcheb2[n-i] = cos((i*PI*1.0/(n + 1)));
        wcheb2[n-i] = PI/(n + 1)*sin(i*PI*1.0/(n + 1))*sin(i*PI/(n + 1));
    }
}




// File remover routine to empty data folders (in MathTools folder)
void file_remover(){

  removedirectoryrecursively("RealDSE");
  removedirectoryrecursively("ComplexDSE/CompDress");
  removedirectoryrecursively("ComplexDSE/Path");
  removedirectoryrecursively("ComplexDSE/Projection");

}



// Generator of p^2 points according to the desired distribution
void mom_gen_0(int np, vdouble & P2){

    vdouble xl(np);
    vdouble wl(np);
    legxw(np,xl,wl);

    for(unsigned int i = 0; i < np; i++){
        // p^2 logarithmic distribution
        P2[i] = sqrt(luv*lir)*pow(sqrt(luv/lir),xl[i]);

        // p^2 in (0,∞) with undefined boundaries
        // double m = 1;
        // P2[i] = pow((1 + xl[i])/(1 - xl[i]),m);
    }
}



// Generator of p^2 points for the numerical integration with jacobian
void mom_gen_1(int ndum, vdouble & P2, vdouble & jac, vdouble & wl){

    vdouble xl(ndum);
    legxw(ndum,xl,wl);

    for(unsigned int i = 0; i < ndum; i++){
        //p^2 logarithmic distribution
        P2[i] = sqrt(luv*lir)*pow(sqrt(luv/lir),xl[i]);
        jac[i] = 1.0/2.0*log(luv/lir)*P2[i];

        // p^2 in (0,∞) with undefined boundaries
        // double m = 1;
        // P2[i] = pow((1 + xl[i])/(1 - xl[i]),m);
        // jac[i] = 2*m*P2[i]/(1 - xl[i]*xl[i]);
    }
}


// Generation of dressings for the real propagator
void dress_gen(int np, vdouble & P2){

    vdouble AP(np);
    vdouble BP(np);

    for(unsigned int i = 0; i < np; i++){

        AP[i] = a0; // A and B take the initial values
        BP[i] = b0;

    // Values are written now in the files

    ofstream aP("RealDSE/As.txt", ios_base::app | ios_base::out);
        aP << P2[i] << "\t" <<  AP[i] << endl;
    ofstream bP("RealDSE/Bs.txt", ios_base::app | ios_base::out);
        bP << P2[i] << "\t" <<  BP[i] << endl;
    ofstream mP("RealDSE/Ms.txt", ios_base::app | ios_base::out);
        mP << P2[i] << "\t" <<  BP[i]/AP[i] << endl;

    }

    // Generate the basis for the interpolation
    splineA.set_points(P2,AP);
    splineB.set_points(P2,BP);
}


void dress_setting(int np, vdouble & P2, vdouble AP, vdouble BP){


    for(unsigned int i = 0; i < np; i++){

    ofstream aP("ABs/As.txt", ios_base::app | ios_base::out);
        aP << P2[i] << "\t" <<  AP[i] << endl;
    ofstream bP("ABs/Bs.txt", ios_base::app | ios_base::out);
        bP << P2[i] << "\t" <<  BP[i] << endl;
    ofstream mP("ABs/Ms.txt", ios_base::app | ios_base::out);
        mP << P2[i] << "\t" <<  BP[i]/AP[i] << endl;

    }

    // Generate the basis for the interpolation
    splineA.set_points(P2,AP);
    splineB.set_points(P2,BP);
}



// Template to read the file of initial contidions
template<typename T, size_t N>
void initialconditions(T (&arr)[N])
{
    ifstream file("initialconditions.txt");
    if(file.is_open())
    {
        for(unsigned int i = 0; i < N; ++i)
        {
            file >> arr[i];
        }
    }
}



// Assignation of the initial conditions to global variables.
void global_variable_reading(){

    double initcond[14]; // Initial conditions vector

    // Evaluation of the template
    initialconditions(initcond);

    np = initcond[0];
    nk = initcond[1];
    nz = initcond[2];
    lir = initcond[3];
    luv = initcond[4];
    mu = pow(initcond[5],2);
    Nc = initcond[6];
    Nf = initcond[7];
    a0 = initcond[8];
    b0 = initcond[9];
    z2 = initcond[10];
    zm = initcond[11];
    M2 = initcond[12]*initcond[12];
    err = initcond[13];
    mq = b0/a0;
    ncp = 4*np;
}




// Maris-Tandy function with Pauli-Vilars regulator
double alphamt(double x){

    double eta = 1.8;
    double lambda = 0.72;
    double gm = 12.0/(11.0*Nc - 2.0*Nf);
    double lt = 1;
    double lqcd = 0.234;

    double PVregulator = (x/(luv) + 1.0);
    double IRterm = PI * pow(eta , 7) * pow(x/(lambda * lambda), 2) * exp(-pow(eta/lambda , 2) * x);
    double UVterm = PI*gm*(1.0 - exp(- x/(lt * lt)))/log(sqrt(exp(2.0)-1.0+pow(1.0 + x/(lqcd * lqcd),2)));

    return (IRterm + UVterm)/PVregulator;
}
