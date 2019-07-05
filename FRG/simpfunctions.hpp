/*

  List of relevant functions

*/



// Generator of Gauss-Legendre weights and abscissas
void legxw(int n, double *xleg, double *wleg){

    double x[n],w[n],error=3.0e-15,p1,p2,p3,pp,z,z1;
    int m = (n+1)/2;

    for(unsigned int i = 1; i <= m; i++){

        z1 = 2.0;
        z = cos(PI*(i-0.25)/(n+0.5));

        while(abs((z-z1)) > error){

            p1 = 1.0;
            p2 = 0.0;

            for(unsigned int j = 1; j <= n; j++){

                p3 = p2;
                p2 = p1;
                p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }

            pp = n*(z*p1-p2)/(z*z-1.0);
            z1 = z;
            z = z1-p1*1.0/pp;
        }

        x[i-1] = -z;
        x[n-i] = z;
        w[i-1] = 2.0/((1.0-z*z)*pp*pp);
        w[n-i] = w[i-1];
    }

    for(unsigned int i = 0; i < n; i++){
        if(abs(x[i]) < error){    xleg[i] = 0;}
        else{xleg[i]=x[i];};
        wleg[i]=w[i];
    }

}



// Generator of gauss-chebishev (2n order) weights and abscissas
void cheb2xw(int n, long double *xcheb2, long double *wcheb2){
    for(unsigned int i = 1; i<= n; i++){
        xcheb2[i-1] = cos((i*PI*1.0/(n+1)));
        wcheb2[i-1] = PI/(n+1)*sin(i*PI*1.0/(n+1))*sin(i*PI/(n+1));
    }
}



// Function to remove previous files
void file_remover(){

    removedirectoryrecursively("Data/Quarkmass");
    removedirectoryrecursively("Data/Sigmadata");
    removedirectoryrecursively("Data/Yukawadata");
    removedirectoryrecursively("Data/Piondata");
    removedirectoryrecursively("Data/Massdata");
    removedirectoryrecursively("Data/k0data/yukawadatak0");
    removedirectoryrecursively("Data/k0data/sigmapropk0");
    removedirectoryrecursively("Data/k0data/quarkpropk0");
    removedirectoryrecursively("Data/k0data/pionpropk0");
    removedirectoryrecursively("Data/k0data/massdatak0");
    removedirectoryrecursively("Data/k0data/dynamicalhadronization");
    removedirectoryrecursively("Data/dynamical");
    remove("Dressings/FP.txt");
    remove("Dressings/FS.txt");
    remove("Dressings/AQ.txt");
    remove("Dressings/HQ.txt");
    remove("Data/values.txt");
}

// Function to remove dressing files
void dressing_remover(){
    remove("Dressings/FP.txt");
    remove("Dressings/FS.txt");
    remove("Dressings/AQ.txt");
    remove("Dressings/HQ.txt");
}



// Function to generate momentum grid with logarithmic distribution
void mom_gen(double *P2){

    double xl2[NL2], wl[NL2];
    legxw(NL2,xl2,wl);

    for(unsigned int i = 0; i < NL2; i++){
        P2[i] = sqrt(luv*lir)*pow(sqrt(luv/lir),xl2[i]);
    }
}


// Template to read initial conditions into an array
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

// Template to read RG-conditions into an array
template<typename T, size_t N>
void initialrgvalues(T (&arr)[N])
{
    ifstream filerg("rgconditions.txt");
    if(filerg.is_open())
    {
        for(unsigned int i = 0; i < N; ++i)
        {
            filerg >> arr[i];
        }
    }
}

// Functions for propagators
double Meson_Propagator(double Z, double p2, double mass2, double reg){
	return pow(Z*p2 + mass2 + reg,-1);
}

double Fermion_Propagator(double Z, double p2, double yukawa, double reg){
	return pow(Nf*p2*(Z + reg)*(Z + reg) + yukawa*yukawa*Rho0,-1);
}

// Function for 1D - Yukawa
double YK(double q2, double p2 , double z){
    double rmom = (q2 + p2 - 2.0*sqrt(q2*p2)*z);
    double tmom = (q2 + p2 + 2.0*sqrt(q2*p2)*z);
	return splineH(rmom/4. + tmom);
}
