// Computing minimum of potential in every step
double Pot_Min(const state_type &x0){
	double c01 = c0;
	state_type y(9);

	for(unsigned int i = 0 ; i < 9 ; i++){
		y[i] = x0[i]/pow(x0[9 + NL2],i);
	}

	Dbrent dbrent;
	Func pot;
	double res = dbrent.minimize(y, 0 , 3*Rho0*x0[9+NL2] , Rho0*x0[9+NL2] , c01,pot);

	return res;
}


// Function to show parameters on the terminal at running time
void terminal_output(const state_type &x0, const double t, double* PP2){

	double rhomin = Pot_Min(x0);

	cout.precision(3);
	cout << fixed<< t << "\t" << sqrt(2*rhomin) << "\t";
	cout << fixed << x0[1]/x0[9 + NL2] << "\t" << x0[2] << "\t" << x0[3] << "\t" << x0[4] << "\t";
	cout << fixed<< x0[9]/x0[9 + 3*NL2]/sqrt(x0[9 + NL2]) << "\t";
	cout << fixed<< x0[9]/x0[9 + 3*NL2]/sqrt(x0[9 + NL2])*sqrt(2*rhomin)/2.0 << "\t";
	cout.precision(7);
	cout << fixed<< x0[9 + NL2] << "\t";
	cout << fixed<< x0[9 + 2*NL2] << "\t";
	cout << fixed<< x0[9 + 3*NL2] << "\t";
	cout.precision(15);
	cout << fixed<< x0[9 + 4*NL2] << "\t";
	cout << "\n";

}

// Function to write terms internal and externally
void file_output(const state_type &x0, const double t, double* PP2);
void fileterms_output(const state_type &x0, const double t, double* PP2)
void step_output(const state_type &x0, const double t, double* PP2){

	terminal_output(x0,t, PP2);
	fileterms_output(x0,t, PP2);

}

// Function necessary from boost
struct write_cout
{
	double* PP2;
	write_cout(double* PPP2) :  PP2 ( PPP2 ) {}
	void operator()( const state_type &x0 , const double t)
	{

	step_output(x0, t, PP2);

	}
};



// Function to print the values to values.txt
void fileterms_output(const state_type &x0, const double t, double* PP2){

	string v;
	v.append("Data/values.txt");
	const char * V = v.c_str();
	ofstream log(V, ios_base::app | ios_base::out);

	log <<  t << "\t";
		for(unsigned int i = 0; i < par; i++)
		{
			log << x0[i] << "\t";
		}

	log << "\n";
}


// Functoin to write terms in file once the system converges
void file_output(double *PP2)
{
    double t;
		double mp, ms, mq, hk, rhomin;
    string line;
    string m,yuk,pion,quark,sigma,dha,hq1,fp1,fs1,zf1;
    m.append("Data/Massdata/massdata.txt");
    yuk.append("Data/Yukawadata/yukawadata.txt");
    pion.append("Data/Piondata/piondata.txt");
    quark.append("Data/Quarkmass/quarkmass.txt");
    sigma.append("Data/Sigmadata/sigmadata.txt");
    dha.append("Data/dynamical/dynamical.txt");
    hq1.append("Dressings/HQ.txt");
    fp1.append("Dressings/FP.txt");
    fs1.append("Dressings/FS.txt");
    zf1.append("Dressings/AQ.txt");

    const char * M = m.c_str();
    const char * YUK = yuk.c_str();
    const char * PION = pion.c_str();
    const char * SIGMA = sigma.c_str();
    const char * QUARK = quark.c_str();
    const char * DH = dha.c_str();
    const char * HQ = hq1.c_str();
    const char * FP = fp1.c_str();
    const char * FS = fs1.c_str();
    const char * ZF = zf1.c_str();

    state_type x0(par);
    ifstream myfile("Data/values.txt");

    while(getline(myfile,line)){
        if(myfile.is_open())
        {
            myfile >> t;
            for(unsigned int i1 = 0; i1 < par; ++i1)
            {
                myfile >> x0[i1];
            }
        }

		rhomin = Pot_Min(x0);
		mp = sqrt(x0[1]/x0[9 + NL2]);
		ms = sqrt(x0[1]/x0[9 + NL2]+ 2.0*rhomin*x0[2]/pow(x0[9 + NL2],2));
		hk = x0[9]/(sqrt(x0[9 + NL2])*x0[9 + 3*NL2]);
		mq = hk*sqrt(2*rhomin)/2.0;


		for(unsigned int i = 0; i < NL2; i++){

		ofstream log(M, ios_base::app | ios_base::out);
		log <<  t << "\t" << PP2[i]
			<< "\t" << sqrt(x0[1]) << "\t"
			<< sqrt(x0[1]+ 2.0*rhomin*x0[2]/pow(x0[9 + NL2 + i],1))
			<< "\t"	<< sqrt(x0[1]/x0[9 + NL2 + i])
			<< "\t" <<  sqrt(x0[1]/x0[9 + NL2 + i]+ 2.0*rhomin*x0[2]/pow(x0[9 + NL2 + i],2))  << "\n";

		ofstream yukawa(YUK, ios_base::app | ios_base::out);
		yukawa <<  t << "\t" << PP2[i] << "\t" << x0[9 + i] << "\t" << x0[9 + i]/(sqrt(x0[9 + NL2 + i])*x0[9 + 2*NL2 +i]) << "\n";

		ofstream pion1(PION, ios_base::app | ios_base::out);
		pion1 <<  t << "\t" << PP2[i] << "\t" << PP2[i] + x0[1] << "\t" << x0[9 + NL2 +i]*PP2[i] + x0[1]  << "\n";

		ofstream sigma1(SIGMA, ios_base::app | ios_base::out);
		sigma1 <<  t << "\t" << PP2[i] << "\t" << PP2[i] + x0[1] + 2.0*rhomin*x0[2] << "\t" << x0[9 + 2*NL2 + i]*PP2[i] + x0[1] + 2.0*rhomin*x0[2] << "\n";

		ofstream quark1(QUARK, ios_base::app | ios_base::out);
		quark1 <<  t << "\t" << PP2[i] << "\t" << x0[9 + 3*NL2 + i] << "\t" << x0[9 + i]*sqrt(2*rhomin)/2.0
		<< "\t" << x0[9 + i]/(sqrt(x0[9 + NL2 + i])*x0[9 + 3*NL2 + i])*sqrt(2*rhomin)/2.0 << "\n";

		ofstream fp(FP, ios_base::app | ios_base::out);
		fp << t << "\t" << PP2[i] << "\t" << x0[9 + NL2 + i] << "\n";

		ofstream fs(FS, ios_base::app | ios_base::out);
		fs << t << "\t" << PP2[i] << "\t" << x0[9 + 2*NL2 + i] << "\n";

		ofstream zf(ZF, ios_base::app | ios_base::out);
		zf << t << "\t" << PP2[i] << "\t" << x0[9 + 3*NL2 + i] << "\n";

		ofstream dymamical1(DH, ios_base::app | ios_base::out);
		dymamical1 << t << "\t" << PP2[i] << "\t" << x0[9 + 4*NL2 + i] << "\n";

		}

}
	cout.precision(3);

	cout << endl << "The Pion mass is " << mp << " MeV." << endl;
	cout << "The Sigma meson mass is " << ms << " MeV." << endl;
	cout << "The Quark mass is " << mq << " MeV." << endl;
	cout << "The Yukawa coupling is " << hk << " MeV." << endl;
	cout << "The Minimum of the potential " << sqrt(2*rhomin) << " MeV." << endl;
}
