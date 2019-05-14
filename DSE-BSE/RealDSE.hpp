/*
	 Real DSE file
*/


void real_prop(const vdouble P2){

	// Assign global variables as initial conditions
	double Z2pre = z2;
	double Zmpre = zm;
  int j1 = 0;

	while( j1 == 0){

		// Computes Z's and compares them to the previous value
		j1 = Zs1(Z2, Zm);

		// Computes the dressings and compares them with previous ones.
		j1 = j1 * realpropagator1(P2);

	}

	cout << "Real propagator solved with Z2 = " << Z2 << " and Zm = " << Zm << " ."<< endl;
}
