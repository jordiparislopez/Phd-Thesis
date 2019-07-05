/*
    Function for initialisation and termination functions
    This file contains all the function involving time.
*/


void initialisation(time_t &stime){

    // stime obtains its initial reference value
    stime = time(NULL);

    // Remove previous data files
    file_remover();
}


void termination(time_t stime, time_t ftime){

    // ftime obtains the final reference time

    ftime = time(NULL);

    // Define precision for the time difference
    cout.precision(0);
  
    cout << endl<<  "Integration process finished in "<< difftime(ftime,stime) << " seconds." << "\n"<< "\n";
}
