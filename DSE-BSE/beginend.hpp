/*
    Function for initialisation and termination functions
    This file contains all the function involving time.
*/


void initialisation(time_t &stime){

    // stime obtains its initial reference value
    stime = time(NULL);

    // Remove previous data files
    file_remover();

    // Call function to read all global variables from the beginning.
    global_variable_reading();

}

void termination(time_t stime, time_t ftime){

    // ftime obtains the final reference time
    ftime = time(NULL);

    // Define precision for the time difference
    cout.precision(6);

    // Print total time.
    cout << endl;
    cout <<  "Total time: "<< difftime(ftime,stime) << " seconds." << "\n";
}
