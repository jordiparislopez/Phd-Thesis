/*
    QMDSE.cc file.
    Main file for the DSE-BSE routine.
*/

// Include headers file
#include "headers.hpp"


int main()
{
  
    // Definition of time variables to measure evaluation time.
    time_t stime, ftime;

    // Initiaisation of the timer.
    initialisation(stime);

    // Calling the DSE function to produce numerical data.
    dse_integration();

    // Calling the BSE function to solve the bound state equations.
    bse_solving();

    // Function to stop the timer.
    termination(stime,ftime);

}
