/*
    QM.cc file.
    Main file for the FRG routine.
*/

// Include headers file
#include "headers.hpp"

int main()
{
    // Definition of time variables to measure evaluation time.
    time_t stime, ftime;

    // Initialisation of the timer
    initialisation(stime);

    // Calling the FRG function
    rg_integration();

    // Function to stop the timer.
    termination(stime,ftime);

    return 0;
}
