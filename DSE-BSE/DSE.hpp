/*
  DSE integration file

  This file determines which routine to use, realDSE, complexDSE or both.
  Comment the unwanted routines
*/

void dse_integration()
{

    // Calls real DSE
    dse_real();

    // Calls complex DSE
    dse_complex();
}
