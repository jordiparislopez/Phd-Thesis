/*
  simpfunctionscomplexDSE.hpp file

  This file contains custom made functions that include complex valued
  parameters.

*/



// Generator of complex parabola, points and jacobian
void comp_mom_gen_0(vcdouble &p2, vcdouble &jac){

     double X = luv;  // UV cutoff to limit the parabola
     double a = 1.0/M2, b = 0, C = -M2/4.0 ; // Parabola values
     double M = sqrt(M2); // Value of the bound state mass
     double CUT = (-b + sqrt(b*b + 4.0*(X - C)*a))/(2.0*a); // Cut points

     // Parametrisation points
     vdouble x1(np),y1(np),t1(np);
     vdouble x2(np),y2(np),t2(np);
     vdouble x3(np),y3(np),t3(np);
     vdouble x4(np),y4(np),t4(np);

     // Quadrature and jacobian points
     vdouble xl(np), wl(np);
     vcdouble jac1(np);
     vcdouble jac2(np);
     vcdouble jac3(np);
     vcdouble jac4(np);

     // dz points
     vcdouble dz1(np);
     vcdouble dz2(np);
     vcdouble dz3(np);
     vcdouble dz4(np);

    legxw(np, xl, wl);

    // Parametrisation C1
    for(unsigned int i = 0; i < np; i++){
      t1[i] = - sqrt(sqrt((X - C)*lir)*pow(sqrt((X - C)/lir),xl[i]));
      x1[i] = t1[i]*t1[i] + C;
      y1[i] = M*t1[i]-b/(2.0*a);
      p2[i] = x1[i] + y1[i]*I;
      dz1[i] = 2.0*t1[i] + I*M;
      jac1[i] = t1[i]/4.0*log((X - C)/lir)*wl[i]*dz1[i];
    }

    // Parametrisation C2
    for(unsigned int i = 0; i < np; i++){
      t2[i] = - sqrt(CUT*lir)*pow(sqrt(CUT/lir),-xl[i]);
      x2[i] = X;
      y2[i] = t2[i];
      p2[i + np] = x2[i] + y2[i]*I;
      dz2[i] = 1.0*I;
      jac2[i] = (-1.0/2.0*t2[i]*log(CUT/lir))*dz2[i]*wl[i];
    }

    // Parametrisation C3
    for(unsigned int i = 0; i < np; i++){
      t3[i] = sqrt(CUT*lir)*pow(sqrt(CUT/lir),xl[i]);
      x3[i] = X;
      y3[i] = t3[i];
      p2[i + 2*np] = x3[i] + y3[i]*I;
      dz3[i] = 1.0*I;
      jac3[i] = (1.0/2.0*t3[i]*log(CUT/lir))*dz3[i]*wl[i];
    }

    // Parametrisation C4
    for(unsigned int i = 0; i < np; i++){
      t4[i] = sqrt(sqrt((X - C)*lir)*pow(sqrt((X - C)/lir),-xl[i]));
      x4[i] = t4[i]*t4[i] + C;
      y4[i] = M*t4[i]-b/(2.0*a);
      p2[i + 3*np] = x4[i] + y4[i]*I;
      dz4[i] = 2.0*t4[i] + M*I;
      jac4[i] = (-t4[i]/4.0*log((X - C)/lir))*dz4[i]*wl[i];
    }



    for(unsigned int i = 0; i < 4*np ; i++){
    if(i < np){jac[i] = jac1[i];}
    if(i > np - 1 && i < 2*np){ jac[i] = jac2[i - np];}
    if(i > 2*np -1 && i < 3*np){jac[i] = jac3[i - 2*np];}
    if(i > 3*np -1 && i < 4*np){jac[i] = jac4[i - 3*np];}
     }

     // Writing path for plot
     for(unsigned int i = 0; i < 4*np;i++){
         ofstream path("ComplexDSE/Path/Path.txt", ios_base::app | ios_base::out);
         path << p2[i] << "\t" << real(p2[i]) << "\t" << imag(p2[i]) << endl;
     }

     // Writing jacobian
     for(unsigned int i = 0; i < 4*np;i++){
         ofstream path1("Complexprop/Path/Jacobian.txt", ios_base::app | ios_base::out);
         path1 << jac[i] << "\t" << real(jac[i]) << "\t" << imag(jac[i]) << endl;
     }
}





//Remove prvious data of complex dressings
void comp_file_remover(){

  removedirectoryrecursively("ComplexDSE/CompDress");

}





// Generation of complex dressings
void comp_dress_gen(vcdouble & P2, vcdouble & jac){


    int npd = P2.size();
    vcdouble AP(npd);
    vcdouble BP(npd);


    for(unsigned int i = 0; i < npd; i++){

        AP[i] = a0;
        BP[i] = b0;

    // Writing values in external file
    ofstream aP("ComplexDSE/CompDress/A.txt", ios_base::app | ios_base::out);
        aP << P2[i] << "\t" <<  AP[i] << endl;
    ofstream bP("ComplexDSE/CompDress/B.txt", ios_base::app | ios_base::out);
        bP << P2[i] << "\t" <<  BP[i] << endl;

    }

    // Generate the 2D dataset for complex interpolation
    complex_parabola_A.set_points(AP,P2,jac);
    complex_parabola_B.set_points(BP,P2,jac);
}





// Set the 2D dataset for complex interpolation
void compdress_setting(const vcdouble & P2, const vcdouble &jac , vcdouble AP, vcdouble BP){

    complex_parabola_A.set_points(AP,P2,jac);
    complex_parabola_B.set_points(BP,P2,jac);
}




// Write complex values in external files
void compdress_print(int ncp, const vcdouble & P2, vcdouble AP, vcdouble BP){

    for(unsigned int i = 0; i < ncp; i++){

      ofstream aP("ComplexDSE/CompDress/Ac.txt", ios_base::app | ios_base::out);
          aP << P2[i] << "\t" <<  AP[i] << endl;
      ofstream bP("ComplexDSE/CompDress/Bc.txt", ios_base::app | ios_base::out);
          bP << P2[i] << "\t" <<  BP[i] << endl;
      ofstream mP("ComplexDSE/CompDress/Mc.txt", ios_base::app | ios_base::out);
          mP << P2[i] << "\t" <<  BP[i]/AP[i] << endl;

    }
}
