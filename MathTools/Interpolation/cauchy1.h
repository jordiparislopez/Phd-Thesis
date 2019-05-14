#ifndef CAUCHY_SPLINE_H1
#define CAUCHY_SPLINE_H1

#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include <complex>


// unnamed namespace only because the implementation is in this
// header file and we don't want to export symbols to the obj files


int factorialint(int a){
    if(a == 0){ return 1.0;}
    else if(a == 1){ return 1.0;}
    else{return a*factorialint(a - 1);}
}

std::complex <double> factorial(int a){
    if(a == 0){ return 1.0;}
    else if(a == 1){ return 1.0;}
    else{return 1.0*a*factorialint(a - 1);}
}


namespace
{

  namespace cauchyint1
  {

    class cauchy1
    {

      private:
      std::vector  < std::complex < double > > m_f ;
      std::vector  < std::complex < double > > m_p2 ;
      std::vector  < std::complex < double > > m_jac ;
      int m_n ;
      double m_X;
      double m_M2;

      public:
      void set_points(std::vector  < std::complex < double > > f ,
                      std::vector  < std::complex < double > > x ,
                      std::vector  < std::complex < double > > v ) ;
      std::complex <double>  operator() (std::complex <double> x) const ;
    } ;


    void cauchy1::set_points(std::vector  < std::complex < double > > f ,
                std::vector  < std::complex < double > > x ,
                std::vector  < std::complex < double > > v )
    {
      m_f = f ;
      m_p2 = x ;
      m_jac = v ;
      m_n = m_p2.size() ;
      m_M2 = -real(m_p2[0]);
      m_X = real(m_p2[m_n/2]);
    }

    std::complex <double> cauchy1::operator() (std::complex <double> x) const
    {
        std::complex <double> sum1 = 0 ;
        std::complex <double> sum2 = 0 ;
        std::complex <double> prod = 0 ;
        std::complex <double> res = 0 ;
	std::complex <double> prod1 = 0.0;

    if(pow(imag(x),2) > 4.0*m_M2*(real(x) + m_M2) ||
       real(x) < -m_M2 ||
       real(x) >  m_X ){
           sum1 = 1;
           sum2 = 0;
       }
    else{

	int mm = 32;
	std::vector <std::complex <double>> Iv(mm);
	std::vector <std::complex <double>> Jv(mm);
	std::vector <std::complex <double>> Fv(mm);

	for(int p = 0; p < mm; p++ ){
		prod = 0.0;
		Iv[p] = 0.0;
		Jv[p] = 0.0;
      		  for(unsigned int i = 0 ; i < m_n ; i++){
        	    prod = m_jac[i]/pow(m_p2[i] - 0.0, p + 1) ;
       		    Iv[p] = Iv[p] + prod*m_f[i] ;
       		    Jv[p] = Jv[p] + prod ;
		}
        }

	Fv[0] = Iv[0]/Jv[0];
	sum2 = Fv[0];

	for(unsigned int m = 1; m < mm; m++){

	std::complex <double> sum1 = 0.0;

		for(unsigned int k = 0; k < m; k++){
		sum1 = sum1 + Fv[k]*Jv[m-k]/factorial(k);
		}

	Fv[m] = factorial(m)/Jv[0]*(Iv[m]-sum1);
	sum2 = sum2 + Fv[m]*pow(x - 0.0,m)/factorial(m);
	}



  }


    return sum2 ;
  }


} // namespace cauchyint

} // namespace

#endif /* TK_SPLINE_H */
