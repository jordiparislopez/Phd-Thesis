#ifndef CAUCHY_SPLINE_H
#define CAUCHY_SPLINE_H

#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include <complex>


// unnamed namespace only because the implementation is in this
// header file and we don't want to export symbols to the obj files
namespace
{

  namespace cauchyint
  {

    class cauchy
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


    void cauchy::set_points(std::vector  < std::complex < double > > f ,
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

    std::complex <double> cauchy::operator() (std::complex <double> x) const
    {
        std::complex <double> sum1 = 0 ;
        std::complex <double> sum2 = 0 ;
        std::complex <double> prod = 0 ;
        std::complex <double> res = 0 ;

    if(pow(imag(x),2) > 4.0*m_M2*(real(x) + m_M2) ||
       real(x) < -m_M2 ||
       real(x) >  m_X ){
           sum1 = 1;
           sum2 = 0;
       }
    else{
        for(unsigned int i = 0 ; i < m_n ; i++){
            prod = m_jac[i]/(m_p2[i] - x) ;
            sum1 = sum1 + prod ;
            sum2 = sum2 + prod*m_f[i] ;
        }
    }

    res = sum2/sum1 ;

    return res ;
  }


} // namespace cauchyint

} // namespace

#endif /* TK_SPLINE_H */
