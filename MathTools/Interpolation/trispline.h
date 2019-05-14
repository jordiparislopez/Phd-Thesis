/*
 * spline.h
 *
 * simple cubic spline interpolation library without external
 * dependencies
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2011, 2014 Tino Kluge (ttk448 at gmail.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 */


#ifndef TTK_BSPLINE_H
#define TTK_BSPLINE_H

#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>


// unnamed namespace only because the implementation is in this
// header file and we don't want to export symbols to the obj files
namespace
{

namespace tritk
{

// band matrix solver
class triband_matrix
{
    private:
    std::vector< std::vector<double> > m_upper;  // upper band
    std::vector< std::vector<double> > m_lower;  // lower band
    public:
    triband_matrix() {};                             // constructor
    triband_matrix(int dim, int n_u, int n_l);       // constructor
    ~triband_matrix() {};                            // destructor
    void resize(int dim, int n_u, int n_l);      // init with dim,n_u,n_l
    int dim() const;                             // matrix dimension
    int num_upper() const
    {
        return m_upper.size()-1;
    }
    int num_lower() const
    {
        return m_lower.size()-1;
    }
    // access operator
    double & operator () (int i, int j);            // write
    double   operator () (int i, int j) const;      // read
    // we can store an additional diogonal (in m_lower)
    double& saved_diag(int i);
    double  saved_diag(int i) const;
    void lu_decompose();
    std::vector<double> r_solve(const std::vector<double>& b) const;
    std::vector<double> l_solve(const std::vector<double>& b) const;
    std::vector<double> lu_solve(const std::vector<double>& b,
                                 bool is_lu_decomposed=false);

};


// spline interpolation
class trispline
{
public:
    enum bd_type {
        first_deriv = 1,
        second_deriv = 2
    };

private:
    std::vector<double> m_x,m_y,m_z,m_w;
    std::vector<double> m_h, m_r, m_F, m_chev;            // x,y coordinates of points
    std::vector<double> m_h1, m_s, m_F1, m_chev1;
    // interpolation parameters
    // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
    std::vector<double> m_a,m_b,m_c;        // spline coefficients
    //std::vector<double> m_b0, m_c0;                     // for left extrapol
    bd_type m_left, m_right;
    double  m_left_value, m_right_value;
    bool    m_force_linear_extrapolation;
    int     m_n;
    int     m_Tord,m_Tord1;
    double  m_IVp, m_IVd, m_IR, m_UV;

    std::vector<double> m_x1, m_y1;            // x,y coordinates of points
    std::vector<double> m_a1, m_b1, m_c1;        // spline coefficients
    double  m_b0, m_c0;                     // for left extrapol
    double m_sum, m_sum1,m_sum0, m_sum01,m_sum00, m_sum02;

    std::vector<double> m_temp;

public:
    // set default boundary condition to be zero curvature at both ends
    trispline(): m_left(second_deriv), m_right(second_deriv),
        m_left_value(0.0), m_right_value(0.0),
        m_force_linear_extrapolation(false)
    {
        ;
    }

    // optional, but if called it has to come be before set_points()
    void set_boundary(bd_type left, double left_value,
                      bd_type right, double right_value,
                      bool force_linear_extrapolation=false);
    void set_points(const std::vector<double>& x,
                    const std::vector<double>& y,
                    const std::vector<double>& w,
                    const std::vector<double>& z,
                    const double UV,
                    const double IR,
                    bool cubic_spline=true);
    double operator() (double x, double y, double z) const;
};



// ---------------------------------------------------------------------
// implementation part, which could be separated into a cpp file
// ---------------------------------------------------------------------


// triband_matrix implementation
// -------------------------

triband_matrix::triband_matrix(int dim, int n_u, int n_l)
{
    resize(dim, n_u, n_l);
}
void triband_matrix::resize(int dim, int n_u, int n_l)
{
    assert(dim>0);
    assert(n_u>=0);
    assert(n_l>=0);
    m_upper.resize(n_u+1);
    m_lower.resize(n_l+1);
    for(size_t i=0; i<m_upper.size(); i++) {
        m_upper[i].resize(dim);
    }
    for(size_t i=0; i<m_lower.size(); i++) {
        m_lower[i].resize(dim);
    }
}
int triband_matrix::dim() const
{
    if(m_upper.size()>0) {
        return m_upper[0].size();
    } else {
        return 0;
    }
}


// defines the new operator (), so that we can access the elements
// by A(i,j), index going from i=0,...,dim()-1
double & triband_matrix::operator () (int i, int j)
{
    int k=j-i;       // what band is the entry
    assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
    assert( (-num_lower()<=k) && (k<=num_upper()) );
    // k=0 -> diogonal, k<0 lower left part, k>0 upper right part
    if(k>=0)   return m_upper[k][i];
    else	    return m_lower[-k][i];
}
double triband_matrix::operator () (int i, int j) const
{
    int k=j-i;       // what band is the entry
    assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
    assert( (-num_lower()<=k) && (k<=num_upper()) );
    // k=0 -> diogonal, k<0 lower left part, k>0 upper right part
    if(k>=0)   return m_upper[k][i];
    else	    return m_lower[-k][i];
}
// second diag (used in LU decomposition), saved in m_lower
double triband_matrix::saved_diag(int i) const
{
    assert( (i>=0) && (i<dim()) );
    return m_lower[0][i];
}
double & triband_matrix::saved_diag(int i)
{
    assert( (i>=0) && (i<dim()) );
    return m_lower[0][i];
}

// LR-Decomposition of a band matrix
void triband_matrix::lu_decompose()
{
    int  i_max,j_max;
    int  j_min;
    double x;

    // preconditioning
    // normalize column i so that a_ii=1
    for(int i=0; i<this->dim(); i++) {
        assert(this->operator()(i,i)!=0.0);
        this->saved_diag(i)=1.0/this->operator()(i,i);
        j_min=std::max(0,i-this->num_lower());
        j_max=std::min(this->dim()-1,i+this->num_upper());
        for(int j=j_min; j<=j_max; j++) {
            this->operator()(i,j) *= this->saved_diag(i);
        }
        this->operator()(i,i)=1.0;          // prevents rounding errors
    }

    // Gauss LR-Decomposition
    for(int k=0; k<this->dim(); k++) {
        i_max=std::min(this->dim()-1,k+this->num_lower());  // num_lower not a mistake!
        for(int i=k+1; i<=i_max; i++) {
            assert(this->operator()(k,k)!=0.0);
            x=-this->operator()(i,k)/this->operator()(k,k);
            this->operator()(i,k)=-x;                         // assembly part of L
            j_max=std::min(this->dim()-1,k+this->num_upper());
            for(int j=k+1; j<=j_max; j++) {
                // assembly part of R
                this->operator()(i,j)=this->operator()(i,j)+x*this->operator()(k,j);
            }
        }
    }
}
// solves Ly=b
std::vector<double> triband_matrix::l_solve(const std::vector<double>& b) const
{
    assert( this->dim()==(int)b.size() );
    std::vector<double> x(this->dim());
    int j_start;
    double sum;
    for(int i=0; i<this->dim(); i++) {
        sum=0;
        j_start=std::max(0,i-this->num_lower());
        for(int j=j_start; j<i; j++) sum += this->operator()(i,j)*x[j];
        x[i]=(b[i]*this->saved_diag(i)) - sum;
    }
    return x;
}
// solves Rx=y
std::vector<double> triband_matrix::r_solve(const std::vector<double>& b) const
{
    assert( this->dim()==(int)b.size() );
    std::vector<double> x(this->dim());
    int j_stop;
    double sum;
    for(int i=this->dim()-1; i>=0; i--) {
        sum=0;
        j_stop=std::min(this->dim()-1,i+this->num_upper());
        for(int j=i+1; j<=j_stop; j++) sum += this->operator()(i,j)*x[j];
        x[i]=( b[i] - sum ) / this->operator()(i,i);
    }
    return x;
}

std::vector<double> triband_matrix::lu_solve(const std::vector<double>& b,
        bool is_lu_decomposed)
{
    assert( this->dim()==(int)b.size() );
    std::vector<double>  x,y;
    if(is_lu_decomposed==false) {
        this->lu_decompose();
    }
    y=this->l_solve(b);
    x=this->r_solve(y);
    return x;
}







// spline implementation
// -----------------------

void trispline::set_boundary(trispline::bd_type left, double left_value,
                          trispline::bd_type right, double right_value,
                          bool force_linear_extrapolation)
{
    assert(m_x.size()==0);          // set_points() must not have happened yet
    m_left=left;
    m_right=right;
    m_left_value=left_value;
    m_right_value=right_value;
    m_force_linear_extrapolation=force_linear_extrapolation;
}

void trispline::set_points(const std::vector<double>& x, const std::vector<double>& y,
                        const std::vector<double>& w, const std::vector<double>& z, const double UV, const double IR, bool cubic_spline)
{
    assert(x.size()*y.size()*w.size()==z.size());

    m_w=w;
    m_x=x;
    m_y=y;
    m_z=z;
    m_IVp=log(sqrt(UV*IR));
    m_IVd=log(sqrt(UV/IR));
    m_IR = IR;
    m_UV = UV;

    int n = x.size();
    int nq = y.size();
    int nz = w.size();

    unsigned int m_Tord = 10;
    unsigned int m_Tord1 = 10;

    m_r.resize(nz);
    m_s.resize(nq);
    m_h1.resize(n);
    m_chev.resize(m_Tord + 1);
    m_chev1.resize(m_Tord1 + 1);

    for(unsigned int k = 0; k < nz ; k++){
        m_sum0  = 0;
        m_sum00 = 0;
        for(unsigned int i = 0; i < n ; i++){
            for(unsigned int j = 0; j < nq ; j++){
                m_sum0  = m_sum0  + z[ i + j*n + k*n*nq];
                m_sum00 = m_sum00 + z[ i + j*n + (nz/2 + 1)*n*nq];
            }
        }
        m_r[k] = m_sum0/m_sum00;
    }

    for(unsigned int iord = 0; iord <= m_Tord ; iord++){
        m_sum = 0;
        for(unsigned int k = 0; k < nz; k++){
            m_sum = m_sum + cos(iord*acos(m_w[k]))*m_r[k];
        }
        m_chev[iord] = (2.0/nz)*m_sum;
    }



    for(unsigned int j = 0; j < nq ; j++){
        m_sum01 = 0;
        m_sum02 = 0;
        for(unsigned int i = 0; i < n ; i++){
            m_sum01 = m_sum01 + z[ i + j*n + (nz/2 + 1)*n*nq];
            m_sum02 = m_sum02 + z[ i + 0*n + (nz/2 + 1)*n*nq];
        }
        m_s[j] = m_sum0/m_sum00;
    }

    for(unsigned int iord = 0; iord <= m_Tord1 ; iord++){
        m_sum1 = 0;
        for(unsigned int j = 0; j < nq; j++){
            m_sum1 = m_sum1 + cos(iord*acos(m_y[j]))*m_s[j];
        }
        m_chev1[iord] = (2.0/nq)*m_sum1;
    }

    for(unsigned int i = 0; i < n ; i++){
        m_h1[i] = z[ i + 0*n + (nz/2 + 1)*n*nq];
    }



    m_a.resize(n);
    m_b.resize(n);
    m_c.resize(n);

            if(cubic_spline==true) { // cubic spline interpolation
            // setting up the matrix and right hand side of the equation system
            // for the parameters b[]
            triband_matrix A(n,1,1);
            std::vector<double>  rhs(n);
            for(int i=1; i<n-1; i++) {
                A(i,i-1)=1.0/3.0*(x[i]-x[i-1]);
                A(i,i)=2.0/3.0*(x[i+1]-x[i-1]);
                A(i,i+1)=1.0/3.0*(x[i+1]-x[i]);
                rhs[i]=(m_h1[(i+1)]-m_h1[i])/(x[i+1]-x[i]) - (m_h1[i]-m_h1[(i-1)])/(x[i]-x[i-1]);
            }
            // boundary conditions
            if(m_left == trispline::second_deriv) {
                // 2*b[0] = f''
                A(0,0)=2.0;
                A(0,1)=0.0;
                rhs[0]=m_left_value;
            } else if(m_left == trispline::first_deriv) {
                // c[0] = f', needs to be re-expressed in terms of b:
                // (2b[0]+b[1])(x[1]-x[0]) = 3 ((y[1]-y[0])/(x[1]-x[0]) - f')
                A(0,0)=2.0*(x[1]-x[0]);
                A(0,1)=1.0*(x[1]-x[0]);
                rhs[0]=3.0*((m_h1[1]-m_h1[0])/(x[1]-x[0])-m_left_value);
            } else {
                assert(false);
            }
            if(m_right == trispline::second_deriv) {
                // 2*b[n-1] = f''
                A(n-1,n-1)=2.0;
                A(n-1,n-2)=0.0;
                rhs[n-1]=m_right_value;
            } else if(m_right == trispline::first_deriv) {
                // c[n-1] = f', needs to be re-expressed in terms of b:
                // (b[n-2]+2b[n-1])(x[n-1]-x[n-2])
                // = 3 (f' - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
                A(n-1,n-1)=2.0*(x[n-1]-x[n-2]);
                A(n-1,n-2)=1.0*(x[n-1]-x[n-2]);
                rhs[n-1]=3.0*(m_right_value-(m_h1[n-1]-m_h1[n-2])/(x[n-1]-x[n-2]));
            } else {
                assert(false);
            }

            // solve the equation system to obtain the parameters b[]
            m_temp = A.lu_solve(rhs);

            for(int i=0; i<n; i++) {
            m_b[i]=m_temp[i];
            }
            // calculate parameters a[] and c[] based on b[]


            for(int i=0; i<n-1; i++) {
                m_a[i]=1.0/3.0*(m_b[i+1]-m_b[i])/(x[i+1]-x[i]);
                m_c[i]=(m_h1[i+1]-m_h1[i])/(x[i+1]-x[i])
                       - 1.0/3.0*(2.0*m_b[i]+m_b[i+1])*(x[i+1]-x[i]);
            }
        } else { // linear interpolation

            for(int i=0; i<n-1; i++) {
                m_a[i]=0.0;
                m_b[i]=0.0;
                m_c[i]=(m_h1[i+1]-m_h1[i])/(m_x[i+1]-m_x[i]);
            }
        }

        // for left extrapolation coefficients
        m_b0 = (m_force_linear_extrapolation==false) ? m_b[0] : 0.0;
        m_c0 = m_c[0];

        // for the right extrapolation coefficients
        // f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
        double h=x[n-1]-x[n-2];
        // m_b[n-1] is determined by the boundary condition
        m_a[n-1]=0.0;
        m_c[n-1]=3.0*m_a[n-2]*h*h+2.0*m_b[n-2]*h+m_c[n-2];   // = f'_{n-2}(x_{n-1})
        if(m_force_linear_extrapolation==true)
            m_b[n-1]=0.0;


}

double trispline::operator() (double x, double y, double z) const
{
    size_t n = m_x.size();

    // find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
    std::vector<double>::const_iterator it;
    it = std::lower_bound(m_x.begin(),m_x.end(),x);
    int idx = std::max( int(it-m_x.begin())-1, 0);

    double h = x-m_x[idx];
    double interpol;
            if(x < m_x[0]) {
                // extrapolation to the left
                interpol = (m_b0*h + m_c0)*h + m_h1[0];
            } else if(x > m_x[n-1]) {
                // extrapolation to the right
                interpol = (m_b[(n-1)]*h + m_c[(n-1)])*h + m_h1[(n-1)];
            } else {
                // interpolation
                interpol = ((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_h1[idx];
            }

    double qq = log(y);
    if(y==0){qq = log(m_IR);}
    double zq = (qq-m_IVp)/m_IVd;

    double sum0 = 0;
    if(zq >  1.0 ){

        for(unsigned int iord = 1; iord <= m_Tord1; iord++){
            sum0 = sum0 + m_chev1[iord]*( 1.0 + iord*iord*(zq - 1.0));
            }
        }

    if(zq <  -1.0 ){

        for(unsigned int iord = 1; iord <= m_Tord1; iord++){
            sum0 = sum0 + m_chev1[iord]*( 1.0 + iord*iord*(zq + 1.0));
            }
        }
    else{

        for(unsigned int iord = 1; iord <= m_Tord1; iord++){
            sum0 = sum0 + m_chev1[iord]*cos(iord*acos(zq));
            }
    }


    double G = interpol*(1.0/2.0*m_chev1[0] + sum0);

    double sum1 = 0;
    for(unsigned int iord = 1; iord <= m_Tord ; iord++){
        sum1 = sum1 + m_chev[iord]*cos(iord*acos(z));
    }

    double H = G*(1.0/2.0*m_chev[0] + sum1);

    return H;
}


} // namespace tk


} // namespace

#endif
