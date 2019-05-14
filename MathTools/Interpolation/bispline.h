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


#ifndef BTK_BSPLINE_H
#define BTK_BSPLINE_H

#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>


// unnamed namespace only because the implementation is in this
// header file and we don't want to export symbols to the obj files
namespace
{

namespace bitk
{

// band matrix solver
class biband_matrix
{
    private:
    std::vector< std::vector<double> > m_upper;  // upper band
    std::vector< std::vector<double> > m_lower;  // lower band
    public:
    biband_matrix() {};                             // constructor
    biband_matrix(int dim, int n_u, int n_l);       // constructor
    ~biband_matrix() {};                            // destructor
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
class bispline
{
public:
    enum bd_type {
        first_deriv = 1,
        second_deriv = 2
    };

private:
    std::vector<double> m_x,m_y,m_z;            // x,y coordinates of points
    // interpolation parameters
    // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
    std::vector<double> m_a,m_b,m_c;        // spline coefficients
    std::vector<double> m_b0, m_c0;                     // for left extrapol
    bd_type m_left, m_right;
    double  m_left_value, m_right_value;
    bool    m_force_linear_extrapolation;
    int     m_n;

    std::vector<double> m_x1,m_y1;            // x,y coordinates of points
    std::vector<double> m_a1,m_b1,m_c1;        // spline coefficients
    double  m_b01, m_c01;                     // for left extrapol


    std::vector<double> m_temp;

public:
    // set default boundary condition to be zero curvature at both ends
    bispline(): m_left(second_deriv), m_right(second_deriv),
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
                    const std::vector<double>& z,
                    int n,
                    bool cubic_spline=true);
    double operator() (double x, double y) const;
};



// ---------------------------------------------------------------------
// implementation part, which could be separated into a cpp file
// ---------------------------------------------------------------------


// biband_matrix implementation
// -------------------------

biband_matrix::biband_matrix(int dim, int n_u, int n_l)
{
    resize(dim, n_u, n_l);
}
void biband_matrix::resize(int dim, int n_u, int n_l)
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
int biband_matrix::dim() const
{
    if(m_upper.size()>0) {
        return m_upper[0].size();
    } else {
        return 0;
    }
}


// defines the new operator (), so that we can access the elements
// by A(i,j), index going from i=0,...,dim()-1
double & biband_matrix::operator () (int i, int j)
{
    int k=j-i;       // what band is the entry
    assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
    assert( (-num_lower()<=k) && (k<=num_upper()) );
    // k=0 -> diogonal, k<0 lower left part, k>0 upper right part
    if(k>=0)   return m_upper[k][i];
    else	    return m_lower[-k][i];
}
double biband_matrix::operator () (int i, int j) const
{
    int k=j-i;       // what band is the entry
    assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
    assert( (-num_lower()<=k) && (k<=num_upper()) );
    // k=0 -> diogonal, k<0 lower left part, k>0 upper right part
    if(k>=0)   return m_upper[k][i];
    else	    return m_lower[-k][i];
}
// second diag (used in LU decomposition), saved in m_lower
double biband_matrix::saved_diag(int i) const
{
    assert( (i>=0) && (i<dim()) );
    return m_lower[0][i];
}
double & biband_matrix::saved_diag(int i)
{
    assert( (i>=0) && (i<dim()) );
    return m_lower[0][i];
}

// LR-Decomposition of a band matrix
void biband_matrix::lu_decompose()
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
std::vector<double> biband_matrix::l_solve(const std::vector<double>& b) const
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
std::vector<double> biband_matrix::r_solve(const std::vector<double>& b) const
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

std::vector<double> biband_matrix::lu_solve(const std::vector<double>& b,
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

void bispline::set_boundary(bispline::bd_type left, double left_value,
                          bispline::bd_type right, double right_value,
                          bool force_linear_extrapolation)
{
    assert(m_x.size()==0);          // set_points() must not have happened yet
    m_left=left;
    m_right=right;
    m_left_value=left_value;
    m_right_value=right_value;
    m_force_linear_extrapolation=force_linear_extrapolation;
}

void bispline::set_points(const std::vector<double>& x, const std::vector<double>& y,
                        const std::vector<double>& z, int n, bool cubic_spline)
{
    assert(x.size()*y.size()==z.size());
    m_x=x;
    m_y=y;
    m_z=z;

    m_a.resize(n*n);
    m_b.resize(n*n);
    m_c.resize(n*n);
    m_b0.resize(n);
    m_c0.resize(n);
    for(unsigned int j = 0; j < n; j++){
        if(cubic_spline==true) { // cubic spline interpolation
            // setting up the matrix and right hand side of the equation system
            // for the parameters b[]
            biband_matrix A(n,1,1);
            std::vector<double>  rhs(n);
            for(int i=1; i<n-1; i++) {
                A(i,i-1)=1.0/3.0*(x[i]-x[i-1]);
                A(i,i)=2.0/3.0*(x[i+1]-x[i-1]);
                A(i,i+1)=1.0/3.0*(x[i+1]-x[i]);
                rhs[i]=(z[(i+1)+n*j]-z[i+n*j])/(x[i+1]-x[i]) - (z[i+n*j]-z[(i-1)+n*j])/(x[i]-x[i-1]);
            }
            // boundary conditions
            if(m_left == bispline::second_deriv) {
                // 2*b[0] = f''
                A(0,0)=2.0;
                A(0,1)=0.0;
                rhs[0]=m_left_value;
            } else if(m_left == bispline::first_deriv) {
                // c[0] = f', needs to be re-expressed in terms of b:
                // (2b[0]+b[1])(x[1]-x[0]) = 3 ((y[1]-y[0])/(x[1]-x[0]) - f')
                A(0,0)=2.0*(x[1]-x[0]);
                A(0,1)=1.0*(x[1]-x[0]);
                rhs[0]=3.0*((z[1+j*n]-z[j*n])/(x[1]-x[0])-m_left_value);
            } else {
                assert(false);
            }
            if(m_right == bispline::second_deriv) {
                // 2*b[n-1] = f''
                A(n-1,n-1)=2.0;
                A(n-1,n-2)=0.0;
                rhs[n-1]=m_right_value;
            } else if(m_right == bispline::first_deriv) {
                // c[n-1] = f', needs to be re-expressed in terms of b:
                // (b[n-2]+2b[n-1])(x[n-1]-x[n-2])
                // = 3 (f' - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
                A(n-1,n-1)=2.0*(x[n-1]-x[n-2]);
                A(n-1,n-2)=1.0*(x[n-1]-x[n-2]);
                rhs[n-1]=3.0*(m_right_value-(z[n-1+n*j]-z[n-2+n*j])/(x[n-1]-x[n-2]));
            } else {
                assert(false);
            }

            // solve the equation system to obtain the parameters b[]
            m_temp = A.lu_solve(rhs);

            for(int i=0; i<n; i++) {
            m_b[i+n*j]=m_temp[i];
            }
            // calculate parameters a[] and c[] based on b[]


            for(int i=0; i<n-1; i++) {
                m_a[i+n*j]=1.0/3.0*(m_b[i+1+n*j]-m_b[i+n*j])/(x[i+1]-x[i]);
                m_c[i+n*j]=(z[i+1+n*j]-z[i+n*j])/(x[i+1]-x[i])
                       - 1.0/3.0*(2.0*m_b[i+n*j]+m_b[i+1+n*j])*(x[i+1]-x[i]);
            }
        } else { // linear interpolation

            for(int i=0; i<n-1; i++) {
                m_a[i+n*j]=0.0;
                m_b[i+n*j]=0.0;
                m_c[i+n*j]=(m_z[i+1+n*j]-m_z[i+n*j])/(m_x[i+1]-m_x[i]);
            }
        }

        // for left extrapolation coefficients
        m_b0[j] = (m_force_linear_extrapolation==false) ? m_b[n*j] : 0.0;
        m_c0[j] = m_c[n*j];

        // for the right extrapolation coefficients
        // f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
        double h=x[n-1]-x[n-2];
        // m_b[n-1] is determined by the boundary condition
        m_a[n-1+n*j]=0.0;
        m_c[n-1+n*j]=3.0*m_a[n-2+n*j]*h*h+2.0*m_b[n-2+n*j]*h+m_c[n-2+n*j];   // = f'_{n-2}(x_{n-1})
        if(m_force_linear_extrapolation==true)
            m_b[n-1+n*j]=0.0;
    }
}


double bispline::operator() (double x, double y) const
{
    size_t n=m_x.size();

    // find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
    std::vector<double>::const_iterator it;
    it=std::lower_bound(m_x.begin(),m_x.end(),x);
    int idx=std::max( int(it-m_x.begin())-1, 0);

    double h=x-m_x[idx];
    double interpol1;
    std::vector<double> interpol(n);

    for(unsigned int j = 0; j < n; j++){
        if(x<m_x[0]) {
            // extrapolation to the left
            interpol[j]=(m_b0[j]*h + m_c0[j])*h + m_z[j*n];
        } else if(x>m_x[n-1]) {
            // extrapolation to the right
            interpol[j]=(m_b[(n-1)+j*n]*h + m_c[(n-1)+j*n])*h + m_z[(n-1)+j*n];
        } else {
            // interpolation
            interpol[j]=((m_a[idx+j*n]*h + m_b[idx+j*n])*h + m_c[idx+j*n])*h + m_z[idx+j*n];
        }
    }

    splinecolumn.set_points(m_y,interpol);

    return splinecolumn(y);
}


} // namespace tk


} // namespace

#endif /* TK_SPLINE_H */
