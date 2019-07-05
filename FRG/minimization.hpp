/*
  Minimisation algorithm taken from the numerical recipes
*/

#include <float.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <sstream>
#include <stdio.h>      /* printf */
#include <math.h>       /* Pow */
#include <stdlib.h>      /* Calloc */
#include <cmath>
#include <iomanip>
#include <complex>
#include <ctime>
#include <string>
#include <limits>       // std::numeric_limits
using namespace std;
typedef std::vector< double > state_type;



struct Func{
    //long double operator()(const long double x){}

    long double func(const long double y,long double y0,long double c,const state_type x){
        return x[0]+
        x[1]*(y-y0)+
        1.0/2.0*x[2]*pow(y-y0,2)+
        1.0/6.0*x[3]*pow(y-y0,3)+
        1.0/24.0*x[4]*pow(y-y0,4)+
        1.0/120.0*x[5]*pow(y-y0,5)+
        1.0/720.0*x[6]*pow(y-y0,6)+
        1.0/5040.0*x[7]*pow(y-y0,7)+
        1.0/40320*x[8]*pow(y-y0,8)-c*sqrt(2*y);
    }
    long double df(const long double y, long double y0,long double c,const state_type x){
        return x[1]+
        x[2]*(y-y0)+
        1.0/2.0*x[3]*pow(y-y0,2)+
        1.0/6.0*x[4]*pow(y-y0,3)+
        1.0/24.0*x[5]*pow(y-y0,4)+
        1.0/120.0*x[6]*pow(y-y0,5)+
        1.0/720.0*x[7]*pow(y-y0,6)+
        1.0/5040.0*x[8]*pow(y-y0,7)-c*1.0/sqrt(2*y);
    }

    long double ddf(const long double y, long double y0,long double c,const state_type x){
        return x[2]+
        x[3]*(y-y0)+
        1.0/2.0*x[4]*pow(y-y0,2)+
        1.0/6.0*x[5]*pow(y-y0,3)+
        1.0/24.0*x[6]*pow(y-y0,4)+
        1.0/120.0*x[7]*pow(y-y0,5)+
        1.0/720.0*x[8]*pow(y-y0,6)+c*1.0/sqrt(2*y)*1.0/(2*y);
    }
};



struct Func2{
    //long double operator()(const long double x){}

    long double func(const long double y,long double y0,long double c,const state_type x){
        return x[0]+
        x[1]*(y-y0)+
        1.0/2.0*x[2]*pow(y-y0,2)+
        1.0/6.0*x[3]*pow(y-y0,3)+
        1.0/24.0*x[4]*pow(y-y0,4)+
        1.0/120.0*x[5]*pow(y-y0,5)+
        1.0/720.0*x[6]*pow(y-y0,6)+
        1.0/5040.0*x[7]*pow(y-y0,7)+
        1.0/40320*x[8]*pow(y-y0,8);
    }
    long double df(const long double y, long double y0,long double c,const state_type x){
        return x[1]+
        x[2]*(y-y0)+
        1.0/2.0*x[3]*pow(y-y0,2)+
        1.0/6.0*x[4]*pow(y-y0,3)+
        1.0/24.0*x[5]*pow(y-y0,4)+
        1.0/120.0*x[6]*pow(y-y0,5)+
        1.0/720.0*x[7]*pow(y-y0,6)+
        1.0/5040.0*x[8]*pow(y-y0,7);
    }

    long double ddf(const long double y, long double y0,long double c,const state_type x){
        return x[2]+
        x[3]*(y-y0)+
        1.0/2.0*x[4]*pow(y-y0,2)+
        1.0/6.0*x[5]*pow(y-y0,3)+
        1.0/24.0*x[6]*pow(y-y0,4)+
        1.0/120.0*x[7]*pow(y-y0,5)+
        1.0/720.0*x[8]*pow(y-y0,6);//+c*1.0/sqrt(2*y)*1.0/(2*y);
    }
};



struct Dbrent {
    long double xmin,fmin;
    long double ax,bx,cx;
    const long double tol;
    Dbrent(const long double toll=3.0e-10) : tol(toll){}
    template <class T>
    long double minimize(const state_type xvec,const long double ax,const long double cx,const long double bx,long double c0, T &funcd)
    {
        const int ITMAX=100;
        const long double ZEPS = numeric_limits<long double>::epsilon()*1.0e-3;
        bool ok1,ok2;
        long double a,b,d=0.0,d1,d2,du,dv,dw,dx,e=0.0;
        long double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
        Func fn,dfn;

// bx is between ax and cx

        if(ax<cx){a=ax; b=cx;};
        if(cx<ax){a=cx; b=ax;};
        x=w=v=bx;
        fw=fv=fx=fn.func(x,bx,c0,xvec);
        dw=dv=dx=dfn.df(x,bx,c0,xvec);

        for(int iter=0;iter<ITMAX;iter++){
            xm=0.5*(a+b);

            tol1=tol*abs(x)+ZEPS;
            tol2=2.0*tol1;

            if(abs(x-xm)<=(tol2-0.5*(b-a))){

                fmin=fx;
                return xmin=x;
            }
            if(abs(e) > tol1){

            d1=2.0*(b-a);
            d2=d1;
            if(dw != dx){d1=(w-x)*dx/(dx-dw);}
            if(dv != dx){d2=(v-x)*dx/(dx-dv);}
            u1=x+d1;
            u2=x+d2;
            ok1=(a-u1)*(u1-b) > 0.0 && dx*d1 <=0.0;
            ok2=(a-u2)*(u2-b) > 0.0 && dx*d2 <=0.0;
            olde=e;
            e=d;
            if(ok1||ok2){
                if(ok1 && ok2){
                    if(abs(d1)<abs(d2)){d=d1;}
                    if(abs(d1)>abs(d2)){d=d2;}
                }
                else if(ok1){d=d1;}
                else{d=d2;}
                if(abs(d)<=abs(0.5*olde)){
                    u=x+d;
                    if(u-a<tol2 || b-u<tol2){
                        d=SIGN(abs(tol1),(xm-x));
                    }
                }
                else{
                    if(dx>=0){e=a-x;}
                    if(dx<0){e=b-x;}
                    d=0.5*e;
                }
            }
            else{
                if(dx>=0){e=a-x;}
                if(dx<0){e=b-x;}
                d=0.5*e;
            }
        }
        else{
            if(dx>=0){e=a-x;}
            if(dx<0){e=b-x;}
            d=0.5*e;
            }

    if(abs(d)>=tol1){
      u=x+d;
      fu=fn.func(u,bx,c0,xvec);
    }

    else{
      u=x+SIGN(abs(tol1),d);
      fu=fn.func(u,bx,c0,xvec);
      if(fu>fx){
        fmin=fx;
        return xmin=x;
      }
    }

    du=dfn.df(u,bx,c0,xvec);
    if(fu<=fx){
      if(u>=x){a=x;}
      else{b=x;}
      mov3(v,fv,dv,w,fw,dw);
      mov3(w,fw,dw,x,fx,dx);
      mov3(x,fx,dx,u,fu,du);
    }
    else{
      if(u<x){a=u;}
      else{b=u;}
      if(fu<=fw||w==x){
        mov3(v,fv,dv,w,fw,dw);
        mov3(w,fw,dw,u,fu,du);
      }
      else if(fu < fv || v==x || v==w){
        mov3(v,fv,dv,u,fu,du);
      }

  }
}
    //cout << "Too many iterations in routine dbrent" << endl;
return xmin;
}
    inline void mov3(long double &a, long double &b, long double &c, const long double d, const long double e,
    const long double f)
    {
    a=d; b=e; c=f;
    }
    long double SIGN(long double a,long double b){
    long double res;

    if(b>=0){res=a;}
    if(b<0){res=-a;}
    return res;
    }
};











struct Dbrent2 {
    long double xmin,fmin;
    long double ax,bx,cx;
    const long double tol;
    Dbrent2(const long double toll=3.0e-10) : tol(toll){}
    template <class T>
    long double minimize(const state_type xvec,const long double ax,const long double cx,const long double bx,long double c0, T &funcd)
    {
        const int ITMAX=100;
        const long double ZEPS = numeric_limits<long double>::epsilon()*1.0e-3;
        bool ok1,ok2;
        long double a,b,d=0.0,d1,d2,du,dv,dw,dx,e=0.0;
        long double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
        Func2 fn,dfn;

// bx is between ax and cx

        if(ax<cx){a=ax; b=cx;};
        if(cx<ax){a=cx; b=ax;};
        x=w=v=bx;
        fw=fv=fx=fn.func(x,bx,c0,xvec);
        dw=dv=dx=dfn.df(x,bx,c0,xvec);

        for(int iter=0;iter<ITMAX;iter++){
            xm=0.5*(a+b);

            tol1=tol*abs(x)+ZEPS;
            tol2=2.0*tol1;

            if(abs(x-xm)<=(tol2-0.5*(b-a))){

                fmin=fx;
                return xmin=x;
            }
            if(abs(e) > tol1){

            d1=2.0*(b-a);
            d2=d1;
            if(dw != dx){d1=(w-x)*dx/(dx-dw);}
            if(dv != dx){d2=(v-x)*dx/(dx-dv);}
            u1=x+d1;
            u2=x+d2;
            ok1=(a-u1)*(u1-b) > 0.0 && dx*d1 <=0.0;
            ok2=(a-u2)*(u2-b) > 0.0 && dx*d2 <=0.0;
            olde=e;
            e=d;
            if(ok1||ok2){
                if(ok1 && ok2){
                    if(abs(d1)<abs(d2)){d=d1;}
                    if(abs(d1)>abs(d2)){d=d2;}
                }
                else if(ok1){d=d1;}
                else{d=d2;}
                if(abs(d)<=abs(0.5*olde)){
                    u=x+d;
                    if(u-a<tol2 || b-u<tol2){
                        d=SIGN(abs(tol1),(xm-x));
                    }
                }
                else{
                    if(dx>=0){e=a-x;}
                    if(dx<0){e=b-x;}
                    d=0.5*e;
                }
            }
            else{
                if(dx>=0){e=a-x;}
                if(dx<0){e=b-x;}
                d=0.5*e;
            }
        }
        else{
            if(dx>=0){e=a-x;}
            if(dx<0){e=b-x;}
            d=0.5*e;
            }

    if(abs(d)>=tol1){
      u=x+d;
      fu=fn.func(u,bx,c0,xvec);
    }

    else{
      u=x+SIGN(abs(tol1),d);
      fu=fn.func(u,bx,c0,xvec);
      if(fu>fx){
        fmin=fx;
        return xmin=x;
      }
    }

    du=dfn.df(u,bx,c0,xvec);
    if(fu<=fx){
      if(u>=x){a=x;}
      else{b=x;}
      mov3(v,fv,dv,w,fw,dw);
      mov3(w,fw,dw,x,fx,dx);
      mov3(x,fx,dx,u,fu,du);
    }
    else{
      if(u<x){a=u;}
      else{b=u;}
      if(fu<=fw||w==x){
        mov3(v,fv,dv,w,fw,dw);
        mov3(w,fw,dw,u,fu,du);
      }
      else if(fu < fv || v==x || v==w){
        mov3(v,fv,dv,u,fu,du);
      }

  }
}
    //cout << "Too many iterations in routine dbrent" << endl;
    }
    inline void mov3(long double &a, long double &b, long double &c, const long double d, const long double e,
    const long double f)
    {
    a=d; b=e; c=f;
    }
    long double SIGN(long double a,long double b){
    long double res;

    if(b>=0){res=a;}
    if(b<0){res=-a;}
    return res;
    }
};
