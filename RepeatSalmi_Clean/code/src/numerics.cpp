#include "salmi/numerics.hpp"
#include "salmi/constants.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace salmi {
std::vector<double> log_grid(double lo,double hi,int n) {
    if (!(lo>0 && hi>lo && n>=2)) throw std::invalid_argument("Invalid log grid");
    std::vector<double> x(n);
    for(int i=0;i<n;++i) x[i]=lo*std::exp(std::log(hi/lo)*i/(n-1));
    x.front()=lo; x.back()=hi; return x;
}
Quadrature gauss_legendre(int n,double lo,double hi) {
    if(n<1 || !(hi>lo)) throw std::invalid_argument("Invalid quadrature");
    Quadrature q; q.x.resize(n); q.w.resize(n);
    for(int i=0;i<(n+1)/2;++i) {
        double z=std::cos(cgs::pi*(i+0.75)/(n+0.5)), deriv=0;
        bool converged=false;
        for(int j=0;j<100;++j) {
            double p0=1,p1=z;
            for(int k=2;k<=n;++k) {double p=((2*k-1)*z*p1-(k-1)*p0)/k;p0=p1;p1=p;}
            if(n==1) {p0=1;p1=z;}
            deriv=n*(z*p1-p0)/(z*z-1);
            double dz=p1/deriv; z-=dz;
            if(std::abs(dz)<2e-15) {converged=true;break;}
        }
        if(!converged) throw std::runtime_error("Gauss root did not converge");
        double w=(hi-lo)/((1-z*z)*deriv*deriv), mid=(hi+lo)/2,half=(hi-lo)/2;
        q.x[i]=mid-half*z;q.x[n-1-i]=mid+half*z;
        q.w[i]=q.w[n-1-i]=w;
    } return q;
}
double integrate(const std::function<double(double)>& f,double lo,double hi,int n) {
    auto q=gauss_legendre(n,lo,hi);double s=0;
    for(int i=0;i<n;++i) s+=q.w[i]*f(q.x[i]);return s;
}
double power_integral(double lo,double hi,double exponent) {
    if(!(lo>0 && hi>lo)) throw std::invalid_argument("Invalid power integral");
    double a=exponent+1, L=std::log(hi/lo);
    if(a==0) return L;
    return std::pow(lo,a)*std::expm1(a*L)/a;
}
double linear_interp(const std::vector<double>& x,const std::vector<double>& y,double v) {
    if(x.size()!=y.size() || x.size()<2 || v<x.front() || v>x.back())
        throw std::out_of_range("Interpolation outside tabulation");
    if(v==x.back()) return y.back();
    auto p=std::upper_bound(x.begin(),x.end(),v); auto i=std::size_t(p-x.begin()-1);
    return y[i]+(y[i+1]-y[i])*(v-x[i])/(x[i+1]-x[i]);
}
double rk4_step(const std::function<double(double,double)>& f,double t,double y,double h) {
    double a=f(t,y),b=f(t+h/2,y+h*a/2),c=f(t+h/2,y+h*b/2),d=f(t+h,y+h*c);
    return y+h*(a+2*b+2*c+d)/6;
}
std::vector<double> tridiagonal(std::vector<double>a,std::vector<double>b,
                               std::vector<double>c,std::vector<double>r) {
    auto n=b.size();
    if(n==0 || a.size()!=n || c.size()!=n || r.size()!=n) throw std::invalid_argument("Tridiagonal size");
    for(std::size_t i=1;i<n;++i) {
        if(std::abs(b[i-1])<1e-300) throw std::runtime_error("Singular tridiagonal");
        double z=a[i]/b[i-1];b[i]-=z*c[i-1];r[i]-=z*r[i-1];
    }
    if(std::abs(b.back())<1e-300) throw std::runtime_error("Singular tridiagonal");
    r.back()/=b.back();
    for(int i=int(n)-2;i>=0;--i) r[i]=(r[i]-c[i]*r[i+1])/b[i];
    return r;
}
}
