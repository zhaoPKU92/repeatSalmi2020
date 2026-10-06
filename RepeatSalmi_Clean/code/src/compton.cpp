#include "salmi/compton.hpp"
#include "salmi/constants.hpp"
#include "salmi/numerics.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace salmi {
namespace {
// Exact azimuth average of F=2+(2-q)d+d^2, d=1/(z+q)-1/z,
// z=A-B cos(phi). Rationalized differences avoid the cancellations in A.22.
template<class Real> double azimuth_average(double Ad,double Bd,double qd) {
    Real a=Ad,b=Bd,q=qd;
    Real m=std::sqrt((a-b)*(a+b)),p=std::sqrt((a+q-b)*(a+q+b));
    Real r=q*(2*a+q)/(p+m),d=-r/(p*m);
    Real d2=r*(a*r*(p*p+3*p*m+m*m)-q*m*m*(2*p+m))/(m*m*m*p*p*p*(m+p));
    return double(2+(2-q)*d+d2);
}
double azimuth_reaction(double a,double b,double q,bool extended=false) {
    double answer=(extended || a-std::abs(b)<1e-4*a)?azimuth_average<long double>(a,b,q):azimuth_average<double>(a,b,q);
    if(!(answer>0) || !std::isfinite(answer))throw std::runtime_error("Invalid analytic azimuth integral");
    return answer;
}
}
double mono_compton_kernel_analytic(double x,double y,double eta,double g,bool extended) {
    if(!(g>=1) || !std::isfinite(g))throw std::invalid_argument("Invalid electron gamma");
    double gs=minimum_electron_gamma(x,y,eta);if(g<gs*(1-1e-13))return 0;g=std::max(g,gs);
    double q=x*y*(1-eta),Q=std::sqrt((x-y)*(x-y)+2*q),p=std::sqrt((g-1)*(g+1));
    double ca=p>0?(g*(y-x)-q)/(p*Q):0;ca=std::max(-1.0,std::min(1.0,ca));
    double ct=(y*eta-x)/Q,st=y*std::sqrt(std::max(0.0,1-eta*eta))/Q;
    double aa=x*(g-p*ca*ct),bb=x*p*std::sqrt(std::max(0.0,1-ca*ca))*st;
    return azimuth_reaction(aa,bb,q,extended)/Q;
}
double klein_nishina_total(double x) {
    if(!(x>0)) throw std::invalid_argument("Klein-Nishina energy <= 0");
    if(x<1e-3) return 1+x*(-2+x*(5.2+x*(-13.3+x*1144.0/35)));
    double l=std::log1p(2*x);
    return 0.75*((1+x)/(x*x*x)*(2*x*(1+x)/(1+2*x)-l)+l/(2*x)-(1+3*x)/std::pow(1+2*x,2));
}
double scaled_bessel_k2(double th,int order) {
    if(!(th>0 && th<=2)) throw std::invalid_argument("Bessel theta outside tested (0,2]");
    // K2(1/th)*exp(1/th), substitute sinh(t/2)=sqrt(th)*v to avoid underflow.
    return integrate([&](double v) {
        double z=th*v*v;
        return 2*std::sqrt(th)*std::exp(-2*v*v)*(1+8*z+8*z*z)/std::sqrt(1+z);
    },0,7,order);
}
double minimum_electron_gamma(double x,double y,double eta) {
    if(!(x>0 && y>0 && eta>=-1 && eta<1) || !std::isfinite(x+y)) throw std::invalid_argument("Invalid kernel kinematics");
    double q=x*y*(1-eta),Q=std::sqrt((x-y)*(x-y)+2*q);
    return std::max(1.0,0.5*(x-y+Q*std::sqrt(1+2/q)));
}
double mono_compton_kernel(double x,double y,double eta,double g,int nphi) {
    if(!(g>=1) || !std::isfinite(g) || nphi<4) throw std::invalid_argument("Invalid electron gamma or azimuth order");
    double gs=minimum_electron_gamma(x,y,eta);
    if(g<gs*(1-1e-13)) return 0;
    g=std::max(g,gs);
    double q=x*y*(1-eta),Q=std::sqrt((x-y)*(x-y)+2*q),p=std::sqrt((g-1)*(g+1));
    double ca=p>0?(g*(y-x)-q)/(p*Q):0;
    ca=std::max(-1.0,std::min(1.0,ca));
    double ct=(y*eta-x)/Q, st=y*std::sqrt(std::max(0.0,1-eta*eta))/Q;
    double sa=std::sqrt(std::max(0.0,1-ca*ca));
    // phi integral by midpoint trapezoid; analytic periodic smooth integrand.
    double sum=0;
    for(int i=0;i<nphi;++i) {
        double cp=std::cos(2*cgs::pi*(i+0.5)/nphi);
        double xi1=x*(g-p*(ca*ct+sa*st*cp)); // incoming electron dot outgoing photon
        double xi=xi1+q;                    // incoming electron dot incoming photon
        double d=-q/(xi*xi1);
        double F=d*d+2*d+xi/xi1+xi1/xi;
        if(!(F>0)) throw std::runtime_error("Nonpositive KN reaction rate");
        sum+=F;
    }return sum/(nphi*Q);
}
double thermal_compton_kernel(double x,double y,double eta,double th,int ngamma,int nphi) {
    if(!(th>0 && th<=2)) throw std::invalid_argument("Kernel theta outside tested domain");
    if(ngamma<8 || nphi<4) throw std::invalid_argument("Insufficient kernel quadrature order");
    double gs=minimum_electron_gamma(x,y,eta),a=(gs-1)/th;
    if(a>700) return 0;
    double z=integrate([&](double u) {
        return std::exp(-u)*mono_compton_kernel(x,y,eta,gs+th*u,nphi);
    },0,40,ngamma);
    return 3/(32*cgs::pi*scaled_bessel_k2(th))*std::exp(-a)*z;
}
double thermal_scattering_rate(double x,double th,int order) {
    auto qa=gauss_legendre(order,-1,1);
    // u=sqrt((gamma-1)/theta) removes square-root endpoint singularity.
    double val=integrate([&](double u) {
        double g=1+th*u*u,p=std::sqrt((g-1)*(g+1)),beta=p/g,s=0;
        for(int i=0;i<order;++i) {
            double v=1-beta*qa.x[i];s+=0.5*qa.w[i]*v*klein_nishina_total(g*x*v);
        }
        return 2*u*g*p*std::exp(-u*u)*s;
    },0,7,order);
    return val/scaled_bessel_k2(th);
}
ThermalKernel::ThermalKernel(double th,int ngamma,int nphi):theta_(th) {
    if(!(th>0 && th<=2) || ngamma<8 || nphi<4)
        throw std::invalid_argument("Invalid cached thermal kernel parameters");
    auto q=gauss_legendre(ngamma,0,40);gamma_u_=q.x;gamma_w_=q.w;
    for(int i=0;i<ngamma;++i)gamma_w_[i]*=std::exp(-gamma_u_[i]);
    for(int i=0;i<nphi;++i)phi_cos_.push_back(std::cos(2*cgs::pi*(i+.5)/nphi));
    prefactor_=3/(32*cgs::pi*scaled_bessel_k2(th));
}
double ThermalKernel::operator()(double x,double y,double eta) const {
    double gs=minimum_electron_gamma(x,y,eta),a=(gs-1)/theta_;
    if(a>700)return 0;
    double q=x*y*(1-eta),Q=std::sqrt((x-y)*(x-y)+2*q);
    double ct=(y*eta-x)/Q,st=y*std::sqrt(std::max(0.0,1-eta*eta))/Q,total=0;
    for(std::size_t i=0;i<gamma_u_.size();++i) {
        double g=gs+theta_*gamma_u_[i],p=std::sqrt((g-1)*(g+1));
        double ca=std::max(-1.0,std::min(1.0,(g*(y-x)-q)/(p*Q)));
        double sa=std::sqrt(std::max(0.0,1-ca*ca)),av=0;
        double aa=x*(g-p*ca*ct),bb=x*p*sa*st;
        av=azimuth_reaction(aa,bb,q);
        total+=gamma_w_[i]*av/Q;
    }
    return prefactor_*std::exp(-a)*total;
}
}
