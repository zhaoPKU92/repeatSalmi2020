#include "salmi/radiation.hpp"
#include "salmi/constants.hpp"
#include "salmi/numerics.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace salmi {
double occupation_intensity_x(double x,double n) {
    if(!(x>0 && n>=0) || !std::isfinite(x+n))throw std::invalid_argument("Invalid intensity conversion");
    return 2*std::pow(cgs::mec2,4)/(std::pow(cgs::h,3)*cgs::c*cgs::c)*x*x*x*n;
}
double planck_kev(double E,double T) {
    if(!(E>0 && T>0)) throw std::invalid_argument("Invalid Planck argument");
    double e=E*cgs::kev,x=e/(cgs::kb*T);
    if(x>700) return 0;
    return 2*e*e*e/(cgs::h*cgs::h*cgs::h*cgs::c*cgs::c)*cgs::kev/std::expm1(x);
}
double free_free_opacity(double E,double T,double rho,double gaunt) {
    if(!(E>0 && T>0 && rho>=0 && gaunt>0)) throw std::invalid_argument("Invalid free-free argument");
    double nu=E*cgs::kev/cgs::h;
    // Rybicki-Lightman form of Kramers law, numerical coefficient in CGS.
    return 3.692e8*gaunt*rho/(cgs::mp*cgs::mp*std::sqrt(T)*nu*nu*nu)*
        (-std::expm1(-E*cgs::kev/(cgs::kb*T)));
}
SCWeights short_characteristic_weights(double d) {
    if(!(d>=0) || !std::isfinite(d)) throw std::invalid_argument("Invalid SC optical path");
    if(d==0) return {1,0,0};
    double e0=-std::expm1(-d),local;
    if(d<1e-3) local=d*(0.5+d*(-1.0/6+d*(1.0/24+d*(-1.0/120+d/720))));
    else local=1-e0/d;
    return {std::exp(-d),e0-local,local};
}
double parabolic_source_integral(double d,double third,double sl,double su,double so) {
    if(!(d>0) || third==0 || third==d || !(sl>=0 && su>=0 && so>=0))
        throw std::invalid_argument("Invalid parabolic source coordinates");
    double q=third/d,u=su-sl,v=so-sl,c=(u-v/q)/(1-q),b=u-c;
    double minimum=std::min(sl,su);
    if(c>0) {double p=-b/(2*c);if(p>0 && p<1)minimum=std::min(minimum,sl+b*p+c*p*p);}
    if(minimum<0) {auto w=short_characteristic_weights(d);return w.local*sl+w.upstream*su;}
    double w0=-std::expm1(-d),w1,w2;
    if(d<.01) {
        w1=0;w2=0;double term=d;
        for(int k=0;k<12;++k){w1+=term/(k+2);w2+=term/(k+3);term*=-d/(k+1);}
    } else {
        double e=std::exp(-d);w1=(1-(d+1)*e)/d;w2=(2-(d*d+2*d+2)*e)/(d*d);
    }
    double value=sl*w0+b*w1+c*w2;
    if(!(value>=0) || !std::isfinite(value)) {auto w=short_characteristic_weights(d);return w.local*sl+w.upstream*su;}
    return value;
}
FormalSolution formal_solve(const std::vector<double>& tau,const std::vector<double>& s,
                           double mu,double bottom,double top) {
    std::size_t n=tau.size();
    if(n<2 || s.size()!=n || !(mu>0 && mu<=1 && bottom>=0 && top>=0))
        throw std::invalid_argument("Invalid formal solution input");
    for(std::size_t i=0;i<n;++i) {
        if(!(s[i]>=0) || !std::isfinite(s[i]) || !std::isfinite(tau[i]) || (i && !(tau[i]>tau[i-1])))
            throw std::invalid_argument("Invalid source or depth grid");
    }
    FormalSolution r;r.inward.resize(n);r.outward.resize(n);r.inward[0]=top;r.outward[n-1]=bottom;
    for(std::size_t i=1;i<n;++i) {
        auto w=short_characteristic_weights((tau[i]-tau[i-1])/mu);
        r.inward[i]=r.inward[i-1]*w.transmission+s[i-1]*w.upstream+s[i]*w.local;
    }
    for(int i=int(n)-2;i>=0;--i) {
        auto w=short_characteristic_weights((tau[i+1]-tau[i])/mu);
        r.outward[i]=r.outward[i+1]*w.transmission+s[i+1]*w.upstream+s[i]*w.local;
    }return r;
}
CoherentSolution coherent_slab(const std::vector<double>& tau,const std::vector<double>& b,
    const std::vector<double>& eps,double bottom,int nmu,double tol,int maxit) {
    auto n=tau.size();
    if(n<2 || b.size()!=n || eps.size()!=n || !(tol>0) || maxit<1) throw std::invalid_argument("Invalid coherent slab");
    auto q=gauss_legendre(nmu,0,1);
    CoherentSolution r;r.source=b;r.mean_intensity.resize(n);r.flux.resize(n);
    std::vector<double> diag(n,0);
    for(std::size_t i=0;i<n;++i) {
        if(!(eps[i]>0 && eps[i]<=1)) throw std::invalid_argument("Destruction probability outside (0,1]");
        for(int a=0;a<nmu;++a) {
            if(i) diag[i]+=0.5*q.w[a]*short_characteristic_weights((tau[i]-tau[i-1])/q.x[a]).local;
            if(i+1<n) diag[i]+=0.5*q.w[a]*short_characteristic_weights((tau[i+1]-tau[i])/q.x[a]).local;
        }
    }
    for(int it=0;it<maxit;++it) {
        std::fill(r.mean_intensity.begin(),r.mean_intensity.end(),0);
        std::fill(r.flux.begin(),r.flux.end(),0);
        for(int a=0;a<nmu;++a) {
            auto f=formal_solve(tau,r.source,q.x[a],bottom);
            for(std::size_t i=0;i<n;++i) {
                r.mean_intensity[i]+=0.5*q.w[a]*(f.outward[i]+f.inward[i]);
                r.flux[i]+=2*cgs::pi*q.w[a]*q.x[a]*(f.outward[i]-f.inward[i]);
            }
        }
        double err=0;
        std::vector<double> next(n);
        for(std::size_t i=0;i<n;++i) {
            double target=eps[i]*b[i]+(1-eps[i])*r.mean_intensity[i];
            err=std::max(err,std::abs(target-r.source[i])/std::max(std::abs(target),1e-100));
            next[i]=r.source[i]+(target-r.source[i])/(1-(1-eps[i])*diag[i]);
            if(!(next[i]>=0) || !std::isfinite(next[i])) throw std::runtime_error("ALI invalid source");
        }
        r.change_history.push_back(err);r.residual=err;r.iterations=it+1;
        // Leave intensities and source from the SAME formal solution on successful exit.
        if(err<tol) return r;
        r.source.swap(next);
    }throw std::runtime_error("Coherent ALI failed to converge");
}
}
