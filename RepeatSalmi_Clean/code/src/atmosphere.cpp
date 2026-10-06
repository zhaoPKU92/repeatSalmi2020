#include "salmi/atmosphere.hpp"
#include "salmi/constants.hpp"
#include <cmath>
#include <stdexcept>
namespace salmi {
GasState GreyAtmosphere::operator()(double tau) const {
    if(!(tau>0 && teff>0 && gravity>0)) throw std::invalid_argument("Invalid grey atmosphere");
    GasState s;
    s.temperature=teff*std::pow(0.75*(tau+2.0/3),0.25);
    s.pressure=gravity*tau/cgs::kappa_e;
    s.ne=s.pressure/(2*cgs::kb*s.temperature); // electrons AND protons
    s.rho=s.ne*cgs::mp;return s;
}
FluxBudget FluxBudget::from_temperatures(double th,double ti) {
    if(!(th>=0 && ti>=0) || !std::isfinite(th+ti)) throw std::invalid_argument("Invalid temperature");
    return {cgs::sigma_sb*std::pow(th,4),cgs::sigma_sb*std::pow(ti,4)};
}
double FluxBudget::number_flux(double g) const {
    if(!(g>1)) throw std::invalid_argument("Mean gamma <= 1");
    return heated/(cgs::mec2*(g-1));
}
double eddington_flux(double g) {
    if(!(g>0)) throw std::invalid_argument("Gravity <= 0");return g*cgs::c/cgs::kappa_e;
}
double goldreich_julian(double b,double p) {
    if(!(b>0 && p>0)) throw std::invalid_argument("Invalid B/P");return b/(cgs::e*p*cgs::c);
}
double incoming_density(double th,double g) {
    return FluxBudget::from_temperatures(th,0).number_flux(g)/cgs::c;
}
double redshift_from_gravity(double g,double mass) {
    if(!(g>0 && mass>0)) throw std::invalid_argument("Invalid gravity/mass");
    double gm=cgs::grav*mass*cgs::solar_mass;
    // u=2GM/Rc^2. g=c^4 u^2/(4GM sqrt(1-u)), monotone on (0,1).
    double lo=0,hi=1;
    for(int i=0;i<100;++i) {
        double u=(lo+hi)/2, v=std::pow(cgs::c,4)*u*u/(4*gm*std::sqrt(1-u));
        if(v>g) hi=u;else lo=u;
    }return 1/std::sqrt(1-(hi+lo)/2)-1;
}
std::vector<double> hydrostatic(const std::vector<double>& m,const std::vector<double>& ram,
    const std::vector<double>& rad,double gravity,double ptop) {
    if(m.size()<2 || ram.size()+1!=m.size() || rad.size()!=ram.size() || !(gravity>0 && ptop>=0))
        throw std::invalid_argument("Hydrostatic input");
    std::vector<double> p(m.size());p[0]=ptop;
    for(std::size_t i=0;i<ram.size();++i) {
        if(!(m[i+1]>m[i])) throw std::invalid_argument("Mass grid not increasing");
        p[i+1]=p[i]+(gravity+ram[i]-rad[i])*(m[i+1]-m[i]);
        if(!(p[i+1]>0)) throw std::runtime_error("Nonpositive hydrostatic pressure");
    }return p;
}
}
