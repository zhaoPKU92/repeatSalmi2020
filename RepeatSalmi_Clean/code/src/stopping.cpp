#include "salmi/stopping.hpp"
#include "salmi/constants.hpp"
#include <algorithm>
#include <cmath>
#include <exception>
#include <stdexcept>
namespace salmi {
namespace {
void check_gamma(double g) {
    if(!(g>=1.1) || !std::isfinite(g)) throw std::invalid_argument("Fast-electron loss requires gamma >= 1.1");
}
}
double phi_coulomb_langmuir(double g,double ne) {
    check_gamma(g);
    if(!(ne>0) || !std::isfinite(ne)) throw std::invalid_argument("Electron density must be positive");
    double b2=(g-1)*(g+1)/(g*g);
    double ep2=cgs::hbar*cgs::hbar*(4*cgs::pi*ne*cgs::e*cgs::e/cgs::me)/(cgs::mec2*cgs::mec2);
    double value=0.75/b2*(std::log((g-1)*b2/(2*ep2))+1+
        0.125*std::pow((g-1)/g,2)-(2*g-1)/(g*g)*std::log(2.0));
    if(!(value>0) || !std::isfinite(value)) throw std::domain_error("Loss approximation outside physical domain");
    return value;
}
double phi_electron_proton(double g) {
    check_gamma(g);double p2=(g-1)*(g+1),p=std::sqrt(p2),b=p/g;
    return 3*cgs::alpha/cgs::pi*g*g*g/(p2+g*g)*
        (std::asinh(p)/b-1.0/3+p2/std::pow(g,6)*
        (2.0/9*g*g-19.0/675*g*p2-0.06*p2*p2/g));
}
double phi_electron_electron(double g) {
    check_gamma(g);
    // Haug's Q_cm times (gamma-1)/sigma_T. The kinetic-energy factor is essential.
    // Primary-source confirmation is tracked in docs/SOURCES.md.
    double shape,energy=(g-1)*cgs::mec2/cgs::kev;
    if(energy<=100) {
        double p=std::sqrt((g-1)*(g+1));
        shape=(g+1)/g*(0.6664+43.935*p-2.272*p*p-3.055*p*p*p)*
              (-std::expm1(-2*cgs::pi*cgs::alpha*g/p));
    }
    else if(energy<=1000) shape=(g+1)/g*(0.5754+1.4492*g-0.0665*g*g);
    else if(energy<=500000) shape=(g+1)*(4.181*g*std::acosh(g)-2.676*g-2.256)/(g*g+1.022*g-3.871);
    else shape=4*(std::log(2*g)-1.0/3);
    return 3*cgs::alpha/(8*cgs::pi)*(g-1)*shape;
}
LossComponents loss_components(double g,double ne) {
    return {phi_coulomb_langmuir(g,ne),phi_electron_proton(g),phi_electron_electron(g)};
}
double loss_rate(double g,double ne,LossMode mode) {
    double v=phi_coulomb_langmuir(g,ne);
    if(mode!=LossMode::CoulombLangmuir) v+=phi_electron_proton(g);
    if(mode==LossMode::FullBremsstrahlung) v+=phi_electron_electron(g);
    return v;
}
Track stopping_track(double g0,const AtmosphereFunction& gas,LossMode mode,double ds,double top,double cut,double relative_tau_step) {
    if(!(g0>cut && cut>=1.1 && ds>0 && ds<=0.2 && top>0 && relative_tau_step>=1e-4 && relative_tau_step<=0.08)
       || !std::isfinite(g0+ds+top+cut+relative_tau_step)) throw std::invalid_argument("Invalid track parameters");
    Track tr;tr.gamma0=g0;tr.cutoff=cut;
    double send=std::log((g0-1)/(cut-1)),s=0,tau=top;
    auto f=[&](double u,double t) {
        double g=std::max(cut,1+(g0-1)*std::exp(-u));
        return (g-1)/loss_rate(g,gas(t).ne,mode);
    };
    for(int iter=0;iter<1000000;++iter) {
        double g=std::max(cut,1+(g0-1)*std::exp(-s));
        if(s==send) g=cut;
        tr.tau.push_back(tau);tr.gamma.push_back(g);tr.phi.push_back(loss_rate(g,gas(tau).ne,mode));
        if(s>=send) return tr;
        // Limits the first steps near the logarithmic ne singularity at tau=0.
        double step=std::min({ds,send-s,relative_tau_step*tau/f(s,tau)});
        double next=rk4_step(f,s,tau,step);
        if(!(next>tau) || !std::isfinite(next)) throw std::runtime_error("Stopping integrator failed");
        tau=next;s=std::min(send,s+step);
        if(send-s<1e-14) s=send;
    }
    throw std::runtime_error("Stopping integrator iteration limit");
}
double Track::kinetic(double t,bool thermalize) const {
    if(t<=tau.front()) return gamma0-1;
    if(t>=tau.back()) return thermalize?0:cutoff-1;
    return linear_interp(tau,gamma,t)-1;
}
double Track::momentum(double t,bool thermalize) const {
    double k=kinetic(t,thermalize);return std::sqrt(k*(k+2));
}
Deposition deposit(const Quadrature& beam,const AtmosphereFunction& gas,LossMode mode,
    const FluxBudget& budget,double mean,const std::vector<double>& edges,double ds,double cut,double relative_tau_step) {
    if(edges.size()<2 || beam.x.empty() || beam.x.size()!=beam.w.size() || !(edges.front()>0))
        throw std::invalid_argument("Invalid deposition input");
    auto n=edges.size()-1;
    for(std::size_t i=0;i<n;++i) if(!(edges[i+1]>edges[i])) throw std::invalid_argument("Depth edges not increasing");
    Deposition d;d.tau_edges=edges;d.average_phi.assign(n,0);d.qplus.assign(n,0);
    d.gram.assign(n,0);d.flux.assign(n+1,budget.total());
    d.quadrature_number=0;d.quadrature_energy=0;double rem=0;
    for(double w:beam.w)if(!(w>=0) || !std::isfinite(w))throw std::invalid_argument("Invalid beam weight");
    // Each trajectory is independent. Store only cell losses, then sum in the
    // original beam order. The summation order is independent of thread count;
    // individual elementary-function evaluations can still differ by roundoff.
    std::vector<double> energy_loss(beam.x.size()*n),momentum_loss(beam.x.size()*n),remaining(beam.x.size());
    std::vector<std::exception_ptr> errors(beam.x.size());
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for(int j=0;j<int(beam.x.size());++j) {
      try {
        double w=beam.w[j];
        auto tr=stopping_track(beam.x[j],gas,mode,ds,edges.front(),cut,relative_tau_step);
        double prev=beam.x[j]-1,pm=std::sqrt(prev*(prev+2));
        for(std::size_t i=0;i<n;++i) {
            double k=tr.kinetic(edges[i+1]),p=std::sqrt(k*(k+2)),dt=edges[i+1]-edges[i];
            energy_loss[j*n+i]=w*(prev-k)/dt;momentum_loss[j*n+i]=w*(pm-p)/dt;prev=k;pm=p;
        }remaining[j]=w*prev;
      }catch(...){errors[j]=std::current_exception();}
    }
    for(std::size_t j=0;j<beam.x.size();++j) {
        if(errors[j])std::rethrow_exception(errors[j]);
        double w=beam.w[j];d.quadrature_number+=w;d.quadrature_energy+=w*(beam.x[j]-1);
        for(std::size_t i=0;i<n;++i){d.average_phi[i]+=energy_loss[j*n+i];d.gram[i]+=momentum_loss[j*n+i];}
        rem+=remaining[j];
    }
    double norm=budget.number_flux(mean)*cgs::mec2;
    double energy=0;
    for(std::size_t i=0;i<n;++i) {
        d.qplus[i]=norm*cgs::kappa_e*d.average_phi[i];
        d.gram[i]*=norm/cgs::c*cgs::kappa_e;
        energy+=d.average_phi[i]*(edges[i+1]-edges[i]);
        d.flux[i+1]=budget.total()-norm*energy;
    }
    d.deposited_fraction=energy/(mean-1);d.remaining_fraction=rem/(mean-1);
    return d;
}
}
