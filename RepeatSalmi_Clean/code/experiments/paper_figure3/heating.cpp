#include "paper.hpp"
#include "salmi/constants.hpp"
#include "salmi/stopping.hpp"
#include "salmi/beam.hpp"
#include "heating_projection.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
namespace salmi {
bool paper_pointwise_heating=false;
bool paper_normalize_point_heating=false;
bool paper_hat_heating=false;
namespace {
void powerlaw_heating(HeatedState& s,double th,double ti,double delta,int panels,double relative_tau_step) {
 // Salmi (16)-(20): number-distribution quadrature, normalized by mean
 // incident kinetic energy, NOT by a representative monoenergetic particle.
 const auto budget=FluxBudget::from_temperatures(th,ti);
 const PowerLaw power(s.beam_gamma_min,s.beam_gamma_max,delta);
 const auto beam=power.quadrature(panels);
 const int D=int(s.mass.size());
 std::vector<double> lm,lt,lp,edge(D+1);
 for(int d=0;d<D;++d){lm.push_back(std::log(s.mass[d]));lt.push_back(std::log(s.temperature[d]));lp.push_back(std::log(s.pressure[d]));}
 edge.front()=s.mass.front();edge.back()=s.mass.back();
 for(int d=1;d<D;++d)edge[d]=.5*(s.mass[d-1]+s.mass[d]);
 AtmosphereFunction gas=[&](double tau){
   double m=tau/cgs::kappa_e,z=std::max(lm.front(),std::log(m));
   // A high-energy track can leave a deliberately truncated atmosphere.
   // Continue only the unmodeled trajectory with the same deep diffusion/
   // hydrostatic extension used by src/heated.cpp; deposited energy inside
   // the grid is still evaluated from its actual entry and exit energies.
   double T,P;
   if(m>s.mass.back()){T=s.temperature.back()*std::pow(m/s.mass.back(),.25);P=s.pressure.back()*m/s.mass.back();}
   else {T=std::exp(linear_interp(lm,lt,z));P=std::exp(linear_interp(lm,lp,z));}
   double ne=P/(2*cgs::kb*T);return GasState{T,P,ne,ne*cgs::mp};
 };
 const double norm=budget.heated/(power.mean_gamma()-1);
 double deposited_energy=0,deposited_momentum=0;
 for(std::size_t j=0;j<beam.x.size();++j){
   const auto tr=stopping_track(beam.x[j],gas,LossMode::CoulombLangmuir,.01,cgs::kappa_e*s.mass.front(),1.1,relative_tau_step);
   const double weight=norm*beam.w[j],remaining=tr.kinetic(cgs::kappa_e*s.mass.back());
   deposited_energy+=weight*(beam.x[j]-1-remaining);
   deposited_momentum+=weight/cgs::c*(std::sqrt(beam.x[j]*beam.x[j]-1)-std::sqrt(remaining*(remaining+2)));
   if(!paper_pointwise_heating){
     for(int d=0;d<D;++d){
       double left=tr.kinetic(cgs::kappa_e*edge[d]),right=tr.kinetic(cgs::kappa_e*edge[d+1]),width=edge[d+1]-edge[d];
       s.heating[d]+=weight*(left-right)/width;
       s.ram[d]+=weight/cgs::c*(std::sqrt(left*(left+2))-std::sqrt(right*(right+2)))/width;
     }
   }else{
     for(int d=0;d<D;++d){
       double tau=cgs::kappa_e*s.mass[d];
       if(tau<tr.tau.back()){
         double gamma=1+tr.kinetic(tau),q=weight*cgs::kappa_e*phi_coulomb_langmuir(gamma,s.rho[d]/cgs::mp);
         s.heating[d]+=q;s.ram[d]+=q/cgs::c*gamma/std::sqrt(gamma*gamma-1);
       }
     }
     double mc=tr.tau.back()/cgs::kappa_e;
     auto it=std::upper_bound(s.mass.begin(),s.mass.end(),mc);
     if(it!=s.mass.begin() && it!=s.mass.end()){
       int hi=int(it-s.mass.begin()),lo=hi-1;double f=(mc-s.mass[lo])/(s.mass[hi]-s.mass[lo]);
       for(int d:{lo,hi}){
         double fraction=d==lo?1-f:f,width=edge[d+1]-edge[d];
         s.heating[d]+=weight*.1*fraction/width;
         s.ram[d]+=weight/cgs::c*std::sqrt(1.1*1.1-1)*fraction/width;
       }
     }
   }
 }
 if(paper_pointwise_heating && paper_normalize_point_heating){
   double energy=0,momentum=0;
   for(int d=1;d<D;++d){double dm=s.mass[d]-s.mass[d-1];energy+=.5*(s.heating[d]+s.heating[d-1])*dm;momentum+=.5*(s.ram[d]+s.ram[d-1])*dm;}
   if(!(energy>0 && momentum>0))throw std::runtime_error("Power-law beam unresolved on heating grid");
   // Preserve transmitted beam energy on truncated grids. Never force all
   // incident energy into the modeled layers when particles leave the bottom.
   for(int d=0;d<D;++d){s.heating[d]*=deposited_energy/energy;s.ram[d]*=deposited_momentum/momentum;}
 }
 s.target_flux[0]=budget.total();
 for(int d=1;d<D;++d)s.target_flux[d]=s.target_flux[d-1]-.5*(s.heating[d]+s.heating[d-1])*(s.mass[d]-s.mass[d-1]);
}
}
void update_beam_heating(HeatedState& s,double th,double ti,double delta,int panels,double relative_tau_step) {
 auto budget=FluxBudget::from_temperatures(th,ti);int D=s.mass.size();
 s.heating.assign(D,0);s.ram.assign(D,0);s.target_flux.assign(D,budget.total());
 if(th==0)return;
 if(s.beam_mono_gamma==0 && !s.beam_bremsstrahlung){
  if(paper_hat_heating)throw std::invalid_argument("Hat projection is currently a monoenergetic diagnostic");
  powerlaw_heating(s,th,ti,delta,panels,relative_tau_step);return;}
 if(!(s.beam_mono_gamma>1.1) || s.beam_bremsstrahlung)throw std::invalid_argument("Figure 3 requires a monoenergetic beam without bremsstrahlung stopping");
 std::vector<double> lm,lt,lp;
 for(int d=0;d<D;++d){lm.push_back(std::log(s.mass[d]));lt.push_back(std::log(s.temperature[d]));lp.push_back(std::log(s.pressure[d]));}
 AtmosphereFunction gas=[&](double tau){
   double m=tau/cgs::kappa_e,z=std::max(lm.front(),std::log(m)); // Protect the top-node round trip from one-ulp underflow.
   double T=std::exp(linear_interp(lm,lt,z)),P=std::exp(linear_interp(lm,lp,z));
   double ne=P/(2*cgs::kb*T);return GasState{T,P,ne,ne*cgs::mp};
 };
 auto tr=stopping_track(s.beam_mono_gamma,gas,LossMode::CoulombLangmuir,.01,
                       cgs::kappa_e*s.mass[0],1.1,relative_tau_step);
 double norm=budget.heated/(s.beam_mono_gamma-1);
 if(paper_hat_heating){
  project_track_hats(tr,s.mass,norm,s.heating,s.ram);
  for(int d=1;d<D;++d)s.target_flux[d]=s.target_flux[d-1]-.5*(s.heating[d]+s.heating[d-1])*(s.mass[d]-s.mass[d-1]);
  return;
 }
 if(!paper_pointwise_heating) {
   // Integrate the SAME continuous stopping trajectory over nodal control
   // volumes. Their widths equal the mass-trapezoid weights. Thus the nodal
   // interpolant conserves the actual energy and momentum deposited by the
   // beam, including the unresolved gamma=1.1 termination, without rescaling
   // the incident beam or fitting a profile. This remapping is our numerical
   // choice; the paper does not specify its coarse-grid interpolation weights.
   std::vector<double> edge(D+1);edge[0]=s.mass.front();edge[D]=s.mass.back();
   for(int d=1;d<D;++d)edge[d]=.5*(s.mass[d-1]+s.mass[d]);
   for(int d=0;d<D;++d) {
     double left=tr.kinetic(cgs::kappa_e*edge[d]),right=tr.kinetic(cgs::kappa_e*edge[d+1]);
     double width=edge[d+1]-edge[d];
     s.heating[d]=norm*(left-right)/width;
     s.ram[d]=norm/cgs::c*(std::sqrt(left*(left+2))-std::sqrt(right*(right+2)))/width;
   }
   // Salmi (21), using the same piecewise-linear Q represented on the
   // radiation grid. The bottom intrinsic flux is recovered to roundoff.
   s.target_flux[0]=budget.total();
   for(int d=1;d<D;++d)s.target_flux[d]=s.target_flux[d-1]-.5*(s.heating[d]+s.heating[d-1])*(s.mass[d]-s.mass[d-1]);
   return;
 }

 for(int d=0;d<D;++d) {
   double tau=cgs::kappa_e*s.mass[d],kin=tr.kinetic(tau);
   s.target_flux[d]=budget.intrinsic+norm*kin;
   if(tau<tr.tau.back()) {
     double gamma=1+kin,ne=s.rho[d]/cgs::mp;
     double phi=phi_coulomb_langmuir(gamma,ne);
     s.heating[d]=norm*cgs::kappa_e*phi;
     s.ram[d]=s.heating[d]/cgs::c*gamma/std::sqrt(gamma*gamma-1);
   }
 }
 // The paper thermalizes the remaining kinetic energy at gamma=1.1.
 // Represent this point contribution on the bracketing nodal hat functions;
 // its mass-trapezoid integral is exactly 0.1*Ndot*mec^2, with no broad tail.
 double mc=tr.tau.back()/cgs::kappa_e;
 auto it=std::upper_bound(s.mass.begin(),s.mass.end(),mc);
 if(it!=s.mass.end() && it!=s.mass.begin()) {
   int hi=it-s.mass.begin(),lo=hi-1;double f=(mc-s.mass[lo])/(s.mass[hi]-s.mass[lo]);
   for(int d:{lo,hi}) {
     double width=d==0?.5*(s.mass[1]-s.mass[0]):d==D-1?.5*(s.mass[d]-s.mass[d-1]):.5*(s.mass[d+1]-s.mass[d-1]);
     double fraction=d==lo?1-f:f;
     s.heating[d]+=norm*.1*fraction/width;
     s.ram[d]+=norm/cgs::c*std::sqrt(1.1*1.1-1)*fraction/width;
   }
 }
 if(paper_normalize_point_heating) {
   // Preserve the sampled Eq.13 profile while making its finite-grid
   // integral equal the prescribed incident energy/momentum. This is an
   // inferred numerical reconstruction, NOT a known private author recipe.
   // Neither normalization uses any plotted author ordinate.
   double energy=0,momentum=0;
   for(int d=1;d<D;++d){double width=s.mass[d]-s.mass[d-1];energy+=.5*(s.heating[d]+s.heating[d-1])*width;momentum+=.5*(s.ram[d]+s.ram[d-1])*width;}
   if(!(energy>0 && momentum>0))throw std::runtime_error("Beam unresolved on heating grid");
   double energy_factor=budget.heated/energy;
   double momentum_factor=norm/cgs::c*std::sqrt(s.beam_mono_gamma*s.beam_mono_gamma-1)/momentum;
   for(int d=0;d<D;++d){s.heating[d]*=energy_factor;s.ram[d]*=momentum_factor;}
   s.target_flux[0]=budget.total();
   for(int d=1;d<D;++d)s.target_flux[d]=s.target_flux[d-1]-.5*(s.heating[d]+s.heating[d-1])*(s.mass[d]-s.mass[d-1]);
 }
}
}
