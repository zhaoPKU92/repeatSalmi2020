#include "paper.hpp"
#include "atlas_numerics.hpp"
#include "salmi/constants.hpp"
#include <algorithm>
#include <cmath>
namespace salmi {
namespace {
double clip(double a,double b){return std::max(-b,std::min(b,a));}
}
double paper_temperature(HeatedState& s,double relax,double maxstep,double weight,double fluxcap) {
 // Nonconvective Avrett--Krook integrating-factor correction, following the
 // public ATLAS tcorr routine. Salmi (21) replaces the constant target flux.
 int D=s.mass.size(),E=paper_response.opacity.size()/D;
 std::vector<double> kap(D),tau(D),lnk(D),c(D,0),gf(D),dtau(D,0),fluxden(D),lf(D),ll(D),sum(D),integ(D,0);
 // ATLAS divides by the SIGNED radiative flux. Inward flux can occur in
 // an unconverged externally heated atmosphere; flooring it to a positive
 // number reverses the Avrett--Krook correction and introduces huge factors.
 // Only a near-zero singularity is regularized, with a relative floor that
 // preserves sign. Such a state still fails the independent flux criterion.
 for(int d=0;d<D;++d){
   kap[d]=cgs::kappa_e*s.rosseland[d];lnk[d]=std::log(kap[d]);
   double floor=1e-12*std::max(std::abs(s.target_flux[d]),1e-30);
   fluxden[d]=std::copysign(std::max(std::abs(s.flux[d]),floor),s.flux[d]);
 }
 auto dk=atlas_derivative(s.mass,kap),dT=atlas_derivative(s.mass,s.temperature);
 tau[0]=kap[0]*s.mass[0];
 for(int d=1;d<D;++d)tau[d]=tau[d-1]+.5*(kap[d]+kap[d-1])*(s.mass[d]-s.mass[d-1]);
 for(int e=0;e<E;++e) {
   std::vector<double> logchi(D);for(int d=0;d<D;++d)logchi[d]=paper_response.opacity[d*E+e];
   auto dc=atlas_derivative(s.mass,logchi);
   for(int d=2;d<D;++d)c[d]+=paper_response.spectral_flux[d*E+e]*dc[d]/logchi[d]/fluxden[d];
 }
 for(int d=2;d<D;++d)c[d]-=dk[d]/kap[d];
 auto intc=atlas_integral(s.mass,c);
 for(int d=0;d<D;++d)gf[d]=std::exp(clip(intc[d],100));
 std::vector<double> gflux(D);for(int d=0;d<D;++d)gflux[d]=gf[d]*(1-s.target_flux[d]/fluxden[d]);
 dtau=atlas_integral(tau,gflux);
 double teff=std::pow(s.target_flux[0]/cgs::sigma_sb,.25);
 for(int d=0;d<D;++d) {
   double shift=clip(dtau[d]/gf[d],.5*tau[d]);
   lf[d]=d<2?0:clip(-shift*dT[d]/kap[d],fluxcap*s.temperature[d]);
   ll[d]=weight*(s.heating[d]-s.cooling[d])/std::max(paper_response.denominator[d],1e-100)*s.temperature[d];
   ll[d]=clip(ll[d],maxstep*s.temperature[d]);
   if(tau[d]>=1) {ll[d]=0;for(int j=std::max(0,d-5);j<d;++j)ll[j]*=.5;}
 }
 // Surface correction and overlap subtraction, ATLAS tcorr nonconvective case.
 for(int d=0;d<D;++d)sum[d]=lf[d]+ll[d];
 integ=atlas_integral(tau,sum);
 double tav=.5*(linear_interp(tau,integ,2)-linear_interp(tau,integ,.1));
 double surf=clip(.25*s.temperature[0]*(1-s.flux[0]/s.target_flux[0]),.04*teff);
 if(surf*tav<=0)tav=0;if(std::abs(tav)>std::abs(surf))tav=surf;surf-=tav;
 // The papers do not give a hot-skin step controller. Keep the published
 // correction direction, but shrink a persistent per-depth trust radius on
 // sign reversals. This changes convergence, not the equilibrium equations.
 static std::vector<double> previous, radius;static std::vector<int> same;
 if(previous.size()!=size_t(D)){previous.assign(D,0);radius.assign(D,maxstep);same.assign(D,0);}
 double maximum=0;
 for(int d=0;d<D;++d) {
   double raw=relax*(lf[d]+ll[d]+surf)/s.temperature[d];
   if(raw*previous[d]<0){radius[d]=std::max(1e-6,.6*radius[d]);same[d]=0;}
   else if(++same[d]>=4){radius[d]=std::min(maxstep,1.15*radius[d]);same[d]=0;}
   double step=clip(raw,radius[d]);previous[d]=raw;
   s.temperature[d]*=1+step;maximum=std::max(maximum,std::abs(step));
 }
 return maximum;
}
}
