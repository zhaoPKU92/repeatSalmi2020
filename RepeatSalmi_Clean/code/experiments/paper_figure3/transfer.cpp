#include "paper.hpp"
#include "salmi/anderson.hpp"
#include "salmi/compton.hpp"
#include "salmi/constants.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif
#include <stdexcept>
namespace salmi {
PaperResponse paper_response;
namespace {
void radiation_seed(const HeatedState& s,const PhotonGrid& g,int inner,double residual) {
 if(paper_response.checkpoint_root.empty() || paper_response.outer_iteration!=0)return;
 std::string path=paper_response.checkpoint_root+"/radiation_seed_0_"+std::to_string(inner);
#ifdef _WIN32
 int result=_mkdir(path.c_str());
#else
 int result=mkdir(path.c_str(),0777);
#endif
 if(result)throw std::runtime_error("Cannot create fresh radiation seed directory");
 {std::ofstream f(path+"/structure.csv");f<<std::setprecision(17)<<"mass,T,P,rho,Qplus,gram,Ftarget\n";
  for(size_t d=0;d<s.mass.size();++d)f<<s.mass[d]<<','<<s.temperature[d]<<','<<s.pressure[d]<<','<<s.rho[d]<<','<<s.heating[d]<<','<<s.ram[d]<<','<<s.target_flux[d]<<'\n';
  f.flush();if(!f)throw std::runtime_error("Radiation seed structure write failed");}
 {std::ofstream f(path+"/occupation.bin",std::ios::binary);int dim[]={int(s.mass.size()),int(g.x.size()),int(g.mu.size())};
  f.write(reinterpret_cast<const char*>(dim),sizeof(dim));f.write(reinterpret_cast<const char*>(s.occupation.data()),s.occupation.size()*sizeof(double));
  f.flush();if(!f)throw std::runtime_error("Radiation seed field write failed");}
 double lo=.5*(std::sqrt(g.dx.front()*g.dx.front()+4*g.x.front()*g.x.front())-g.dx.front())*cgs::mec2/cgs::kev;
 double hi=.5*(std::sqrt(g.dx.back()*g.dx.back()+4*g.x.back()*g.x.back())+g.dx.back())*cgs::mec2/cgs::kev;
 {std::ofstream f(path+"/status.json");f<<std::setprecision(17)<<"{\"iteration_checkpoint\":true,\"discrete_converged\":false,\"partial_radiation\":true,\"depth_nodes\":"<<s.mass.size()<<",\"energy_bins\":"<<g.x.size()<<",\"directions_per_hemisphere\":"<<g.mu.size()/2<<",\"energy_min_keV\":"<<lo<<",\"energy_max_keV\":"<<hi<<",\"transfer_residual\":"<<residual<<",\"inner_iteration\":"<<inner<<",\"scope\":\"Unconverged radiation seed at fixed temperature; not an atmosphere result. Use parent invocation for physics.\"}\n";f.flush();if(!f)throw std::runtime_error("Radiation seed metadata write failed");}
 {std::ofstream f(path+"/complete.json");f<<"{\"state_files_closed\":true,\"partial_radiation\":true}\n";f.flush();if(!f)throw std::runtime_error("Radiation seed commit failed");}
}
}

TransferDiagnostics atmosphere_transfer(HeatedState& s,const PhotonGrid& g,
 const std::vector<SparseScattering>& kernels,double /*bottom_flux*/,int maxit,double tol) {
 const int D=s.mass.size(), E=g.x.size(), A=g.mu.size(), N=g.size();
 const size_t M=size_t(D)*N;
 const bool full_departure=std::getenv("SALMI_PAPER_FULL_DEPARTURE")!=nullptr;
 const char* departure_option=std::getenv("SALMI_RT_RECONCILE_DEPARTURE");
 const bool reconcile_departure=!departure_option || std::string(departure_option)!="0";
 std::vector<double> B(D*E),dB(D*E),ff(D*E),cs(D*E),chi(M),source_delta(M),
   source(M),scatter_fraction(M),rowI(M),beta(M),dep(M),nextdep(M),next(M),bottom(N);
 const double scaleI=occupation_intensity_x(1,1);
 // Optional absolute INTENSITY normalization. This never clips the field.
 // Use the fixed total effective temperature, not the hot local gas temperature.
 std::vector<double> intensity_floor_n(E,0);
 if(paper_response.rt_planck_intensity_floor_fraction>0) {
   if(s.target_flux.empty() || !(s.target_flux.front()>0))
     throw std::runtime_error("Intensity normalization requires positive target flux");
   const double teff=std::pow(s.target_flux.front()/cgs::sigma_sb,.25);
   const double up=2.821439372122079,xp=up*cgs::kb*teff/cgs::mec2;
   for(int e=0;e<E;++e)intensity_floor_n[e]=paper_response.rt_planck_intensity_floor_fraction*std::pow(xp/g.x[e],3)/std::expm1(up);
   std::cout<<"RT intensity normalization fraction="<<paper_response.rt_planck_intensity_floor_fraction<<" reference_Teff="<<teff<<std::endl;
 }
 s.rosseland.assign(D,0);paper_response.denominator.assign(D,0);
 paper_response.opacity.resize(D*E);paper_response.spectral_flux.assign(D*E,0);
 paper_response.ali_guards=0;
 #pragma omp parallel for schedule(dynamic)
 for(int d=0;d<D;++d) {
   double rb=0,rw=0,theta=cgs::kb*s.temperature[d]/cgs::mec2;
   std::vector<double> unitI(N),zero(E,0);
   for(int e=0;e<E;++e) {
     int j=d*E+e;double u=g.x[e]/theta;
     B[j]=u>700?0:1/std::expm1(u);dB[j]=B[j]*(1+B[j])*u;
     ff[j]=s.absorption(g.x[e]*cgs::mec2/cgs::kev,d)/cgs::kappa_e;
     cs[j]=thermal_scattering_rate(g.x[e],theta,32);
     paper_response.opacity[j]=cgs::kappa_e*(ff[j]+cs[j]);
     double wt=std::pow(g.x[e],3)*g.dx[e]*dB[j]; rb+=wt;rw+=wt/(ff[j]+cs[j]);
     for(int a=0;a<A;++a)unitI[e*A+a]=1/std::pow(g.x[e],3);
   }
   s.rosseland[d]=rb/rw;
   // B.14: x^2 integral R dx1/x1^2 dmu1. Convert the occupation gain to I_x.
   auto ri=kernels[d].linearized_gain(unitI.data(),zero.data());
   for(int i=0;i<N;++i)rowI[size_t(d)*N+i]=ri[i]*std::pow(g.x[i/A],3);
 }
 if(s.occupation.size()!=M) {
   s.occupation.resize(M);
   for(int d=0;d<D;++d)for(int i=0;i<N;++i)s.occupation[size_t(d)*N+i]=B[d*E+i/A];
   for(int e=0;e<E;++e)for(int a=0;a<A/2;++a)s.occupation[e*A+a]=0;
 }
 for(int d=0;d<D;++d)for(int i=0;i<N;++i)dep[size_t(d)*N+i]=s.occupation[size_t(d)*N+i]-B[d*E+i/A];
 AndersonMixing extrapolation(paper_response.rt_history);
 TransferDiagnostics result{1,0,0};
 for(int it=0;it<maxit;++it) {
   #pragma omp parallel for schedule(dynamic)
   for(int d=0;d<D;++d) {
     auto c=kernels[d].terms(&s.occupation[size_t(d)*N]);
     for(int i=0;i<N;++i) {
       size_t p=size_t(d)*N+i;int j=d*E+i/A;
       // SPW12 (12),(14): full stimulated LOSS opacity, source (1+n)*gain.
       // Algebraically exact departure evaluation avoids subtracting two LTE terms.
       chi[p]=ff[j]+c.loss[i];scatter_fraction[p]=c.loss[i]/chi[p];
       source_delta[p]=c.gain[i]; // Temporary gain; convert below after opacities at all depths are known.
       source[p]=(ff[j]*B[j]+(1+s.occupation[p])*c.gain[i])/chi[p];
     }
   }
   #pragma omp parallel for schedule(dynamic)
   for(int d=0;d<D;++d) {
     // A stable SC sweep can request the departure source at this node
     // or an adjacent upwind node. The absorption test is a conservative
     // superset of all branches that use it; no used contribution is omitted.
     std::vector<double> mask(E,full_departure?1.:0.);
     if(!full_departure)for(int e=0;e<E;++e)
       for(int c=std::max(0,d-1);c<=std::min(D-1,d+1);++c)
         for(int a=0;a<A;++a)if(ff[c*E+e]>.05*chi[size_t(c)*N+e*A+a])mask[e]=1;
     auto dg=kernels[d].linearized_gain(&dep[size_t(d)*N],&B[d*E],mask.data());
     for(int i=0;i<N;++i) {
       size_t p=size_t(d)*N+i;
       source_delta[p]=mask[i/A]?(dg[i]+source_delta[p]*dep[p])/chi[p]:0;
     }
   }
   #pragma omp parallel for schedule(static)
   for(int e=0;e<E;++e)for(int a=0;a<A/2;++a) {
     int i=e*A+a,ip=e*A+A-1-a;double mu=-g.mu[a];
     nextdep[i]=-B[e];next[i]=0;beta[i]=0;
     for(int d=1;d<D;++d) {
       size_t p=size_t(d)*N+i,u=p-N;
       double dt=(s.mass[d]-s.mass[d-1])*cgs::kappa_e*.5*(chi[p]+chi[u])/mu;
       auto w=short_characteristic_weights(dt);beta[p]=w.local;
       if(dt>5 && ff[d*E+e]>.05*chi[p]) {
         nextdep[p]=w.transmission*nextdep[u]+(w.transmission+w.upstream)*(B[(d-1)*E+e]-B[d*E+e])+w.upstream*source_delta[u]+w.local*source_delta[p];
         next[p]=B[d*E+e]+nextdep[p];
       } else {
         next[p]=w.transmission*next[u]+w.upstream*source[u]+w.local*source[p];
         nextdep[p]=next[p]-B[d*E+e];
       }
     }
     size_t p=size_t(D-1)*N+ip,im=size_t(D-1)*N+i;
     double dt=(s.mass[D-1]-s.mass[D-2])*cgs::kappa_e*.5*(chi[p]+chi[p-N])/mu;
     // SPW12 B.10. The bottom flux is solved, not renormalized to the target.
     bottom[ip]=2*(B[(D-1)*E+e]-B[(D-2)*E+e])/dt;
     nextdep[p]=nextdep[im]+bottom[ip];next[p]=next[im]+bottom[ip];beta[p]=1;
     for(int d=D-2;d>=0;--d) {
       size_t q=size_t(d)*N+ip,u=q+N;
       double t=(s.mass[d+1]-s.mass[d])*cgs::kappa_e*.5*(chi[q]+chi[u])/mu;
       auto w=short_characteristic_weights(t);beta[q]=w.local;
       if(t>5 && ff[d*E+e]>.05*chi[q]) {
         nextdep[q]=w.transmission*nextdep[u]+(w.transmission+w.upstream)*(B[(d+1)*E+e]-B[d*E+e])+w.upstream*source_delta[u]+w.local*source_delta[q];
         next[q]=B[d*E+e]+nextdep[q];
       } else {
         next[q]=w.transmission*next[u]+w.upstream*source[u]+w.local*source[q];
         nextdep[q]=next[q]-B[d*E+e];
       }
     }
   }
   double error=0,worst_mean=0;int guards=0,worst_depth=0,worst_energy=0;
   for(int d=0;d<D;++d)for(int e=0;e<E;++e) {
     double oldj=0,newj=0;
     for(int a=0;a<A;++a) {
       size_t p=size_t(d)*N+e*A+a;
       oldj+=g.dmu[a]*s.occupation[p];newj+=g.dmu[a]*next[p];
     }
     double relative=std::abs(newj-oldj)/std::max({std::abs(newj),std::abs(oldj),paper_response.rt_occupation_floor,2*intensity_floor_n[e]});
     if(relative>error){error=relative;worst_depth=d;worst_energy=e;worst_mean=.5*newj;}
     for(int a=0;a<A;++a) {
       size_t p=size_t(d)*N+e*A+a,r=size_t(d)*N+e*A+A-1-a;
       double lambda=.25*(beta[p]+beta[r])*rowI[p];
       double denominator=1-scatter_fraction[p]*lambda;
       if(denominator<.1){denominator=.1;++guards;}
       bool absolute=s.occupation[p]<.1*B[d*E+e] || next[p]<.1*B[d*E+e];
       double raw=absolute?next[p]-s.occupation[p]:nextdep[p]-dep[p];
       double correction=raw/denominator*s.transfer_relaxation;
       double candidate=absolute?s.occupation[p]+correction:B[d*E+e]+dep[p]+correction;
       if(candidate<0){correction=.5*raw;candidate=s.occupation[p]+correction;++guards;}
       nextdep[p]=absolute?candidate-B[d*E+e]:dep[p]+correction;
       next[p]=candidate;
       if(next[p]<0 || !std::isfinite(next[p])) {
         std::ostringstream msg;msg<<std::setprecision(17)<<"Invalid linear-SC intensity: d="<<d<<" e="<<e<<" a="<<a
           <<" m="<<s.mass[d]<<" E_keV="<<g.x[e]*cgs::mec2/cgs::kev<<" T="<<s.temperature[d]
           <<" old="<<s.occupation[p]<<" candidate="<<next[p]<<" B="<<B[d*E+e]<<" dep="<<dep[p]
           <<" raw="<<raw<<" correction="<<correction<<" denominator="<<denominator<<" absolute="<<absolute
           <<" source="<<source[p]<<" source_delta="<<source_delta[p]<<" chi="<<chi[p]<<" ff="<<ff[d*E+e];
         throw std::runtime_error(msg.str());
       }
     }
   }
   for(int e=0;e<E;++e)for(int a=0;a<A/2;++a) {
     next[e*A+a]=0;nextdep[e*A+a]=-B[e];
     size_t i=size_t(D-1)*N+e*A+a,p=size_t(D-1)*N+e*A+A-1-a;
     nextdep[p]=nextdep[i]+bottom[e*A+A-1-a];next[p]=next[i]+bottom[e*A+A-1-a];
   }
   // Optional acceleration of the SAME published full-opacity ALI map.
   // Acceptance still tests its ordinary formal-solution residual, not the
   // extrapolated step size. The terminating update is always plain ALI.
   if(paper_response.rt_extrapolate && it>=2 && error>=tol) {
     std::vector<double> state(2*M),correction(2*M),scale(M);
     for(size_t p=0;p<M;++p) {
       state[p]=s.occupation[p];state[M+p]=dep[p];
       correction[p]=next[p]-s.occupation[p];correction[M+p]=nextdep[p]-dep[p];
       scale[p]=std::max({std::abs(next[p]),std::abs(s.occupation[p]),paper_response.rt_occupation_floor,intensity_floor_n[(p%N)/A]});
     }
     auto candidate=extrapolation.step(state,correction,scale,1,int(M));
     double max_inconsistency=0;
     for(size_t p=0;p<M;++p){
       next[p]=candidate[p];nextdep[p]=candidate[M+p];
       const double planck=B[(p/N)*E+(p%N)/A];
       if(next[p]<=.5*planck || next[p]>=2*planck) {
         const double exact_departure=next[p]-planck;
         max_inconsistency=std::max(max_inconsistency,std::abs(nextdep[p]-exact_departure)/std::max({std::abs(next[p]),planck,1e-300}));
         // Far from LTE, subtraction is well conditioned. Near LTE preserve
         // the independent departure to resolve the tiny diffusion flux.
         if(reconcile_departure)nextdep[p]=exact_departure;
       }
     }
     if((it+1)%25==0 && max_inconsistency>1e-6)
       std::cout<<"  LTE departure inconsistency="<<max_inconsistency<<" reconciled="<<reconcile_departure<<std::endl;
     for(int e=0;e<E;++e)for(int a=0;a<A/2;++a) {
       next[e*A+a]=0;nextdep[e*A+a]=-B[e];
       size_t i=size_t(D-1)*N+e*A+a,p=size_t(D-1)*N+e*A+A-1-a;
       nextdep[p]=nextdep[i]+bottom[e*A+A-1-a];next[p]=next[i]+bottom[e*A+A-1-a];
     }
   }
   dep.swap(nextdep);s.occupation.swap(next);result={error,it+1,guards};
   paper_response.ali_guards+=guards;
   bool stop_requested=!paper_response.checkpoint_root.empty() && bool(std::ifstream(paper_response.checkpoint_root+"/STOP_REQUEST"));
   if((it+1)%25==0 || stop_requested)std::cout<<"  ALI "<<it+1<<" residual="<<error<<" at m="<<s.mass[worst_depth]<<" E_keV="<<g.x[worst_energy]*cgs::mec2/cgs::kev<<" mean_n="<<worst_mean<<" guards="<<guards<<std::endl;
   if((it+1)%200==0 || stop_requested)radiation_seed(s,g,it+1,error);
   if(stop_requested)break;
   if(error<tol)break;
 }
 s.flux.assign(D,0);s.cooling.assign(D,0);s.radiative_acceleration.assign(D,0);s.local_derivative.assign(D,0);s.temperature_response.assign(D,0);
 #pragma omp parallel for schedule(dynamic)
 for(int d=0;d<D;++d) {
   auto c=kernels[d].terms(&s.occupation[size_t(d)*N]);
   std::vector<double> mask(E,full_departure?1.:0.);
   if(!full_departure)for(int e=0;e<E;++e)for(int a=0;a<A;++a)
     if(std::abs(dep[size_t(d)*N+e*A+a])<.1*B[d*E+e])mask[e]=1;
   auto dg=kernels[d].linearized_gain(&dep[size_t(d)*N],&B[d*E],mask.data());
   double radiation=0,free_free_response=0;
   for(int e=0;e<E;++e) {
     double ld=0;bool absolute_flux=false;
     for(int a=0;a<A;++a){double v=s.occupation[size_t(d)*N+e*A+a];if(v<.5*B[d*E+e] || v>2*B[d*E+e])absolute_flux=true;}
     for(int a=0;a<A;++a) {
       size_t p=size_t(d)*N+e*A+a;int j=d*E+e;
       double w=2*cgs::pi*scaleI*g.weight[e*A+a]*g.x[e];
       double net=std::abs(dep[p])<.1*B[j]?
           dg[e*A+a]-(ff[j]+c.loss[e*A+a]-c.gain[e*A+a])*dep[p]:
           c.net[e*A+a]+ff[j]*(B[j]-s.occupation[p]);
       double flux=w*g.mu[a]*(absolute_flux?s.occupation[p]:dep[p]);s.flux[d]+=flux;paper_response.spectral_flux[j]+=flux;
       s.cooling[d]+=cgs::kappa_e*w*net;s.radiative_acceleration[d]-=cgs::kappa_e*w*g.mu[a]*net/cgs::c;
       radiation+=w*s.occupation[p];s.local_derivative[d]+=cgs::kappa_e*w*ff[j]*dB[j];
       if(paper_response.full_thermal_response) {
         // At fixed gas pressure rho is proportional to1/T. With gff=1,
         // dln(kappa_ff)/dlnT = -1.5 - u/(exp(u)-1).
         double u=g.x[e]*cgs::mec2/(cgs::kb*s.temperature[d]);
         free_free_response+=cgs::kappa_e*w*ff[j]*(dB[j]+(1.5+u*B[j])*dep[p]);
       }
       ld+=.5*g.dmu[a]*beta[p];
     }
     int j=d*E+e;double alpha=cs[j]/(ff[j]+cs[j]);
     // SPW12 (28), with dB/dlnT; physical Q requires the common 4*pi factor.
     paper_response.denominator[d]+=4*cgs::pi*scaleI*std::pow(g.x[e],3)*g.dx[e]*cgs::kappa_e*ff[j]*dB[j]*(1-ld)/std::max(1-alpha*ld,1e-100);
   }
   double theta=cgs::kb*s.temperature[d]/cgs::mec2;
   s.local_derivative[d]+=cgs::kappa_e*4*theta*(1+8*theta)*radiation;
   s.temperature_response[d]=paper_response.full_thermal_response?
       free_free_response+compton_cooling_temperature_response(kernels[d],&s.occupation[size_t(d)*N],paper_response.kernel_temperature_step):s.local_derivative[d];
   if(!std::isfinite(s.temperature_response[d]))throw std::runtime_error("Nonfinite thermal response");
 }
 double minimum=1e300;
 for(int e=0;e<E;++e) {
   double ta=0,tt=0;
   for(int d=1;d<D;++d) {
     double dm=s.mass[d]-s.mass[d-1];
     ta+=dm*cgs::kappa_e*.5*(ff[d*E+e]+ff[(d-1)*E+e]);
     tt+=dm*.5*(paper_response.opacity[d*E+e]+paper_response.opacity[(d-1)*E+e]);
   }
   minimum=std::min(minimum,std::sqrt(ta*tt));
 }
 paper_response.bottom_thermalization=minimum;
 return result;
}
}
