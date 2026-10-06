#include "balanced_temperature.hpp"
#include "lapack_solver.hpp"
#include "salmi/constants.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <numeric>
#include <chrono>
namespace salmi {
namespace {
using Vec=std::vector<double>;using Mat=std::vector<Vec>;
struct State {Vec m,T,k,F,target,Q,cool,L,response;};
State state(const HeatedState&s){return {s.mass,s.temperature,s.rosseland,s.flux,s.target_flux,s.heating,s.cooling,s.local_derivative,s.temperature_response.empty()?s.local_derivative:s.temperature_response};}
double norm(const Vec&v){double q=0;for(double x:v)q=std::hypot(q,x);return q;}
Vec optical_depth(const State&s){int n=s.m.size();Vec t(n);t[0]=cgs::kappa_e*s.k[0]*s.m[0];for(int i=1;i<n;++i)t[i]=t[i-1]+cgs::kappa_e*.5*(s.k[i]+s.k[i-1])*(s.m[i]-s.m[i-1]);return t;}
// One-sided Jacobi SVD; column rotations avoid squaring the condition number.
// Retain the same relative singular-value cutoff as the NumPy prototype.
Vec least_squares_jacobi(const Mat&a,const Vec&rhs){
 int m=a.size(),n=a[0].size();Mat b(n,Vec(m)),v(n,Vec(n));
 for(int j=0;j<n;++j){v[j][j]=1;for(int i=0;i<m;++i)b[j][i]=a[i][j];}
 bool settled=false;
 for(int sweep=0;sweep<80;++sweep){int rotations=0;
  for(int p=0;p<n;++p)for(int q=p+1;q<n;++q){
   double aa=0,bb=0,ab=0;for(int i=0;i<m;++i){aa+=b[p][i]*b[p][i];bb+=b[q][i]*b[q][i];ab+=b[p][i]*b[q][i];}
   if(aa==0||bb==0||std::abs(ab)<=1e-12*std::sqrt(aa)*std::sqrt(bb))continue;
   double z=(bb-aa)/(2*ab),t=std::copysign(1./(std::abs(z)+std::hypot(1.,z)),z),c=1/std::sqrt(1+t*t),sn=c*t;
   for(int i=0;i<m;++i){double x=b[p][i],y=b[q][i];b[p][i]=c*x-sn*y;b[q][i]=sn*x+c*y;}
   for(int i=0;i<n;++i){double x=v[p][i],y=v[q][i];v[p][i]=c*x-sn*y;v[q][i]=sn*x+c*y;}++rotations;
  }
  if(!rotations){settled=true;break;}
 }
 if(!settled)throw std::runtime_error("Temperature Jacobian SVD did not converge");
 Vec square(n);double maximum=0;for(int j=0;j<n;++j){for(double x:b[j])square[j]+=x*x;maximum=std::max(maximum,square[j]);}
 Vec x(n,0);for(int j=0;j<n;++j)if(square[j]>maximum*1e-12){double c=0;for(int i=0;i<m;++i)c+=b[j][i]*rhs[i];c/=square[j];for(int i=0;i<n;++i)x[i]+=v[j][i]*c;}
 for(double z:x)if(!std::isfinite(z))throw std::runtime_error("Nonfinite temperature least-squares step");return x;
}
Vec least_squares(const Mat&a,const Vec&rhs){
 auto start=std::chrono::steady_clock::now();const char* dll=std::getenv("SALMI_LAPACK_DLL");bool lapack=dll && *dll;
 Vec x=lapack?lapack_temperature_solve(a,rhs):least_squares_jacobi(a,rhs);
 if(const char* path=std::getenv("SALMI_TEMPERATURE_MATRIX_DUMP")){
  if(*path){std::ofstream f(path,std::ios::binary);int dims[]={int(a.size()),int(a.front().size())};
   f.write(reinterpret_cast<char*>(dims),sizeof(dims));for(const auto& row:a)f.write(reinterpret_cast<const char*>(row.data()),row.size()*sizeof(double));
   f.write(reinterpret_cast<const char*>(rhs.data()),rhs.size()*sizeof(double));f.write(reinterpret_cast<const char*>(x.data()),x.size()*sizeof(double));
   if(!f)throw std::runtime_error("Cannot save temperature matrix diagnostic");}
 }
 if(std::getenv("SALMI_TEMPERATURE_TIMING"))std::cout<<"TEMPERATURE_LINEAR_SOLVE backend="<<(lapack?"DGELSD":"Jacobi")<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
 return x;
}
struct Controller {
 bool initialized=false,hot=false;int last_probe_iteration=-1000000;HeatedState accepted_full;State accepted;Vec x,r,width;Mat J;std::vector<int> energy,gradient;double radius=0,merit=0;int failures=0;
 void initialize(const State&s,double cap,bool hot_response){
  initialized=true;hot=hot_response;accepted=s;radius=cap;failures=0;int n=s.m.size();width.resize(n);auto tau=optical_depth(s);
  energy.clear();gradient.clear();
  for(int d=0;d<n;++d){width[d]=d==0?s.m[1]-s.m[0]:d==n-1?s.m[n-1]-s.m[n-2]:.5*(s.m[d+1]-s.m[d-1]);if(tau[d]<1)energy.push_back(d);if(d&&tau[d]>1)gradient.push_back(d);}
  x.resize(n);for(int d=0;d<n;++d)x[d]=std::log(s.T[d]);r=residual(s);merit=norm(r);J=jacobian(s);
 }
 bool refresh_regions(const State& now){
  // Acceptance uses the CURRENT Rosseland depth. Follow the same regions in
  // the temperature residual; a cell crossing tau_R=1 must not be omitted.
  // Keep learned Jacobian rows with unchanged physical identities. Changing
  // row counts without this remapping would corrupt the Broyden update.
  const int n=now.m.size();auto tau=optical_depth(now);
  std::vector<int> next_energy,next_gradient;
  for(int d=0;d<n;++d){if(tau[d]<1)next_energy.push_back(d);if(d&&tau[d]>1)next_gradient.push_back(d);}
  if(next_energy==energy && next_gradient==gradient)return false;
  auto old_energy=energy,old_gradient=gradient;Mat old_J=std::move(J);
  energy=std::move(next_energy);gradient=std::move(next_gradient);
  J=jacobian(accepted);
  for(int d=0;d<n;++d)J[d]=std::move(old_J[d]);
  auto preserve=[&](const std::vector<int>& before,const std::vector<int>& after,int old_offset,int new_offset){
   for(size_t i=0;i<after.size();++i){auto p=std::lower_bound(before.begin(),before.end(),after[i]);
    if(p!=before.end() && *p==after[i])J[new_offset+i]=std::move(old_J[old_offset+(p-before.begin())]);
   }
  };
  preserve(old_energy,energy,n,n);
  preserve(old_gradient,gradient,n+old_energy.size(),n+energy.size());
  // Reevaluate the accepted state in the SAME residual coordinates as the
  // current trial before comparing merits or taking a secant difference.
  r=residual(accepted);merit=norm(r);
  if(!paper_response.checkpoint_root.empty()){
   std::string path=paper_response.checkpoint_root+"/temperature_regions.csv";bool first=!std::ifstream(path);std::ofstream f(path,std::ios::app);
   if(first)f<<"iteration,old_upper_count,new_upper_count,old_deep_count,new_deep_count\n";
   f<<paper_response.outer_iteration<<','<<old_energy.size()<<','<<energy.size()<<','<<old_gradient.size()<<','<<gradient.size()<<'\n';
   if(!f)throw std::runtime_error("Cannot save temperature region audit");
  }
  std::cout<<"BALANCED_REGIONS upper="<<old_energy.size()<<"->"<<energy.size()<<" deep="<<old_gradient.size()<<"->"<<gradient.size()<<std::endl;
  return true;
 }
 double energy_weight(const State&s,int d)const{return paper_response.balanced_energy_scale*std::max(1000.,s.L[d]*width[d]/s.target[0]*10000.);}
 Vec residual(const State&s)const{
  int n=s.m.size();Vec a;auto tau=optical_depth(s);for(int d=0;d<n;++d)a.push_back((s.F[d]-s.target[d])/s.target[d]*1000.);
  for(int d:energy)a.push_back((s.cool[d]-s.Q[d])/s.L[d]*energy_weight(s,d));
  for(int d:gradient)a.push_back(((s.F[d]-s.target[d])-(s.F[d-1]-s.target[d-1]))/(tau[d]-tau[d-1])/s.target[0]*10000.);return a;
 }
 Mat jacobian(const State&s)const{
  int n=s.m.size();Vec response(n),factor(n,1),dt(n-1);Mat flux(n,Vec(n)),mean(n,Vec(n));double teff=std::pow(s.target[0]/cgs::sigma_sb,.25);
  for(int d=0;d<n;++d){response[d]=4*cgs::sigma_sb/cgs::pi*std::pow(s.T[d],4);
   if(hot&&s.Q[d]>0&&s.T[d]>3*teff){response[d]=std::min(response[d],s.L[d]/(4*cgs::pi*cgs::kappa_e*s.k[d]));double theta=cgs::kb*s.T[d]/cgs::mec2;factor[d]=-4*theta*(1+8*theta)/s.k[d];}
   if(d)dt[d-1]=cgs::kappa_e*.5*(s.k[d]+s.k[d-1])*(s.m[d]-s.m[d-1]);
  }
  const double mu[]={.23861918608319690863,.66120938646626451366,.93246951420315202781};
  const double wt[]={.46791393457269104739,.36076157304813860757,.17132449237917034504};
  auto weights=[](double t){double tr=std::exp(-t),one=-std::expm1(-t),local=t<1e-4?t/2-t*t/6+t*t*t/24-t*t*t*t/120:1-one/t;return Vec{tr,one-local,local};};
  for(int a=0;a<3;++a){Mat in(n,Vec(n)),out(n,Vec(n));
   for(int d=1;d<n;++d){auto w=weights(dt[d-1]/mu[a]);for(int j=0;j<n;++j)in[d][j]=w[0]*in[d-1][j];in[d][d-1]+=w[1]*response[d-1];in[d][d]+=w[2]*response[d];}
   out[n-1]=in[n-1];out[n-1][n-1]+=2*mu[a]*response[n-1]/dt[n-2];out[n-1][n-2]-=2*mu[a]*response[n-2]/dt[n-2];
   for(int d=n-2;d>=0;--d){auto w=weights(dt[d]/mu[a]);for(int j=0;j<n;++j)out[d][j]=w[0]*out[d+1][j];out[d][d+1]+=w[1]*response[d+1];out[d][d]+=w[2]*response[d];}
   for(int d=0;d<n;++d)for(int j=0;j<n;++j){flux[d][j]+=2*cgs::pi*wt[a]*mu[a]*(out[d][j]-in[d][j]);mean[d][j]+=.5*wt[a]*(out[d][j]+in[d][j]);}
  }
  Mat a;for(int d=0;d<n;++d){Vec row(n);for(int j=0;j<n;++j)row[j]=flux[d][j]/s.target[d]*1000.;a.push_back(row);}
  for(int d:energy){Vec row(n);for(int j=0;j<n;++j)row[j]=((j==d?s.response[d]/s.L[d]:0.)-mean[d][j]*factor[d]/response[d])*energy_weight(s,d);a.push_back(row);}
  for(int d:gradient){Vec row(n);for(int j=0;j<n;++j)row[j]=(flux[d][j]-flux[d-1][j])/dt[d-1]/s.target[0]*10000.;a.push_back(row);}return a;
 }
 double update(HeatedState&s,double cap,bool hot_response,const TemperatureProbe& probe){
  State now=state(s);int n=now.m.size();bool took=true;double measured;
  if(!initialized){initialize(now,cap,hot_response);measured=merit;if(paper_response.fd_columns)accepted_full=s;}
  else{
   refresh_regions(now);
   Vec rn=residual(now),dx(n);double square=0;for(int d=0;d<n;++d){dx[d]=std::log(now.T[d])-x[d];square+=dx[d]*dx[d];}measured=norm(rn);
   if(square>1e-20)for(size_t i=0;i<J.size();++i){double error=rn[i]-r[i];for(int d=0;d<n;++d)error-=J[i][d]*dx[d];for(int d=0;d<n;++d)J[i][d]+=error*dx[d]/square;}
   took=measured<merit;
   if(took){accepted=now;if(paper_response.fd_columns)accepted_full=s;for(int d=0;d<n;++d)x[d]=std::log(now.T[d]);r.swap(rn);merit=measured;failures=0;radius=std::min(cap,1.3*radius);}
   else{++failures;radius*=.5;if(failures>=8||radius<1e-7){J=jacobian(accepted);radius=std::min(cap,.001);failures=0;}}
  }
  if(paper_response.fd_columns && (last_probe_iteration<0 || paper_response.outer_iteration-last_probe_iteration>=paper_response.fd_every)) {
   if(!probe)throw std::runtime_error("Numerical Jacobian probes requested without physical evaluator");
   // Reconverge the baseline as well as both perturbations. Otherwise a loose
   // construction radiation field would contaminate every numerical column.
   auto base=probe(accepted_full,-1,0);accepted_full=base;accepted=state(base);s=base;
   refresh_regions(accepted);
   for(int d=0;d<n;++d)x[d]=std::log(accepted.T[d]);r=residual(accepted);merit=norm(r);
   std::vector<int> selected=energy;
   std::stable_sort(selected.begin(),selected.end(),[&](int a,int b){return std::abs((accepted.cool[a]-accepted.Q[a])/accepted.L[a])*energy_weight(accepted,a)>std::abs((accepted.cool[b]-accepted.Q[b])/accepted.L[b])*energy_weight(accepted,b);});
   for(int k=0;k<std::min(paper_response.fd_columns,int(selected.size()));++k){int d=selected[k];double h=paper_response.fd_step;
    auto minus=probe(base,d,-h),plus=probe(base,d,h);auto rm=residual(state(minus)),rp=residual(state(plus));Vec measured_column(r.size()),difference(r.size()),predicted(r.size());
    for(size_t i=0;i<r.size();++i){measured_column[i]=(rp[i]-rm[i])/(2*h);predicted[i]=J[i][d];difference[i]=predicted[i]-measured_column[i];if(!std::isfinite(measured_column[i]))throw std::runtime_error("Nonfinite numerical Jacobian column");}
    double size=norm(measured_column);if(!(size>0))throw std::runtime_error("Zero numerical Jacobian column");
    if(!paper_response.checkpoint_root.empty()){std::string path=paper_response.checkpoint_root+"/fd_columns.csv";bool first=!std::ifstream(path);std::ofstream f(path,std::ios::app);if(first)f<<"iteration,depth,mass,logT_step,predicted_column_norm,measured_column_norm,relative_column_error\n";f<<std::setprecision(17)<<paper_response.outer_iteration<<','<<d<<','<<accepted.m[d]<<','<<h<<','<<norm(predicted)<<','<<size<<','<<norm(difference)/size<<'\n';f.flush();if(!f)throw std::runtime_error("Cannot save Jacobian audit");}
    std::cout<<"FD_COLUMN depth="<<d<<" relative_prediction_error="<<norm(difference)/size<<std::endl;
    for(size_t i=0;i<r.size();++i)J[i][d]=measured_column[i];
   }
   last_probe_iteration=paper_response.outer_iteration;
  }
  Vec rhs=r;for(double&v:rhs)v=-v;auto step=least_squares(J,rhs);double maximum=0;for(double v:step)maximum=std::max(maximum,std::abs(v));double scale=std::min(1.,radius/std::max(maximum,1e-100)),change=0;
  Vec applied=step;for(double& v:applied)v*=scale;
  if(paper_response.compare_clipped_step){Vec clipped=step;for(double& v:clipped)v=std::max(-radius,std::min(radius,v));
   auto model_norm=[&](const Vec& dx){Vec predicted=r;for(size_t i=0;i<J.size();++i)for(int d=0;d<n;++d)predicted[i]+=J[i][d]*dx[d];return norm(predicted);};
   if(model_norm(clipped)<model_norm(applied))applied.swap(clipped);
  }
  double applied_max=0;
  for(int d=0;d<n;++d){applied_max=std::max(applied_max,std::abs(applied[d]));double proposed=std::exp(x[d]+applied[d]);change=std::max(change,std::abs(proposed/s.temperature[d]-1));s.temperature[d]=proposed;}
  if(!paper_response.checkpoint_root.empty()){
   std::string path=paper_response.checkpoint_root+"/temperature_solver.csv";bool first=paper_response.outer_iteration==0;std::ofstream f(path,first?std::ios::out:std::ios::app);if(!f)throw std::runtime_error("Cannot write temperature solver diagnostics");
   if(first)f<<"evaluation,trial_accepted,measured_joint_norm,accepted_joint_norm,trust_radius,proposed_max_log_step\n";
   f<<std::setprecision(17)<<paper_response.outer_iteration<<','<<(took?1:0)<<','<<measured<<','<<merit<<','<<radius<<','<<applied_max<<'\n';
  }
  return change;
 }
};
Controller controller;
}
double balanced_temperature(HeatedState&s,double cap,bool hot_response,const TemperatureProbe& probe){return controller.update(s,cap,hot_response,probe);}
}
