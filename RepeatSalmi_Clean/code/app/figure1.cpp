// Figure 1 only: prescribed atmospheres and individual electron trajectories.
#include "salmi/stopping.hpp"
#include "salmi/constants.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
using namespace salmi;
int main(int argc,char**argv){try{
 if(argc<4 || argc>5)throw std::runtime_error("Usage: figure1 OUTPUT STEP TAU_TOP [MASS_NE_PROFILE]");
 std::string out=argv[1];double step=std::stod(argv[2]),top=std::stod(argv[3]);
 auto csv=[&](std::ofstream& f,std::string name,std::string header){f.open(out+"/"+name);if(!f)throw std::runtime_error("Create output directory first");f<<std::setprecision(17)<<header<<'\n';};
 GreyAtmosphere grey;grey.teff=4.64e6;grey.gravity=1e14;
 std::ofstream gasfile;csv(gasfile,"grey.csv","tau,mass,T,P,ne,rho");
 for(double tau:log_grid(top,100,600)){auto a=grey(tau);gasfile<<tau<<','<<tau/cgs::kappa_e<<','<<a.temperature<<','<<a.pressure<<','<<a.ne<<','<<a.rho<<'\n';}
 std::ofstream losses;csv(losses,"losses.csv","gamma,ne,CL,ep,ee");
 for(double ne:{1e18,1e22,1e26})for(double g:log_grid(1.1,2e5,400)){auto a=loss_components(g,ne);losses<<g<<','<<ne<<','<<a.coulomb_langmuir<<','<<a.ep<<','<<a.ee<<'\n';}
 std::ofstream ranges;csv(ranges,"ranges.csv","case,gamma0,tau_stop,gamma_stop,unresolved_kinetic_fraction");
 auto track=[&](std::string name,double g0,AtmosphereFunction gas,LossMode mode){
  auto t=stopping_track(g0,gas,mode,step,top,1.1,.02);
  std::ofstream f;csv(f,name+".csv","tau,gamma,phi,ne,T,CL,ep,ee");
  for(std::size_t i=0;i<t.tau.size();++i){auto a=gas(t.tau[i]);auto l=loss_components(t.gamma[i],a.ne);f<<t.tau[i]<<','<<t.gamma[i]<<','<<t.phi[i]<<','<<a.ne<<','<<a.temperature<<','<<l.coulomb_langmuir<<','<<l.ep<<','<<l.ee<<'\n';}
  ranges<<name<<','<<g0<<','<<t.tau.back()<<','<<t.gamma.back()<<','<<.1/(g0-1)<<'\n';
 };
 for(int b=0;b<2;++b)for(int g:{2,10,50,500})track("g"+std::to_string(g)+"_b"+std::to_string(b),g,grey,b?LossMode::FullBremsstrahlung:LossMode::CoulombLangmuir);
 if(argc==5){
  std::ifstream f(argv[4]);std::string line;std::getline(f,line);std::vector<double> lm,ln;
  while(std::getline(f,line)){std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);double m,ne;if(!(row>>m>>ne)||m<=0||ne<=0)throw std::runtime_error("Invalid mass,ne profile");lm.push_back(std::log(m));ln.push_back(std::log(ne));}
  if(lm.size()<3)throw std::runtime_error("Missing density profile");
  if(!std::is_sorted(lm.begin(),lm.end())||std::adjacent_find(lm.begin(),lm.end())!=lm.end())throw std::runtime_error("Profile mass must strictly increase");
  AtmosphereFunction gas=[&](double tau){double x=std::log(tau/cgs::kappa_e);int j=std::upper_bound(lm.begin(),lm.end(),x)-lm.begin();j=std::max(1,std::min(int(lm.size())-1,j));double ne=std::exp(ln[j-1]+(ln[j]-ln[j-1])*(x-lm[j-1])/(lm[j]-lm[j-1]));return GasState{NAN,NAN,ne,ne*cgs::mp};};
  track("g500_heated_reference_density",500,gas,LossMode::FullBremsstrahlung);
 }
 std::cout<<"Figure 1 trajectories computed from initial gamma; step="<<step<<" top="<<top<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}}
