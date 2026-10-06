#include "coupled_probe.hpp"
#include "salmi/constants.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <fstream>
#include <iostream>
#include <stdexcept>
namespace salmi {
HeatedState coupled_temperature_probe(const HeatedState& initial,const PhotonGrid& g,KernelLibrary& library,const CoupledProbeOptions& o) {
 struct Restore {PaperResponse saved;Restore():saved(paper_response){}~Restore(){paper_response=std::move(saved);}} restore;
 paper_response.outer_iteration=-1; // Keep STOP_REQUEST; suppress construction radiation seeds.
 HeatedState s=initial;const int D=int(s.mass.size());
 std::vector<SparseScattering> kernels;kernels.reserve(D);
 for(double T:s.temperature){kernels.push_back(library.at(T,o.coefficient_cache_mb?3:o.compact_kernel));if(o.compact_kernel>=2)library.release_dense_tables();}
 if(o.coefficient_cache_mb){std::vector<int> order(D);std::iota(order.begin(),order.end(),0);std::stable_sort(order.begin(),order.end(),[&](int a,int b){return s.temperature[a]<s.temperature[b];});std::size_t used=0,budget=std::size_t(o.coefficient_cache_mb)*1024*1024;
  for(int d:order){auto& k=kernels[d];std::size_t bytes=k.start.back()*sizeof(float);if(bytes>budget-used)continue;std::vector<float> v(k.start.back());for(int i=0;i<g.size();++i)for(std::size_t z=k.start[i];z<k.start[i+1];++z)v[z]=float(k.symmetric_value(i,z));k.symmetric_kernel=std::move(v);used+=bytes;}}
 library.release_dense_tables();
 for(int pass=0;pass<12;++pass){
  double change=1;
  for(int h=0;h<30 && change>1e-10;++h){update_hydrostatic(s,o.gravity,1.);update_beam_heating(s,o.heated_temperature,o.intrinsic_temperature,o.delta,o.beam_panels,o.beam_tau_step);auto before=s.pressure;update_hydrostatic(s,o.gravity,1.);change=0;for(int d=0;d<D;++d)change=std::max(change,std::abs(s.pressure[d]/before[d]-1));}
  if(change>1e-8)throw std::runtime_error("Probe beam/pressure iteration failed");
  auto tr=atmosphere_transfer(s,g,kernels,cgs::sigma_sb*std::pow(o.intrinsic_temperature,4),o.max_radiation_iterations,o.radiation_tolerance);
  if(tr.residual>o.radiation_tolerance)throw std::runtime_error("Probe radiation not converged; refusing noisy Jacobian column");
  double P=(o.gravity+s.ram[0]-s.radiative_acceleration[0])*s.mass[0],error=0;
  for(int d=0;d<D;++d){if(d)P+=(s.mass[d]-s.mass[d-1])*(o.gravity+.5*(s.ram[d]+s.ram[d-1]-s.radiative_acceleration[d]-s.radiative_acceleration[d-1]));error=std::max(error,std::abs(s.pressure[d]/P-1));}
  std::cout<<"FD_PROBE pass="<<pass<<" RT="<<tr.residual<<" HSE="<<error<<std::endl;
  if(error<=1e-7)return s;
 }
 throw std::runtime_error("Probe coupled hydrostatic/radiation iteration failed");
}
}
