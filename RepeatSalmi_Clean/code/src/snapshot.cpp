#include "salmi/snapshot.hpp"
#include <cerrno>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif
namespace salmi {
namespace {
void folder(const std::string& path,bool allow_existing) {
#ifdef _WIN32
 int r=_mkdir(path.c_str());
#else
 int r=mkdir(path.c_str(),0777);
#endif
 if(r && !(allow_existing && errno==EEXIST))throw std::runtime_error("Cannot create fresh snapshot directory: "+path);
}
void closed(std::ofstream& f) {f.flush();if(!f)throw std::runtime_error("Snapshot write failed");f.close();if(!f)throw std::runtime_error("Snapshot close failed");}
}
void write_iteration_snapshot(const HeatedState& s,const PhotonGrid& g,const std::string& out,
 int iteration,double emin,double emax,const std::string& gaunt,double residual) {
 int D=int(s.mass.size()),E=int(g.x.size()),A=int(g.mu.size());
 if(D<3 || iteration<0 || !(emin>0 && emax>emin && residual>=0 && std::isfinite(residual)) || s.occupation.size()!=std::size_t(D)*g.size())throw std::invalid_argument("Invalid snapshot dimensions or metadata");
 for(const auto* v:{&s.mass,&s.temperature,&s.pressure,&s.rho,&s.heating,&s.ram,&s.target_flux,&s.flux,&s.cooling,&s.radiative_acceleration,&s.rosseland,&s.local_derivative}) {
  if(v->size()!=std::size_t(D))throw std::invalid_argument("Inconsistent snapshot columns");
  for(double x:*v)if(!std::isfinite(x))throw std::invalid_argument("Nonfinite snapshot column");
 }
 for(double x:s.occupation)if(!(x>=0 && std::isfinite(x)))throw std::invalid_argument("Invalid snapshot radiation");
 std::string parent=out+"/checkpoints",dir=parent+"/iteration_"+std::to_string(iteration);folder(parent,true);folder(dir,false);
 std::ofstream st(dir+"/structure.csv");st<<std::setprecision(17)<<"mass,T,P,rho,Qplus,gram,Ftarget,Fcomputed,cooling,grad,rosseland,local_derivative\n";
 for(int d=0;d<D;++d)st<<s.mass[d]<<','<<s.temperature[d]<<','<<s.pressure[d]<<','<<s.rho[d]<<','<<s.heating[d]<<','<<s.ram[d]<<','<<s.target_flux[d]<<','<<s.flux[d]<<','<<s.cooling[d]<<','<<s.radiative_acceleration[d]<<','<<s.rosseland[d]<<','<<s.local_derivative[d]<<'\n';
 closed(st);
 std::ofstream field(dir+"/occupation.bin",std::ios::binary);int dimensions[]={D,E,A};
 field.write(reinterpret_cast<const char*>(dimensions),sizeof(dimensions));field.write(reinterpret_cast<const char*>(s.occupation.data()),s.occupation.size()*sizeof(double));closed(field);
 std::ofstream meta(dir+"/status.json");meta<<std::setprecision(17)<<"{\"iteration_checkpoint\":true,\"discrete_converged\":false,\"iteration\":"<<iteration<<",\"depth_nodes\":"<<D<<",\"energy_bins\":"<<E<<",\"directions_per_hemisphere\":"<<A/2<<",\"energy_min_keV\":"<<emin<<",\"energy_max_keV\":"<<emax<<",\"gaunt_table\":"<<std::quoted(gaunt)<<",\"scattering\":\""<<(s.coherent_scattering?"coherent_thomson":"Compton")<<"\",\"transfer_residual\":"<<residual<<",\"scope\":\"Completed radiation iterate only; use the parent invocation for physical parameters. Not certified atmosphere equilibrium.\"}\n";closed(meta);
 // Commit marker is written last, after all state files have been closed.
 std::ofstream marker(dir+"/complete.json");marker<<"{\"state_files_closed\":true,\"iteration\":"<<iteration<<"}\n";closed(marker);
}
}
