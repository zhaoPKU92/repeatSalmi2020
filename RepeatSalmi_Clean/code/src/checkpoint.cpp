#include "salmi/heated.hpp"
#include "salmi/constants.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <regex>
#include <sstream>
#include <stdexcept>
namespace salmi {
namespace {
struct Bracket {int lo,hi;double fraction;};
Bracket bracket(const std::vector<double>& coordinates,double x) {
    if(x<=coordinates.front())return {0,0,0};
    int end=int(coordinates.size())-1;if(x>=coordinates.back())return {end,end,0};
    int hi=int(std::lower_bound(coordinates.begin(),coordinates.end(),x)-coordinates.begin()),lo=hi-1;
    return {lo,hi,(x-coordinates[lo])/(coordinates[hi]-coordinates[lo])};
}
double blend(double a,double b,double f) {
    if(f<=0)return a;if(f>=1)return b;
    return a>0 && b>0?std::exp((1-f)*std::log(a)+f*std::log(b)):(1-f)*a+f*b;
}
}
void initialize_radiation_checkpoint(HeatedState& s,const PhotonGrid& grid,const std::string& directory) {
    std::ifstream manifest(directory+"/status.json");
    if(!manifest)throw std::runtime_error("Cannot read radiation initialization manifest");
    std::string metadata((std::istreambuf_iterator<char>(manifest)),std::istreambuf_iterator<char>());
    if(metadata.find("\"iteration_checkpoint\":true")!=std::string::npos) {
        std::ifstream committed(directory+"/complete.json");
        if(!committed)throw std::runtime_error("Iteration checkpoint was not committed");
    }
    auto parameter=[&](const std::string& name) {
        std::regex pattern("\""+name+"\"\\s*:\\s*([-+0-9.eE]+)");std::smatch match;
        if(!std::regex_search(metadata,match,pattern))throw std::runtime_error("Missing radiation grid metadata: "+name);
        return std::stod(match[1]);
    };
    std::ifstream field(directory+"/occupation.bin",std::ios::binary);int dim[3]={};
    if(!field.read(reinterpret_cast<char*>(dim),sizeof(dim)))throw std::runtime_error("Cannot read radiation initialization data");
    if(dim[0]!=int(parameter("depth_nodes")) || dim[1]!=int(parameter("energy_bins")) || dim[2]!=2*int(parameter("directions_per_hemisphere")))throw std::runtime_error("Radiation initialization dimensions do not match manifest");
    auto old=PhotonGrid::logarithmic(parameter("energy_min_keV")*cgs::kev/cgs::mec2,parameter("energy_max_keV")*cgs::kev/cgs::mec2,dim[1],dim[2]/2);
    std::vector<double> n(static_cast<std::size_t>(dim[0])*old.size());
    std::vector<double> mass,energy;
    if(!field.read(reinterpret_cast<char*>(n.data()),n.size()*sizeof(double)))throw std::runtime_error("Truncated radiation initialization");
    for(double value:n)if(!(value>=0) || !std::isfinite(value))throw std::runtime_error("Invalid radiation initialization value");
    std::ifstream state(directory+"/structure.csv");std::string line;std::getline(state,line);
    while(std::getline(state,line)) {
        std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);double m;
        if(!(row>>m) || !(m>0))throw std::runtime_error("Invalid radiation initialization mass");
        mass.push_back(std::log(m));
    }
    if(int(mass.size())!=dim[0])throw std::runtime_error("Radiation initialization depth mismatch");
    for(double x:old.x)energy.push_back(std::log(x));
    s.occupation.resize(s.mass.size()*grid.size());
    for(int d=0;d<int(s.mass.size());++d)for(int e=0;e<int(grid.x.size());++e)for(int a=0;a<int(grid.mu.size());++a) {
        auto m=bracket(mass,std::log(s.mass[d])),en=bracket(energy,std::log(grid.x[e])),mu=bracket(old.mu,grid.mu[a]);
        auto angular=[&](int dd,int ee) {
            double l=n[std::size_t(dd)*old.size()+old.state(ee,mu.lo)],r=n[std::size_t(dd)*old.size()+old.state(ee,mu.hi)];
            return (1-mu.fraction)*l+mu.fraction*r;
        };
        double lo=blend(angular(m.lo,en.lo),angular(m.lo,en.hi),en.fraction);
        double hi=blend(angular(m.hi,en.lo),angular(m.hi,en.hi),en.fraction);
        double value=blend(lo,hi,m.fraction);
        // Newly added energy bins need a decaying seed, not a constant
        // occupation plateau. This only initializes the iteration; all bins
        // participate in the subsequent full physical transfer solve.
        if(grid.x[e]>old.x.back())value*=std::exp(-(grid.x[e]-old.x.back())/(cgs::kb*s.temperature[d]/cgs::mec2));
        else if(grid.x[e]<old.x.front())value*=old.x.front()/grid.x[e];
        s.occupation[std::size_t(d)*grid.size()+grid.state(e,a)]=d==0 && grid.mu[a]<0?0:value;
    }
}
}
