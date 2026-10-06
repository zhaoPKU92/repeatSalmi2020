#include "salmi/radiation.hpp"
#include "salmi/constants.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
namespace salmi {
GauntTable::GauntTable(const std::string& path) {
    std::ifstream in(path);if(!in)throw std::runtime_error("Cannot read Gaunt table: "+path);
    std::stringstream tokens;std::string line;
    while(std::getline(in,line))tokens<<line.substr(0,line.find('#'))<<'\n';
    int magic=0;tokens>>magic>>nx_>>ny_>>xmin_>>ymin_>>step_;
    if((magic!=20140210 && magic!=20140510)||nx_<4||ny_<4||nx_>1000||ny_>1000||!(step_>0))
        throw std::runtime_error("Invalid Gaunt table header");
    logarithmic_=magic==20140510;
    if(logarithmic_) {int z=0;tokens>>z;if(z!=1)throw std::runtime_error("Only hydrogen Gaunt tables are supported");}
    values_.resize(nx_*ny_);
    for(double& v:values_) {
        if(!(tokens>>v)||!(v>0)||!std::isfinite(v))throw std::runtime_error("Invalid or truncated Gaunt table");
        if(logarithmic_)v=std::log(v);
    }
}
double GauntTable::at(double gam2,double u) const {
    if(!(gam2>0 && u>0))throw std::invalid_argument("Invalid Gaunt coordinates");
    double x=(std::log10(gam2)-xmin_)/step_,y=(std::log10(u)-ymin_)/step_;
    if(!std::isfinite(x+y)||x < -1e-12||y < -1e-12||x>nx_-1+1e-12||y>ny_-1+1e-12)
        throw std::out_of_range("Gaunt table extrapolation is not permitted");
    int ix=std::max(0,std::min(nx_-4,int(std::floor(x))-1));
    int iy=std::max(0,std::min(ny_-4,int(std::floor(y))-1));
    double result=0;
    for(int j=0;j<4;++j)for(int i=0;i<4;++i) {
        double w=1;
        for(int k=0;k<4;++k) {if(k!=i)w*=(x-ix-k)/double(i-k);if(k!=j)w*=(y-iy-k)/double(j-k);}
        result+=w*values_[(iy+j)*nx_+ix+i];
    }
    result=logarithmic_?std::exp(result):result;
    if(!(result>0)||!std::isfinite(result))throw std::runtime_error("Invalid interpolated Gaunt factor");
    return result;
}
double GauntTable::hydrogen(double E,double T) const {
    // Rydberg energy in eV, converted using the same cgs eV as the opacity.
    return at(13.605693122994*.001*cgs::kev/(cgs::kb*T),E*cgs::kev/(cgs::kb*T));
}
}
