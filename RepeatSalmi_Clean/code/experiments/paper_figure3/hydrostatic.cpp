#include "salmi/heated.hpp"
#include "salmi/constants.hpp"
#include "salmi/stopping.hpp"
#include "salmi/beam.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
namespace salmi {
void update_hydrostatic(HeatedState& s,double gravity,double relaxation) {
    int D=int(s.mass.size());if(D<2 || !(gravity>0))throw std::invalid_argument("Invalid hydrostatic model");
    if(s.pressure.size()!=s.mass.size())s.pressure.assign(D,0);
    double surface_force=gravity;
    if(s.ram.size()==s.mass.size())surface_force+=s.ram.front();
    if(s.radiative_acceleration.size()==s.mass.size())surface_force-=s.radiative_acceleration.front();
    s.rho.resize(D);double previous=surface_force*s.mass[0];
    for(int d=0;d<D;++d) {
        if(d) {
            double force=gravity;
            if(s.ram.size()==s.mass.size())force+=.5*(s.ram[d]+s.ram[d-1]);
            if(s.radiative_acceleration.size()==s.mass.size())force-=.5*(s.radiative_acceleration[d]+s.radiative_acceleration[d-1]);
            if(!(force>0))throw std::runtime_error("Hydrostatic force becomes nonpositive");
            previous+=force*(s.mass[d]-s.mass[d-1]);
        }
        s.pressure[d]=s.pressure[d]>0?(1-relaxation)*s.pressure[d]+relaxation*previous:previous;
        s.rho[d]=s.pressure[d]*cgs::mp/(2*cgs::kb*s.temperature[d]);
    }
}
}
