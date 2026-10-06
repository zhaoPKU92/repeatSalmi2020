#pragma once
#include "salmi/heated.hpp"
namespace salmi {
// A versioned restart snapshot of one completed radiation solve. This does
// not certify temperature equilibrium; parent invocation defines the model.
void write_iteration_snapshot(const HeatedState&,const PhotonGrid&,
    const std::string& output_directory,int iteration,double energy_min_keV,
    double energy_max_keV,const std::string& gaunt_path,double radiation_residual);
}
