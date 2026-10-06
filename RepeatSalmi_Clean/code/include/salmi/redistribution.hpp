#pragma once
#include "salmi/numerics.hpp"
#include <vector>
namespace salmi {
// 0: historical midpoint operator; 1: energy-transfer fitted bin quadrature.
// Selected explicitly with SALMI_COMPTON_BIN_INTEGRATION=moment_v1. The
// integrator remains opt-in: its atmosphere validation covers the reported
// nonheated cases, rather than all externally heated models.
int redistribution_integration_version();
// Discrete occupation numbers n=CI_x/x^3. x is in electron rest-mass units.
// State order: energy index * number of directions + direction index.
struct PhotonGrid {
    std::vector<double> x,dx,mu,dmu,weight;
    int size() const {return int(x.size()*mu.size());}
    int state(int e,int a) const {return e*int(mu.size())+a;}
    static PhotonGrid logarithmic(double xmin,double xmax,int ne,int nmu_half);
};
struct CollisionTerms {
    std::vector<double> gain,loss,net;
};
struct Redistribution {
    PhotonGrid grid;
    double theta;
    // P[target,source]: discrete spontaneous probability per Thomson depth.
    // The midpoint and energy-transfer-fitted quadratures give different
    // approximations to the continuum operator; see docs/energy_bin_integration.md.
    std::vector<double> probability,rate;
    static Redistribution assemble(const PhotonGrid& grid,double theta,
        int nazimuth=16,int ngamma=32,int electron_azimuth=24);
    static Redistribution assemble_bins(const PhotonGrid& grid,double theta,
        int nazimuth=16,int ngamma=32,int electron_azimuth=24,
        int energy_order=8,int rate_order=64);
    CollisionTerms collision(const std::vector<double>& occupation,bool induced=true) const;
    std::vector<double> bose_einstein(double chemical_potential=0) const;
};
struct ComptonSlab {
    std::vector<double> tau,occupation,flux,photon_flux,history;
    int iterations;
    double residual;
};
// Prescribed-temperature, purely scattering slab. Boundary data are occupation
// numbers for each energy/direction; only incoming directions are used.
// This is the transfer building block, not a self-consistent atmosphere model.
ComptonSlab solve_compton_slab(const Redistribution& kernel,
    const std::vector<double>& tau,const std::vector<double>& bottom,
    const std::vector<double>& top,bool induced=true,double tolerance=1e-7,
    int max_iterations=2000);
}
