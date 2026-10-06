#pragma once
#include "salmi/atmosphere.hpp"
#include "salmi/numerics.hpp"
#include <string>
namespace salmi {
enum class LossMode { CoulombLangmuir, CoulombLangmuirEP, FullBremsstrahlung };
struct LossComponents {double coulomb_langmuir,ep,ee;};
double phi_coulomb_langmuir(double gamma,double ne); // Salmi (13)
double phi_electron_proton(double gamma);           // Salmi (11)
double phi_electron_electron(double gamma);         // Haug (2004), gamma >= 1.1
LossComponents loss_components(double gamma,double ne);
double loss_rate(double gamma,double ne,LossMode mode);
struct Track {
    std::vector<double> tau,gamma,phi;
    double gamma0,cutoff;
    // Kinetic energy in units mec^2. Residual at cutoff is thermalized locally.
    double kinetic(double tau_query,bool thermalize=true) const;
    double momentum(double tau_query,bool thermalize=true) const;
};
// Integrate tau(s), s=log[(gamma0-1)/(gamma-1)]. RK4 with relative tau step cap.
Track stopping_track(double gamma0,const AtmosphereFunction& gas,LossMode mode,
    double ds=0.01,double tau_top=1e-10,double cutoff=1.1,double relative_tau_step=0.08);
struct Deposition {
    std::vector<double> tau_edges,average_phi,qplus,gram,flux;
    double deposited_fraction,remaining_fraction,quadrature_number,quadrature_energy;
};
// Conservative finite-volume deposition; each cell gets kinetic energy and momentum loss.
Deposition deposit(const Quadrature& beam,const AtmosphereFunction& gas,LossMode mode,
    const FluxBudget& budget,double mean_gamma,const std::vector<double>& tau_edges,
    double ds=0.01,double cutoff=1.1,double relative_tau_step=0.08);
}
