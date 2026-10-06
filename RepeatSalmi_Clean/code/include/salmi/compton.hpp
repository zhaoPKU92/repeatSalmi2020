#pragma once
#include <vector>
namespace salmi {
double mono_compton_kernel_analytic(double x,double y,double eta,double gamma,bool extended_precision=false);
// All energies x are h nu / (me c^2); theta=kT/(me c^2).
double klein_nishina_total(double x);
double scaled_bessel_k2(double theta,int order=64);
double minimum_electron_gamma(double x_out,double x_in,double eta);
// Exact monoenergetic kernel, evaluated from the positive azimuth integral SPW12 A.19.
double mono_compton_kernel(double x_out,double x_in,double eta,double gamma,int nphi=24);
double thermal_compton_kernel(double x_out,double x_in,double eta,double theta,
                             int ngamma=32,int nphi=24);
// Independent electron-angle average of relative rate * total Klein-Nishina cross section.
double thermal_scattering_rate(double x,double theta,int order=48);
// Reuses quadrature nodes and normalization when assembling many kernel values.
class ThermalKernel {
    double theta_,prefactor_;
    std::vector<double> gamma_u_,gamma_w_,phi_cos_;
public:
    explicit ThermalKernel(double theta,int ngamma=32,int nphi=24);
    double operator()(double x_out,double x_in,double eta) const;
};
}
