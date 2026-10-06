#pragma once
#include <functional>
#include <vector>
namespace salmi {
struct GasState { double temperature, pressure, ne, rho; };
// Eddington grey reference, T^4=3 Teff^4 (tau+2/3)/4; ignores ram/radiation force.
struct GreyAtmosphere {
    double teff=4.64e6, gravity=1e14;
    GasState operator()(double tau) const;
};
using AtmosphereFunction=std::function<GasState(double)>;
struct FluxBudget {
    double heated, intrinsic;
    static FluxBudget from_temperatures(double th, double ti);
    double total() const {return heated+intrinsic;}
    double number_flux(double mean_gamma) const;
};
double eddington_flux(double gravity);
double goldreich_julian(double field, double period);
double incoming_density(double th, double mean_gamma);
double redshift_from_gravity(double gravity, double mass_solar=1.4);
// dm positive inward. Cell accelerations -> pressures on cell faces.
std::vector<double> hydrostatic(const std::vector<double>& m,
    const std::vector<double>& ram,const std::vector<double>& rad,double gravity,double p_top);
}
