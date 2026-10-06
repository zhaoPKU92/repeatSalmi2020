#pragma once
#include "salmi/heated.hpp"
namespace salmi {
extern bool paper_pointwise_heating;
extern bool paper_normalize_point_heating;
extern bool paper_hat_heating;

struct PaperResponse {
    std::string checkpoint_root; int outer_iteration=-1;
    bool rt_extrapolate=false; int rt_history=4; double rt_occupation_floor=1e-30;
    double rt_planck_intensity_floor_fraction=0;
    double balanced_energy_scale=1;
    bool compare_clipped_step=false; int fd_columns=0,fd_every=8; double fd_step=.02,fd_tolerance=1e-7;
    bool full_thermal_response=false; double kernel_temperature_step=.1;
    std::vector<double> denominator, spectral_flux, opacity;
    int ali_guards=0; double bottom_thermalization=0;
};
extern PaperResponse paper_response;
double compton_cooling_temperature_response(const SparseScattering&,const double*,double);
double paper_temperature(HeatedState&,double,double,double,double);
}
