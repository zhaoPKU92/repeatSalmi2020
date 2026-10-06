#pragma once
#include <vector>
#include <string>
namespace salmi {
// Optional published hydrogen Gaunt tables (van Hoof 2014/2015). The
// nonrelativistic table interpolates g; the merged table interpolates log(g),
// following the author's respective reference interpolation programs.
class GauntTable {
    int nx_=0,ny_=0;
    double xmin_=0,ymin_=0,step_=0;
    bool logarithmic_=false;
    std::vector<double> values_;
public:
    explicit GauntTable(const std::string& path);
    double at(double gamma_squared,double u) const;
    double hydrogen(double energy_kev,double temperature) const;
};
double planck_kev(double energy_kev,double temperature);
// Specific intensity per dimensionless x, with occupation n=1/(exp(x/theta)-1).
double occupation_intensity_x(double x,double occupation);
// Kramers free-free mass opacity, stimulated-emission correction included; gff explicit.
double free_free_opacity(double energy_kev,double temperature,double rho,double gaunt=1);
struct SCWeights { double transmission, upstream, local; };
SCWeights short_characteristic_weights(double optical_path);
// Integral of a quadratic source through local/upstream/third point. The third
// signed optical coordinate may be downstream (<0) or farther upstream (>d).
// A negative interpolated source falls back to positive linear interpolation.
double parabolic_source_integral(double d,double third_coordinate,
    double local_source,double upstream_source,double third_source);
struct FormalSolution {std::vector<double> outward,inward;};
// tau increases inward, mu > 0, source sampled on nodes, exact linear SC.
FormalSolution formal_solve(const std::vector<double>& tau,const std::vector<double>& source,
                           double mu,double bottom_intensity,double top_intensity=0);
struct CoherentSolution {
    std::vector<double> source,mean_intensity,flux,change_history;
    double residual;int iterations;
};
// Isotropic coherent-scattering benchmark, not Salmi's Compton atmosphere.
CoherentSolution coherent_slab(const std::vector<double>& tau,const std::vector<double>& planck,
    const std::vector<double>& destruction,double bottom_intensity,int nmu=8,
    double tolerance=1e-9,int max_iterations=20000);
}
