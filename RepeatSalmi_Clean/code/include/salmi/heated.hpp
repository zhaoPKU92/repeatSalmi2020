#pragma once
#include "salmi/redistribution.hpp"
#include "salmi/radiation.hpp"
#include <map>
#include <memory>
#include <string>
namespace salmi {
struct ScatteringPairTemplate {
    std::vector<std::size_t> start;
    std::vector<int> column;
    std::vector<double> low,high;
    std::size_t used=0;
    std::size_t bytes() const {return start.size()*sizeof(std::size_t)+column.size()*sizeof(int)+(low.size()+high.size())*sizeof(double);}
};
struct SparseScattering {
    PhotonGrid grid;
    double temperature;
    std::vector<std::size_t> start;
    std::vector<int> column;
    std::vector<double> gain_coefficient,loss_coefficient;
    // Optional memory-saving representation: one shared, symmetric kernel
    // value per pair. Detailed balance factors remain double precision.
    std::vector<float> symmetric_kernel;
    std::vector<double> boltzmann_pair;
    // Storage 3 retains one bracketing template across all depths in the
    // temperature interval. Actual-temperature coefficients are reconstructed
    // with the identical cutoff and float rounding used by storage 2.
    // Storage 4 uses the same shared layout and stores those rounded values
    // once per depth in symmetric_kernel; no rate or cutoff is changed.
    std::shared_ptr<const ScatteringPairTemplate> shared_template;
    double temperature_fraction=0;
    bool paired=false;
    double symmetric_value(int row,std::size_t entry) const;
    double gain(int i,std::size_t z) const {
        if(symmetric_kernel.empty() && !shared_template)return gain_coefficient[z];
        int A=int(grid.mu.size()),E=int(grid.x.size()),j=shared_template?shared_template->column[z]:column[z],e=i/A,f=j/A;
        double base=symmetric_value(i,z)*(grid.x[f]/grid.x[e])*grid.dx[f]*grid.dmu[j%A];
        return e>f?base*boltzmann_pair[e*E+f]:base;
    }
    double loss(int row,std::size_t entry) const;
    std::vector<double> linearized_gain(const double* departure,const double* planck,
                                       const double* target_mask=nullptr) const;
    CollisionTerms terms(const double* occupation) const;
};
// Temperature interpolation of downward rates, followed by detailed balance
// at the actual local temperature. Cache metadata checks the logarithmic grid
// dimensions, endpoints, temperature spacing and quadrature orders.
class KernelLibrary {
    PhotonGrid grid_;
    std::string directory_;
    double log_step_;
    std::map<int,std::shared_ptr<Redistribution>> tables_;
    std::map<int,std::shared_ptr<ScatteringPairTemplate>> pairs_;
    std::size_t pair_clock_=0;
    const Redistribution& table(int index);
public:
    KernelLibrary(const PhotonGrid& grid,const std::string& cache,double log_step=.25);
    SparseScattering at(double temperature,int compact=0);
    // Returned matrices own all their coefficients. Dense assembly tables
    // can be released while those matrices continue to perform transport.
    std::size_t dense_cache_bytes() const {
        std::size_t bytes=0;
        for(const auto& table:tables_)bytes+=(table.second->probability.size()+table.second->rate.size())*sizeof(double);
        return bytes;
    }
    void release_dense_tables() {tables_.clear();}
};
struct HeatedState {
    double beam_gamma_min=10,beam_gamma_max=2e5,beam_mono_gamma=0;
    bool beam_bremsstrahlung=false;
    std::shared_ptr<const GauntTable> gaunt_table;
    double absorption(double energy_kev,int depth) const;
    bool parabolic=true;
    bool conservative_cells=false;
    // Coherent Thomson benchmark: exact elastic opacity and no Compton T derivative.
    bool coherent_scattering=false;
    double transfer_relaxation=.6;
    bool anderson=false;
    std::vector<double> mass,temperature,pressure,rho,heating,ram,target_flux;
    std::vector<double> occupation,flux,cooling,radiative_acceleration,rosseland,local_derivative;
    std::vector<double> temperature_response; // Optional signed d(cooling)/dlnT for the temperature preconditioner.
    std::vector<double> face_flux,cell_width;
    std::vector<double> target_face_flux;
    std::vector<double> last_v_limiter;
};
struct TransferDiagnostics {double residual;int iterations;int effective_opacity_fallbacks;};
// Interpolate a saved computed radiation field onto a new grid as an initial
// iterate. Final transfer and equilibrium tolerances still apply unchanged.
void initialize_radiation_checkpoint(HeatedState& state,const PhotonGrid& grid,
                                    const std::string& directory);
// Vacuum surface, prescribed bottom flux via diffusion boundary. The source
// field is solved for a fixed T,rho profile; temperature correction is separate.
TransferDiagnostics atmosphere_transfer(HeatedState& state,const PhotonGrid& grid,
    const std::vector<SparseScattering>& scattering,double bottom_flux,
    int max_iterations=500,double tolerance=1e-6);
TransferDiagnostics conservative_transfer(HeatedState& state,const PhotonGrid& grid,
    const std::vector<SparseScattering>& scattering,double bottom_flux,
    int max_iterations=500,double tolerance=1e-6);
void update_hydrostatic(HeatedState& state,double gravity,double relaxation=1);
void update_beam_heating(HeatedState& state,double th,double ti,double delta,
                        int beam_panels=64,double relative_tau_step=0.08);
double correct_temperature(HeatedState& state,double relaxation=.3,double max_log_step=.08,
                           double local_weight=1,double max_flux_log_step=0);
std::vector<double> surface_intensity(const HeatedState& state,const PhotonGrid& grid,
    const std::vector<double>& energies_kev,const std::vector<double>& mus,
    bool include_scattering=true,int azimuth_order=12,int energy_order=24);
}
