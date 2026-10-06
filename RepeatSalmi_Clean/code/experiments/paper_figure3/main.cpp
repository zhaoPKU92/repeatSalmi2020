#include "paper.hpp"
#include "balanced_temperature.hpp"
#include "coupled_probe.hpp"
#include "salmi/anderson.hpp"
#include "salmi/atmosphere.hpp"
#include "salmi/constants.hpp"
#include "salmi/radiation.hpp"
#include "salmi/snapshot.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <set>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace salmi;
int main(int argc,char** argv) {
 try {
    std::map<std::string,std::string> args;
    const std::set<std::string> allowed={"--rt-planck-intensity-floor","--output","--cache","--delta","--energy","--depth","--angles",
        "--iterations","--inner","--initial-profile","--rays","--kernel-step","--relaxation","--formal-order","--solver","--transfer-relax",
        "--energy-min","--energy-max","--mass-min","--mass-max","--beam-panels","--require-equilibrium","--compact-kernel","--threads","--anderson","--initial-field","--transfer-tolerance","--temperature-step","--gaunt-table",
        "--heated-temperature","--intrinsic-temperature","--logg","--gamma-min","--gamma-max","--mono-gamma","--beam-brems","--startup-iterations","--heating-ramp-iterations","--flux-step","--schedule-offset","--temperature-tolerance","--mass-grid","--beam-relative-tau-step","--checkpoint-every","--heating-remap","--rt-extrapolate","--rt-history","--rt-occupation-floor","--temperature-method","--restore-hydrostatic","--balanced-hot-response","--balanced-staged-radiation","--coefficient-cache-mb","--balanced-energy-scale","--balanced-full-thermal-response","--balanced-step-control","--balanced-fd-columns","--balanced-fd-every","--balanced-fd-step","--balanced-fd-tolerance"};
    for(int i=1;i<argc;i+=2) {
        if(std::string(argv[i])=="--help") {
            std::cout<<"Modular heated hydrogen atmosphere (CGS, K, keV). Create output/cache directories first.\n"
                <<"--output DIR --cache DIR --delta 0|1|2|3 --energy N --depth N --angles N\n"
                <<"--iterations N --inner N --initial-profile CSV --initial-field MODEL_DIR --rays 0|1\n"
                <<"--kernel-step FLOAT --relaxation FLOAT --temperature-step FLOAT\n"
                <<"--temperature-method paper|balanced|external --balanced-hot-response 0|1 --balanced-staged-radiation 0|1 --balanced-energy-scale FLOAT --balanced-full-thermal-response 0|1\n"
                <<"--balanced-step-control scaled|clip_compare --balanced-fd-columns 0..8 --balanced-fd-every N --balanced-fd-step .02 --balanced-fd-tolerance 1e-7\n"
                <<"balanced: native C++ joint energy/flux solve; temperature-step is the log-T trust radius.\n"
                <<"--restore-hydrostatic 0|1 imports same-grid computed pressure and radiative acceleration.\n"
                <<"--solver feautrier|short_characteristic --formal-order linear|parabolic --transfer-relax FLOAT\n"
                <<"--transfer-tolerance FLOAT --anderson 0|1 --compact-kernel 0|1|2|3|4 --threads N\n"
                <<"--coefficient-cache-mb N (compact4 only; 0 unlimited, positive bounds cached coefficients)\n"
                <<"--temperature-tolerance FLOAT (default .001; may only be tightened)\n"
                <<"--checkpoint-every N (default 1; 0 disables completed-iteration restart snapshots)\n"
                <<"--energy-min KEV --energy-max KEV --mass-min CGS --mass-max CGS\n"
                <<"--mass-grid CSV optional single mass column, strictly increasing;\n"
                <<"its count and endpoints must match --depth, --mass-min and --mass-max.\n"
                <<"--beam-panels N --gaunt-table FILE --require-equilibrium 0|1\n"
                <<"--beam-relative-tau-step FLOAT (1e-4 .. .08; default .08)\n"
                <<"--heated-temperature K --intrinsic-temperature K --logg FLOAT\n"
                <<"--gamma-min FLOAT --gamma-max FLOAT --mono-gamma FLOAT --beam-brems 0|1\n"
                <<"--startup-iterations N --heating-ramp-iterations N --flux-step FLOAT\n"
                <<"--schedule-offset N continues that schedule from a saved computed profile;\n"
                <<"the radiation field is reconverged, not assumed to be an exact process snapshot.\n"
                <<"Defaults reproduce Figure 3 physical parameters; delta=0 is the non-heated reference.\n"
                <<"Heated startup/ramp defaults:30/70 iterations. When enabled: Th starts at Ti, flux-only startup,\n"
                <<"then linear Th and local-energy-weight ramps; the final flux-step cap defaults to .01.\n"
                <<"Convergence on one grid is not a certificate of agreement with the paper.\n";
            return 0;
        }
        if(!allowed.count(argv[i]))throw std::invalid_argument("Unknown option: "+std::string(argv[i]));
        if(i+1>=argc)throw std::invalid_argument("Missing option value");
        args[argv[i]]=argv[i+1];
    }
    auto get=[&](const std::string& k,const std::string& def){auto it=args.find(k);return it==args.end()?def:it->second;};
    std::string out=get("--output","output/heated"),cache=get("--cache","output/kernel_cache"),profile=get("--initial-profile","");
    const std::string origin=profile.empty()?"grey":profile.find("figure4_initial")!=std::string::npos?"digitized_paper_temperature":"computed_profile_checkpoint";
    int ne=std::stoi(get("--energy","300")),nd=std::stoi(get("--depth","98")),na=std::stoi(get("--angles","3"));
    int outer=std::stoi(get("--iterations","600")),inner=std::stoi(get("--inner","2000")),delta=std::stoi(get("--delta","1"));
    const bool external_temperature=get("--temperature-method","paper")=="external";
    const bool native_balanced=get("--temperature-method","paper")=="balanced";
    const bool coupled_temperature=external_temperature||native_balanced;
    if(get("--temperature-method","paper")!="paper" && !coupled_temperature)throw std::invalid_argument("Invalid temperature method");
    const bool balanced_staged_radiation=get("--balanced-staged-radiation","0")=="1";
    if(get("--balanced-staged-radiation","0")!="0" && get("--balanced-staged-radiation","0")!="1")throw std::invalid_argument("Invalid staged-radiation switch");
    if(get("--balanced-full-thermal-response","0")!="0" && get("--balanced-full-thermal-response","0")!="1")throw std::invalid_argument("Invalid full thermal-response switch");
    paper_response.full_thermal_response=get("--balanced-full-thermal-response","0")=="1";
    paper_response.balanced_energy_scale=std::stod(get("--balanced-energy-scale","1"));
    if(!(paper_response.balanced_energy_scale>0 && paper_response.balanced_energy_scale<=100))throw std::invalid_argument("Balanced energy scale must be positive and at most100");
    std::string step_control=get("--balanced-step-control","scaled");
    if(step_control!="scaled" && step_control!="clip_compare")throw std::invalid_argument("Invalid balanced step control");
    paper_response.compare_clipped_step=step_control=="clip_compare";
    paper_response.fd_columns=std::stoi(get("--balanced-fd-columns","0"));paper_response.fd_every=std::stoi(get("--balanced-fd-every","8"));
    paper_response.fd_step=std::stod(get("--balanced-fd-step",".02"));paper_response.fd_tolerance=std::stod(get("--balanced-fd-tolerance","1e-7"));
    if(paper_response.fd_columns<0 || paper_response.fd_columns>8 || paper_response.fd_every<1 || !(paper_response.fd_step>=.001 && paper_response.fd_step<=.05) || !(paper_response.fd_tolerance>0 && paper_response.fd_tolerance<=1e-7))throw std::invalid_argument("Invalid numerical Jacobian probe controls");
    if(paper_response.fd_columns && !native_balanced)throw std::invalid_argument("Jacobian probes require native balanced temperature method");
    const bool balanced_hot_response=get("--balanced-hot-response","1")=="1";
    if(get("--balanced-hot-response","1")!="1" && get("--balanced-hot-response","1")!="0")throw std::invalid_argument("Invalid hot-response switch");
    bool rays=std::stoi(get("--rays","0"))!=0;
    bool conservative=false;
    if(get("--solver","short_characteristic")!="short_characteristic")throw std::invalid_argument("Paper experiment requires linear SC");
    if(get("--solver","short_characteristic")!="feautrier" && get("--solver","short_characteristic")!="short_characteristic")throw std::invalid_argument("Unknown transfer solver");
    int beam_panels=std::stoi(get("--beam-panels","64"));
    int checkpoint_every=std::stoi(get("--checkpoint-every","1"));
    if(checkpoint_every<0)throw std::invalid_argument("Checkpoint interval must be nonnegative");
    double beam_relative_tau_step=std::stod(get("--beam-relative-tau-step",".08"));
    if(!(beam_relative_tau_step>=1e-4 && beam_relative_tau_step<=.08))throw std::invalid_argument("Beam relative tau step must be in [1e-4,.08]");
    bool require_equilibrium=std::stoi(get("--require-equilibrium","0"))!=0;
    int compact_kernel=std::stoi(get("--compact-kernel","4"));
    if(compact_kernel<0 || compact_kernel>4)throw std::invalid_argument("Kernel storage must be 0, 1, 2, 3 or 4");
    int coefficient_cache_mb=std::stoi(get("--coefficient-cache-mb","0"));
    if(coefficient_cache_mb<0 || coefficient_cache_mb>8192 || (coefficient_cache_mb && compact_kernel!=4))throw std::invalid_argument("Coefficient cache budget requires compact-kernel4 and must be0..8192 MB");
    int threads=std::stoi(get("--threads","8"));if(threads<1)throw std::invalid_argument("Threads must be positive");
#ifdef _OPENMP
    omp_set_dynamic(0);omp_set_num_threads(threads);
#else
    if(threads!=1)throw std::runtime_error("Rebuild with OpenMP to use --threads > 1");
#endif
    if(delta<0 || delta>3 || nd<10 || outer<0 || inner<1)throw std::invalid_argument("Invalid model parameters");
    std::string beam_remap=get("--heating-remap","normalized_point");
    if(beam_remap!="conservative" && beam_remap!="pointwise" && beam_remap!="normalized_point" && beam_remap!="hat")throw std::invalid_argument("Unknown beam remapping");
    paper_hat_heating=beam_remap=="hat";
    paper_pointwise_heating=beam_remap!="conservative";
    paper_normalize_point_heating=beam_remap=="normalized_point";
    if(get("--rt-extrapolate","0")!="0" && get("--rt-extrapolate","0")!="1")throw std::invalid_argument("Invalid radiation extrapolation switch");
    paper_response.rt_extrapolate=get("--rt-extrapolate","0")=="1";
    paper_response.rt_history=std::stoi(get("--rt-history","4"));
    paper_response.rt_planck_intensity_floor_fraction=std::stod(get("--rt-planck-intensity-floor","0"));
    if(!std::isfinite(paper_response.rt_planck_intensity_floor_fraction) || paper_response.rt_planck_intensity_floor_fraction<0 || paper_response.rt_planck_intensity_floor_fraction>.1)throw std::invalid_argument("Invalid RT intensity floor fraction");
    paper_response.rt_occupation_floor=std::stod(get("--rt-occupation-floor","1e-30"));
    if(!(paper_response.rt_occupation_floor>0) || !std::isfinite(paper_response.rt_occupation_floor))throw std::invalid_argument("Invalid RT occupation floor");
    if(paper_response.rt_history<1 || paper_response.rt_history>20)throw std::invalid_argument("RT history must be 1 to20");
    double step=std::stod(get("--kernel-step","0.1")),relax=std::stod(get("--relaxation","0.5"));
    paper_response.kernel_temperature_step=step;
    if(paper_response.full_thermal_response && compact_kernel<3)throw std::invalid_argument("Full thermal response requires compact kernel3 or4");
    double transfer_tolerance=std::stod(get("--transfer-tolerance","1e-6")),temperature_step=std::stod(get("--temperature-step",".08"));
    double temperature_tolerance=std::stod(get("--temperature-tolerance",".001"));
    if(!(temperature_tolerance>0 && temperature_tolerance<=.001))throw std::invalid_argument("Temperature tolerance must be positive and at most .001");
    if(!(transfer_tolerance>0 && transfer_tolerance<=1e-3 && temperature_step>0 && temperature_step<=.5))throw std::invalid_argument("Invalid iteration tolerance or temperature step");
    if(!(relax>0 && relax<=1))throw std::invalid_argument("Relaxation must be in (0,1]");
    if(get("--formal-order","linear")!="linear" && get("--formal-order","linear")!="parabolic")
        throw std::invalid_argument("Formal order must be linear or parabolic");
    double emin=std::stod(get("--energy-min",".00041356677")),emax=std::stod(get("--energy-max","413.56677"));
    double mmin=std::stod(get("--mass-min","1e-6")),mmax=std::stod(get("--mass-max","1e6"));
    if(beam_panels<4)throw std::invalid_argument("Too few beam panels");
    auto g=PhotonGrid::logarithmic(emin*cgs::kev/cgs::mec2,emax*cgs::kev/cgs::mec2,ne,na);
    KernelLibrary library(g,cache,step);HeatedState s;s.parabolic=false;
    std::string gaunt_path=get("--gaunt-table","");std::replace(gaunt_path.begin(),gaunt_path.end(),'\\','/');
    if(gaunt_path.find('"')!=std::string::npos)throw std::invalid_argument("Invalid Gaunt filename");
    if(!gaunt_path.empty())throw std::invalid_argument("This experiment requires constant Gaunt factor 1");
    s.conservative_cells=false; // Paper pointwise energy balance; preserve microscopic cooling.
    s.anderson=false;
    if(get("--anderson","0")!="0" || get("--formal-order","linear")!="linear")throw std::invalid_argument("Use the paper linear-SC ALI method");
    s.transfer_relaxation=std::stod(get("--transfer-relax","1"));
    if(!(s.transfer_relaxation>0 && s.transfer_relaxation<=1))throw std::invalid_argument("Transfer relaxation must be in (0,1]");
    std::string mass_grid_path=get("--mass-grid","");
    std::replace(mass_grid_path.begin(),mass_grid_path.end(),'\\','/');
    if(mass_grid_path.find('"')!=std::string::npos)throw std::invalid_argument("Invalid mass grid filename");
    if(mass_grid_path.empty())s.mass=log_grid(mmin,mmax,nd);
    else {
        if(!(mmin>0 && mmax>mmin) || !std::isfinite(mmin+mmax))throw std::invalid_argument("Invalid mass bounds");
        std::ifstream grid_file(mass_grid_path);std::string line;
        if(!grid_file || !std::getline(grid_file,line) || line!="mass")throw std::invalid_argument("Mass grid must have a single mass header");
        while(std::getline(grid_file,line)) {
            std::istringstream row(line);double mass;std::string extra;
            if(!(row>>mass) || (row>>extra) || !(mass>0) || !std::isfinite(mass) ||
               (!s.mass.empty() && mass<=s.mass.back()))throw std::invalid_argument("Mass nodes must be finite, positive, strictly increasing single values");
            s.mass.push_back(mass);
        }
        if(int(s.mass.size())!=nd)throw std::invalid_argument("Mass grid count does not match --depth");
        if(std::abs(s.mass.front()/mmin-1)>1e-12 || std::abs(s.mass.back()/mmax-1)>1e-12)
            throw std::invalid_argument("Mass grid endpoints do not match mass bounds");
    }
    s.temperature.resize(nd);
    s.beam_gamma_min=std::stod(get("--gamma-min","10"));s.beam_gamma_max=std::stod(get("--gamma-max","200000"));
    s.beam_mono_gamma=std::stod(get("--mono-gamma",delta?"50":"0"));s.beam_bremsstrahlung=std::stoi(get("--beam-brems","0"))!=0;
    if(!(s.beam_gamma_min>1.1 && s.beam_gamma_max>s.beam_gamma_min && (s.beam_mono_gamma==0 || s.beam_mono_gamma>1.1)))throw std::invalid_argument("Invalid beam distribution");
    double th_final=std::stod(get("--heated-temperature",delta?"5000000":"0"));
    double ti=args.count("--intrinsic-temperature")?std::stod(get("--intrinsic-temperature","")):(delta?th_final*std::pow(.01,.25):5000000.);
    double logg=std::stod(get("--logg","14.3")),gravity=std::pow(10.,logg);
    int startup=std::stoi(get("--startup-iterations",delta?"30":"0")),ramp=std::stoi(get("--heating-ramp-iterations",delta?"70":"0"));
    int schedule_offset=std::stoi(get("--schedule-offset","0"));
    double final_flux_step=std::stod(get("--flux-step",delta?".01":".04"));
    bool continuation=startup>0 || ramp>0;
    if(!(th_final>=0 && ti>0 && gravity>0 && std::isfinite(th_final+ti+gravity) && startup>=0 && ramp>=0 && final_flux_step>0))throw std::invalid_argument("Invalid atmosphere or continuation parameters");
    if(schedule_offset<0 || (schedule_offset>0 && (profile.empty() || !continuation)))
        throw std::invalid_argument("Schedule offset requires a computed initial profile and an enabled continuation schedule");
    if(continuation && (th_final<ti || ramp<1 || double(outer)+schedule_offset<double(startup)+ramp))throw std::invalid_argument("Continuation requires Th>=Ti, a positive ramp and sufficient outer iterations");
    auto progress=[&](int it){return continuation?std::max(0.,std::min(1.,(double(it)+schedule_offset-startup)/ramp)):1.;};
    if(native_balanced && progress(0)!=1)throw std::invalid_argument("Native balanced solver requires a fixed flux budget: use a computed final-budget seed, or disable startup/ramp explicitly");
    double th=continuation?ti+(th_final-ti)*progress(0):th_final;
    auto flux=FluxBudget::from_temperatures(th,ti);
    std::vector<double> pm,pt;
    if(!profile.empty()) {
        std::ifstream in(profile);if(!in)throw std::runtime_error("Cannot read initial T profile");
        std::string line;std::getline(in,line);
        while(std::getline(in,line)){std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);double m,T;if(row>>m>>T){pm.push_back(std::log(m));pt.push_back(std::log(T));}}
        if(pm.size()<3)throw std::runtime_error("Initial T profile has insufficient data");
    }
    for(int d=0;d<nd;++d) {
        double lm=std::log(s.mass[d]);
        if(pm.empty())s.temperature[d]=std::pow((std::pow(th,4)+std::pow(ti,4))*.75*(cgs::kappa_e*s.mass[d]+2./3),.25);
        else if(lm<pm.front())s.temperature[d]=std::exp(pt.front());
        else if(lm>pm.back())s.temperature[d]=std::exp(pt.back()+.25*(lm-pm.back()));
        else s.temperature[d]=std::exp(linear_interp(pm,pt,lm));
    }
    update_hydrostatic(s,gravity);update_beam_heating(s,th,ti,delta,beam_panels,beam_relative_tau_step);update_hydrostatic(s,gravity);
    // Optional same-grid hydrostatic warm state. This is not a full process restart.
    if(get("--restore-hydrostatic","0")=="1") {
        if(profile.empty())throw std::invalid_argument("Hydrostatic restore requires a profile");
        std::ifstream input(profile);std::string line;std::getline(input,line);int d=0;
        s.radiative_acceleration.resize(nd);
        while(std::getline(input,line)) {
            std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);
            double m,T,P,rho,Q,ram,target,F,cooling,grad;
            if(!(row>>m>>T>>P>>rho>>Q>>ram>>target>>F>>cooling>>grad) || d>=nd || std::abs(m/s.mass[d]-1)>1e-12 || !(P>0) || !std::isfinite(grad))throw std::runtime_error("Invalid same-grid hydrostatic checkpoint");
            s.pressure[d]=P;s.rho[d]=P*cgs::mp/(2*cgs::kb*s.temperature[d]);s.radiative_acceleration[d]=grad;++d;
        }
        if(d!=nd)throw std::runtime_error("Incomplete hydrostatic checkpoint");
    }
    // Figure 3: allow fresh grey starts and arbitrary computed initial T profiles.
    // Pressure and density were rebuilt above from the new physical inputs.
    AndersonMixing temperature_mixing(4);int thermal_only_count=0;
    if(!get("--initial-field","").empty())initialize_radiation_checkpoint(s,g,get("--initial-field",""));
    std::ofstream history(out+"/iteration.csv");if(!history)throw std::runtime_error("Create output directory first");
    history<<std::setprecision(17)<<"iteration,transfer_iterations,transfer_residual,max_flux_error,max_cell_energy_error,temperature_step,Tmax,opacity_fallbacks,temperature_residual,hydrostatic_residual,Th,local_energy_weight,global_nodal_energy,deep_flux_derivative,stage_transfer_tolerance,min_flux_over_target,negative_flux_nodes,near_zero_flux_nodes\n";
    bool converged=false,stationary=false,requested_stop=false;int stationary_count=0,completed_iteration=0;double previous_flux_error=1e300;double ferr=1,cerr=1,change=0;TransferDiagnostics tr={1,0,0};
    double temperature_residual=get("--initial-field","").empty()?1:0,hydrostatic_residual=1;int anderson_fallbacks=0;
    for(int it=0;it<=outer;++it) {
        double weight=progress(it);
        th=continuation?ti+(th_final-ti)*weight:th_final;
        flux=FluxBudget::from_temperatures(th,ti);
        std::cout<<"Atmosphere delta="<<delta<<" iteration="<<it<<" Tmax="<<*std::max_element(s.temperature.begin(),s.temperature.end())<<std::endl;
        if(coupled_temperature) {
            // Evaluate the beam on the proposed temperature's density, not the
            // preceding temperature iterate's density. Iterate ram force and
            // pressure at fixed T; radiation force is refreshed by the RT solve.
            double pressure_change=1;
            for(int coupled=0;coupled<20 && pressure_change>1e-9;++coupled) {
                update_hydrostatic(s,gravity,1.);
                update_beam_heating(s,th,ti,delta,beam_panels,beam_relative_tau_step);
                auto old_pressure=s.pressure;update_hydrostatic(s,gravity,1.);
                pressure_change=0;for(int d=0;d<nd;++d)pressure_change=std::max(pressure_change,std::abs(s.pressure[d]/old_pressure[d]-1));
            }
            if(pressure_change>1e-6)throw std::runtime_error("Fixed-temperature beam/hydrostatic iteration did not converge");
        } else {
            if(it)update_beam_heating(s,th,ti,delta,beam_panels,beam_relative_tau_step);
            update_hydrostatic(s,gravity,.5);
        }
        std::vector<SparseScattering> kernels;kernels.reserve(nd);
        std::size_t dense_bytes=0;
        for(double T:s.temperature) {
            kernels.push_back(library.at(T,coefficient_cache_mb?3:compact_kernel));
            // Paired templates own their copied values. Once a layer is
            // assembled, subsequent layers in the same temperature bin use
            // that template directly and no longer need the dense tables.
            if(compact_kernel>=2) {
                dense_bytes+=library.dense_cache_bytes();library.release_dense_tables();
            }
        }
        if(coefficient_cache_mb) {
            // Storage3 and4 reconstruct exactly the same rounded coefficient.
            // Prioritize cooler layers, where departure-source evaluations use
            // the coefficients as well as the full collision operator.
            std::vector<int> order(nd);std::iota(order.begin(),order.end(),0);
            std::stable_sort(order.begin(),order.end(),[&](int a,int b){return s.temperature[a]<s.temperature[b];});
            std::size_t used=0,budget=std::size_t(coefficient_cache_mb)*1024*1024;int cached=0;
            for(int d:order) {
                auto& kernel=kernels[d];std::size_t bytes=kernel.start.back()*sizeof(float);
                if(bytes>budget-used)continue;
                std::vector<float> values(kernel.start.back());
                for(int i=0;i<g.size();++i)for(std::size_t z=kernel.start[i];z<kernel.start[i+1];++z)values[z]=float(kernel.symmetric_value(i,z));
                kernel.symmetric_kernel=std::move(values);used+=bytes;++cached;
            }
            std::cout<<"  coefficient cache bytes="<<used<<" budget="<<budget<<" cached layers="<<cached<<"/"<<nd<<std::endl;
        }
        dense_bytes+=library.dense_cache_bytes();library.release_dense_tables();
        std::cout<<"  released dense assembly cache bytes="<<dense_bytes<<std::endl;
        // SPW12 B.11 uses 1e-4 during atmosphere construction. Tighten as
        // equilibrium is approached; final acceptance retains the requested tolerance.
        // Both local thermal balance and the global flux profile must be
        // close before expensive tight radiation solves are useful. A small
        // local thermal correction alone can coexist with a large deep flux
        // deficit after remeshing. Final acceptance below is unchanged.
        const double construction_error=std::max(temperature_residual,ferr);
        double stage_tolerance=construction_error>.01?1e-4:construction_error>.001?1e-5:transfer_tolerance;
        // A saved field is only a warm start, not evidence of an equilibrated
        // atmosphere. Start staged construction loosely even when one exists.
        if(it==0 && native_balanced && balanced_staged_radiation && outer>0)stage_tolerance=1e-4;
        stage_tolerance=std::max(transfer_tolerance,stage_tolerance);
        if(it==outer || external_temperature || (native_balanced && !balanced_staged_radiation))stage_tolerance=transfer_tolerance;
        // Staging affects construction only. Final acceptance below still
        // requires transfer_tolerance; a zero-update verification is always tight.
        paper_response.checkpoint_root=out;paper_response.outer_iteration=it;
        tr=atmosphere_transfer(s,g,kernels,flux.intrinsic,inner,stage_tolerance);
        if(tr.residual>stage_tolerance)std::cout<<"RT NOT CONVERGED: saving diagnostic state before stopping."<<std::endl;
        ferr=0;cerr=0;
        std::vector<double> tauR(nd);tauR[0]=cgs::kappa_e*s.rosseland[0]*s.mass[0];
        for(int d=1;d<nd;++d)tauR[d]=tauR[d-1]+.5*cgs::kappa_e*(s.rosseland[d]+s.rosseland[d-1])*(s.mass[d]-s.mass[d-1]);
        double global_nodal_energy=0,deep_flux_derivative=0,min_flux_ratio=1e300;int negative_flux_nodes=0,near_zero_flux_nodes=0;

        for(int d=0;d<nd;++d) {
            ferr=std::max(ferr,std::abs(s.flux[d]/s.target_flux[d]-1));
            min_flux_ratio=std::min(min_flux_ratio,s.flux[d]/s.target_flux[d]);negative_flux_nodes+=(s.flux[d]<0);near_zero_flux_nodes+=(std::abs(s.flux[d])<1e-12*std::abs(s.target_flux[d]));
            double width=(conservative || s.conservative_cells)?s.cell_width[d]:d==0?s.mass[1]-s.mass[0]:d==nd-1?s.mass[nd-1]-s.mass[nd-2]:.5*(s.mass[d+1]-s.mass[d-1]);
            double local_energy=std::abs(s.cooling[d]-s.heating[d])*width/flux.total();
            global_nodal_energy=std::max(global_nodal_energy,local_energy);
            if(tauR[d]<1)cerr=std::max(cerr,local_energy);
            if(d && tauR[d]>1)deep_flux_derivative=std::max(deep_flux_derivative,std::abs((s.flux[d]-s.target_flux[d])-(s.flux[d-1]-s.target_flux[d-1]))/((tauR[d]-tauR[d-1])*flux.total()));
        }
        if(conservative || s.conservative_cells)for(int d=0;d<=nd;++d)ferr=std::max(ferr,std::abs(s.face_flux[d]/s.target_face_flux[d]-1));
        temperature_residual=0;hydrostatic_residual=0;
        double desired_pressure=(gravity+s.ram[0]-s.radiative_acceleration[0])*s.mass[0];
        for(int d=0;d<nd;++d) {
            if(tauR[d]<1)temperature_residual=std::max(temperature_residual,std::abs(s.cooling[d]-s.heating[d])/std::max(s.local_derivative[d],1e-100));
            if(d)desired_pressure+=(s.mass[d]-s.mass[d-1])*(gravity+.5*(s.ram[d]+s.ram[d-1]-s.radiative_acceleration[d]-s.radiative_acceleration[d-1]));
            hydrostatic_residual=std::max(hydrostatic_residual,std::abs(s.pressure[d]/desired_pressure-1));
        }
        converged=weight==1 && tr.residual<=transfer_tolerance && ferr<.001 && cerr<.0001 && deep_flux_derivative<.0001 && temperature_residual<temperature_tolerance && hydrostatic_residual<.001;
        completed_iteration=it;
        if(it>0 && weight==1 && change<1e-5 && std::abs(ferr-previous_flux_error)<2e-6 && temperature_residual<.002 && tr.residual<=transfer_tolerance)++stationary_count;
        else stationary_count=0;
        stationary=stationary_count>=30 && !converged;previous_flux_error=ferr;
        history<<it<<','<<tr.iterations<<','<<tr.residual<<','<<ferr<<','<<cerr<<','<<change<<','<<*std::max_element(s.temperature.begin(),s.temperature.end())<<','<<tr.effective_opacity_fallbacks<<','<<temperature_residual<<','<<hydrostatic_residual<<','<<th<<','<<weight<<','<<global_nodal_energy<<','<<deep_flux_derivative<<','<<stage_tolerance<<','<<min_flux_ratio<<','<<negative_flux_nodes<<','<<near_zero_flux_nodes<<'\n';history.flush();
        std::cout<<"  transfer="<<tr.residual<<" flux="<<ferr<<" upper energy="<<cerr<<" local T="<<temperature_residual<<" deep dF="<<deep_flux_derivative<<" HSE="<<hydrostatic_residual<<std::endl;
        std::ofstream state(out+"/structure.csv");state<<std::setprecision(17)<<"mass,T,P,rho,Qplus,gram,Ftarget,Fcomputed,cooling,grad,rosseland,local_derivative,temperature_response\n";
        for(int d=0;d<nd;++d)state<<s.mass[d]<<','<<s.temperature[d]<<','<<s.pressure[d]<<','<<s.rho[d]<<','<<s.heating[d]<<','<<s.ram[d]<<','<<s.target_flux[d]<<','<<s.flux[d]<<','<<s.cooling[d]<<','<<s.radiative_acceleration[d]<<','<<s.rosseland[d]<<','<<s.local_derivative[d]<<','<<s.temperature_response[d]<<'\n';
        if(conservative || s.conservative_cells) {
            std::ofstream balance(out+"/cell_balance.csv");balance<<std::setprecision(17)<<"mass,width,Fleft,Fright,cooling,Qplus\n";
            for(int d=0;d<nd;++d)balance<<s.mass[d]<<','<<s.cell_width[d]<<','<<s.face_flux[d]<<','<<s.face_flux[d+1]<<','<<s.cooling[d]<<','<<s.heating[d]<<'\n';
        }
        requested_stop=bool(std::ifstream(out+"/STOP_REQUEST"));
        if(requested_stop)std::cout<<"STOP_REQUEST observed; saving current diagnostic state."<<std::endl;
        if(checkpoint_every && (it%checkpoint_every==0 || converged || stationary || requested_stop || it==outer))
            write_iteration_snapshot(s,g,out,it,emin,emax,get("--gaunt-table",""),tr.residual);
        if(converged || stationary || requested_stop || it==outer || tr.residual>stage_tolerance)break;
        if(external_temperature) {
            state.close();std::cout<<"TEMPERATURE_READY "<<it<<std::endl;
            std::string request;if(!(std::cin>>request) || request=="STOP")break;
            if(request!="TEMPERATURE")throw std::runtime_error("Invalid external temperature request");
            std::vector<double> proposed(nd);change=0;
            for(int d=0;d<nd;++d) {
                if(!(std::cin>>proposed[d]) || !(proposed[d]>0) || !std::isfinite(proposed[d]))throw std::runtime_error("Invalid external temperature");
                change=std::max(change,std::abs(proposed[d]/s.temperature[d]-1));
            }
            s.temperature.swap(proposed);continue;
        }
        if(native_balanced){
            // Free current coefficient vectors before probes; templates remain cacheable.
            kernels.clear();
            TemperatureProbe probe=[&](const HeatedState& base,int depth,double shift){
                HeatedState input=base;if(depth>=0)input.temperature[depth]*=std::exp(shift);
                std::cout<<"FD_PROBE start iteration="<<it<<" depth="<<depth<<" logT_shift="<<shift<<std::endl;
                CoupledProbeOptions o{gravity,th,ti,beam_relative_tau_step,paper_response.fd_tolerance,delta,beam_panels,compact_kernel,coefficient_cache_mb,inner};
                auto result=coupled_temperature_probe(input,g,library,o);
                std::string name=out+"/fd_"+std::to_string(it)+"_"+std::to_string(depth)+(shift<0?"_minus":shift>0?"_plus":"_base")+".csv";
                std::ofstream f(name);f<<std::setprecision(17)<<"mass,T,P,rho,Qplus,gram,Ftarget,Fcomputed,cooling,grad,rosseland,local_derivative,temperature_response\n";
                for(int d=0;d<nd;++d)f<<result.mass[d]<<','<<result.temperature[d]<<','<<result.pressure[d]<<','<<result.rho[d]<<','<<result.heating[d]<<','<<result.ram[d]<<','<<result.target_flux[d]<<','<<result.flux[d]<<','<<result.cooling[d]<<','<<result.radiative_acceleration[d]<<','<<result.rosseland[d]<<','<<result.local_derivative[d]<<','<<result.temperature_response[d]<<'\n';
                f.flush();if(!f)throw std::runtime_error("Cannot save numerical Jacobian probe");return result;
            };
            change=balanced_temperature(s,temperature_step,balanced_hot_response,probe);continue;
        }
        change=paper_temperature(s,relax,temperature_step,weight,temperature_step*(1-weight)+final_flux_step*weight);

    }
    std::ofstream spectrum(out+"/intensity_grid.csv");spectrum<<std::setprecision(17)<<"energy_keV,mu,I_E\n";
    double limiter_flux_fraction=0;
    if(conservative) {
        std::ofstream limiter(out+"/positivity_limiter.csv");limiter<<std::setprecision(17)<<"depth_index,mass,energy_keV,mu,delta_n_odd,I_E_change_bound,F_direction_change_bound\n";
        for(int d=0;d<nd;++d) {
            double bound=0;
            for(int e=0;e<ne;++e)for(int a=0;a<na;++a) {
                double change=s.last_v_limiter[d*ne*na+e*na+a];if(change==0)continue;
                int i=g.state(e,na+a);double I=occupation_intensity_x(g.x[e],change)*cgs::kev/cgs::mec2;
                double F=4*cgs::pi*occupation_intensity_x(1,1)*g.weight[i]*g.x[e]*g.mu[na+a]*change;
                bound+=F;limiter<<d<<','<<s.mass[d]<<','<<g.x[e]*cgs::mec2/cgs::kev<<','<<g.mu[na+a]<<','<<change<<','<<I<<','<<F<<'\n';
            }
            limiter_flux_fraction=std::max(limiter_flux_fraction,bound/s.target_flux[d]);
        }
    }
    // Full field checkpoint permits independent ray and collision diagnostics
    // without repeating the converged atmospheric iteration.
    std::ofstream field(out+"/occupation.bin",std::ios::binary);
    const int dimensions[]={nd,ne,2*na};
    field.write(reinterpret_cast<const char*>(dimensions),sizeof(dimensions));
    field.write(reinterpret_cast<const char*>(s.occupation.data()),s.occupation.size()*sizeof(double));
    if(!field)throw std::runtime_error("Cannot save radiation checkpoint");
    for(int e=0;e<ne;++e)for(int a=na;a<2*na;++a) {
        int i=g.state(e,a);double E=g.x[e]*cgs::mec2/cgs::kev,I=occupation_intensity_x(g.x[e],s.occupation[i])*cgs::kev/cgs::mec2;
        spectrum<<E<<','<<g.mu[a]<<','<<I<<'\n';
    }
    std::ofstream manifest(out+"/status.json");manifest<<std::setprecision(17)<<"{\"delta\":"<<delta<<",\"Th\":"<<th<<",\"Ti\":"<<ti
        <<",\"logg\":"<<logg<<",\"gamma_min\":"<<s.beam_gamma_min<<",\"gamma_max\":"<<s.beam_gamma_max
        <<",\"mono_gamma\":"<<s.beam_mono_gamma<<",\"bremsstrahlung_beam\":"<<(s.beam_bremsstrahlung?"true":"false")
        <<",\"startup_iterations\":"<<startup<<",\"heating_ramp_iterations\":"<<ramp<<",\"schedule_offset\":"<<schedule_offset<<",\"final_flux_step\":"<<final_flux_step
        <<",\"gaunt\":"<<(gaunt_path.empty()?"1":"null")<<",\"gaunt_table\":\""<<gaunt_path<<"\",\"energy_bins\":"<<ne
        <<",\"depth_nodes\":"<<nd<<",\"directions_per_hemisphere\":"<<na<<",\"kernel_log_temperature_step\":"<<step
        <<",\"energy_min_keV\":"<<emin<<",\"energy_max_keV\":"<<emax<<",\"mass_min\":"<<mmin<<",\"mass_max\":"<<mmax<<",\"beam_panels\":"<<beam_panels
        <<",\"mass_grid_file\":\""<<mass_grid_path<<"\""
        <<",\"beam_relative_tau_step\":"<<beam_relative_tau_step
        <<",\"beam_remap\":\""<<beam_remap<<"\""
        <<",\"rt_extrapolation\":\""<<(paper_response.rt_extrapolate?"Anderson_on_published_ALI":"none")<<"\""
        <<",\"rt_history\":"<<paper_response.rt_history
        <<",\"rt_occupation_floor\":"<<paper_response.rt_occupation_floor
        <<",\"rt_planck_intensity_floor_fraction\":"<<paper_response.rt_planck_intensity_floor_fraction
        <<",\"pair_cache_budget_mb\":"<<(std::getenv("SALMI_PAIR_CACHE_MB")?std::stoull(std::getenv("SALMI_PAIR_CACHE_MB")):512)
        <<",\"compact_kernel\":"<<compact_kernel<<",\"coefficient_cache_budget_mb\":"<<coefficient_cache_mb
        <<",\"dense_cache_release_policy\":\""<<(compact_kernel>=2?"after_each_paired_layer":"after_all_layers")<<"\""
        <<",\"balanced_full_thermal_response\":"<<(paper_response.full_thermal_response?"true":"false")
        <<",\"balanced_step_control\":\""<<step_control<<"\",\"balanced_fd_columns\":"<<paper_response.fd_columns<<",\"balanced_fd_every\":"<<paper_response.fd_every<<",\"balanced_fd_step\":"<<paper_response.fd_step<<",\"balanced_fd_tolerance\":"<<paper_response.fd_tolerance
        <<",\"compton_integration\":\""<<(redistribution_integration_version()?"moment_v1":"midpoint")<<"\""
        <<",\"intrinsic_to_heated_flux_ratio\":"<<(th>0?std::pow(ti/th,4):0)
        <<",\"balanced_energy_scale\":"<<paper_response.balanced_energy_scale
        <<",\"balanced_staged_radiation\":"<<(balanced_staged_radiation?"true":"false")
        <<",\"threads\":"<<threads
        <<",\"requested_stop\":"<<(requested_stop?"true":"false")<<",\"stationary_non_equilibrium\":"<<(stationary?"true":"false")
        <<",\"completed_temperature_iteration\":"<<completed_iteration
        <<",\"anderson\":"<<(s.anderson?"true":"false")
        <<",\"radiation_initialization\":\""<<(get("--initial-field","").empty()?"Planck":"computed_checkpoint")<<"\""
        <<",\"transfer_tolerance\":"<<transfer_tolerance<<",\"temperature_tolerance\":"<<temperature_tolerance<<",\"temperature_relaxation\":"<<relax<<",\"max_log_temperature_step\":"<<temperature_step<<",\"anderson_fallbacks\":"<<anderson_fallbacks
        <<",\"max_positivity_flux_change_bound\":"<<limiter_flux_fraction
        <<",\"solver\":\""<<(conservative?"feautrier":"short_characteristic")<<"\",\"formal_order\":\""<<(s.parabolic?"hybrid_parabolic":"linear")<<"\",\"temperature_iterations_requested\":"<<outer<<",\"initialization\":\""<<origin
        <<"\",\"discrete_converged\":"<<(converged?"true":"false")<<",\"transfer_residual\":"<<tr.residual
        <<",\"energy_discretization\":\"pointwise_paper_equation22\",\"temperature_update\":\""<<(external_temperature?"external_residual_solver":native_balanced?"native_balanced_Broyden_SVD":"SPW12_eq28_ATLAS_AK_surface_signed_flux")<<"\",\"max_flux_error\":"<<ferr<<",\"max_cell_energy_error\":"<<cerr<<",\"temperature_residual\":"<<temperature_residual<<",\"hydrostatic_residual\":"<<hydrostatic_residual<<",\"figure3_certified\":false}\n";manifest.close();
    std::ofstream extra(out+"/paper_diagnostics.json");
    extra<<std::setprecision(17)<<"{\"bottom_min_sqrt_tau_abs_tau_total\":"<<paper_response.bottom_thermalization<<",\"ali_guard_count\":"<<paper_response.ali_guards<<"}\n";
    if(stationary)std::cout<<"Stationary iteration, but equilibrium tolerances NOT met. Requires a numerical/discretization check."<<std::endl;
    std::cout<<"Discrete equilibrium: "<<(converged?"PASS":"NOT CONVERGED")<<"; Figure 3 not certified.\n";
    return require_equilibrium && !converged?3:0;
 }catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 2;}
}
