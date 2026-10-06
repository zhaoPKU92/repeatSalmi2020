#pragma once
#include "paper.hpp"
namespace salmi {
struct CoupledProbeOptions {
 double gravity,heated_temperature,intrinsic_temperature,beam_tau_step,radiation_tolerance;
 int delta,beam_panels,compact_kernel,coefficient_cache_mb,max_radiation_iterations;
};
// Re-solve radiation, beam heating and hydrostatic balance at fixed T.
// Used for numerical Jacobian probes, never a substituted physical operator.
HeatedState coupled_temperature_probe(const HeatedState&,const PhotonGrid&,KernelLibrary&,const CoupledProbeOptions&);
}
