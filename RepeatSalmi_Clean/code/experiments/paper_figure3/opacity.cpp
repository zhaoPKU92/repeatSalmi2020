#include "salmi/heated.hpp"
namespace salmi {
double HeatedState::absorption(double E,int d) const {
    // Explicit user choice for this reconstruction: constant free-free Gaunt factor.
    return free_free_opacity(E,temperature[d],rho[d],1.0);
}
}
