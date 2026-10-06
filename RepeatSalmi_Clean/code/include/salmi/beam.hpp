#pragma once
#include "salmi/numerics.hpp"
namespace salmi {
class PowerLaw {
public:
    double minimum,maximum,delta,normalization;
    PowerLaw(double minimum,double maximum,double delta);
    double density(double gamma) const;
    double mean_gamma() const;
    // Integrate f(gamma) dgamma in log(gamma), including the Jacobian gamma.
    Quadrature quadrature(int panels, int order=4) const;
};
}
