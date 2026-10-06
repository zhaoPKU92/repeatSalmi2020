#pragma once
#include "paper.hpp"
#include <functional>
namespace salmi {
// Joint local flux, microscopic upper energy, and deep flux-gradient correction.
// Gray transfer is a numerical Jacobian preconditioner only.
using TemperatureProbe=std::function<HeatedState(const HeatedState&,int,double)>;
double balanced_temperature(HeatedState&,double maximum_log_step,bool hot_response,const TemperatureProbe& probe=TemperatureProbe{});
}
