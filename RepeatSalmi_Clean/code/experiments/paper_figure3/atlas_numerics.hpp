#pragma once
#include <vector>
namespace salmi {
std::vector<double> atlas_derivative(const std::vector<double>& x,const std::vector<double>& f);
std::vector<double> atlas_integral(const std::vector<double>& x,const std::vector<double>& f,double start=0);
}
