#pragma once
#include <functional>
#include <vector>
namespace salmi {
struct Quadrature { std::vector<double> x, w; };
std::vector<double> log_grid(double lo, double hi, int n);
Quadrature gauss_legendre(int n, double lo, double hi);
double integrate(const std::function<double(double)>& f, double lo, double hi, int n=64);
double power_integral(double lo, double hi, double exponent);
double linear_interp(const std::vector<double>& x, const std::vector<double>& y, double v);
double rk4_step(const std::function<double(double,double)>& f, double t, double y, double h);
std::vector<double> tridiagonal(std::vector<double> a, std::vector<double> b,
                               std::vector<double> c, std::vector<double> r);
}
