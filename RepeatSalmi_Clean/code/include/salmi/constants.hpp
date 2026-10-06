#pragma once
// CGS throughout. Energies accepted by radiation routines are explicitly in keV.
namespace salmi { namespace cgs {
constexpr double pi = 3.1415926535897932384626433832795;
constexpr double c = 2.99792458e10;
constexpr double h = 6.62607015e-27;
constexpr double hbar = h/(2*pi);
constexpr double kb = 1.380649e-16;
constexpr double me = 9.1093837015e-28;
constexpr double mp = 1.67262192369e-24;
constexpr double e = 4.803204712570263e-10;
constexpr double alpha = 7.2973525693e-3;
constexpr double sigma_t = 6.6524587321e-25;
constexpr double sigma_sb = 5.670374419e-5;
constexpr double kev = 1.602176634e-9;
constexpr double mec2 = me*c*c;
constexpr double kappa_e = sigma_t/mp;
constexpr double grav = 6.67430e-8;
constexpr double solar_mass = 1.98847e33;
} }
