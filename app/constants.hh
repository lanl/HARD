
#pragma once

namespace hard {

namespace constants {

namespace cgs {

//
// Adopted from 2024 Review of Particle Physics by Particle Data Group
// (see https://pdg.lbl.gov/2024/reviews/rpp2024-rev-phys-constants.pdf)
//
constexpr double speed_of_light{2.99792458e10};
constexpr double stefan_boltzmann_constant{5.670374419e-5};
constexpr double boltzmann_constant{1.380649e-16};

// a = 4 * stefan_boltzmann_constant / speed of light
constexpr double radiation_constant{7.56573325e-15};

constexpr double proton_mass{1.67262192369e-24};

} // namespace cgs

} // namespace constants

} // namespace hard
