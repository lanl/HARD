#pragma once

#include "../../constants.hh"
#include "../../types.hh"
#include <cmath>
#include <cstddef>

namespace hard::tasks::initial_data {

namespace rad_shock {

struct rad_rankine_hugoniot {
  // Values for left and right states are taken
  // from section 4.2 in Moens et al., A&A 657,A81 (22)
  static constexpr double rL = 1.0e-2;
  static constexpr double uL = 1.0e9;
  static constexpr double vL = 0.0;
  static constexpr double wL = 0.0;
  static constexpr double TL = 1.0e4;
  static constexpr double EL = 5.0e15; // total gas energy
  static constexpr double EradL = 75.6; // radiation energy

  static constexpr double rR = 6.85847e-2;
  static constexpr double uR = rL * uL / rR;
  static constexpr double vR = 0.0;
  static constexpr double wR = 0.0;
  // static constexpr double TR = 153.34928811;
  static constexpr double TR = 42285513.562735;
  // static constexpr double TR = 4.239e7;
  static constexpr double ER = 1.93e15;
  static constexpr double EradR = 2.44e16;

  static constexpr double x0 = 5.0e4;
  static constexpr double y0 = 0.0;
}; // struct rad_rankine_hugoniot

} // namespace rad_shock

/*----------------------------------------------------------------------------*
  Radiation modified Rankine Hugoniot.
 *----------------------------------------------------------------------------*/

template<typename T, std::size_t D>
auto
rad_RH(typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<wo, ro> r_a,
  typename field<vec<D>>::template accessor<wo, ro> ru_a,
  field<double>::accessor<wo, ro> rE_a,
  field<double>::accessor<wo, ro> Erad_a,
  const double gamma,
  single<double>::accessor<ro> particle_mass_a) {
  auto r = m.template mdcolex<is::cells>(r_a);
  auto ru = m.template mdcolex<is::cells>(ru_a);
  auto rE = m.template mdcolex<is::cells>(rE_a);
  auto Erad = m.template mdcolex<is::cells>(Erad_a);
  auto const particle_mass = *particle_mass_a;
  const double mult = 1.0 / (gamma - 1.0);

  const double kb = hard::constants::cgs::boltzmann_constant;
  const double a = hard::constants::cgs::radiation_constant;

  if constexpr(D == 1) {
    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      const auto x = m.template center<ax::x>(i);

      if(x < T::x0) {
        r(i) = T::rL;
        ru(i).x = T::rL * T::uL;
        rE(i) = mult * kb * T::TL * T::rL / (particle_mass) +
                (0.5 * T::rL * T::uL * T::uL);
        Erad(i) = a * T::TL * T::TL * T::TL * T::TL;
      }
      else {
        r(i) = T::rR;
        ru(i).x = T::rR * T::uR;
        // rE(i) = T::ER;
        rE(i) = mult * kb * T::TR * T::rR / (particle_mass) +
                (0.5 * T::rR * T::uR * T::uR);
        Erad(i) = a * T::TR * T::TR * T::TR * T::TR;
        // Erad(i) = T::EradR;
      } // if
      if(rE(i) <= 0) {
        std::cout << "TotalE is negative for i = " << i << std::endl;
      }
    } // for
  }
  else if constexpr(D == 2) {
    for(auto j : m.template cells<ax::y, dm::quantities>()) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        const auto x = m.template center<ax::x>(i);

        if(x < T::x0) {
          r(i, j) = T::rL;
          ru(i, j).x = T::rL * T::uL;
          ru(i, j).y = T::rL * T::vL;
          // rE(i, j) = T::EL;
          // Erad(i, j) = T::EradL;
          rE(i, j) = mult * kb * T::TL * T::rL / (particle_mass) +
                     (0.5 * T::rL * T::uL * T::uL);
          Erad(i, j) = a * T::TL * T::TL * T::TL * T::TL;
        }
        else {
          r(i, j) = T::rR;
          ru(i, j).x = T::rR * T::uR;
          ru(i, j).y = T::rR * T::vR;
          // rE(i, j) = T::ER;
          // Erad(i, j) = T::EradR;
          rE(i, j) = mult * kb * T::TR * T::rR / (particle_mass) +
                     (0.5 * T::rR * T::uR * T::uR);
          Erad(i, j) = a * T::TR * T::TR * T::TR * T::TR;
        } // if
      } // for
    } // for
  }
  else /* D == 3 */ {
    for(auto k : m.template cells<ax::z, dm::quantities>()) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          const auto x = m.template center<ax::x>(i);

          if(x < T::x0) {
            r(i, j, k) = T::rL;
            ru(i, j, k).x = T::rL * T::uL;
            ru(i, j, k).y = T::rL * T::vL;
            ru(i, j, k).z = T::rL * T::wL;
            // rE(i, j, k) = T::EL;
            // Erad(i, j, k) = T::EradL;
            rE(i, j, k) = mult * kb * T::TL * T::rL / (particle_mass) +
                          (0.5 * T::rL * T::uL * T::uL);
            Erad(i, j, k) = a * T::TL * T::TL * T::TL * T::TL;
          }
          else {
            r(i, j, k) = T::rR;
            ru(i, j, k).x = T::rR * T::uR;
            ru(i, j, k).y = T::rR * T::vR;
            ru(i, j, k).z = T::rR * T::wR;
            // rE(i, j, k) = T::ER;
            // Erad(i, j, k) = T::EradR;
            rE(i, j, k) = mult * kb * T::TR * T::rR / (particle_mass) +
                          (0.5 * T::rR * T::uR * T::uR);
            Erad(i, j, k) = a * T::TR * T::TR * T::TR * T::TR;
          } // if
        } // for
      } // for
    } // for
  } // if
} // rad_RH

} // namespace hard::tasks::initial_data
