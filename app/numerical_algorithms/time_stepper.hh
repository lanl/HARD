
#pragma once

namespace flastro::time_stepper {
//
// IMEX-SSP2(2,2,2)  (Pareschi & Russo 2004, arxiv:1009.2757)
//
// gamma = 1 - 1/sqrt(2)
constexpr double time_stepper_gamma{0.29289321881345254};
constexpr double one_minus_2_time_stepper_gamma{0.4142135623730949};

enum rk_stage { First, Second, Update };

} // namespace flastro::time_stepper
