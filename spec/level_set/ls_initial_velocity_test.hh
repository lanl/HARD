#ifndef SPEC_LS_INITIAL_VELOCITY_TEST
#define SPEC_LS_INITIAL_VELOCITY_TEST
#include <cmath>

// This is the verification purpose only
// The setup is adapted from Tariq et al. 2025

namespace spec::level_set {

struct TimeReversal1D {
  double T{2.0};
  inline void sample(double x,
    double,
    double,
    double t,
    double & u,
    double & v,
    double & w) const {
    u = 2.0 * std::sin(2.0 * M_PI * x) * std::cos(M_PI * t / T);
    v = 0.0;
    w = 0.0;
  }
};

struct TimeReversal2D {
  double T{2.0};
  inline void sample(double x,
    double y,
    double t,
    double & u,
    double & v,
    double & w) const {
    u = -2.0 * std::cos(M_PI * y) * std::sin(M_PI * x) * std::sin(M_PI * x) *
        std::sin(M_PI * y) * std::cos(M_PI * t / T);
    v = 2.0 * std::cos(M_PI * x) * std::sin(M_PI * x) * std::sin(M_PI * y) *
        std::sin(M_PI * y) * std::cos(M_PI * t / T);
    w = 0.0;
  }
};

struct TimeReversal3D {
  double T{1.0};
  inline void sample(double x,
    double y,
    double z,
    double t,
    double & u,
    double & v,
    double & w) const {
    const double sx = std::sin(M_PI * x), sy = std::sin(M_PI * y),
                 sz = std::sin(M_PI * z);
    u = 2.0 * sx * sx * std::sin(2 * M_PI * y) * std::sin(2 * M_PI * z) *
        std::cos(M_PI * t / T);
    v = -2.0 * std::sin(2 * M_PI * x) * sy * sy * std::sin(2 * M_PI * z) *
        std::cos(M_PI * t / T);
    w = -2.0 * std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y) * sz * sz *
        std::cos(M_PI * t / T);
  }
};
} // namespace spec::level_set

#endif
