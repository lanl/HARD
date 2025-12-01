#ifndef SPEC_LS_TASK
#define SPEC_LS_TASK

#include "spec/level_set/ls_fields.hh"
#include "spec/level_set/ls_velocity.hh"
#include <flecsi/field/field.hh>
#include <flecsi/runtime/scheduler.hh>

namespace spec::level_set {

struct StepParams {
  double t, dt;
  double dx, dy, dz;
  int nx, ny, nz;
  bool periodic{true};
};

// === Task: advance once with SSPRK(3,3) + Algo1 ===
template<class Vel>
void
advance_task(hard::mesh::accessor<flecsi::ro> topo,
  flecsi::field<double>::accessor<flecsi::rw> psi_a,
  flecsi::field<int>::accessor<flecsi::rw> lab_a,
  StepParams p,
  Vel vel) {

  using Exec = Kokkos::DefaultExecutionSpace;
  const int N = p.nx * p.ny * p.nz;

  auto psi_span = psi_a.span(); // contiguous (device) storage
  auto lab_span = lab_a.span();

  // Work arrays (stage outputs & L)
  // (You can also allocate these as FleCSI scratch fields; kept simple here.)
  Kokkos::View<double *, Exec> psi1("psi1", N), psi2("psi2", N),
    psi3("psi3", N), L("L", N);

  hard::ls::GridView gv{p.nx, p.ny, p.nz, p.dx, p.dy, p.dz, p.periodic};

  // Stage 1
  {
    UpwindHJ<Exec> hj{Kokkos::View<const double *, Exec>(psi_span.data(), N),
      Kokkos::View<const int *, Exec>(lab_span.data(), N),
      L,
      gv,
      p.t,
      /*T=*/(p.nz > 1 ? 1.0 : 2.0)};
    hj(vel);
    Kokkos::parallel_for(
      "rk1", Kokkos::RangePolicy<Exec>(0, N), KOKKOS_LAMBDA(int i) {
        psi1(i) = psi_span[i] + p.dt * L(i);
      });
  }
  // Stage 2
  {
    UpwindHJ<Exec> hj{Kokkos::View<const double *, Exec>(psi1.data(), N),
      Kokkos::View<const int *, Exec>(lab_span.data(), N),
      L,
      gv,
      p.t + p.dt,
      (p.nz > 1 ? 1.0 : 2.0)};
    hj(vel);
    Kokkos::parallel_for(
      "rk2", Kokkos::RangePolicy<Exec>(0, N), KOKKOS_LAMBDA(int i) {
        psi2(i) = 0.75 * psi_span[i] + 0.25 * (psi1(i) + p.dt * L(i));
      });
  }
  // Stage 3
  {
    UpwindHJ<Exec> hj{Kokkos::View<const double *, Exec>(psi2.data(), N),
      Kokkos::View<const int *, Exec>(lab_span.data(), N),
      L,
      gv,
      p.t + 0.5 * p.dt,
      (p.nz > 1 ? 1.0 : 2.0)};
    hj(vel);
    Kokkos::parallel_for(
      "rk3", Kokkos::RangePolicy<Exec>(0, N), KOKKOS_LAMBDA(int i) {
        psi3(i) =
          (1.0 / 3.0) * psi_span[i] + (2.0 / 3.0) * (psi2(i) + p.dt * L(i));
      });
  }

  // Algorithm 1 (4-neigh in 2D / 6-neigh in 3D)
  {
    Algo1Fix<Exec> a1{psi3,
      Kokkos::View<const int *, Exec>(lab_span.data(), N),
      Kokkos::View<double *, Exec>(psi_span.data(), N),
      Kokkos::View<int *, Exec>(lab_span.data(), N),
      gv};
    a1();
  }
}

} // namespace spec::level_set
