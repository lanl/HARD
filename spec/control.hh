#ifndef SPEC_CONTROL_HH
#define SPEC_CONTROL_HH

#include "types.hh"
#include <flecsi/execution.hh>
#include <flecsi/flog.hh>
#include <flecsi/run/control.hh>

#include <fstream>

namespace spec {

/// Control Points.
enum class cp {
  ///
  initialize,
  ///
  advance,
  ///
  analyze,
  ///
  finalize
};

inline const char *
operator*(cp control_point) {
  switch(control_point) {
    case cp::initialize:
      return "initialize";
    case cp::advance:
      return "advance";
    case cp::analyze:
      return "analyze";
    case cp::finalize:
      return "finalize";
  }
  flog_fatal("invalid control point");
}

/// MUSCL-Hancock Solver Control Model.
template<template<std::size_t> typename S, std::size_t D>
struct control_policy : flecsi::run::control_base {

  static constexpr std::size_t dimension = D;

#ifdef HARD_ENABLE_LEGION_TRACING
  flecsi::exec::trace tracing;
  std::optional<flecsi::exec::trace::guard> guard;
#endif

  using control_points_enum = cp;

  static void init_dt(
    flecsi::field<double, flecsi::data::single>::template accessor<flecsi::rw>
      t,
    double t_) {
    t = t_;
  }

  control_policy(double t0,
    double tf,
    std::size_t max_steps,
    double cfl,
    double max_dt,
    std::size_t log_frequency,
    std::size_t output_frequency)
    : t0_(t0), tf_(tf), t_(t0), max_steps_(max_steps), cfl_(cfl),
      max_dt_(max_dt), log_frequency_(log_frequency),
      output_frequency_(output_frequency) {}

  S<D> & state() {
    return state_;
  }

  std::size_t step() const {
    return step_;
  }

  std::size_t output_frequency() const {
    return output_frequency_;
  }

  std::size_t max_steps() const {
    return max_steps_;
  }

  auto time() const {
    return t_;
  }

  auto max_time() const {
    return tf_;
  }

  static void compute_dt(
    typename single<double>::template accessor<flecsi::wo> t,
    typename single<double>::template accessor<flecsi::wo> dt,
    flecsi::future<double> dtmin,
    double tf,
    double max_dt,
    double cfl) {
    dt = cfl * dtmin.get();
    dt = t + dt > tf ? tf - t : dt;
    dt = std::min(*dt, max_dt);
    t += dt;
  }

  static std::tuple<double, double> compute_dt_mpi(
    typename single<double>::template accessor<flecsi::wo> t,
    typename single<double>::template accessor<flecsi::wo> dt,
    flecsi::future<double> dtmin,
    double tf,
    double max_dt,
    double cfl) {
    dt = cfl * dtmin.get();
    dt = t + dt > tf ? tf - t : dt;
    dt = std::min(*dt, max_dt);
    t += dt;
    return std::make_tuple(t, dt);
  }

  static bool cycle_control(control_policy & cp) {
#ifdef HARD_BENCHMARK_MODE

    // Time each cycle
    if(cp.step_ == 0) {
      // initialize timer
      cp.start_timer_ = std::chrono::system_clock::now();
    }
    else {
      std::chrono::time_point<std::chrono::system_clock> stop_timer =
        std::chrono::system_clock::now();
      double runtime = (stop_timer - cp.start_timer_).count() / 1e9;
      cp.runtimes_ = cp.runtimes_ + std::to_string(runtime) + ";";
      // updates for next iteration
      cp.start_timer_ = stop_timer;
      cp.total_runtime_ += runtime;
      if(cp.step_ == cp.max_steps_ && flecsi::process() == 0) {
#if defined(FLECSI_ENABLE_HPX)
        auto threads = hpx::get_os_thread_count();
#else
        auto threads = 1; // omp_get_num_threads();
#endif
        std::ofstream runtime_file;
        runtime_file.open("result/runtimes.txt", std::ios_base::app);
        runtime_file << flecsi::processes() << ";" << threads << ";"
                     << cp.max_steps_ << ";" << cp.total_runtime_ << ";"
                     << cp.runtimes_ << std::endl;
        runtime_file.close();
      }
    }

#endif

    bool exec_cycle = cp.step_ < cp.max_steps_;

    auto & s = cp.state();

#if FLECSI_BACKEND == FLECSI_BACKEND_legion
    flecsi::execute<compute_dt>(
      s.t(*s.gt), s.dt(*s.gt), s.dtmin_, cp.tf_, cp.max_dt_, cp.cfl_);

    if((cp.step_ % cp.log_frequency_) == 0 || cp.step_ == cp.max_steps_) {
      flog(info) << "step: " << cp.step_ << "/" << cp.max_steps_ << std::endl;
      flecsi::flog::flush();
    } // if

#else
    auto [t, dt] = flecsi::execute<compute_dt_mpi>(
      s.t(*s.gt), s.dt(*s.gt), s.dtmin_, cp.tf_, cp.max_dt_, cp.cfl_)
                     .get();
    cp.t_ = t;

    if((cp.step_ % cp.log_frequency_) == 0 || cp.step_ == cp.max_steps_ ||
       cp.t_ == cp.tf_) {
      flog(info) << "step: " << cp.step_ << " time: " << cp.t_ << " dt: " << dt
                 << std::endl;
      flecsi::flog::flush();
    } // if

    exec_cycle = exec_cycle && t <= cp.tf_ && dt != 0.0;

#endif

    ++cp.step_;
    return exec_cycle;
  } // cycle_control

  using evolve = cycle<cycle_control, point<cp::advance>, point<cp::analyze>>;

  using control_points =
    list<point<cp::initialize>, evolve, point<cp::finalize>>;

private:
  std::size_t step_{0};
  double t0_;
  double tf_;
  double t_;
  std::size_t max_steps_;
  double cfl_;
  double max_dt_;
  std::size_t log_frequency_;
  std::size_t output_frequency_;
  S<D> state_;
#ifdef HARD_BENCHMARK_MODE

  std::chrono::time_point<std::chrono::system_clock> start_timer_;
  double total_runtime_{0.0};
  std::string runtimes_{""};

#endif
}; // struct control_policy

} // namespace spec

#endif // SPEC_CONTROL_HH
