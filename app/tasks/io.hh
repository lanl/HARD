#ifndef HARD_TASKS_IO_HH
#define HARD_TASKS_IO_HH

#include "../state.hh"

#include <fstream>
#include <iomanip>
#include <spec/io.hh>

namespace hard::tasks::io {

template<std::size_t D>
void inline raw(spec::io::name const & base,
  single<double>::accessor<ro> time,
  multi<typename mesh<D>::template accessor<ro>> mm,
  multi<field<double>::accessor<ro, ro>> r_ma,
  multi<field<double>::accessor<ro, ro>> p_ma,
  multi<field<double>::accessor<ro, ro>> c_ma,
  multi<typename field<vec<D>>::template accessor<ro, ro>> u_ma,
  multi<typename field<vec<D>>::template accessor<ro, ro>> ru_ma,
  multi<field<double>::accessor<ro, ro>> rE_ma,
  multi<field<double>::accessor<ro, ro>> radE_ma) {

  for(uint32_t i{0}; i < mm.depth(); ++i) {
    const auto m = mm.accessors()[i];

    auto r_a = r_ma.accessors()[i];
    auto p_a = p_ma.accessors()[i];
    auto c_a = c_ma.accessors()[i];
    auto u_a = u_ma.accessors()[i];
    auto ru_a = ru_ma.accessors()[i];
    auto rE_a = rE_ma.accessors()[i];
    auto radE_a = radE_ma.accessors()[i];

    auto r = m.template mdcolex<is::cells>(r_a);
    auto p = m.template mdcolex<is::cells>(p_a);
    auto c = m.template mdcolex<is::cells>(c_a);
    auto u = m.template mdcolex<is::cells>(u_a);
    auto ru = m.template mdcolex<is::cells>(ru_a);
    auto rE = m.template mdcolex<is::cells>(rE_a);
    auto radE = m.template mdcolex<is::cells>(radE_a);

    if constexpr(D == 1) {
      std::ofstream file(
        base.str() + "-1D-" + std::to_string(flecsi::process()) + ".raw");

      file << "#" << m.template size<ax::x, dm::quantities>() << std::endl;
      file << "#" << m.template size<ax::x, dm::global>() << std::endl;

      {
        auto ccoords = m.color_indeces();
        file << "#" << ccoords[ax::x] << std::endl;
        auto ccolors = m.axis_colors();
        file << "#" << ccolors[ax::x] << std::endl;
      } // scope

      // File format, columns:
      // time, cellid_x, coord_x, density, pressure, velocity,
      // fluid total energy density, radiation energy density
      file << "#time\tcellid_x\tcoord_"
              "x\tdensity\tpressure\tvelocity\ttotalE\tRadE\tsoundspeed"
           << std::endl;
      file << std::fixed; // just to get a uniform file format
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        file << std::setprecision(6) << std::scientific << time << "\t" << i
             << "\t" << m.template center<ax::x>(i) << std::setprecision(12)
             << "\t" << r(i) << "\t" << p(i) << "\t" << u(i).x << "\t" << rE(i)
             << "\t" << radE(i) << "\t" << c(i) << std::endl;
      } // for
      file << std::endl; // empty line needed for gnuplot's splot
    }
    else if constexpr(D == 2) {
      std::ofstream file(
        base.str() + "-2D-" + std::to_string(flecsi::process()) + ".raw");

      file << "#" << m.template size<ax::x, dm::quantities>() << " "
           << m.template size<ax::y, dm::quantities>() << std::endl;
      file << "#" << m.template size<ax::x, dm::global>() << " "
           << m.template size<ax::y, dm::global>() << std::endl;

      {
        auto ccoords = m.color_indeces();
        file << "#" << ccoords[ax::x] << " " << ccoords[ax::y] << std::endl;
        auto ccolors = m.axis_colors();
        file << "#" << ccolors[ax::x] << " " << ccolors[ax::y] << std::endl;
      } // scope

      // File format, columns:
      // time, cellid_x, cellid_y, coord_x, coord_y, density, pressure,
      // velocity, fluid total energy density, radiation energy density
      file << "#time\tcellidx\tcellidy\tx\ty\tdensity\tpressure\tvx\tvy\ttotalE"
              "\tRadE"
           << std::endl;
      file << std::fixed; // just to get a uniform file format
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          file << std::setprecision(6) << std::scientific << time << "\t" << i
               << "\t" << j << "\t" << m.template center<ax::x>(i) << "\t"
               << m.template center<ax::y>(j) << std::setprecision(9) << "\t"
               << r(i, j) << "\t" << p(i, j) << "\t" << u(i, j).x << "\t"
               << u(i, j).y << "\t" << rE(i, j) << "\t" << radE(i, j)
               << std::endl;
        } // for
        file << std::endl; // empty line needed for gnuplot's splot
      } // for
    }
    else /* D == 3 */ {
      std::ofstream file(
        base.str() + "-" + std::to_string(flecsi::process()) + ".raw");

      file << m.template size<ax::x, dm::quantities>() << " "
           << m.template size<ax::y, dm::quantities>() << " "
           << m.template size<ax::z, dm::quantities>() << std::endl;
      file << m.template size<ax::x, dm::global>() << " "
           << m.template size<ax::y, dm::global>() << " "
           << m.template size<ax::z, dm::global>() << std::endl;

      {
        auto ccoords = m.color_indeces();
        file << ccoords[ax::x] << " " << ccoords[ax::y] << " " << ccoords[ax::z]
             << std::endl;
        auto ccolors = m.axis_colors();
        file << ccolors[ax::x] << " " << ccolors[ax::y] << " " << ccolors[ax::z]
             << std::endl;
      } // scope

      // Density
      for(auto k : m.template cells<ax::z, dm::quantities>()) {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            file << r(i, j, k) << std::endl;
          } // for
        } // for
      } // for

      // Momentum
      for(auto k : m.template cells<ax::z, dm::quantities>()) {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            file << ru(i, j, k).x << " " << ru(i, j, k).y << " "
                 << ru(i, j, k).z << std::endl;
          } // for
        } // for
      } // for

      // Total Energy
      for(auto k : m.template cells<ax::z, dm::quantities>()) {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto i : m.template cells<ax::x, dm::quantities>()) {
            file << rE(i, j, k) << std::endl;
          } // for
        } // for
      } // for
    } // for
  } // if
} // raw

} // namespace hard::tasks::io

#endif // HARD_TASKS_IO_HH
