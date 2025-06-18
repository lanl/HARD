#ifndef HARD_TASKS_IO_HH
#define HARD_TASKS_IO_HH

#include "../state.hh"

#include <fstream>
#include <iomanip>
#include <spec/io.hh>

namespace hard::tasks::io {

template<std::size_t D>
void inline csv(flecsi::exec::cpu s,
  spec::io::name const & base,
  single<double>::accessor<ro> time,
  multi<typename mesh<D>::template accessor<ro>> mm,
  std::vector<multi<field<double>::accessor<ro, ro>>> field_ma,
  std::vector<multi<typename field<vec<D>>::template accessor<ro, ro>>>
    field_vector_ma,
  std::vector<std::string> field_name,
  std::vector<std::string> field_vector_name) {

  for(uint32_t i{0}; i < mm.depth(); ++i) {
    const auto m = mm.accessors()[i];
    std::vector<flecsi::util::mdcolex<const double, D>> fields;
    std::vector<flecsi::util::mdcolex<const vec<D>, D>> fields_vectors;

    std::transform(field_ma.begin(),
      field_ma.end(),
      std::back_inserter(fields),
      [&m, &i](const multi<field<double>::accessor<ro, ro>> & ma)
        -> flecsi::util::mdcolex<const double, D> {
        return m.template mdcolex<is::cells>(ma.accessors()[i]);
      });
    std::transform(field_vector_ma.begin(),
      field_vector_ma.end(),
      std::back_inserter(fields_vectors),
      [&m, &i](const multi<typename field<vec<D>>::accessor<ro, ro>> & ma)
        -> flecsi::util::mdcolex<const vec<D>, D> {
        return m.template mdcolex<is::cells>(ma.accessors()[i]);
      });

    if constexpr(D == 1) {
      std::ofstream file("output-1D-" + std::to_string(s.launch().index) + "-" +
                         base.str() + ".csv");

      file << "time,cellid_x,coord_x";
      for(auto n : field_name)
        file << "," << n;
      for(auto n : field_vector_name)
        file << "," << n << "_x";
      file << std::endl;
      file << std::fixed; // just to get a uniform file format

      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        file << std::setprecision(6) << std::scientific << time << "," << i
             << "," << m.template center<ax::x>(i) << std::setprecision(12);
        for(auto f_ma : fields)
          file << "," << f_ma(i);
        for(auto f_ma : fields_vectors)
          file << "," << f_ma(i)[0];
        file << std::endl;
      }
    }
    else if constexpr(D == 2) {

      std::ofstream file("output-2D-" + std::to_string(s.launch().index) + "-" +
                         base.str() + ".csv");

      file << "time,cellid_x,cellid_y,coord_x,coord_y";
      for(auto n : field_name)
        file << "," << n;
      for(auto n : field_vector_name)
        file << "," << n << "_x," << n << "_y";
      file << std::endl;
      file << std::fixed; // just to get a uniform file format

      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {

          file << std::setprecision(6) << std::scientific << time << "," << i
               << "," << j << "," << m.template center<ax::x>(i) << ","
               << m.template center<ax::y>(j) << std::setprecision(12);
          for(auto f_ma : fields)
            file << "," << f_ma(i, j);
          for(auto f_ma : fields_vectors)
            file << "," << f_ma(i, j)[0] << "," << f_ma(i, j)[1];
          file << std::endl;
        }
      }
    }
    else /* D == 3 */ {

      std::ofstream file("output-3D-" + std::to_string(s.launch().index) + "-" +
                         base.str() + ".csv");

      file << "time,cellid_x,cellid_y,cellid_z,coord_x,coord_y,coord_z";
      for(auto n : field_name)
        file << "," << n;
      for(auto n : field_vector_name)
        file << "," << n << "_x," << n << "_y," << n << "_z";
      file << std::endl;
      file << std::fixed; // just to get a uniform file format

      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        for(auto j : m.template cells<ax::y, dm::quantities>()) {
          for(auto k : m.template cells<ax::z, dm::quantities>()) {

            file << std::setprecision(6) << std::scientific << time << "," << i
                 << "," << j << "," << k << "," << m.template center<ax::x>(i)
                 << "," << m.template center<ax::y>(j) << ","
                 << m.template center<ax::z>(k) << std::setprecision(12);
            for(auto f_ma : fields)
              file << "," << f_ma(i, j, k);
            for(auto f_ma : fields_vectors)
              file << "," << f_ma(i, j, k)[0] << "," << f_ma(i, j, k)[1] << ","
                   << f_ma(i, j, k)[2];
            file << std::endl;
          }
        }
      }
    } // if
  } // for
} // csv

} // namespace hard::tasks::io

#endif // HARD_TASKS_IO_HH
