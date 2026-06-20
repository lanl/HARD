#ifndef HARD_TASKS_IO_CSV_HH
#define HARD_TASKS_IO_CSV_HH

#include "common.hh"
#include "state.hh"
#include <../modules/spec/io.hh>
#include <fstream>
#include <iomanip>

namespace hard::tasks::io {

template<std::size_t D>
void inline csv(flecsi::exec::cpu s,
  spec::io::name const & base,
  single<double>::accessor<ro> time,
  multi<typename mesh<D>::template accessor<ro>> mm,
  std::vector<std::tuple<multi<field<double>::accessor<ro, ro>>,
    const std::string &>> field_ma,
  std::vector<
    std::tuple<typename vector_field_double<D>::template multi_accessor<ro, ro>,
      const std::string &>> field_vector_ma) {

  for(uint32_t i{0}; i < mm.depth(); ++i) {
    const auto m = mm.accessors()[i];
    const auto [fields, fields_vectors, field_name, field_vector_name] =
      field_extractor<D>(m, i, field_ma, field_vector_ma);

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

      for(auto idx : m.template cells<ax::x, dm::quantities>()) {
        file << std::setprecision(16) << std::scientific << time << "," << idx
             << "," << m.template center<ax::x>(idx) << std::setprecision(16);
        for(const auto & f_ma : fields)
          file << "," << f_ma(idx);
        for(const auto & f_ma_vec : fields_vectors)
          file << "," << f_ma_vec[0](idx);
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

      for(auto idx_i : m.template cells<ax::x, dm::quantities>()) {
        for(auto idx_j : m.template cells<ax::y, dm::quantities>()) {

          file << std::setprecision(16) << std::scientific << time << ","
               << idx_i << "," << idx_j << ","
               << m.template center<ax::x>(idx_i) << ","
               << m.template center<ax::y>(idx_j) << std::setprecision(16);
          for(const auto & f_ma : fields)
            file << "," << f_ma(idx_i, idx_j);
          for(const auto & f_ma_vec : fields_vectors)
            file << "," << f_ma_vec[0](idx_i, idx_j) << ","
                 << f_ma_vec[1](idx_i, idx_j);
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

      for(auto idx_i : m.template cells<ax::x, dm::quantities>()) {
        for(auto idx_j : m.template cells<ax::y, dm::quantities>()) {
          for(auto idx_k : m.template cells<ax::z, dm::quantities>()) {

            file << std::setprecision(16) << std::scientific << time << ","
                 << idx_i << "," << idx_j << "," << idx_k << ","
                 << m.template center<ax::x>(idx_i) << ","
                 << m.template center<ax::y>(idx_j) << ","
                 << m.template center<ax::z>(idx_k) << std::setprecision(16);
            for(const auto & f_ma : fields)
              file << "," << f_ma(idx_i, idx_j, idx_k);
            for(const auto & f_ma_vec : fields_vectors)
              file << "," << f_ma_vec[0](idx_i, idx_j, idx_k) << ","
                   << f_ma_vec[1](idx_i, idx_j, idx_k) << ","
                   << f_ma_vec[2](idx_i, idx_j, idx_k);
            file << std::endl;
          }
        }
      }
    } // if
  } // for
} // csv

} // namespace hard::tasks::io

#endif // HARD_TASKS_IO_CSV_HH
