#ifndef HARD_TASKS_IO_VTI_HH
#define HARD_TASKS_IO_VTI_HH

#include "common.hh"
#include "state.hh"
#include <../modules/spec/io.hh>
#include <format>
#include <fstream>
#include <sstream>

namespace hard::tasks::io {

template<std::size_t D, dm::domain L = dm::quantities>
void inline vtk(flecsi::exec::cpu s,
  spec::io::name const & base,
  single<double>::accessor<ro> time,
  multi<typename mesh<D>::template accessor<ro>> mm,
  std::vector<std::tuple<multi<field<double>::accessor<ro, ro>>, std::string>>
    field_ma,
  std::vector<
    std::tuple<multi<typename field<vec<D>>::template accessor<ro, ro>>,
      std::string>> field_vector_ma) noexcept {

  const auto rank = s.launch().index;
  const auto num_ranks = s.launch().size;

  for(uint32_t i{0}; i < mm.depth(); ++i) {
    const auto m = mm.accessors()[i];
    const auto [fields, fields_vectors, field_name, field_vector_name] =
      field_extractor<D>(m, i, field_ma, field_vector_ma);

    // Color indices (position in the decomposition)
    const auto color_idx = m.color_indeces();
    const auto axis_colors = m.axis_colors();

    // Get mesh information and compute extents
    std::array<std::size_t, 6> extent = {0, 0, 0, 1, 0, 1};
    std::array<std::size_t, 6> whole_extent = {0, 0, 0, 1, 0, 1};
    std::array<double, 3> spacing = {0.0, 0.0, 0.0};
    std::array<double, 3> origin = {0.0, 0.0, 0.0};
    std::array<std::size_t, 3> local_size = {0, 1, 1};
    std::array<std::size_t, 3> global_size = {0, 1, 1};

    // X direction (always present)
    spacing[0] = m.template delta<ax::x>();
    local_size[0] = m.template size<ax::x, L>();
    global_size[0] = m.template size<ax::x, dm::global>();
    const auto first_x = *m.template cells<ax::x, L>().begin();
    const auto x_min = m.template global_id<ax::x>(first_x);
    extent[0] = x_min;
    extent[1] = x_min + local_size[0];
    whole_extent[1] = global_size[0];
    origin[0] = m.template center<ax::x>(first_x) - 0.5 * spacing[0];

    if constexpr(D >= 2) {
      spacing[1] = m.template delta<ax::y>();
      local_size[1] = m.template size<ax::y, L>();
      global_size[1] = m.template size<ax::y, dm::global>();
      const auto first_y = *m.template cells<ax::y, L>().begin();
      const auto y_min = m.template global_id<ax::y>(first_y);
      extent[2] = y_min;
      extent[3] = y_min + local_size[1];
      whole_extent[3] = global_size[1];
      origin[1] = m.template center<ax::y>(first_y) - 0.5 * spacing[1];
    }

    if constexpr(D == 3) {
      spacing[2] = m.template delta<ax::z>();
      local_size[2] = m.template size<ax::z, L>();
      global_size[2] = m.template size<ax::z, dm::global>();
      const auto first_z = *m.template cells<ax::z, L>().begin();
      const auto z_min = m.template global_id<ax::z>(first_z);
      extent[4] = z_min;
      extent[5] = z_min + local_size[2];
      whole_extent[5] = global_size[2];
      origin[2] = m.template center<ax::z>(first_z) - 0.5 * spacing[2];
    }

    // Write piece file (.vti)
    std::string piece_filename = "output-" + std::to_string(D) + "D-" +
                                 std::to_string(rank) + "-" + base.str() +
                                 ".vti";
    std::ofstream piece_file(piece_filename);

    piece_file << "<?xml version=\"1.0\"?>\n";
    piece_file << "<VTKFile type=\"ImageData\" version=\"1.0\" "
                  "byte_order=\"LittleEndian\">\n";
    piece_file << "  <ImageData WholeExtent=\"" << whole_extent[0] << " "
               << whole_extent[1] << " " << whole_extent[2] << " "
               << whole_extent[3] << " " << whole_extent[4] << " "
               << whole_extent[5] << "\" Origin=\"" << std::scientific
               << std::setprecision(16) << origin[0] << " " << origin[1] << " "
               << origin[2] << "\" Spacing=\"" << spacing[0] << " "
               << spacing[1] << " " << spacing[2] << "\">\n";
    piece_file << "    <Piece Extent=\"" << extent[0] << " " << extent[1] << " "
               << extent[2] << " " << extent[3] << " " << extent[4] << " "
               << extent[5] << "\">\n";

    // FieldData for time
    piece_file << "      <FieldData>\n";
    piece_file << "        <DataArray type=\"Float64\" Name=\"TimeValue\" "
                  "NumberOfTuples=\"1\" format=\"ascii\">\n";
    piece_file << "          " << std::scientific << std::setprecision(16)
               << time << "\n";
    piece_file << "        </DataArray>\n";
    piece_file << "      </FieldData>\n";

    // CellData
    const std::size_t num_cells = local_size[0] * local_size[1] * local_size[2];
    piece_file << "      <CellData>\n";

    // Write scalar fields
    for(std::size_t f = 0; f < fields.size(); ++f) {
      piece_file << "        <DataArray type=\"Float64\" Name=\""
                 << field_name[f] << "\" format=\"ascii\">\n";
      piece_file << "          ";

      if constexpr(D == 1) {
        for(auto i : m.template cells<ax::x, L>()) {
          piece_file << std::scientific << std::setprecision(16) << fields[f](i)
                     << " ";
        }
      }
      else if constexpr(D == 2) {
        for(auto j : m.template cells<ax::y, L>()) {
          for(auto i : m.template cells<ax::x, L>()) {
            piece_file << std::scientific << std::setprecision(16)
                       << fields[f](i, j) << " ";
          }
        }
      }
      else /* D == 3 */ {
        for(auto k : m.template cells<ax::z, L>()) {
          for(auto j : m.template cells<ax::y, L>()) {
            for(auto i : m.template cells<ax::x, L>()) {
              piece_file << std::scientific << std::setprecision(16)
                         << fields[f](i, j, k) << " ";
            }
          }
        }
      }

      piece_file << "\n        </DataArray>\n";
    }

    // Write vector fields
    for(std::size_t f = 0; f < fields_vectors.size(); ++f) {
      piece_file << "        <DataArray type=\"Float64\" Name=\""
                 << field_vector_name[f]
                 << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
      piece_file << "          ";

      if constexpr(D == 1) {
        for(auto i : m.template cells<ax::x, L>()) {
          piece_file << std::scientific << std::setprecision(16)
                     << fields_vectors[f](i)[0] << " 0.0 0.0 ";
        }
      }
      else if constexpr(D == 2) {
        for(auto j : m.template cells<ax::y, L>()) {
          for(auto i : m.template cells<ax::x, L>()) {
            piece_file << std::scientific << std::setprecision(16)
                       << fields_vectors[f](i, j)[0] << " "
                       << fields_vectors[f](i, j)[1] << " 0.0 ";
          }
        }
      }
      else /* D == 3 */ {
        for(auto k : m.template cells<ax::z, L>()) {
          for(auto j : m.template cells<ax::y, L>()) {
            for(auto i : m.template cells<ax::x, L>()) {
              piece_file << std::scientific << std::setprecision(16)
                         << fields_vectors[f](i, j, k)[0] << " "
                         << fields_vectors[f](i, j, k)[1] << " "
                         << fields_vectors[f](i, j, k)[2] << " ";
            }
          }
        }
      }

      piece_file << "\n        </DataArray>\n";
    }

    piece_file << "      </CellData>\n";
    piece_file << "    </Piece>\n";
    piece_file << "  </ImageData>\n";
    piece_file << "</VTKFile>\n";
    piece_file.close();

    // Rank 0 writes the master .pvti file
    if(rank == 0) {
      std::string master_filename =
        "output-" + std::to_string(D) + "D-" + base.str() + ".pvti";
      std::ofstream master_file(master_filename);

      master_file << "<?xml version=\"1.0\"?>\n";
      master_file << "<VTKFile type=\"PImageData\" version=\"1.0\" "
                     "byte_order=\"LittleEndian\">\n";
      master_file << "  <PImageData WholeExtent=\"" << whole_extent[0] << " "
                  << whole_extent[1] << " " << whole_extent[2] << " "
                  << whole_extent[3] << " " << whole_extent[4] << " "
                  << whole_extent[5] << "\" Origin=\"" << std::scientific
                  << std::setprecision(16) << origin[0] << " " << origin[1]
                  << " " << origin[2] << "\" Spacing=\"" << spacing[0] << " "
                  << spacing[1] << " " << spacing[2] << "\">\n";

      // FieldData
      master_file << "    <PFieldData>\n";
      master_file << "      <PDataArray type=\"Float64\" Name=\"TimeValue\" "
                     "NumberOfTuples=\"1\"/>\n";
      master_file << "    </PFieldData>\n";

      // PCellData
      master_file << "    <PCellData>\n";
      for(const auto & name : field_name) {
        master_file << "      <PDataArray type=\"Float64\" Name=\"" << name
                    << "\"/>\n";
      }
      for(const auto & name : field_vector_name) {
        master_file << "      <PDataArray type=\"Float64\" Name=\"" << name
                    << "\" NumberOfComponents=\"3\"/>\n";
      }
      master_file << "    </PCellData>\n";

      // Piece references
      // Compute extent for each rank assuming uniform decomposition
      std::array<std::size_t, 3> per_rank = {
        global_size[0] / axis_colors[0], 1, 1};
      if constexpr(D >= 2)
        per_rank[1] = global_size[1] / axis_colors[1];
      if constexpr(D == 3)
        per_rank[2] = global_size[2] / axis_colors[2];

      for(std::size_t r = 0; r < num_ranks; ++r) {
        // Compute color indices for rank r
        std::size_t ix, iy, iz;
        if constexpr(D == 1) {
          ix = r;
          iy = 0;
          iz = 0;
        }
        else if constexpr(D == 2) {
          ix = r % axis_colors[0];
          iy = r / axis_colors[0];
          iz = 0;
        }
        else /* D == 3 */ {
          ix = r % axis_colors[0];
          iy = (r / axis_colors[0]) % axis_colors[1];
          iz = r / (axis_colors[0] * axis_colors[1]);
        }

        // Compute extent for this piece
        std::array<std::size_t, 6> piece_extent = {ix * per_rank[0],
          (ix + 1) * per_rank[0],
          iy * per_rank[1],
          (iy + 1) * per_rank[1],
          iz * per_rank[2],
          (iz + 1) * per_rank[2]};

        std::string piece_name = "output-" + std::to_string(D) + "D-" +
                                 std::to_string(r) + "-" + base.str() + ".vti";

        master_file << "    <Piece Extent=\"" << piece_extent[0] << " "
                    << piece_extent[1] << " " << piece_extent[2] << " "
                    << piece_extent[3] << " " << piece_extent[4] << " "
                    << piece_extent[5] << "\" Source=\"" << piece_name
                    << "\"/>\n";
      }

      master_file << "  </PImageData>\n";
      master_file << "</VTKFile>\n";
      master_file.close();
    }
  } // for mm.depth()
} // vtk

} // namespace hard::tasks::io

#endif // HARD_TASKS_IO_VTI_HH
