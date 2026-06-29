#ifndef HARD_TASKS_IO_XDMF_HH
#define HARD_TASKS_IO_XDMF_HH

#include "common.hh"
#include "hdf5.hh"
#include "state.hh"
#include <../modules/spec/io.hh>
#include <format>
#include <fstream>
#include <sstream>

namespace hard::tasks::io {

template<std::size_t D, dm::domain L = dm::quantities>
void inline xdmf(flecsi::exec::cpu s,
  spec::io::name const & base,
  single<double>::accessor<ro> time,
  multi<typename mesh<D>::template accessor<ro>> mm,
  std::vector<std::tuple<multi<field<double>::accessor<ro, ro>>, std::string>>
    field_ma,
  std::vector<
    std::tuple<multi<typename field<vec<D>>::template accessor<ro, ro>>,
      std::string>> field_vector_ma) noexcept {

  const auto rank = s.launch().index;

  // Use MPI_COMM_WORLD for parallel HDF5 operations
  MPI_Comm comm = MPI_COMM_WORLD;

  for(uint32_t i{0}; i < mm.depth(); ++i) {
    const auto m = mm.accessors()[i];
    const auto [fields, fields_vectors, field_name, field_vector_name] =
      field_extractor<D>(m, i, field_ma, field_vector_ma);

    // Get mesh information
    std::array<double, 3> spacing = {0.0, 0.0, 0.0};
    std::array<double, 3> origin = {0.0, 0.0, 0.0};
    std::array<std::size_t, 3> local_size = {0, 1, 1};
    std::array<std::size_t, 3> global_size = {0, 1, 1};
    std::array<hsize_t, 3> offset = {0, 0, 0};

    // X direction (always present)
    spacing[0] = m.template delta<ax::x>();
    local_size[0] = m.template size<ax::x, L>();
    global_size[0] = m.template size<ax::x, dm::global>();
    const auto first_x = *m.template cells<ax::x, L>().begin();
    offset[0] = m.template global_id<ax::x>(first_x);
    origin[0] = m.template center<ax::x>(first_x) - 0.5 * spacing[0];

    if constexpr(D >= 2) {
      spacing[1] = m.template delta<ax::y>();
      local_size[1] = m.template size<ax::y, L>();
      global_size[1] = m.template size<ax::y, dm::global>();
      const auto first_y = *m.template cells<ax::y, L>().begin();
      offset[1] = m.template global_id<ax::y>(first_y);
      origin[1] = m.template center<ax::y>(first_y) - 0.5 * spacing[1];
    }

    if constexpr(D == 3) {
      spacing[2] = m.template delta<ax::z>();
      local_size[2] = m.template size<ax::z, L>();
      global_size[2] = m.template size<ax::z, dm::global>();
      const auto first_z = *m.template cells<ax::z, L>().begin();
      offset[2] = m.template global_id<ax::z>(first_z);
      origin[2] = m.template center<ax::z>(first_z) - 0.5 * spacing[2];
    }

    // Create HDF5 filename
    std::string h5_filename = std::format("output-{}D-{}.h5", D, base.str());

    // HDF5 writing scope - file will close when exiting this block
    {
      // Open HDF5 file collectively
      hdf5::file h5file(h5_filename, comm);

      // Create groups
      hdf5::group mesh_group(h5file, "/mesh");
      hdf5::group fields_group(h5file, "/fields");

      // Write mesh metadata (all ranks write the same data collectively)
      hdf5::write_1d_collective(
        mesh_group, "origin", {origin[0], origin[1], origin[2]});
      hdf5::write_1d_collective(
        mesh_group, "spacing", {spacing[0], spacing[1], spacing[2]});
      hdf5::write_1d_collective(mesh_group,
        "dimensions",
        {static_cast<double>(global_size[0]),
          static_cast<double>(global_size[1]),
          static_cast<double>(global_size[2])});

      // Write scalar fields
      for(std::size_t f = 0; f < fields.size(); ++f) {
        // Collect local data
        std::vector<double> local_data;
        local_data.reserve(local_size[0] * local_size[1] * local_size[2]);

        if constexpr(D == 1) {
          for(auto i : m.template cells<ax::x, L>()) {
            local_data.push_back(fields[f](i));
          }
        }
        else if constexpr(D == 2) {
          for(auto j : m.template cells<ax::y, L>()) {
            for(auto i : m.template cells<ax::x, L>()) {
              local_data.push_back(fields[f](i, j));
            }
          }
        }
        else /* D == 3 */ {
          for(auto k : m.template cells<ax::z, L>()) {
            for(auto j : m.template cells<ax::y, L>()) {
              for(auto i : m.template cells<ax::x, L>()) {
                local_data.push_back(fields[f](i, j, k));
              }
            }
          }
        }

        // Create dataset with global dimensions
        std::vector<hsize_t> global_dims = {
          global_size[2], global_size[1], global_size[0]};
        hdf5::dataspace space(global_dims);
        hdf5::dataset dset(
          fields_group, field_name[f], H5T_NATIVE_DOUBLE, space);

        // Write hyperslab
        std::array<hsize_t, 3> count = {
          local_size[2], local_size[1], local_size[0]};
        std::array<hsize_t, 3> file_offset = {offset[2], offset[1], offset[0]};
        hdf5::write_3d_hyperslab(dset, file_offset, count, local_data);
      }

      // Write vector fields
      for(std::size_t f = 0; f < fields_vectors.size(); ++f) {
        // Collect local data (interleaved components)
        std::vector<double> local_data;
        local_data.reserve(local_size[0] * local_size[1] * local_size[2] * D);

        if constexpr(D == 1) {
          for(auto i : m.template cells<ax::x, L>()) {
            local_data.push_back(fields_vectors[f](i)[0]);
          }
        }
        else if constexpr(D == 2) {
          for(auto j : m.template cells<ax::y, L>()) {
            for(auto i : m.template cells<ax::x, L>()) {
              local_data.push_back(fields_vectors[f](i, j)[0]);
              local_data.push_back(fields_vectors[f](i, j)[1]);
            }
          }
        }
        else /* D == 3 */ {
          for(auto k : m.template cells<ax::z, L>()) {
            for(auto j : m.template cells<ax::y, L>()) {
              for(auto i : m.template cells<ax::x, L>()) {
                local_data.push_back(fields_vectors[f](i, j, k)[0]);
                local_data.push_back(fields_vectors[f](i, j, k)[1]);
                local_data.push_back(fields_vectors[f](i, j, k)[2]);
              }
            }
          }
        }

        // Create dataset with global dimensions [nz, ny, nx, D]
        std::vector<hsize_t> global_dims = {
          global_size[2], global_size[1], global_size[0], D};
        hdf5::dataspace space(global_dims);
        hdf5::dataset dset(
          fields_group, field_vector_name[f], H5T_NATIVE_DOUBLE, space);

        // For vector fields, we need to write a 4D hyperslab
        // Select hyperslab in file
        hdf5::dataspace filespace(dset);
        std::array<hsize_t, 4> vec_offset = {
          offset[2], offset[1], offset[0], 0};
        std::array<hsize_t, 4> vec_count = {
          local_size[2], local_size[1], local_size[0], D};
        H5Sselect_hyperslab(filespace,
          H5S_SELECT_SET,
          vec_offset.data(),
          nullptr,
          vec_count.data(),
          nullptr);

        // Create memory dataspace
        std::vector<hsize_t> mem_dims = {
          local_size[2], local_size[1], local_size[0], D};
        hdf5::dataspace memspace(mem_dims);

        // Collective write
        hdf5::plist xfer_plist(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
        H5Dwrite(dset,
          H5T_NATIVE_DOUBLE,
          memspace,
          filespace,
          xfer_plist,
          local_data.data());
      }
    } // Close HDF5 scope - file and groups are now closed

    // Rank 0 writes XDMF metadata file after HDF5 file is closed
    if(rank == 0) {
      std::string xmf_filename =
        std::format("output-{}D-{}.xmf", D, base.str());
      std::ofstream xmf(xmf_filename);

      xmf << "<?xml version=\"1.0\"?>\n";
      xmf << "<Xdmf Version=\"3.0\">\n";
      xmf << "  <Domain>\n";
      xmf << "    <Grid Name=\"mesh\" GridType=\"Uniform\">\n";
      xmf << "      <Time Value=\"" << std::scientific << std::setprecision(16)
          << time << "\"/>\n";

      // Topology (dimensions are nz+1, ny+1, nx+1 for nodes)
      if constexpr(D == 1) {
        xmf << "      <Topology TopologyType=\"1DCoRectMesh\" Dimensions=\""
            << global_size[0] + 1 << "\"/>\n";
      }
      else if constexpr(D == 2) {
        xmf << "      <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\""
            << global_size[1] + 1 << " " << global_size[0] + 1 << "\"/>\n";
      }
      else {
        xmf << "      <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\""
            << global_size[2] + 1 << " " << global_size[1] + 1 << " "
            << global_size[0] + 1 << "\"/>\n";
      }

      // Geometry
      xmf << "      <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
      xmf << "        <DataItem Dimensions=\"3\" NumberType=\"Float\" "
             "Precision=\"8\" Format=\"HDF\">\n";
      xmf << "          " << h5_filename << ":/mesh/origin\n";
      xmf << "        </DataItem>\n";
      xmf << "        <DataItem Dimensions=\"3\" NumberType=\"Float\" "
             "Precision=\"8\" Format=\"HDF\">\n";
      xmf << "          " << h5_filename << ":/mesh/spacing\n";
      xmf << "        </DataItem>\n";
      xmf << "      </Geometry>\n";

      // Scalar fields
      for(const auto & fname : field_name) {
        xmf << "      <Attribute Name=\"" << fname
            << "\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
        xmf << "        <DataItem Dimensions=\"";
        if constexpr(D == 1) {
          xmf << global_size[0];
        }
        else if constexpr(D == 2) {
          xmf << global_size[1] << " " << global_size[0];
        }
        else {
          xmf << global_size[2] << " " << global_size[1] << " "
              << global_size[0];
        }
        xmf << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
        xmf << "          " << h5_filename << ":/fields/" << fname << "\n";
        xmf << "        </DataItem>\n";
        xmf << "      </Attribute>\n";
      }

      // Vector fields
      for(const auto & fname : field_vector_name) {
        xmf << "      <Attribute Name=\"" << fname
            << "\" AttributeType=\"Vector\" Center=\"Cell\">\n";
        xmf << "        <DataItem Dimensions=\"";
        if constexpr(D == 1) {
          xmf << global_size[0] << " " << D;
        }
        else if constexpr(D == 2) {
          xmf << global_size[1] << " " << global_size[0] << " " << D;
        }
        else {
          xmf << global_size[2] << " " << global_size[1] << " "
              << global_size[0] << " " << D;
        }
        xmf << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
        xmf << "          " << h5_filename << ":/fields/" << fname << "\n";
        xmf << "        </DataItem>\n";
        xmf << "      </Attribute>\n";
      }

      xmf << "    </Grid>\n";
      xmf << "  </Domain>\n";
      xmf << "</Xdmf>\n";
    }

  } // for mm.depth()
} // xdmf

} // namespace hard::tasks::io

#endif // HARD_TASKS_IO_XDMF_HH
