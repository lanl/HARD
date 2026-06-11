#ifndef HARD_SPEC_TASKS_HDF5_HH
#define HARD_SPEC_TASKS_HDF5_HH

#include <array>
#include <cstddef>
#include <hdf5.h>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace hard::tasks::io::hdf5 {

// Exception type for HDF5 errors
struct hdf5_error : std::runtime_error {
  using runtime_error::runtime_error;
};

// Check HDF5 return value and throw on error
inline void
check_h5(hid_t ret, const char * msg) {
  if(ret < 0) {
    throw hdf5_error(msg);
  }
}

// RAII wrapper for HDF5 file
struct file {
  file() = default;
  file(const std::string & filename, MPI_Comm comm) {
    // Create file access property list for parallel I/O
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    check_h5(fapl_id, "Failed to create file access property list");

    herr_t ret = H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL);
    check_h5(ret, "Failed to set MPI-IO file access");

    // Create file
    id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    check_h5(id, "Failed to create HDF5 file");

    H5Pclose(fapl_id);
  }

  ~file() {
    if(id >= 0) {
      H5Fclose(id);
    }
  }

  file(const file &) = delete;
  file & operator=(const file &) = delete;
  file(file && other) noexcept : id(other.id) {
    other.id = -1;
  }
  file & operator=(file && other) noexcept {
    if(this != &other) {
      if(id >= 0)
        H5Fclose(id);
      id = other.id;
      other.id = -1;
    }
    return *this;
  }

  operator hid_t() const {
    return id;
  }

private:
  hid_t id = -1;
};

// RAII wrapper for HDF5 group
struct group {
  group() = default;
  group(hid_t loc_id, const std::string & name) {
    // Check if group exists
    if(H5Lexists(loc_id, name.c_str(), H5P_DEFAULT) > 0) {
      id = H5Gopen2(loc_id, name.c_str(), H5P_DEFAULT);
    }
    else {
      // Create group with intermediate groups
      hid_t gcpl_id = H5Pcreate(H5P_LINK_CREATE);
      H5Pset_create_intermediate_group(gcpl_id, 1);
      id = H5Gcreate2(loc_id, name.c_str(), gcpl_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Pclose(gcpl_id);
    }
    check_h5(id, "Failed to create/open group");
  }

  ~group() {
    if(id >= 0) {
      H5Gclose(id);
    }
  }

  group(const group &) = delete;
  group & operator=(const group &) = delete;
  group(group && other) noexcept : id(other.id) {
    other.id = -1;
  }
  group & operator=(group && other) noexcept {
    if(this != &other) {
      if(id >= 0)
        H5Gclose(id);
      id = other.id;
      other.id = -1;
    }
    return *this;
  }

  operator hid_t() const {
    return id;
  }

private:
  hid_t id = -1;
};

// RAII wrapper for HDF5 dataset
struct dataset {
  dataset() = default;
  dataset(hid_t loc_id,
    const std::string & name,
    hid_t type_id,
    hid_t space_id) {
    id = H5Dcreate2(loc_id,
      name.c_str(),
      type_id,
      space_id,
      H5P_DEFAULT,
      H5P_DEFAULT,
      H5P_DEFAULT);
    check_h5(id, "Failed to create dataset");
  }

  ~dataset() {
    if(id >= 0) {
      H5Dclose(id);
    }
  }

  dataset(const dataset &) = delete;
  dataset & operator=(const dataset &) = delete;
  dataset(dataset && other) noexcept : id(other.id) {
    other.id = -1;
  }
  dataset & operator=(dataset && other) noexcept {
    if(this != &other) {
      if(id >= 0)
        H5Dclose(id);
      id = other.id;
      other.id = -1;
    }
    return *this;
  }

  operator hid_t() const {
    return id;
  }

  hid_t get_space() const {
    return H5Dget_space(id);
  }

private:
  hid_t id = -1;
};

// RAII wrapper for HDF5 dataspace
struct dataspace {
  dataspace() = default;
  explicit dataspace(const std::vector<hsize_t> & dims) {
    id = H5Screate_simple(dims.size(), dims.data(), nullptr);
    check_h5(id, "Failed to create dataspace");
  }

  explicit dataspace(hid_t dset_id) {
    id = H5Dget_space(dset_id);
    check_h5(id, "Failed to get dataspace");
  }

  ~dataspace() {
    if(id >= 0) {
      H5Sclose(id);
    }
  }

  dataspace(const dataspace &) = delete;
  dataspace & operator=(const dataspace &) = delete;
  dataspace(dataspace && other) noexcept : id(other.id) {
    other.id = -1;
  }
  dataspace & operator=(dataspace && other) noexcept {
    if(this != &other) {
      if(id >= 0)
        H5Sclose(id);
      id = other.id;
      other.id = -1;
    }
    return *this;
  }

  operator hid_t() const {
    return id;
  }

private:
  hid_t id = -1;
};

// RAII wrapper for HDF5 property list
struct plist {
  plist() = default;
  explicit plist(hid_t cls_id) {
    id = H5Pcreate(cls_id);
    check_h5(id, "Failed to create property list");
  }

  ~plist() {
    if(id >= 0) {
      H5Pclose(id);
    }
  }

  plist(const plist &) = delete;
  plist & operator=(const plist &) = delete;
  plist(plist && other) noexcept : id(other.id) {
    other.id = -1;
  }
  plist & operator=(plist && other) noexcept {
    if(this != &other) {
      if(id >= 0)
        H5Pclose(id);
      id = other.id;
      other.id = -1;
    }
    return *this;
  }

  operator hid_t() const {
    return id;
  }

private:
  hid_t id = -1;
};

// Helper to write a 1D array collectively
inline void
write_1d_collective(hid_t loc_id,
  const std::string & name,
  const std::vector<double> & data) {
  // Create dataspace
  std::vector<hsize_t> dims = {data.size()};
  dataspace space(dims);

  // Create dataset
  dataset dset(loc_id, name, H5T_NATIVE_DOUBLE, space);

  // Create transfer property list for collective I/O
  plist xfer_plist(H5P_DATASET_XFER);
  herr_t ret = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
  check_h5(ret, "Failed to set collective transfer");

  // Write data
  ret = H5Dwrite(
    dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, xfer_plist, data.data());
  check_h5(ret, "Failed to write 1D data");
}

// Helper to write a 3D hyperslab collectively
inline void
write_3d_hyperslab(hid_t dset_id,
  const std::array<hsize_t, 3> & offset,
  const std::array<hsize_t, 3> & count,
  const std::vector<double> & data) {
  // Get file dataspace and select hyperslab
  dataspace filespace(dset_id);
  herr_t ret = H5Sselect_hyperslab(
    filespace, H5S_SELECT_SET, offset.data(), nullptr, count.data(), nullptr);
  check_h5(ret, "Failed to select hyperslab in file");

  // Create memory dataspace
  std::vector<hsize_t> mem_dims = {count[0], count[1], count[2]};
  dataspace memspace(mem_dims);

  // Create transfer property list for collective I/O
  plist xfer_plist(H5P_DATASET_XFER);
  ret = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
  check_h5(ret, "Failed to set collective transfer");

  // Write data
  ret = H5Dwrite(
    dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, xfer_plist, data.data());
  check_h5(ret, "Failed to write hyperslab");
}

} // namespace hard::tasks::io::hdf5

#endif // HARD_SPEC_TASKS_HDF5_HH
