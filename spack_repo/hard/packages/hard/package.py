# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack_repo.builtin.build_systems.cuda import CudaPackage
from spack.package import *

class Hard(CMakePackage, CudaPackage):
    """A FleCSI-based radiation-hydrodynamics solver suite
    for the study of astrophysical phenomena"""

    git = "https://github.com/lanl/hard"
    maintainers("JulienLoiseau")

    version("main", branch="main")

    variant("hdf5", default=True, description="Enable HDF5-based I/O (XDMF output)")
    variant("radiation", default=True, description="Enable support for radiation physics")
    variant("tests", default=False, description="Enable unit tests")
    variant("verification", default=False, description="Enable physics verification tests")
    variant("format", default=False, description="Enable format target")

    depends_on("flecsi@2.4: +flog")
    depends_on("flecsi +cuda", when="+cuda")
    depends_on("flecsolve")
    depends_on("flecsolve+cuda", when="+cuda")
    depends_on("hdf5@1.10: +mpi", when="+hdf5")
    depends_on("py-pybind11")
    depends_on("c", type="build")
    depends_on("cxx", type="build")

    depends_on("singularity-eos@1.11.0: ~closure~fortran+hdf5 +spiner build_extra=sesame")
    depends_on("singularity-eos~closure~fortran~eospac+kokkos+kokkos-kernels+cuda", when="+cuda")
    depends_on("ports-of-call@2.0.1:")

    depends_on("cmake@3.27:")
    depends_on("llvm@20:20", type="build", when="+format")
    depends_on("python", when="+tests")
    depends_on("py-numpy", when="+tests")
    depends_on("py-scipy", when="+tests")
    depends_on("py-exactpack", when="+tests")
    depends_on("py-matplotlib", when="+tests")

    requires("%clang@17:", when="+cuda", msg="CUDA version only supports Clang compiler")

    # Propagate cuda_arch requirement to dependencies
    for _flag in CudaPackage.cuda_arch_values:
        requires(f"+cuda cuda_arch={_flag}", when=f"^kokkos +cuda cuda_arch={_flag}")
        depends_on(f"kokkos cuda_arch={_flag}", when=f"+cuda cuda_arch={_flag}")
        depends_on(f"singularity-eos cuda_arch={_flag}", when=f"+cuda cuda_arch={_flag}")
        depends_on(f"flecsolve cuda_arch={_flag}", when=f"+cuda cuda_arch={_flag}")
        depends_on(f"flecsi cuda_arch={_flag}", when=f"+cuda cuda_arch={_flag}")

    def cmake_args(self):
        options = [
            self.define_from_variant("ENABLE_UNIT_TESTS", "tests"),
            self.define_from_variant("ENABLE_VERIFICATION", "verification"),
            self.define_from_variant("ENABLE_HDF5", "hdf5"),
            self.define_from_variant("ENABLE_FORMAT", "format"),
        ]

        return options
