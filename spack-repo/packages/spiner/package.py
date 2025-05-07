# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class Spiner(CMakePackage):
    """Spiner:
    Performance portable routines for generic, tabulated, multi-dimensional data"""

    homepage = "https://github.com/lanl/spiner"
    url = "https://github.com/lanl/spiner/archive/refs/tags/1.4.0.tar.gz"
    git = "https://github.com/lanl/spiner.git"

    maintainers = ["rbberger"]

    version("1.6.3", commit="c379f5d2cce9625a60b149192ec43c8cff7b74b1")

    # When overriding/overloading varaints, the last variant is always used, except for
    # "when" clauses. Therefore, call the whens FIRST then the non-whens.
    # https://spack.readthedocs.io/en/latest/packaging_guide.html#overriding-variants
    variant(
        "kokkos",
        default=False,
        description="Enable kokkos",
    )

    variant("hdf5", default=False, description="Enable hdf5")
    variant("mpi", default=False, description="Support parallel hdf5")

    variant("python", default=False, description="Python, Numpy & Matplotlib Support")

    variant("test", default=False, description="Build tests")

    depends_on("cmake@3.12:", when="@:1.5.1")
    depends_on("cmake@3.19:", when="@1.6.0:")
    depends_on("catch2@3.0.1:", when="@1.6.3: +test")
    depends_on("catch2@2.13.4:2.13.9", when="@:1.6.2 +test")
    depends_on("ports-of-call@1.2.0:", when="@:1.5.1")
    depends_on("ports-of-call@1.5.1:", when="@1.6.0:")
    depends_on("ports-of-call@main", when="@main")

    # Currently the raw cuda backend of ports-of-call is not supported.
    depends_on("ports-of-call portability_strategy=Kokkos", when="@:1.5.1 +kokkos")
    depends_on("ports-of-call portability_strategy=None", when="@:1.5.1 ~kokkos")
    depends_on("kokkos@3.3.00:", when="+kokkos")
    requires("^kokkos+cuda_lambda+cuda_constexpr", when="+kokkos ^kokkos+cuda")

    depends_on("hdf5+hl~mpi", when="+hdf5~mpi")
    depends_on("hdf5+hl+mpi", when="+hdf5+mpi")

    depends_on("python", when="+python")
    depends_on("py-numpy", when="+python")
    depends_on("py-matplotlib", when="+python")

    conflicts("+mpi", when="~hdf5")

    def cmake_args(self):
        if self.spec.satisfies("@1.6.0:"):
            use_kokkos_option = "SPINER_TEST_USE_KOKKOS"
        else:
            use_kokkos_option = "SPINER_USE_KOKKOS"

        args = [
            self.define("BUILD_TESTING", self.run_tests),
            self.define("SPINER_BUILD_TESTS", self.run_tests),
            self.define("SPINER_TEST_USE_KOKKOS", self.run_tests and self.spec.satisfies("+kokkos")),
            self.define_from_variant(use_kokkos_option, "kokkos"),
            self.define_from_variant("SPINER_USE_HDF", "hdf5"),
        ]
        if self.spec.satisfies("^kokkos+cuda"):
            args.append(
                self.define("CMAKE_CUDA_ARCHITECTURES", self.spec["kokkos"].variants["cuda_arch"].value)
            )
        if self.spec.satisfies("^kokkos+rocm"):
            args.append(self.define("CMAKE_CXX_COMPILER", self.spec["hip"].hipcc))
            args.append(self.define("CMAKE_C_COMPILER", self.spec["hip"].hipcc))
        if self.spec.satisfies("^kokkos+cuda"):
            args.append(self.define("CMAKE_CXX_COMPILER", self.spec["kokkos"].kokkos_cxx))
        return args
