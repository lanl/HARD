# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack_repo.builtin.build_systems.cuda import CudaPackage
from spack_repo.builtin.build_systems.rocm import ROCmPackage
from spack.package import *

class Flecsolve(CMakePackage, CudaPackage, ROCmPackage):
    """Solvers package built on top of FleCSI"""

    homepage = "https://github.com/lanl/flecsolve.git"
    git = "https://github.com/lanl/flecsolve.git"

    version("main", commit="661877df93b594e2447790fc313cca4446837fad")

    variant("tests", default=False, description="Enable unit tests")
    variant("standard", default=False, description="Standard setup for flecsolve")

    depends_on("flecsi@2.4:")
    depends_on("eigen")
    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build")

    depends_on("cuda", when="+cuda")

    def cmake_args(self):
        args = [
            self.define_from_variant("FLECSOLVE_ENABLE_UNIT_TESTS", "tests"),
            self.define_from_variant("FLECSOLVE_ENABLE_AMP", "standard"),
        ]

        if "+cuda" in self.spec:
            cuda_prefix = self.spec["cuda"].prefix
            args.append(self.define("CMAKE_CXX_FLAGS", f"--cuda-path={cuda_prefix}"))

        return args
