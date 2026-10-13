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

    version("main", branch="main")
    version("1.0.1", commit="f5e324da8bb876f71930eb1406e2cb0b132ddede", tag="v1.0.1")

    variant("tests", default=False, description="Enable unit tests")
    variant("standard", default=False, description="Standard setup for flecsolve")

    depends_on("flecsi@2.4:")
    depends_on("eigen")
    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build", when="+standard")

    def cmake_args(self):
        args = [
            self.define_from_variant("FLECSOLVE_ENABLE_UNIT_TESTS", "tests"),
            self.define_from_variant("FLECSOLVE_ENABLE_AMP", "standard"),
        ]
        return args
