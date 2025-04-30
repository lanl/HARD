# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *

class Hard(CMakePackage):
    """A FleCSI-based radiation-hydrodynamics solver suite
    for the study of astrophysical phenomena"""

    git = "https://github.com/lanl/hard"
    maintainers("JulienLoiseau")

    version("main", branch="main")

    variant("catalyst", default=False, description="Enable catalyst for paraview interface")
    variant("radiation", default=True, description="Enable support for radiation physics")
    variant("tests", default=False, description="Enable unit tests")
    variant("format", default=False, description="Enable format target")

    depends_on("flecsi@2.3.0 +flog")
    depends_on("libcatalyst", when="+catalyst")
    depends_on("yaml-cpp@0.8:")
    depends_on("singularity-eos@1.9.2 +hdf5 +spiner +eospac build_extra=sesame")
    depends_on("llvm@13.0.0", type="build", when="+format")
    depends_on("python", when="+tests")
    depends_on("py-numpy", when="+tests")
    depends_on("py-yamlreader", when="+tests")
    depends_on("py-scipy", when="+tests")
    depends_on("py-exactpack", when="+tests")
    depends_on("py-matplotlib", when="+tests")

    def cmake_args(self):
        options = [
            self.define_from_variant("ENABLE_UNIT_TESTS", "tests"),
            self.define_from_variant("ENABLE_CATALYST", "catalyst"),
            self.define_from_variant("ENABLE_RADIATION", "radiation"),
        ]

        return options
