# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
from spack.spec import ConflictsInSpecError

class Hard(CMakePackage):

    """ A FleCSI-based radiation-hydrodynamics solver suite for the study of astrphysical phenomena """

    git = "https://github.com/lanl/hard"

    version("develop", branch="main")

    variant("catalyst", default=False, description="Enable catalyst for paraview interface")
    variant("radiation", default=True, description="Enable support for radiation physics")
    variant("tests", default=False, description="Enable unit tests")
    
    depends_on("flecsi@2.3.0:")
    depends_on("libcatalyst", when="+catalyst")
    depends_on("yaml-cpp@0.8.0")

    def cmake_args(self):
        options = [
            self.define_from_variant("ENABLE_UNIT_TESTS", "tests"),
            self.define_from_variant("ENABLE_CATALYST", "catalyst"),
        ]

        if(self.spec.satisfies("+radiation")):
            options.append("DISABLE_RADIATION=OFF");
        else:
            options.append("DISABLE_RADIATION=ON"); 
