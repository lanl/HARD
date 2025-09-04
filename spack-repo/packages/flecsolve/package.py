# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *

class Flecsolve(CMakePackage):
    """Solvers package built on top of FleCSI"""

    homepage="https://github.com/lanl/flecsolve.git"
    git = "https://github.com/lanl/flecsolve.git"
    # TO DO update to stable when available
    version("main", branch="flecsi-2.4")

    variant("tests", default=False, description="Enable unit tests")
    variant("standard", default=False, description="Standard setup for flecsolve")

    depends_on('flecsi@2.4:')
    #depends_on('amp+hypre', when="+standard") # Might want to enable at some point.
    #depends_on('stacktrace+shared', when="+standard")
    #depends_on('lapackwrappers@main', when="+standard")

    def cmake_args(self):
        args = [
            self.define_from_variant("FLECSOLVE_ENABLE_UNIT_TESTS", "tests"),
            self.define_from_variant("FLECSOLVE_ENABLE_AMP", "standard"),
        ]
        return args


