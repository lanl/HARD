# dependency package for singulary-opac

import os
from spack import *
from spack.package import *

class SingularityOpacDeps(CMakePackage):
    homepage    = "https://github.com/lanl/singularity-opac"
    git         = "git@github.com:lanl/singularity-opac.git"

    version("main", branch="main")

    #depends_on("cmake")
    #depends_on("hdf5")
    #depends_on("kokkos")

    #phases=["install"]
    
    def install(self, spec, prefix):
        install_tree(join_path(self.stage.source_path, "include"), prefix.include)

