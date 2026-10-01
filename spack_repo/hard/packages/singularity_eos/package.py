from spack.package import *
from spack_repo.builtin.packages.singularity_eos.package import SingularityEos as BuiltinSingularityEos

class SingularityEos(BuiltinSingularityEos):
    """Extended singularity-eos package with additional versions."""

    version("1.12.1", commit="760ac3f8e106addc13dad8a47b9d4ad75e44ea48")
    version("1.11.1", commit="7365053a5bd59839ac47e6133426620540aca7e3")

