from spack.package import *
from spack_repo.builtin.packages.flecsi.package import Flecsi

class Flecsi(Flecsi):
    """
    Additional named versions for FleCSI
    """
    version("2.4.1", commit="f57a634e3f1f136e8932ad81f267d3b69657ae15")
