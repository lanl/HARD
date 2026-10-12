from spack.package import *
from spack_repo.builtin.packages.flecsi.package import Flecsi

class Flecsi(Flecsi):
    """
    Additional named versions for FleCSI
    """
    version("2.develop", branch="2")
