from spack.package import *
from spack_repo.builtin.packages.flecsi.package import Flecsi

class Flecsi(Flecsi):
    """
    Additional named versions for FleCSI
    """
    # not a real version
    version("2.5", commit="1ce3933784e1337a44cb54b4065115a804a42856")
