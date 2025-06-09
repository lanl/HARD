from spack.package import *
from spack.pkg.builtin.flecsi import Flecsi

class Flecsi(Flecsi):
    """
    Additional named versions for FleCSI
    """
    version("2.4.0", commit="ff213702c291c303add0ed94eaf58c138af67f39")
