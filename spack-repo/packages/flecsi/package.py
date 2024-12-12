from spack.package import *
from spack.pkg.builtin.flecsi import Flecsi

class Flecsi(Flecsi):
    """
    Additional named versions for FleCSI.
    """
    version("2.3-cdss", commit="5b814a9efbe9022defd0a631a0b903f265e98d12")
