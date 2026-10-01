from spack.package import *
from spack_repo.builtin.packages.flecsi.package import Flecsi

class Flecsi(Flecsi):
    """
    Additional named versions for FleCSI
    """

    version("2.develop", commit="1ce3933784e1337a44cb54b4065115a804a42856")
    version("2.4.2", commit="59545852cc377fe3c66f919e5910853c9bb0232a")
    version("2.4.1", commit="f57a634e3f1f136e8932ad81f267d3b69657ae15")
