
# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class PyExactpack(PythonPackage):
    """ExactPack Python package"""

    git = "https://github.com/lanl/exactpack"

    version("1.7.10", commit="dc331d9ac450c0ebc94ada908a275ff26949e0e5")

    depends_on("python@3.6:", type=("build", "run"))
    depends_on("py-setuptools", type="build")
