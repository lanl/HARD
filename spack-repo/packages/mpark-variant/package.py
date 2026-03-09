# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
from spack_repo.builtin.packages.mpark_variant.package import MparkVariant

class MparkVariant(MparkVariant):
    patch("gpu_compatibility.patch")
