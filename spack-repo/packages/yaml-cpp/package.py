from spack import *
from spack.pkg.builtin.yaml_cpp import YamlCpp

class YamlCpp(YamlCpp):
    """
    Additional named versions for Yaml-Cpp.
    """
    version("0.8.0", sha256="fbe74bbdcee21d656715688706da3c8becfd946d92cd44705cc6098bb23b3a16")
