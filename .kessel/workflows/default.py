from kessel.workflows import environment
from kessel.workflows.base.spack import BuildEnvironment
from kessel.workflows.base.cmake import CMake

class Default(BuildEnvironment, CMake):
    steps = ["env", "configure", "build", "test", "install"]

    project_spec = environment("hard+tests+verification")

    def ci_message(self, args):
        return super().ci_message(args, post_alloc_init="source .gitlab/kessel.sh")

    def build(self, args):
        """Build (with FLOG)"""
        cmake_args = [
            self.define("ENABLE_FLOG", True),
            self.define("ENABLE_UNIT_TESTS", True),
            self.define("ENABLE_DEVELOPER_WARNINGS", True)
        ]
        super().build(args, cmake_args)
