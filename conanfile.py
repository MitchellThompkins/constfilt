from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, cmake_layout
from conan.tools.files import copy
import os


class ConstfiltConan(ConanFile):
    name = "constfilt"
    description = (
        "Header-only C++17 constexpr library for compile-time "
        "IIR digital filter design"
    )
    license = "Apache-2.0"
    url = "https://github.com/MitchellThompkins/constfilt"
    homepage = "https://github.com/MitchellThompkins/constfilt"
    topics = (
        "constexpr",
        "iir",
        "filter",
        "butterworth",
        "dsp",
        "header-only",
        "compile-time",
    )
    package_type = "header-library"
    no_copy_source = True

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["CONSTFILT_VERSION"] = self.version
        tc.variables["CONSTFILT_BUILD_TESTS"] = False
        tc.generate()

    def package(self):
        copy(
            self,
            "LICENSE",
            src=self.source_folder,
            dst=os.path.join(self.package_folder, "licenses"),
        )
        cmake = CMake(self)
        cmake.configure()
        cmake.install()

    def package_info(self):
        self.cpp_info.bindirs = []
        self.cpp_info.libdirs = []
        self.cpp_info.set_property("cmake_file_name", "constfilt")
        self.cpp_info.set_property("cmake_target_name", "constfilt::constfilt")
