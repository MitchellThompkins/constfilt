from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, cmake_layout
from conan.tools.files import copy, get, rmdir
from conan.tools.build import check_min_cppstd
import os

required_conan_version = ">=2.1"


class ConstfiltConan(ConanFile):
    name = "constfilt"
    description = (
        "Header-only C++17 constexpr library for compile-time "
        "IIR digital filter design"
    )
    license = "Apache-2.0"
    url = "https://github.com/conan-io/conan-center-index"
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
    settings = "os", "arch", "compiler", "build_type"
    no_copy_source = True
    implements = ["auto_header_only"]

    def layout(self):
        cmake_layout(self, src_folder="src")

    def validate(self):
        check_min_cppstd(self, 17)

    def source(self):
        get(self, **self.conan_data["sources"][self.version], strip_root=True)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["CONSTFILT_VERSION"] = self.version
        tc.variables["CONSTFILT_BUILD_TESTS"] = False
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        copy(
            self,
            "LICENSE",
            src=self.source_folder,
            dst=os.path.join(self.package_folder, "licenses"),
        )
        cmake = CMake(self)
        cmake.install()
        rmdir(self, os.path.join(self.package_folder, "lib"))

    def package_info(self):
        self.cpp_info.bindirs = []
        self.cpp_info.libdirs = []
        self.cpp_info.set_property("cmake_file_name", "constfilt")
        self.cpp_info.set_property("cmake_target_name", "constfilt::constfilt")
