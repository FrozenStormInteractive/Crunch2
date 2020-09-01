import os
from conans import ConanFile, CMake, tools

class Crunch2Conan(ConanFile):
    name = "crunch2"
    description = "Advanced DXTc texture compression and transcoding library"
    homepage = "https://github.com/FrozenStormInteractive/Crunch2"
    url = "https://github.com/FrozenStormInteractive/Crunch2"
    license = "Zlib"
    topics = ("conan", "crunch", "texture", "compression", "decompression", "transcoding")
    settings = "os", "compiler", "arch", "build_type"
    exports_sources = ["CMakeLists.txt", "license.txt", "crnlib/*", "crunch/*", "inc/*", "3rdparty/*"]
    generators = "cmake"
    options = {
        "fPIC": [True, False],
        "shared": [True, False],
    }
    default_options = {
        "fPIC": True,
        "shared": False,
    }
    
    _cmake = None

    @property
    def _build_subfolder(self):
        return "build_subfolder"

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def configure(self):
        if self.options.shared:
            del self.options.fPIC

    def _configure_cmake(self):
        if self._cmake:
            return self._cmake
        self._cmake = CMake(self)
        self._cmake.definitions["CRN_BUILD_EXAMPLES"] = False
        self._cmake.definitions["CRN_BUILD_SHARED_LIBS"] = self.options.shared
        self._cmake.configure(build_folder=self._build_subfolder)
        return self._cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        self.copy("license.txt", src=self.source_folder, dst="licenses")
        cmake = self._configure_cmake()
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = tools.collect_libs(self)

        bindir = os.path.join(self.package_folder, "bin")
        self.output.info("Appending PATH environment variable: {}".format(bindir))
        self.env_info.PATH.append(bindir)
