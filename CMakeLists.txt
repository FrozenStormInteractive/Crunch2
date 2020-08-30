cmake_minimum_required(VERSION 3.1)
project(crunch)

include(GNUInstallDirs)

option(CRN_BUILD_SHARED_LIBS "Build crnlib as shared library." ${BUILD_SHARED_LIBS})
option(CRN_BUILD_EXAMPLES "Build examples." OFF)

set(BUILD_SHARED_LIBS ${CRN_BUILD_SHARED_LIBS})

add_subdirectory(3rdparty/lzma)
add_subdirectory(crnlib)
add_subdirectory(crunch)

add_subdirectory(3rdparty/jpeg-compressor)

if (CRN_BUILD_EXAMPLES)
	add_subdirectory(examples)
    add_subdirectory(3rdparty/stb)
endif(CRN_BUILD_EXAMPLES)
