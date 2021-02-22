cmake_minimum_required(VERSION 3.12)
project(Crunch2 VERSION "1.2.0")

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/conanbuildinfo.cmake")
    include("${CMAKE_CURRENT_SOURCE_DIR}/conanbuildinfo.cmake")
    conan_basic_setup()
elseif(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/conanbuildinfo.cmake")
    include("${CMAKE_CURRENT_BINARY_DIR}/conanbuildinfo.cmake")
    conan_basic_setup()
elseif(DEFINED CONANBUILD_DIR AND EXISTS "${CONANBUILD_DIR}/conanbuildinfo.cmake")
    include("${CONANBUILD_DIR}/conanbuildinfo.cmake")
    conan_basic_setup()
endif()

include(GNUInstallDirs)

option(CRN_BUILD_SHARED_LIBS "Build crnlib as shared library." ${BUILD_SHARED_LIBS})
option(CRN_BUILD_EXAMPLES "Build examples." OFF)

set(BUILD_SHARED_LIBS ${CRN_BUILD_SHARED_LIBS})

add_subdirectory(crnlib)
add_subdirectory(crunch)

add_subdirectory(3rdparty/lzma)
add_subdirectory(3rdparty/jpeg-compressor)
add_subdirectory(3rdparty/miniz)
add_subdirectory(3rdparty/stb)

if (CRN_BUILD_EXAMPLES)
    if(NOT WIN32)
        message(WARNING "Examples aren't supported this platform.")
    endif()
	add_subdirectory(examples)
endif(CRN_BUILD_EXAMPLES)
