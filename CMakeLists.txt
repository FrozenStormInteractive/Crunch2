cmake_minimum_required(VERSION 3.1)
project(crunch)

option(BUILD_SHARED_LIBS "Build crunchlib as shared library and link executables against it." OFF)
option(CRN_BUILD_EXAMPLES "Build examples." OFF)

include(GNUInstallDirs)

add_subdirectory(crnlib)
add_subdirectory(crunch)

if (CRN_BUILD_EXAMPLES)
	add_subdirectory(examples)
endif(CRN_BUILD_EXAMPLES)
