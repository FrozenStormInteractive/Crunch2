cmake_minimum_required(VERSION 3.0)
project(miniz)

set(BUILD_SHARED_LIBS OFF)

add_library(miniz EXCLUDE_FROM_ALL miniz.c miniz.h)
target_include_directories(miniz PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(miniz PROPERTIES FOLDER "3rdparty")
