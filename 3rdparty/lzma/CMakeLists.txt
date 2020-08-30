cmake_minimum_required(VERSION 3.1)
project(lzma_sdk C)

file(GLOB lzma_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
file(GLOB lzma_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
if(MSVC)
    list(APPEND lzma_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/7zVersion.rc)
endif()

add_library(lzma STATIC ${lzma_SOURCES} ${lzma_HEADERS})
target_include_directories(lzma PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
