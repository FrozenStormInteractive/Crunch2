cmake_minimum_required(VERSION 3.0)
project(jpeg-compressor)

set(BUILD_SHARED_LIBS OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(JPGD_SRC_LIST jpgd.cpp)
set(JPGD_HDR_LIST jpgd.h jpgd_idct.h)
add_library(jpgd ${JPGD_SRC_LIST} ${JPGD_HDR_LIST})
target_include_directories(jpgd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(JPGE_SRC_LIST jpge.cpp)
set(JPGE_HDR_LIST jpge.h jpge.h)
add_library(jpge EXCLUDE_FROM_ALL ${JPGE_SRC_LIST} ${JPGE_HDR_LIST})
target_include_directories(jpge PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(jpgd jpge PROPERTIES
    FOLDER "3rdparty"
    POSITION_INDEPENDENT_CODE ON
)
