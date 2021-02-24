cmake_minimum_required(VERSION 3.0)
project(jpeg-compressor)

find_package(jpeg-compressor QUIET)

if(jpeg-compressor_FOUND)

    if(NOT TARGET jpeg-compressor)
        if(TARGET jpeg-compressor::jpeg-compressor)
            add_library(jpeg-compressor ALIAS jpeg-compressor::jpeg-compressor)
        else()
            add_library(jpeg-compressor INTERFACE)
            if(jpeg-compressor_INCLUDE_DIRS)
                target_include_directories(jpeg-compressor INTERFACE ${jpeg-compressor_INCLUDE_DIRS})
            endif()
            if(jpeg-compressor_LIBRARIES)
                target_link_libraries(jpeg-compressor INTERFACE ${jpeg-compressor_LIBRARIES})
            endif()
            if(jpeg-compressor_COMPILE_DEFINITIONS)
                target_compile_definitions(jpeg-compressor INTERFACE ${jpeg-compressor_COMPILE_DEFINITIONS})
            endif()
            if(jpeg-compressor_COMPILE_OPTIONS_LIST)
                target_compile_options(jpeg-compressor INTERFACE ${jpeg-compressor_COMPILE_OPTIONS_LIST})
            endif()
        endif()
    endif()

else()

    message(STATUS "jpeg-compressor not found")
    message(STATUS "Building jpeg-compressor from 3rdparty sources")

    set(JPGC_SRC_LIST jpgd.cpp jpge.cpp)
    set(JPGC_HDR_LIST jpgd.h jpgd_idct.h jpge.h)
    add_library(jpeg-compressor OBJECT EXCLUDE_FROM_ALL ${JPGC_SRC_LIST} ${JPGC_HDR_LIST})
    target_include_directories(jpeg-compressor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

    set_target_properties(jpeg-compressor PROPERTIES
        FOLDER "3rdparty"
        POSITION_INDEPENDENT_CODE ON
    )

endif()

