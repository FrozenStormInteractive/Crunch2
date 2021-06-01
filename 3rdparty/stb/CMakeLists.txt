cmake_minimum_required(VERSION 3.0)
project(stb)

find_package(stb QUIET)

if(stb_FOUND)

    if(NOT TARGET stb)
        add_library(stb INTERFACE)
        if(stb_INCLUDE_DIRS)
            target_include_directories(stb INTERFACE ${stb_INCLUDE_DIRS})
        endif()
        if(stb_LIBRARIES)
            target_link_libraries(stb INTERFACE ${stb_LIBRARIES})
        endif()
        if(stb_COMPILE_DEFINITIONS)
            target_compile_definitions(stb INTERFACE ${stb_COMPILE_DEFINITIONS})
        endif()
        if(stb_COMPILE_OPTIONS_LIST)
            target_compile_options(stb INTERFACE ${stb_COMPILE_OPTIONS_LIST})
        endif()
    endif()

else()

    message(STATUS "stb not found")
    message(STATUS "Building stb from 3rdparty sources")

    add_library(stb INTERFACE)
    target_include_directories(stb INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

endif()
