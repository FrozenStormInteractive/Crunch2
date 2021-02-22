cmake_minimum_required(VERSION 3.0)
project(miniz)

find_package(miniz QUIET)

if(miniz_FOUND)

    if(NOT TARGET miniz)
        if(TARGET miniz::miniz)
            add_library(miniz ALIAS miniz::miniz)
        else()
            add_library(miniz INTERFACE)
            if(miniz_INCLUDE_DIRS)
                target_include_directories(miniz INTERFACE ${miniz_INCLUDE_DIRS})
            endif()
            if(miniz_LIBRARIES)
                target_link_libraries(miniz INTERFACE ${miniz_LIBRARIES})
            endif()
            if(miniz_COMPILE_DEFINITIONS)
                target_compile_definitions(miniz INTERFACE ${miniz_COMPILE_DEFINITIONS})
            endif()
            if(miniz_COMPILE_OPTIONS_LIST)
                target_compile_options(miniz INTERFACE ${miniz_COMPILE_OPTIONS_LIST})
            endif()
        endif()
    endif()

else()

    message(STATUS "Miniz not found")
    message(STATUS "Building miniz from 3rdparty sources")

    add_library(miniz STATIC EXCLUDE_FROM_ALL miniz.c miniz.h)
    target_include_directories(miniz PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

    set_target_properties(miniz PROPERTIES
        FOLDER "3rdparty"
        POSITION_INDEPENDENT_CODE ON
    )

endif()
