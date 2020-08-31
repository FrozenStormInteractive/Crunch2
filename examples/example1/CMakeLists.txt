set(EXAMPLE1_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/example1.cpp
)
add_executable(example1 ${EXAMPLE1_SRCS})
target_link_libraries(example1 crn stb)

set_target_properties(example1 PROPERTIES FOLDER "examples")
