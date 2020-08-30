set(EXAMPLE1_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/example1.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/stb_image.h
)
add_executable(example1 ${EXAMPLE1_SRCS})
target_link_libraries(example1 crn)
