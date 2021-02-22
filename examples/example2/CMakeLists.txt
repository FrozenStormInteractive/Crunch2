set(EXAMPLE2_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/example2.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/timer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/timer.h
)

add_executable(example2 ${EXAMPLE2_SRCS})
target_link_libraries(example2 crn)
set_target_properties(example2 PROPERTIES CXX_STANDARD 11 FOLDER "examples")
