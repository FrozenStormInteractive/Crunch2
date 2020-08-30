set(CRUNCH_SRCS
	${CMAKE_CURRENT_SOURCE_DIR}/corpus_gen.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/corpus_gen.h
	${CMAKE_CURRENT_SOURCE_DIR}/corpus_test.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/corpus_test.h
	${CMAKE_CURRENT_SOURCE_DIR}/crunch.cpp
)

add_executable(crunch ${CRUNCH_SRCS})
target_link_libraries(crunch crn)

install(TARGETS crunch)
