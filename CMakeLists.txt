cmake_minimum_required(VERSION 3.7)
project(arpes)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -Werror=return-type -Wno-deprecated-declarations -Wno-psabi -m64 -march=native -fmax-errors=2")
set(CMAKE_CXX_FLAGS_DEBUG "-Wextra -O0 -ggdb3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

file(GLOB_RECURSE CPP_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_executable(arpes "${CPP_FILES}")

if (NOT TARGET es_util)
	add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/es_util")
endif()
target_link_libraries(arpes es_util)

if (NOT TARGET es_la)
	add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/es_la")
endif()
target_link_libraries(arpes es_la)

if (NOT TARGET es_fe)
	add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/es_fe")
endif()
target_link_libraries(arpes es_fe)

