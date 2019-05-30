cmake_minimum_required(VERSION 3.13)
project(arpes)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -Wextra -Werror=return-type -m64 -march=native -fmax-errors=5")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Wno-unused-parameter -Wno-deprecated-declarations -O3 -DNDEBUG")

set(MKL_LIB "$ENV{MKLROOT}/lib/intel64")

file(GLOB_RECURSE CPP_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_executable(arpes "${CPP_FILES}")
target_link_directories(arpes PRIVATE "${MKL_LIB}")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/es_util")
target_link_libraries(arpes es_util)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/es_la")
target_link_libraries(arpes es_la)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/es_fe")
target_link_libraries(arpes es_fe)
