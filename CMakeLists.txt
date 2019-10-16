cmake_minimum_required(VERSION 3.13)
project(arpes)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#set(MKL_LIB "$ENV{MKLROOT}/lib/intel64")

add_executable(arpes "src/arpes.cpp")
target_compile_features(arpes PUBLIC cxx_std_17)
target_compile_options(arpes PUBLIC -Wall -Wpedantic -Wextra -Werror=return-type -m64 -march=native 
					   $<$<CONFIG:DEBUG>:-O0 -g> $<$<CONFIG:RELEASE>:-Wno-unused-parameter -Wno-deprecated-declarations -O3 -DNDEBUG>)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	target_compile_options(arpes PUBLIC -Wno-psabi)
endif()

#target_link_directories(arpes PRIVATE "${MKL_LIB}")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/eslib")
target_link_libraries(arpes eslib)
