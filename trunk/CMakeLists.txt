cmake_minimum_required(VERSION 2.8)

include_directories (include)

FIND_PACKAGE(OpenMP)

if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
message("OpenMP not supported by the compiler! \nTo get optimal performance for the \"chunk\" algorithm, use a C++ compiler with OpenMP support (e.g., GCC).\nTo use a different compiler, pass it to cmake in the variable CMAKE_CXX_COMPILER: \n  cmake . -DCMAKE_CXX_COMPILER=g++")
endif()

add_executable (simple_example src/simple_example.cpp)
add_executable (self_test src/self_test.cpp)
add_executable (phat src/phat.cpp)
