cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE "Release")

INCLUDE(CheckCXXSourceCompiles)

include_directories (include)
include_directories (/opt/local/include)

FIND_PACKAGE(OpenMP)

set (CMAKE_REQUIRED_FLAGS ${OpenMP_CXX_FLAGS})
CHECK_CXX_SOURCE_COMPILES("
#include <omp.h>
int main() {
#if (_OPENMP >= 200805)
  return 0;
#else
  breaks_on_purpose
#endif
}
" OPENMP_VERSION)

if(OPENMP_VERSION)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
message("

OpenMP 3.0 not supported by the compiler (${CMAKE_CXX_COMPILER})! 
To get optimal performance for the \"chunk\" algorithm, use a C++ compiler with OpenMP 3.0 support (e.g., GCC >=4.4).
To use a different compiler, pass it to cmake in the variable CMAKE_CXX_COMPILER: 
 cmake . -DCMAKE_CXX_COMPILER=g++-4.7

")
endif()

add_executable (simple_example src/simple_example.cpp)
add_executable (self_test src/self_test.cpp)
add_executable (phat src/phat.cpp)
add_executable (reduction_test reduction_test.cpp)
add_executable (benchmark src/benchmark.cpp)


set(CMAKE_MODULE_PATH "modules")
set(CGAL_DIR "/opt/local/lib/cmake")

FIND_PACKAGE(CGAL)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

FIND_PACKAGE(Boost)

add_executable (alpha_3 addons/alpha_3.cpp)
