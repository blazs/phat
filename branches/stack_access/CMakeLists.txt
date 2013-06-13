cmake_minimum_required(VERSION 2.8)

project("Phat")

set(CMAKE_BUILD_TYPE "Release")

INCLUDE(CheckCXXSourceCompiles)

include_directories (include)

FIND_PACKAGE(OpenMP)

set (CMAKE_REQUIRED_FLAGS ${OpenMP_CXX_FLAGS})
CHECK_CXX_SOURCE_COMPILES("
#include <omp.h>
int main() {
#if (_OPENMP >= 200805 || _MSC_VER >= 1500)
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
 
FILE(GLOB_RECURSE all_headers "include/phat/*.h")

FILE(GLOB general_includes "include/phat/common/*.h")
FILE(GLOB random_access_includes "include/phat/random_access/*.h")
FILE(GLOB random_access_representations_includes "include/phat/random_access/representations/*.h")
FILE(GLOB random_access_reducers_includes "include/phat/random_access/reducers/*.h")
FILE(GLOB stack_access_includes "include/phat/stack_access/*.h")
FILE(GLOB stack_access_representations_includes "include/phat/stack_access/representations/*.h")
FILE(GLOB stack_access_reducers_includes "include/phat/stack_access/reducers/*.h")
FILE(GLOB auto_reducing_includes "include/phat/auto_reducing/*.h")
FILE(GLOB auto_reducing_representations_includes "include/phat/auto_reducing/representations/*.h")
FILE(GLOB auto_reducing_reducers_includes "include/phat/auto_reducing/reducers/*.h")

add_executable (random_access_example src/random_access_example.cpp ${all_headers})
add_executable (stack_access_example src/stack_access_example.cpp ${all_headers})
add_executable (auto_reducing_example src/auto_reducing_example.cpp ${all_headers})
add_executable (self_test src/self_test.cpp ${all_headers})
add_executable (phat src/phat.cpp ${all_headers})
add_executable (info src/info.cpp ${all_headers})
add_executable (benchmark src/benchmark.cpp ${all_headers})

source_group(Header\ Files\\common FILES ${general_includes})
source_group(Header\ Files\\random_access FILES ${random_access_includes})
source_group(Header\ Files\\random_access\\representations FILES ${random_access_representations_includes})
source_group(Header\ Files\\random_access\\reducers FILES ${random_access_reducers_includes})
source_group(Header\ Files\\stack_access FILES ${stack_access_includes})
source_group(Header\ Files\\stack_access\\representations FILES ${stack_access_representations_includes})
source_group(Header\ Files\\stack_access\\reducers FILES ${stack_access_reducers_includes})
source_group(Header\ Files\\auto_reducing FILES ${auto_reducing_includes})
source_group(Header\ Files\\auto_reducing\\representations FILES ${auto_reducing_representations_includes})
source_group(Header\ Files\\auto_reducing\\reducers FILES ${auto_reducing_reducers_includes})