# CMakeLists files in this project can
# refer to the root source directory of the project as ${SUPERMESH_SOURCE_DIR} and
# to the root binary directory of the project as ${SUPERMESH_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (SUPERMESH C CXX Fortran)
enable_language (Fortran)
enable_testing ()

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")

string(TOLOWER "${CMAKE_BUILD_TYPE}" lower_cmake_build_type)
if(lower_cmake_build_type STREQUAL debug)
  message(STATUS "Enabling debugging")
  set(libsupermesh_DEBUG 1)
elseif(lower_cmake_build_type STREQUAL release)
  message(STATUS "Disabling debugging")
  set(libsupermesh_DEBUG 0)
elseif(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Disabling debugging")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type. Valid values are \"Debug\" and \"Release\"." FORCE)
  set(libsupermesh_DEBUG 0)
else()
   message(SEND_ERROR "Unrecognised build type")
endif()
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

find_package(MPI REQUIRED)
set(link_libraries ${MPI_Fortran_LIBRARIES})

option(ENABLE_JUDY "Enable use of the Judy library" OFF)
if(ENABLE_JUDY)
  find_library(JUDY_LIBRARY NAMES Judy)
  find_path(JUDY_INCLUDE NAMES Judy.h)
  find_package_handle_standard_args(Judy DEFAULT_MSG REQUIRED_VARS JUDY_LIBRARY JUDY_INCLUDE)
  set(link_libraries ${link_libraries} ${JUDY_LIBRARY})
  include_directories(${JUDY_INCLUDE})
else()
  set(JUDY_LIBRARY "")
  set(JUDY_INCLUDE "")
endif()

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (CMAKE_C_FLAGS)
# CMAKE_C_FLAGS has been defined by the USER. Do not modify his settings.
  SET (local_IGNORE_ALL_C_FLAGS "TRUE")
  MESSAGE ("CMAKE_C_FLAGS have been already set.")
endif (CMAKE_C_FLAGS)

if (CMAKE_CXX_FLAGS)
# CMAKE_CXX_FLAGS has been defined by the USER. Do not modify his settings.
  SET (local_IGNORE_ALL_CXX_FLAGS "TRUE")
  MESSAGE ("CMAKE_CXX_FLAGS have been already set.")
endif (CMAKE_CXX_FLAGS)

if (CMAKE_Fortran_FLAGS)
# CMAKE_Fortran_FLAGS has been defined by the USER. Do not modify his settings.
  SET (local_IGNORE_ALL_Fortran_FLAGS "TRUE")
  MESSAGE ("CMAKE_Fortran_FLAGS have been already set.")
endif (CMAKE_Fortran_FLAGS)


if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # gfortran
  if (NOT local_IGNORE_ALL_Fortran_FLAGS MATCHES "TRUE")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fno-f2c -pipe ")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -DNDEBUG ")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Wall -fcheck=all -ffpe-trap=zero,overflow,underflow -frecord-marker=4 -fno-realloc-lhs -finit-integer=-66666 -finit-real=nan -fimplicit-none -fbacktrace ")
    if (BUILD_SHARED_LIBS)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC ")
    endif ()
  endif ()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # ifort (untested)
  if (NOT local_IGNORE_ALL_Fortran_FLAGS MATCHES "TRUE")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ")
    set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O2 -xhost -r8")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g  -check noarg_temp_created -fpe0 -warn alignments -warn declarations -warn general -warn interfaces -warn truncated_source -warn uncalled -warn uninitialized -warn usage -common_args -warn unused -fp-stack-check -check bounds -check uninit -check format -r8")
  endif ()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  # cray
  if (NOT local_IGNORE_ALL_Fortran_FLAGS MATCHES "TRUE")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -s real64 -fpic ")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -s real64 -fpic ")
  endif ()
endif ()

if (CMAKE_C_COMPILER_ID MATCHES "GNU")
  # gcc
  if (NOT local_IGNORE_ALL_C_FLAGS MATCHES "TRUE")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe -D_GNU_SOURCE ")
    set (CMAKE_C_FLAGS_RELEASE "-O2 -march=native -ftree-vectorize -DNDEBUG ")
    set (CMAKE_C_FLAGS_DEBUG   "-O0 -g -Wall ")
    if (BUILD_SHARED_LIBS)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC ")
    endif ()
  endif ()
elseif (CMAKE_C_COMPILER_ID MATCHES "Intel")
  # icc (untested)
  if (NOT local_IGNORE_ALL_C_FLAGS MATCHES "TRUE")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ")
    set (CMAKE_C_FLAGS_RELEASE "-O2 -xhost ")
    set (CMAKE_C_FLAGS_DEBUG   "-O0 -g")
  endif ()
elseif (CMAKE_C_COMPILER_ID MATCHES "Cray")
  # icc (untested)
  if (NOT local_IGNORE_ALL_C_FLAGS MATCHES "TRUE")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ")
    set (CMAKE_C_FLAGS_RELEASE "-O2 -fpic ")
    set (CMAKE_C_FLAGS_DEBUG   "-O0 -g -fpic ")
  endif ()
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # g++
  if (NOT local_IGNORE_ALL_CXX_FLAGS MATCHES "TRUE")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe ")
    set (CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -DNDEBUG ")
    set (CMAKE_CXX_FLAGS_DEBUG   "-O0 -g -Wall ")
    if (BUILD_SHARED_LIBS)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC ")
    endif ()
  endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  # icpc (untested)
  if (NOT local_IGNORE_ALL_CXX_FLAGS MATCHES "TRUE")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
    set (CMAKE_CXX_FLAGS_RELEASE "-O2 -xhost ")
    set (CMAKE_CXX_FLAGS_DEBUG   "-O0 -g")
  endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Cray")
  # icpc (untested)
  if (NOT local_IGNORE_ALL_CXX_FLAGS MATCHES "TRUE")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
    set (CMAKE_CXX_FLAGS_RELEASE "-O2 -fpic ")
    set (CMAKE_CXX_FLAGS_DEBUG   "-O0 -g -fpic ")
  endif ()
endif ()

SET( SUBPROJECTS_C_FLAGS ${CMAKE_C_FLAGS})
SET( SUBPROJECTS_CXX_FLAGS ${CMAKE_CXX_FLAGS})
SET( SUBPROJECTS_Fortran_FLAGS ${CMAKE_Fortran_FLAGS})

set( SUBPROJECTS_debugging "NO")
string( TOLOWER "${CMAKE_BUILD_TYPE}" lower_build_type)
if (lower_build_type MATCHES release)
  if (NOT local_IGNORE_ALL_C_FLAGS MATCHES "TRUE")
    set( SUBPROJECTS_C_FLAGS "${SUBPROJECTS_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
  endif ()
  
  if (NOT local_IGNORE_ALL_CXX_FLAGS MATCHES "TRUE")
    set( SUBPROJECTS_CXX_FLAGS "${SUBPROJECTS_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  endif ()

  if (NOT local_IGNORE_ALL_Fortran_FLAGS MATCHES "TRUE")
    set( SUBPROJECTS_Fortran_FLAGS "${SUBPROJECTS_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
  endif ()
elseif (lower_build_type MATCHES debug)
  if (NOT local_IGNORE_ALL_C_FLAGS MATCHES "TRUE")
    set( SUBPROJECTS_C_FLAGS "${SUBPROJECTS_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
  endif ()

  if (NOT local_IGNORE_ALL_CXX_FLAGS MATCHES "TRUE")
    set( SUBPROJECTS_CXX_FLAGS "${SUBPROJECTS_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  endif ()

  if (NOT local_IGNORE_ALL_Fortran_FLAGS MATCHES "TRUE")
    set( SUBPROJECTS_Fortran_FLAGS "${SUBPROJECTS_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
  endif ()

  set( SUBPROJECTS_debugging "YES")
endif (lower_build_type MATCHES release)

if (SUBPROJECTS_debugging MATCHES YES)
  set(local_debugging "--enable-debug")
endif (SUBPROJECTS_debugging MATCHES YES)

option (USE_PROFILE "Use profile flags" OFF)
if (USE_PROFILE)
  SET( CMAKE_CXX_FLAGS             "${CMAKE_CXX_FLAGS} -pg -g")
  SET( CMAKE_C_FLAGS               "${CMAKE_C_FLAGS} -pg -g" )
  SET( CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -pg -g")
  SET( CMAKE_EXE_LINKER_FLAGS      "${CMAKE_EXE_LINKER_FLAGS} -pg -g" )
  SET( CMAKE_SHARED_LINKER_FLAGS   "${CMAKE_SHARED_LINKER_FLAGS} -pg -g" )
  SET( CMAKE_MODULE_LINKER_FLAGS   "${CMAKE_MODULE_LINKER_FLAGS} -pg -g" )
endif (USE_PROFILE)

option (USE_PROFPROFILE "Use perftools profile flags" OFF)
if (USE_PROFPROFILE)
  SET( CMAKE_CXX_FLAGS             "${CMAKE_CXX_FLAGS} -lprofiler -g")
  SET( CMAKE_C_FLAGS               "${CMAKE_C_FLAGS} -lprofiler -g" )
  SET( CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -lprofiler -g")
  SET( CMAKE_EXE_LINKER_FLAGS      "${CMAKE_EXE_LINKER_FLAGS} -lprofiler -g" )
  SET( CMAKE_SHARED_LINKER_FLAGS   "${CMAKE_SHARED_LINKER_FLAGS} -lprofiler -g" )
  SET( CMAKE_MODULE_LINKER_FLAGS   "${CMAKE_MODULE_LINKER_FLAGS} -lprofiler -g" )
endif (USE_PROFPROFILE)

message ("MAIN CMAKE_C_FLAGS:${CMAKE_C_FLAGS}")
message ("MAIN CMAKE_CXX_FLAGS:${CMAKE_CXX_FLAGS}")
message ("MAIN CMAKE_Fortran_FLAGS:${CMAKE_Fortran_FLAGS}")

message ("MAIN SUBPROJECTS_C_FLAGS:${SUBPROJECTS_C_FLAGS}")
message ("MAIN SUBPROJECTS_CXX_FLAGS:${SUBPROJECTS_CXX_FLAGS}")
message ("MAIN SUBPROJECTS_Fortran_FLAGS:${SUBPROJECTS_Fortran_FLAGS}")

message ("MAIN PROJECT_SOURCE_DIR:${PROJECT_SOURCE_DIR}")
message ("MAIN PROJECT_BINARY_DIR:${PROJECT_BINARY_DIR}")

include_directories(${MPI_Fortran_INCLUDE_PATH})

# The version number.
set (libsupermesh_VERSION_MAJOR 0)
set (libsupermesh_VERSION_MINOR 9)
set (libsupermesh_VERSION_MINOR_CHANGE 0)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/config/libsupermesh.pc.in"
  "${PROJECT_BINARY_DIR}/include/libsupermesh.pc"
  )

option (USE_CODE_PROFILE
        "Add profile timers and print statements" OFF)
if (USE_CODE_PROFILE)
  set (libsupermesh_PROFILE 1)
else()
  set (libsupermesh_PROFILE 0)
endif (USE_CODE_PROFILE)

option (USE_OVERLAP_COMPUTE_COMMS
        "Use overlap compute/comms. Make sure that your MPI impementation supports overlapping" OFF)
if (USE_OVERLAP_COMPUTE_COMMS)
  set (OVERLAP_COMPUTE_COMMS 1)
else()
  set (OVERLAP_COMPUTE_COMMS 0)
endif (USE_OVERLAP_COMPUTE_COMMS)

configure_file (
  "${PROJECT_SOURCE_DIR}/include/libsupermesh_configuration.h.in"
  "${PROJECT_BINARY_DIR}/include/libsupermesh_configuration.h"
  )
# add the binary tree/include to the search path for include files
# so that we will find the files
include_directories("${PROJECT_BINARY_DIR}/include")

# should we use our own LIBRARY?
option (USE_MYLIB 
        "Use library implementation" OFF) 

include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories(${PROJECT_SOURCE_DIR}/spatialindex-1.8.5/include/)
include_directories(${PROJECT_SOURCE_DIR}/spatialindex-1.8.5/src/)

add_subdirectory(spatialindex-1.8.5 EXCLUDE_FROM_ALL)
set(spatialindex_object_files
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tprtree/Index.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tprtree/TPRTree.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tprtree/Leaf.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tprtree/Statistics.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tprtree/Node.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/storagemanager/MemoryStorageManager.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/storagemanager/RandomEvictionsBuffer.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/storagemanager/Buffer.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/storagemanager/DiskStorageManager.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/mvrtree/Index.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/mvrtree/Leaf.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/mvrtree/MVRTree.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/mvrtree/Statistics.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/mvrtree/Node.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/SpatialIndexImpl.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/Point.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/MovingRegion.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/MovingPoint.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/TimeRegion.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/Region.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/LineSegment.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/TimePoint.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tools/Tools.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/rtree/Index.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/rtree/Leaf.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/rtree/RTree.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/rtree/Statistics.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/rtree/BulkLoader.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/rtree/Node.cc${CMAKE_CXX_OUTPUT_EXTENSION}
)
check_function_exists(srand48 HAVE_SRAND48)
if(NOT HAVE_SRAND48)
  set(spatialindex_object_files ${spatialindex_object_files}
    ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tools/rand48.cc${CMAKE_CXX_OUTPUT_EXTENSION})
endif()
foreach(spatialindex_object_file ${spatialindex_object_files})
  set_source_files_properties(${spatialindex_object_file} PROPERTIES GENERATED 1)
endforeach()
  
set (CPP_TEST_SOURCE_FILES "")
set (CPP_TEST_SOURCE_FILES
  ${PROJECT_SOURCE_DIR}/src/tests/test_main.cpp
)


# Create a library called "supermesh" which includes the source files.
# The extension is already found. Any number of sources could be listed here.
FILE(GLOB CPP_SRCFILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
FILE(GLOB C_SRCFILES ${PROJECT_SOURCE_DIR}/src/*.c)
FILE(GLOB FORTRAN_SRCFILES ${PROJECT_SOURCE_DIR}/src/*.F90)
FILE(GLOB FORTRAN_TEST_INTERVAL_INTERSECTOR_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_interval_intersector.F90)
FILE(GLOB FORTRAN_TEST_TRI_INTERSECTOR_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_tri_intersector.F90)
FILE(GLOB FORTRAN_TEST_TET_INTERSECTOR_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_tet_intersector.F90)
FILE(GLOB FORTRAN_TEST_ELEMENT_INTERSECTOR_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_element_intersector.F90)
FILE(GLOB FORTRAN_TEST_INTERSECTION_FINDER_2D_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_intersection_finder_2d.F90)
FILE(GLOB FORTRAN_TEST_INTERSECTION_FINDER_3D_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_intersection_finder_3d.F90)
FILE(GLOB FORTRAN_TEST_INTERSECTION_FINDER_COMPLETENESS_1D_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_intersection_finder_completeness_1d.F90)
FILE(GLOB FORTRAN_TEST_INTERSECTION_FINDER_COMPLETENESS_2D_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_intersection_finder_completeness_2d.F90)
FILE(GLOB FORTRAN_TEST_INTERSECTION_FINDER_COMPLETENESS_3D_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_intersection_finder_completeness_3d.F90)
FILE(GLOB FORTRAN_TEST_PARALLEL_PARTITION_COMPLETE_AB_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_parallel_partition_complete_ab.F90)
FILE(GLOB FORTRAN_TEST_PARALLEL_PARTITION_COMPLETE_AB_SAME_ALGO_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_parallel_partition_complete_ab_same_algo.F90)
FILE(GLOB FORTRAN_TEST_PARALLEL_UNIVERSAL_NODE_NUMBERING_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_parallel_universal_node_numbering.F90)
FILE(GLOB FORTRAN_BENCHMARK_PARALLEL_COMPLETE_FILES ${PROJECT_SOURCE_DIR}/tests/benchmark_parallel_complete.F90)
FILE(GLOB FORTRAN_BENCHMARK_PARALLEL_COMPLETE_3D_FILES ${PROJECT_SOURCE_DIR}/tests/benchmark_parallel_complete_3d.F90)
FILE(GLOB FORTRAN_BENCHMARK_SERIAL_FILES ${PROJECT_SOURCE_DIR}/tests/benchmark_serial.F90)
FILE(GLOB FORTRAN_BENCHMARK_SERIAL_3D_FILES ${PROJECT_SOURCE_DIR}/tests/benchmark_serial_3d.F90)
FILE(GLOB FORTRAN_BENCHMARK_SERIAL_SAME_ALGO_FILES ${PROJECT_SOURCE_DIR}/tests/benchmark_serial_same_algo.F90)
FILE(GLOB FORTRAN_BENCHMARK_SERIAL_SAME_ALGO_3D_FILES ${PROJECT_SOURCE_DIR}/tests/benchmark_serial_same_algo_3d.F90)
FILE(GLOB FORTRAN_BENCHMARK_PARALLEL_P1_INNER_PRODUCT_FILES ${PROJECT_SOURCE_DIR}/tests/benchmark_parallel_p1_inner_product.F90)
FILE(GLOB FORTRAN_BENCHMARK_PARALLEL_P1_INNER_PRODUCT_3D_FILES ${PROJECT_SOURCE_DIR}/tests/benchmark_parallel_p1_inner_product_3d.F90)
FILE(GLOB FORTRAN_BENCHMARK_PARALLEL_P2_INNER_PRODUCT_FILES ${PROJECT_SOURCE_DIR}/tests/benchmark_parallel_p2_inner_product.F90)

FILE(GLOB TEST_DATA_FILES ${PROJECT_SOURCE_DIR}/src/tests/data/*)

add_library (supermesh ${FORTRAN_SRCFILES} ${CPP_SRCFILES} ${C_SRCFILES} ${spatialindex_object_files})
add_dependencies(supermesh spatialindex)
  
set (UNITTESTS
  test_interval_intersector
  test_tri_intersector
  test_tet_intersector
  test_element_intersector
  test_intersection_finder_2d
  test_intersection_finder_3d
  test_intersection_finder_completeness_1d
  test_intersection_finder_completeness_2d
  test_intersection_finder_completeness_3d
  test_parallel_partition_complete_ab
  test_parallel_partition_complete_ab_same_algo
  test_parallel_universal_node_numbering
)

set (BENCHMARKTESTS
  benchmark_parallel_complete
  benchmark_parallel_complete_3d
  benchmark_serial
  benchmark_serial_3d
  benchmark_serial_same_algo
  benchmark_serial_same_algo_3d
  benchmark_parallel_p1_inner_product
  benchmark_parallel_p1_inner_product_3d
  benchmark_parallel_p2_inner_product
)

add_executable (test_interval_intersector ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_INTERVAL_INTERSECTOR_FILES})
add_executable (test_tri_intersector ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_TRI_INTERSECTOR_FILES})
add_executable (test_tet_intersector ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_TET_INTERSECTOR_FILES})
add_executable (test_element_intersector ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_ELEMENT_INTERSECTOR_FILES})
add_executable (test_intersection_finder_2d ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_INTERSECTION_FINDER_2D_FILES})
add_executable (test_intersection_finder_3d ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_INTERSECTION_FINDER_3D_FILES})
add_executable (test_intersection_finder_completeness_1d ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_INTERSECTION_FINDER_COMPLETENESS_1D_FILES})
add_executable (test_intersection_finder_completeness_2d ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_INTERSECTION_FINDER_COMPLETENESS_2D_FILES})
add_executable (test_intersection_finder_completeness_3d ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_INTERSECTION_FINDER_COMPLETENESS_3D_FILES})
add_executable (test_parallel_partition_complete_ab ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_PARALLEL_PARTITION_COMPLETE_AB_FILES})
add_executable (test_parallel_partition_complete_ab_same_algo ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_PARALLEL_PARTITION_COMPLETE_AB_SAME_ALGO_FILES})
add_executable (test_parallel_universal_node_numbering ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_PARALLEL_UNIVERSAL_NODE_NUMBERING_FILES})
add_executable (benchmark_parallel_complete EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_PARALLEL_COMPLETE_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_parallel_complete_3d EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_PARALLEL_COMPLETE_3D_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_serial EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_SERIAL_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_serial_3d EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_SERIAL_3D_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_serial_same_algo EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_SERIAL_SAME_ALGO_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_serial_same_algo_3d EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_SERIAL_SAME_ALGO_3D_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_parallel_p1_inner_product EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_PARALLEL_P1_INNER_PRODUCT_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_parallel_p1_inner_product_3d EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_PARALLEL_P1_INNER_PRODUCT_3D_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_parallel_p2_inner_product EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_PARALLEL_P2_INNER_PRODUCT_FILES} ${CPP_TEST_SOURCE_FILES})

foreach (test ${UNITTESTS})
 SET_PROPERTY(TARGET "${test}" PROPERTY COMPILE_DEFINITIONS "TESTNAME=${test}")
 if (test MATCHES parallel)
   ADD_TEST("${test}" ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${test})
 else()
   ADD_TEST("${test}" ${test})
 endif (test MATCHES parallel)
 SET_TESTS_PROPERTIES("${test}" PROPERTIES FAIL_REGULAR_EXPRESSION "Fail:")
endforeach()
foreach (benchmark ${BENCHMARKTESTS})
 SET_PROPERTY(TARGET "${benchmark}" PROPERTY COMPILE_DEFINITIONS "TESTNAME=${benchmark}")
endforeach()

add_custom_target(benchmark DEPENDS benchmark_parallel_complete benchmark_parallel_complete_3d benchmark_serial benchmark_serial_3d benchmark_serial_same_algo benchmark_serial_same_algo_3d benchmark_parallel_p1_inner_product benchmark_parallel_p1_inner_product_3d benchmark_parallel_p2_inner_product)

add_custom_command(TARGET test_tri_intersector
  POST_BUILD
  COMMAND mkdir -p unittests/data
  COMMAND cp test_tri_intersector unittests/
  COMMAND mkdir -p data
  COMMAND cp -r "${PROJECT_SOURCE_DIR}/src/tests/data/" ./
  COMMAND cp -r "${PROJECT_SOURCE_DIR}/src/tests/data/" unittests/
  COMMAND cp -r "${PROJECT_SOURCE_DIR}/tests/data/" ./
  COMMAND unzip -o data/2D_09.zip -d data/
  COMMAND unzip -o data/2D_05.zip -d data/
  COMMAND unzip -o data/2D_9.zip -d data/
  COMMAND unzip -o data/2D_5.zip -d data/
  COMMAND unzip -o data/3D_9.zip -d data/
  COMMAND unzip -o data/3D_5.zip -d data/
  COMMAND unzip -o data/2D_0_01.zip -d data/
#  COMMAND cd unittests && ./test_tri_intersector && cd ..
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Copying test data files"
  VERBATIM
)

add_custom_target(distclean 
  COMMAND cd "${PROJECT_BINARY_DIR}"
  COMMAND make clean
  COMMAND rm -rf unittests/
  COMMAND rm -rf data/
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

set(test_link_libraries supermesh ${link_libraries} ${MPI_CXX_LIBRARIES})
target_link_libraries(supermesh ${link_libraries})
foreach(test ${UNITTESTS})
  target_link_libraries(${test} ${test_link_libraries})
endforeach()
foreach(benchmark ${BENCHMARKTESTS})
  target_link_libraries(${benchmark} ${test_link_libraries})
endforeach()

set_target_properties(supermesh            PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
set_target_properties(supermesh            PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${LAPACK_LINKER_FLAGS}")
target_link_libraries(supermesh ${MPI_Fortran_LIBRARIES} rt)

if(MPI_Fortran_COMPILE_FLAGS)
  set_target_properties(supermesh            PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
  foreach (test ${UNITTESTS})
    set_target_properties(${test} PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
  endforeach()
  foreach (benchmark ${BENCHMARKTESTS})
    set_target_properties(${benchmark} PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
  endforeach()
endif()

if(MPI_Fortran_LINK_FLAGS)
  set_target_properties(supermesh            PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
  foreach (test ${UNITTESTS})
    set_target_properties(${test} PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
  endforeach()
  foreach (benchmark ${BENCHMARKTESTS})
    set_target_properties(${benchmark} PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
  endforeach()
endif()

install(TARGETS supermesh DESTINATION lib)
install(FILES "${PROJECT_BINARY_DIR}/libsupermesh_graphs.mod" DESTINATION include)
install(FILES "${PROJECT_BINARY_DIR}/libsupermesh_integer_hash_table.mod" DESTINATION include)
install(FILES "${PROJECT_BINARY_DIR}/libsupermesh_integer_set.mod" DESTINATION include)
install(FILES "${PROJECT_BINARY_DIR}/libsupermesh_intersection_finder.mod" DESTINATION include)
install(FILES "${PROJECT_BINARY_DIR}/libsupermesh_interval_intersection.mod" DESTINATION include)
install(FILES "${PROJECT_BINARY_DIR}/libsupermesh_parallel_supermesh.mod" DESTINATION include)
install(FILES "${PROJECT_BINARY_DIR}/libsupermesh_read_triangle.mod" DESTINATION include)
install(FILES "${PROJECT_BINARY_DIR}/libsupermesh_supermesh.mod" DESTINATION include)
install(FILES "${PROJECT_BINARY_DIR}/libsupermesh_tet_intersection.mod" DESTINATION include)
install(FILES "${PROJECT_BINARY_DIR}/libsupermesh_tri_intersection.mod" DESTINATION include)
install(FILES "${PROJECT_BINARY_DIR}/include/libsupermesh.pc" DESTINATION lib/pkgconfig)
