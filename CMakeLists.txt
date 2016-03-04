cmake_minimum_required(VERSION 2.8)
project(libsupermesh C CXX Fortran)
enable_language(Fortran)
enable_testing()
set(libsupermesh_VERSION_MAJOR 0)
set(libsupermesh_VERSION_MINOR 9)
set(libsupermesh_VERSION_SUBMINOR 0)
set(libsupermesh_VERSION_RELEASE 0)
exec_program("git" ${CMAKE_SOURCE_DIR} ARGS "rev-parse HEAD" OUTPUT_VARIABLE libsupermesh_GIT_REVISION RETURN_VALUE git_error)
if(NOT "${git_error}" STREQUAL 0)
  set(libsupermesh_GIT_REVISION unknown)
endif()
message(STATUS "Git revision: ${libsupermesh_GIT_REVISION}")
if(libsupermesh_VERSION_RELEASE STREQUAL 0)
  set(libsupermesh_VERSION ${libsupermesh_VERSION_MAJOR}.${libsupermesh_VERSION_MINOR}.${libsupermesh_VERSION_SUBMINOR}+)
else()
  set(libsupermesh_VERSION ${libsupermesh_VERSION_MAJOR}.${libsupermesh_VERSION_MINOR}.${libsupermesh_VERSION_SUBMINOR})
endif()

if(${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
  message(SEND_ERROR "In-source builds are not permitted")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake-modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/include_local)

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
include_directories(${MPI_Fortran_INCLUDE_PATH})

option(ENABLE_JUDY "Enable use of the Judy library" OFF)
if(ENABLE_JUDY)
  find_library(JUDY_LIBRARY NAMES Judy)
  find_path(JUDY_INCLUDE NAMES Judy.h)
  find_package_handle_standard_args(Judy DEFAULT_MSG REQUIRED_VARS JUDY_LIBRARY JUDY_INCLUDE)
  list(APPEND link_libraries ${JUDY_LIBRARY})
  include_directories(${JUDY_INCLUDE})
else()
  set(JUDY_LIBRARY "")
  set(JUDY_INCLUDE "")
endif()

option(ENABLE_TIMERS "Enable internal timers" OFF)
if(ENABLE_TIMERS)
  set(libsupermesh_TIMERS 1)
else()
  set(libsupermesh_TIMERS 0)
endif()

option(OVERLAP_COMPUTE_COMMS "Overlap computation and communication. Efficiency of this depends upon the specific MPI implementation." OFF)
if(OVERLAP_COMPUTE_COMMS)
  set(OVERLAP_COMPUTE_COMMS 1)
else()
  set(OVERLAP_COMPUTE_COMMS 0)
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
    set (CMAKE_C_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -DNDEBUG ")
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

message(STATUS "CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}")
message(STATUS "CMAKE_Fortran_FLAGS_RELEASE=${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG}")
message(STATUS "CMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}")

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
  list(APPEND spatialindex_object_files
    ${PROJECT_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tools/rand48.cc${CMAKE_CXX_OUTPUT_EXTENSION})
endif()
foreach(spatialindex_object_file ${spatialindex_object_files})
  set_source_files_properties(${spatialindex_object_file} PROPERTIES GENERATED 1)
endforeach()

file(GLOB source_files ${PROJECT_SOURCE_DIR}/src/*.F90 ${PROJECT_SOURCE_DIR}/src/*.c ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_library(supermesh ${source_files} ${spatialindex_object_files})
add_dependencies(supermesh spatialindex)

set(test_main ${PROJECT_SOURCE_DIR}/src/tests/test_main.cpp)
file(GLOB unittest_files ${PROJECT_SOURCE_DIR}/src/tests/*.F90)
set(unittests "")
foreach(unittest_file ${unittest_files})
  get_filename_component(unittest ${unittest_file} NAME_WE)
  add_executable(${unittest} ${unittest_file} ${test_main})
  list(APPEND unittests ${unittest})
  set_property(TARGET ${unittest} PROPERTY COMPILE_DEFINITIONS "TESTNAME=${unittest}")
  if(unittest MATCHES parallel)
    add_test(${unittest} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${unittest})
  else()
    add_test(${unittest} ${unittest})
  endif()
  set_tests_properties(${unittest} PROPERTIES FAIL_REGULAR_EXPRESSION "Fail:")
endforeach()

file(GLOB benchmark_files ${PROJECT_SOURCE_DIR}/tests/*.F90)
set(benchmarks "")
foreach(benchmark_file ${benchmark_files})
  get_filename_component(benchmark ${benchmark_file} NAME_WE)
  add_executable(${benchmark} EXCLUDE_FROM_ALL ${benchmark_file} ${test_main})
  list(APPEND benchmarks ${benchmark})
  set_property(TARGET ${benchmark} PROPERTY COMPILE_DEFINITIONS "TESTNAME=${benchmark}")
endforeach()
add_custom_target(benchmark DEPENDS ${benchmarks})

set(test_link_libraries supermesh ${link_libraries} ${MPI_CXX_LIBRARIES})
target_link_libraries(supermesh ${link_libraries})
foreach(test ${unittests})
  target_link_libraries(${test} ${test_link_libraries})
endforeach()
foreach(benchmark ${benchmarks})
  target_link_libraries(${benchmark} ${test_link_libraries})
endforeach()

set_target_properties(supermesh PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
set_target_properties(supermesh PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
foreach(unittest ${unittests})
  set_target_properties(${unittest} PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
  set_target_properties(${unittest} PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
endforeach()
foreach(benchmark ${benchmarks})
  set_target_properties(${benchmark} PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
  set_target_properties(${benchmark} PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
endforeach()

configure_file (
  "${PROJECT_SOURCE_DIR}/config/libsupermesh.pc.in"
  "${PROJECT_BINARY_DIR}/config/libsupermesh.pc"
  )
configure_file (
  "${PROJECT_SOURCE_DIR}/include/libsupermesh.h.in"
  "${PROJECT_BINARY_DIR}/include/libsupermesh.h"
  )
configure_file (
  "${PROJECT_SOURCE_DIR}/include/libsupermesh_configuration.h.in"
  "${PROJECT_BINARY_DIR}/include_local/libsupermesh_configuration.h"
  )

install(TARGETS supermesh DESTINATION lib)
install(DIRECTORY ${CMAKE_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ${PROJECT_BINARY_DIR}/config/libsupermesh.pc DESTINATION lib/pkgconfig)

file(GLOB data_files ${PROJECT_SOURCE_DIR}/src/tests/data/*.ele
                     ${PROJECT_SOURCE_DIR}/src/tests/data/*.node)
file(COPY ${data_files} DESTINATION ${PROJECT_BINARY_DIR}/data)
file(GLOB data_files ${PROJECT_SOURCE_DIR}/src/tests/data/*.tar.bz2)
foreach(data_file ${data_files})
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvjf ${data_file} WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/data)
endforeach()