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

if(DEFINED ENV{ARCHER})
# https://github.com/pism/pism/issues/316
  set (PETSC_EXECUTABLE_RUNS "YES")
#Change this to libcraypetsc_cray_real.so for Cray compiler
  set (MY_PETSC_LIBRARY "libcraypetsc_gnu_real.so")
  set (MY_LIBRARY_SUBDIR "/lib/")
  set (MY_PETSC_LIBRARY_DIR $ENV{CRAY_PETSC_PREFIX_DIR})
  set (PETSC_LIBRARY_SINGLE "${MY_PETSC_LIBRARY_DIR}${MY_LIBRARY_SUBDIR}${MY_PETSC_LIBRARY}")
endif()

FIND_PACKAGE(PETSc)
find_package(MPI REQUIRED)
find_package(Threads)
if(DEFINED ENV{ARCHER})

else()
  find_package(LAPACK REQUIRED)
endif()

#if ("${PETSC_VERSION}" VERSION_LESS 3.4)
#  message(FATAL_ERROR "PETSc version 3.4 is required, but not found. (Found version:" ${PETSC_VERSION} ".")
#elseif ("${PETSC_VERSION}" VERSION_GREATER 3.6)
#  message(FATAL_ERROR "PETSc version 3.4 is required, but not found.. (Found version:" ${PETSC_VERSION} ".")
#endif ()

if (PETSC_FOUND)
  set (PETSC_EXECUTABLE_RUNS YES                  CACHE BOOL "Disable checking if this setup works" FORCE)
  set (PETSC_FOUND           YES                  CACHE BOOL "PETSc was found (manually)" FORCE)
  set (PETSC_INCLUDES        ${PETSC_INCLUDES}    CACHE STRING "Semicolon-delimited list of PETSc include directories" FORCE)
  set (PETSC_LIBRARIES       ${PETSC_LIBRARIES}   CACHE STRING "Semicolon-delimited list of PETSc libraries" FORCE)
  set (PETSC_COMPILER        ${PETSC_COMPILER}    CACHE FILEPATH "PETSc compiler; helpful to find a compatible MPI" FORCE)
  set (PETSC_DEFINITIONS     ${PETSC_DEFINITIONS} CACHE STRING "PETSc definitions" FORCE)
  set (PETSC_MPIEXEC         ${PETSC_MPIEXEC}     CACHE FILEPATH "Executable for running PETSc MPI programs" FORCE)
  set (PETSC_VERSION         ${PETSC_VERSION}     CACHE STRING "PETSc version: MAJOR.MINOR.SUBMINOR" FORCE)
  mark_as_advanced (PETSC_INCLUDES PETSC_LIBRARIES
    PETSC_COMPILER PETSC_DEFINITIONS
    PETSC_MPIEXEC PETSC_EXECUTABLE_RUNS PETSC_VERSION)
endif ()

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

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
    set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize ")
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
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe ")
    set (CMAKE_C_FLAGS_RELEASE "-O2 -march=native -ftree-vectorize ")
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
    set (CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize ")
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


add_subdirectory(spatialindex-1.8.0 EXCLUDE_FROM_ALL)
#add_subdirectory(libjudy)

include(ExternalProject)
ExternalProject_Add(libjudy
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/libjudy
  PREFIX     ${PROJECT_BINARY_DIR}/libjudy
  BINARY_DIR ${PROJECT_BINARY_DIR}/libjudy
#  CONFIGURE_COMMAND CXX=mpicxx CC=mpicc FC=mpif90 ${CMAKE_CURRENT_SOURCE_DIR}/configure --prefix=${PROJECT_BINARY_DIR} --enable-shared --enable-static --enable-debug
# errors  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/configure "CFLAGS=-O3 -march=native" --prefix=${PROJECT_BINARY_DIR} --enable-shared --enable-static
  CONFIGURE_COMMAND ./configure "CFLAGS=${SUBPROJECTS_C_FLAGS}" "CXXFLAGS=${SUBPROJECTS_CXX_FLAGS}" "FFLAGS=${SUBPROJECTS_Fortran_FLAGS}" --prefix=${PROJECT_BINARY_DIR} --enable-shared --enable-static ${local_debugging}
#  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}
  BUILD_COMMAND make
#  BUILD_DIR ${PROJECT_BINARY_DIR}
#  BUILD_IN_SOURCE 1
  INSTALL_COMMAND echo "installing"
)
# copy SRC directory from source to build
ExternalProject_Add_Step(libjudy copy_libjudy
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/libjudy ${PROJECT_BINARY_DIR}/libjudy
  COMMENT copying libjudy from source to build
  DEPENDEES download
)

include_directories(${MPI_Fortran_INCLUDE_PATH} ${PETSC_INCLUDES})
ADD_DEFINITIONS(${PETSC_DEFINITIONS})

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
  set (libsupermesh_OVERLAP_COMPUTE_COMMS 1)
else()
  set (libsupermesh_OVERLAP_COMPUTE_COMMS 0)
endif (USE_OVERLAP_COMPUTE_COMMS)

if (lower_build_type MATCHES debug)
  set (libsupermesh_DEBUG 1)
else ()
  set (libsupermesh_DEBUG 0)
endif (lower_build_type MATCHES debug)
if (PETSC_FOUND)
  set (libsupermesh_HAVE_PETSC 1)
else ()
  set (libsupermesh_HAVE_PETSC 0)
endif (PETSC_FOUND)

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
include_directories ("${PROJECT_SOURCE_DIR}/spatialindex-1.8.0/include/")
include_directories ("${PROJECT_BINARY_DIR}/libjudy/src/")

set (CPP_SPATIALINDEX_O_FILES "")
set (CPP_SPATIALINDEX_O_FILES
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/tprtree/Index.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/tprtree/TPRTree.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/tprtree/Leaf.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/tprtree/Statistics.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/tprtree/Node.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/storagemanager/MemoryStorageManager.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/storagemanager/RandomEvictionsBuffer.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/storagemanager/Buffer.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/storagemanager/DiskStorageManager.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/mvrtree/Index.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/mvrtree/Leaf.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/mvrtree/MVRTree.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/mvrtree/Statistics.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/mvrtree/Node.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/spatialindex/SpatialIndexImpl.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/spatialindex/Point.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/spatialindex/MovingRegion.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/spatialindex/MovingPoint.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/spatialindex/TimeRegion.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/spatialindex/Region.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/spatialindex/LineSegment.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/spatialindex/TimePoint.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/tools/Tools.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/rtree/Index.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/rtree/Leaf.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/rtree/RTree.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/rtree/Statistics.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/rtree/BulkLoader.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex.dir/rtree/Node.cc.o
)
foreach (test ${CPP_SPATIALINDEX_O_FILES})
 SET_SOURCE_FILES_PROPERTIES(${test} PROPERTIES GENERATED 1)
endforeach()

set (CPP_LIBJUDY_O_FILES "")
set (CPP_LIBJUDY_O_FILES
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/Judy1FreeArray.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/Judy1Set.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/Judy1Count.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/libcount_la-Judy1ByCount.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/Judy1Unset.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/Judy1Test.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/Judy1MallocIF.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/Judy1Tables.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/Judy1Decascade.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/Judy1Cascade.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/Judy1CreateBranch.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/Judy1InsertBranch.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/Judy1First.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/libnext_la-Judy1Next.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/libprev_la-Judy1Prev.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/libnext_la-Judy1NextEmpty.o
  ${PROJECT_BINARY_DIR}/libjudy/src/Judy1/libprev_la-Judy1PrevEmpty.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/JudyLFreeArray.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/JudyLIns.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/JudyLCount.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/JudyLGet.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/JudyLDel.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/libcount_la-JudyLByCount.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/JudyLTables.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/JudyLMallocIF.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/JudyLDecascade.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/JudyLCascade.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/JudyLInsertBranch.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/JudyLCreateBranch.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/JudyLFirst.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/libnext_la-JudyLNext.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/libprev_la-JudyLPrev.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/libnext_la-JudyLNextEmpty.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyL/libprev_la-JudyLPrevEmpty.o
  ${PROJECT_BINARY_DIR}/libjudy/src/JudyCommon/JudyMalloc.o
)
foreach (test ${CPP_LIBJUDY_O_FILES})
 SET_SOURCE_FILES_PROPERTIES(${test} PROPERTIES GENERATED 1)
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

add_library (supermesh ${FORTRAN_SRCFILES} ${CPP_SRCFILES} ${C_SRCFILES} ${CPP_SPATIALINDEX_O_FILES} ${CPP_LIBJUDY_O_FILES})

add_dependencies(supermesh libjudy spatialindex)

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
  COMMAND make clean
  COMMAND cd "${PROJECT_SOURCE_DIR}/libjudy" && rm -rf tmp/ && rm -rf src/libjudy-stamp/libjudy-done && rm -rf src/libjudy-stamp/ 
  COMMAND rm -rf unittests/
  COMMAND rm -rf data/
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#  COMMAND cd "${PROJECT_SOURCE_DIR}/libjud" && make distcleanlibjudy
)

if(UNIX AND NOT APPLE)
  target_link_libraries(supermesh ${MPI_Fortran_LIBRARIES} rt)

  foreach (test ${UNITTESTS})
    target_link_libraries(${test} supermesh ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PETSC_LIBRARIES} ${LAPACK_LIBRARIES} rt)
  endforeach()
  foreach (benchmark ${BENCHMARKTESTS})
    target_link_libraries(${benchmark} supermesh ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PETSC_LIBRARIES} ${LAPACK_LIBRARIES} rt)
  endforeach()
else ()
  target_link_libraries(supermesh ${MPI_Fortran_LIBRARIES})

  foreach (test ${UNITTESTS})
    target_link_libraries(${test} supermesh ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PETSC_LIBRARIES} ${LAPACK_LIBRARIES})
  endforeach()
  foreach (benchmark ${BENCHMARKTESTS})
    target_link_libraries(${benchmark} supermesh ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PETSC_LIBRARIES} ${LAPACK_LIBRARIES})
  endforeach()
endif()

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
