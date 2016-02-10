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

FIND_PACKAGE(PETSc REQUIRED)
find_package(MPI REQUIRED)
find_package(Threads)
if(DEFINED ENV{ARCHER})

else()
  find_package(LAPACK REQUIRED)
endif()
#set(Python_ADDITIONAL_VERSIONS 2.7)
#find_package(PythonInterp 2.7 REQUIRED)
#find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)

if ("${PETSC_VERSION}" VERSION_LESS 3.4)
  message(FATAL_ERROR "PETSc version 3.4 is required, but not found. (Found version:" ${PETSC_VERSION} ".")
elseif ("${PETSC_VERSION}" VERSION_GREATER 3.6)
  message(FATAL_ERROR "PETSc version 3.4 is required, but not found.. (Found version:" ${PETSC_VERSION} ".")
endif ()

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
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0 -fno-f2c -march=native -O3 -ftree-vectorize -pipe -g ")
    set (CMAKE_Fortran_FLAGS_RELEASE "-fno-f2c -march=native -O3 -ftree-vectorize -pipe -g ")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fcheck=all -ffpe-trap=zero,overflow,underflow -frecord-marker=4 -fno-realloc-lhs -finit-integer=-66666 -finit-real=nan -fno-f2c -Og -g -ffree-line-length-0 -fimplicit-none -fbacktrace -Wextra -Wall -Wimplicit-interface -Wno-surprising -pipe")
  endif ()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # ifort (untested)
  if (NOT local_IGNORE_ALL_Fortran_FLAGS MATCHES "TRUE")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ")
    set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O2 -xhost -r8")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g  -check noarg_temp_created -fpe0 -warn alignments -warn declarations -warn general -warn interfaces -warn truncated_source -warn uncalled -warn uninitialized -warn usage -common_args -warn unused -fp-stack-check -check bounds -check uninit -check format -r8")
  endif ()
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  if (NOT local_IGNORE_ALL_Fortran_FLAGS MATCHES "TRUE")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ")
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O2 -m32 ")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
  endif ()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  # cray
  if (NOT local_IGNORE_ALL_Fortran_FLAGS MATCHES "TRUE")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -s real64 -fpic ")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -s real64 -fpic ")
  endif ()
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  if (NOT local_IGNORE_ALL_Fortran_FLAGS MATCHES "TRUE")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2 ")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  endif ()
endif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

if (CMAKE_C_COMPILER_ID MATCHES "GNU")
  # gcc
  if (NOT local_IGNORE_ALL_C_FLAGS MATCHES "TRUE")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -ftree-vectorize ")
    set (CMAKE_C_FLAGS_RELEASE "-march=native -O2 -pipe -g -ftree-vectorize ") 
    set (CMAKE_C_FLAGS_DEBUG   "-Og -g -Wextra -Wall -pipe")
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
endif (CMAKE_C_COMPILER_ID MATCHES "GNU")

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # g++
  if (NOT local_IGNORE_ALL_CXX_FLAGS MATCHES "TRUE")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -g -ftree-vectorize")
    set (CMAKE_CXX_FLAGS_RELEASE "-march=native -O2 -pipe -g -ftree-vectorize") 
    set (CMAKE_CXX_FLAGS_DEBUG   "-Og -g -Wextra -Wall -pipe")
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
endif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")

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

set( SUBPROJECTS_profilling      "")

option (USE_PROFILE "Use profile flags" OFF)
if (USE_PROFILE)
  SET( CMAKE_CXX_FLAGS             "${CMAKE_CXX_FLAGS} -pg -g")
  SET( CMAKE_C_FLAGS               "${CMAKE_C_FLAGS} -pg -g" )
  SET( CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -pg -g")
  SET( CMAKE_EXE_LINKER_FLAGS      "${CMAKE_EXE_LINKER_FLAGS} -pg -g" )
  SET( CMAKE_SHARED_LINKER_FLAGS   "${CMAKE_SHARED_LINKER_FLAGS} -pg -g" )
  SET( CMAKE_MODULE_LINKER_FLAGS   "${CMAKE_MODULE_LINKER_FLAGS} -pg -g" )
  set( SUBPROJECTS_profilling      "-pg -g" )
endif (USE_PROFILE)

option (USE_PROFPROFILE "Use perftools profile flags" OFF)
if (USE_PROFPROFILE)
  SET( CMAKE_CXX_FLAGS             "${CMAKE_CXX_FLAGS} -lprofiler -g")
  SET( CMAKE_C_FLAGS               "${CMAKE_C_FLAGS} -lprofiler -g" )
  SET( CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -lprofiler -g")
  SET( CMAKE_EXE_LINKER_FLAGS      "${CMAKE_EXE_LINKER_FLAGS} -lprofiler -g" )
  SET( CMAKE_SHARED_LINKER_FLAGS   "${CMAKE_SHARED_LINKER_FLAGS} -lprofiler -g" )
  SET( CMAKE_MODULE_LINKER_FLAGS   "${CMAKE_MODULE_LINKER_FLAGS} -lprofiler -g" )
  set( SUBPROJECTS_profilling      "-lprofiler -g" )
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
  CONFIGURE_COMMAND ./configure "CFLAGS=${SUBPROJECTS_C_FLAGS}" "CXXFLAGS=${SUBPROJECTS_CXX_FLAGS}" "FFLAGS=${SUBPROJECTS_Fortran_FLAGS}" --prefix=${PROJECT_BINARY_DIR} --enable-shared --enable-static ${SUBPROJECTS_profilling} ${local_debugging}
#  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}
  BUILD_COMMAND make
#  BUILD_DIR ${PROJECT_BINARY_DIR}
#  BUILD_IN_SOURCE 1
  INSTALL_COMMAND echo "installing"
)
# copy SRC directory from source to build
ExternalProject_Add_Step(libjudy copy_libjudy
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/libjudy ${CMAKE_BINARY_DIR}/libjudy
  COMMENT copying libjudy from source to build
  DEPENDEES download
)

include_directories(${MPI_Fortran_INCLUDE_PATH} ${PETSC_INCLUDES})
ADD_DEFINITIONS(${PETSC_DEFINITIONS})

# Set installation prefix automatically to custom path if not explicitly specified on the command line
#if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
#    set (CMAKE_INSTALL_PREFIX "/home/ipanourg/fluidity_git/fluidity" CACHE PATH "default install path" FORCE )
#endif()

# The version number.
set (SuperMesh_VERSION_MAJOR 0)
set (SuperMesh_VERSION_MINOR 0)
set (SuperMesh_VERSION_MINOR_CHANGE 1)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/include/SuperMeshConfig.h.in"
  "${PROJECT_BINARY_DIR}/include/SuperMeshConfig.h"
  )

configure_file (
  "${PROJECT_SOURCE_DIR}/config/SuperMesh.pc.in"
  "${PROJECT_BINARY_DIR}/include/SuperMesh.pc"
  )

set (SuperMesh_HAVE_MPI 1)

option (USE_CODE_PROFILE
        "Add profile timers and print statements" OFF)
if (USE_CODE_PROFILE)
  set (SuperMesh_PROFILE 1)
else()
  set (SuperMesh_PROFILE 0)
endif (USE_CODE_PROFILE)

option (USE_OVERLAP_COMPUTE_COMMS
        "Use overlap compute/comms. Make sure that your MPI impementation supports overlapping" OFF)
if (USE_OVERLAP_COMPUTE_COMMS)
  set (SuperMesh_OVERLAP_COMPUTE_COMMS 1)
else()
  set (SuperMesh_OVERLAP_COMPUTE_COMMS 0)
endif (USE_OVERLAP_COMPUTE_COMMS)

if (lower_build_type MATCHES debug)
  set (SuperMesh_DEBUG 1)
else ()
  set (SuperMesh_DEBUG 0)
endif (lower_build_type MATCHES debug)
if (PETSC_FOUND)
  set (SuperMesh_HAVE_PETSC 1)
else ()
  set (SuperMesh_HAVE_PETSC 0)
endif (PETSC_FOUND)

configure_file (
  "${PROJECT_SOURCE_DIR}/include/confdefs.h.in"
  "${PROJECT_BINARY_DIR}/include/confdefs.h"
  )
# add the binary tree/include to the search path for include files
# so that we will find the files
include_directories("${PROJECT_BINARY_DIR}/include")

# should we use our own LIBRARY?
option (USE_MYLIB 
        "Use library implementation" OFF) 

# add the library?
#
if (USE_MYLIB)
#  include_directories ("${PROJECT_SOURCE_DIR}/MathFunctions")
#  add_subdirectory (MathFunctions)
#  set (EXTRA_LIBS ${EXTRA_LIBS} MathFunctions)
endif (USE_MYLIB)

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
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex_c.dir/capi/ObjVisitor.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex_c.dir/capi/Utility.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex_c.dir/capi/BoundsQuery.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex_c.dir/capi/Error.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex_c.dir/capi/Index.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex_c.dir/capi/CountVisitor.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex_c.dir/capi/CustomStorage.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex_c.dir/capi/IdVisitor.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex_c.dir/capi/LeafQuery.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex_c.dir/capi/sidx_api.cc.o
  ${PROJECT_BINARY_DIR}/spatialindex-1.8.0/src/CMakeFiles/spatialindex_c.dir/capi/DataStream.cc.o
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


# Create a library called "SuperMesh" which includes the source files.
# The extension is already found. Any number of sources could be listed here.
#add_library (SuperMesh STATIC ${PROJECT_SOURCE_DIR}/src/Element_Intersection.cpp ${PROJECT_SOURCE_DIR}/src/Intersection_finder.F90 ${PROJECT_SOURCE_DIR}/src/Linked_Lists.F90)
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
FILE(GLOB FORTRAN_TEST_PARALLEL_PARTITION_A_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_parallel_partition_a.F90)
FILE(GLOB FORTRAN_TEST_PARALLEL_PARTITION_B_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_parallel_partition_b.F90)
FILE(GLOB FORTRAN_TEST_PARALLEL_PARTITION_AB_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_parallel_partition_ab.F90)
FILE(GLOB FORTRAN_TEST_PARALLEL_PARTITION_AB_OPTIMAL_TRANSFER_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_parallel_partition_ab_optimal_transfer.F90)
FILE(GLOB FORTRAN_TEST_PARALLEL_PARTITION_AB_OPTIMAL_TRANSFER_DATA_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_parallel_partition_ab_optimal_transfer_data.F90)
FILE(GLOB FORTRAN_TEST_PARALLEL_PARTITION_AB_OPTIMAL_TRANSFER_DATA_RTREE_ONCE_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_parallel_partition_ab_optimal_transfer_data_rtree_once.F90)
FILE(GLOB FORTRAN_TEST_PARALLEL_PARTITION_AB_OPTIMAL_TRANSFER_DATA_RTREE_ONCE_PAR_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_parallel_partition_ab_optimal_transfer_data_rtree_once_par.F90)
FILE(GLOB FORTRAN_TEST_PARALLEL_PARTITION_COMPLETE_AB_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_parallel_partition_complete_ab.F90)
FILE(GLOB FORTRAN_TEST_PARALLEL_PARTITION_COMPLETE_AB_SAME_ALGO_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_parallel_partition_complete_ab_same_algo.F90)
FILE(GLOB FORTRAN_TEST_PARALLEL_UNIVERSAL_NODE_NUMBERING_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_parallel_universal_node_numbering.F90)
FILE(GLOB FORTRAN_BENCHMARK_PARALLEL_PARTITION_AB_FILES ${PROJECT_SOURCE_DIR}/tests/benchmark_parallel_partition_ab.F90)
FILE(GLOB FORTRAN_BENCHMARK_PARALLEL_PARTITION_AB_SERIAL_FILES ${PROJECT_SOURCE_DIR}/tests/benchmark_parallel_partition_ab_serial.F90)
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

foreach (file ${FORTRAN_SRCFILES})
  if (${file} MATCHES "Tetrahedron_intersection")
    SET_SOURCE_FILES_PROPERTIES(${file} PROPERTIES COMPILE_FLAGS -O2)
  endif()
endforeach()

add_library (SuperMesh STATIC ${FORTRAN_SRCFILES} ${CPP_SRCFILES} ${C_SRCFILES} ${CPP_SPATIALINDEX_O_FILES} ${CPP_LIBJUDY_O_FILES})

#add_dependencies(SuperMesh spatialindex spatialindex_c libjudy)
add_dependencies(SuperMesh libjudy spatialindex spatialindex_c)

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
  test_parallel_partition_a
  test_parallel_partition_b
  test_parallel_partition_ab
  test_parallel_partition_ab_optimal_transfer
  test_parallel_partition_ab_optimal_transfer_data
  test_parallel_partition_ab_optimal_transfer_data_rtree_once
  test_parallel_partition_ab_optimal_transfer_data_rtree_once_par
  test_parallel_partition_complete_ab
  test_parallel_partition_complete_ab_same_algo
  test_parallel_universal_node_numbering
)

set (BENCHMARKTESTS
  benchmark_parallel_partition_ab
  benchmark_parallel_partition_ab_serial
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

#add_executable (test_tri_intersector EXCLUDE_FROM_ALL ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_TRI_INTERSECTOR_FILES})
add_executable (test_interval_intersector ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_INTERVAL_INTERSECTOR_FILES})
add_executable (test_tri_intersector ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_TRI_INTERSECTOR_FILES})
add_executable (test_tet_intersector ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_TET_INTERSECTOR_FILES})
add_executable (test_element_intersector ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_ELEMENT_INTERSECTOR_FILES})
add_executable (test_intersection_finder_2d ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_INTERSECTION_FINDER_2D_FILES})
add_executable (test_intersection_finder_3d ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_INTERSECTION_FINDER_3D_FILES})
add_executable (test_intersection_finder_completeness_1d ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_INTERSECTION_FINDER_COMPLETENESS_1D_FILES})
add_executable (test_intersection_finder_completeness_2d ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_INTERSECTION_FINDER_COMPLETENESS_2D_FILES})
add_executable (test_intersection_finder_completeness_3d ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_INTERSECTION_FINDER_COMPLETENESS_3D_FILES})
add_executable (test_parallel_partition_a ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_PARALLEL_PARTITION_A_FILES})
add_executable (test_parallel_partition_b ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_PARALLEL_PARTITION_B_FILES})
add_executable (test_parallel_partition_ab ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_PARALLEL_PARTITION_AB_FILES})
add_executable (test_parallel_partition_ab_optimal_transfer ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_PARALLEL_PARTITION_AB_OPTIMAL_TRANSFER_FILES})
add_executable (test_parallel_partition_ab_optimal_transfer_data ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_PARALLEL_PARTITION_AB_OPTIMAL_TRANSFER_DATA_FILES})
add_executable (test_parallel_partition_ab_optimal_transfer_data_rtree_once ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_PARALLEL_PARTITION_AB_OPTIMAL_TRANSFER_DATA_RTREE_ONCE_FILES})
add_executable (test_parallel_partition_ab_optimal_transfer_data_rtree_once_par ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_PARALLEL_PARTITION_AB_OPTIMAL_TRANSFER_DATA_RTREE_ONCE_PAR_FILES})
add_executable (test_parallel_partition_complete_ab ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_PARALLEL_PARTITION_COMPLETE_AB_FILES})
add_executable (test_parallel_partition_complete_ab_same_algo ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_PARALLEL_PARTITION_COMPLETE_AB_SAME_ALGO_FILES})
add_executable (test_parallel_universal_node_numbering ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_PARALLEL_UNIVERSAL_NODE_NUMBERING_FILES})
add_executable (benchmark_parallel_partition_ab EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_PARALLEL_PARTITION_AB_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_parallel_partition_ab_serial EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_PARALLEL_PARTITION_AB_SERIAL_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_parallel_complete EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_PARALLEL_COMPLETE_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_parallel_complete_3d EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_PARALLEL_COMPLETE_3D_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_serial EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_SERIAL_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_serial_3d EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_SERIAL_3D_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_serial_same_algo EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_SERIAL_SAME_ALGO_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_serial_same_algo_3d EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_SERIAL_SAME_ALGO_3D_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_parallel_p1_inner_product EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_PARALLEL_P1_INNER_PRODUCT_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_parallel_p1_inner_product_3d EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_PARALLEL_P1_INNER_PRODUCT_3D_FILES} ${CPP_TEST_SOURCE_FILES})
add_executable (benchmark_parallel_p2_inner_product EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_PARALLEL_P2_INNER_PRODUCT_FILES} ${CPP_TEST_SOURCE_FILES})

#add_executable (benchmark_parallel EXCLUDE_FROM_ALL ${FORTRAN_BENCHMARK_PARALLEL_INTERSECTOR_FILES} ${CPP_TEST_SOURCE_FILES})

foreach (test ${UNITTESTS})
 SET_PROPERTY(TARGET "${test}" PROPERTY COMPILE_DEFINITIONS "TESTNAME=${test}_;HAVE_MPI")
 if (test MATCHES parallel)
   ADD_TEST("test_${test}" ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${test})
 else()
   ADD_TEST("test_${test}" ${test})
 endif (test MATCHES parallel)
 SET_TESTS_PROPERTIES("test_${test}" PROPERTIES FAIL_REGULAR_EXPRESSION "Fail:")
endforeach()
foreach (benchmark ${BENCHMARKTESTS})
 SET_PROPERTY(TARGET "${benchmark}" PROPERTY COMPILE_DEFINITIONS "TESTNAME=${benchmark}_;HAVE_MPI")
endforeach()

add_custom_target(benchmark DEPENDS benchmark_parallel_partition_ab benchmark_parallel_partition_ab_serial benchmark_parallel_complete benchmark_parallel_complete_3d benchmark_serial benchmark_serial_3d benchmark_serial_same_algo benchmark_serial_same_algo_3d benchmark_parallel_p1_inner_product benchmark_parallel_p1_inner_product_3d benchmark_parallel_p2_inner_product)

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
  target_link_libraries(SuperMesh ${MPI_Fortran_LIBRARIES} rt)

  foreach (test ${UNITTESTS})
    target_link_libraries(${test} SuperMesh ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PETSC_LIBRARIES} ${LAPACK_LIBRARIES} rt)
  endforeach()
  foreach (benchmark ${BENCHMARKTESTS})
    target_link_libraries(${benchmark} SuperMesh ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PETSC_LIBRARIES} ${LAPACK_LIBRARIES} rt)
  endforeach()
else ()
  target_link_libraries(SuperMesh ${MPI_Fortran_LIBRARIES})

  foreach (test ${UNITTESTS})
    target_link_libraries(${test} SuperMesh ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PETSC_LIBRARIES} ${LAPACK_LIBRARIES})
  endforeach()
  foreach (benchmark ${BENCHMARKTESTS})
    target_link_libraries(${benchmark} SuperMesh ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PETSC_LIBRARIES} ${LAPACK_LIBRARIES})
  endforeach()
endif()

#get_property(supermesh_sources_files TARGET SuperMesh PROPERTY SOURCES)
#message(STATUS "SuperMesh SOURCES:" ${supermesh_sources_files} ".")

#if (MPI_Fortran_FOUND)
  set_target_properties(SuperMesh            PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
  set_target_properties(SuperMesh            PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${LAPACK_LINKER_FLAGS}")
  target_link_libraries(SuperMesh ${MPI_Fortran_LIBRARIES} rt)
#endif()

if(MPI_Fortran_COMPILE_FLAGS)
  set_target_properties(SuperMesh            PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
  foreach (test ${UNITTESTS})
    set_target_properties(${test} PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
  endforeach()
  foreach (benchmark ${BENCHMARKTESTS})
    set_target_properties(${benchmark} PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
  endforeach()
endif()

if(MPI_Fortran_LINK_FLAGS)
  set_target_properties(SuperMesh            PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
  foreach (test ${UNITTESTS})
    set_target_properties(${test} PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
  endforeach()
  foreach (benchmark ${BENCHMARKTESTS})
    set_target_properties(${benchmark} PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
  endforeach()
endif()

#if(CMAKE_THREAD_LIBS_INIT)
#  set_target_properties(SuperMesh            PROPERTIES LINK_FLAGS_DEBUG "${CMAKE_THREAD_LIBS_INIT}")
####  set_target_properties(test_tri_intersector PROPERTIES LINK_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
#  target_link_libraries(SuperMesh            ${CMAKE_THREAD_LIBS_INIT})
####  target_link_libraries(test_tri_intersector ${CMAKE_THREAD_LIBS_INIT})
#endif()



install (FILES "${PROJECT_BINARY_DIR}/libSuperMesh.a"                       DESTINATION lib)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_construction.mod"        DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_fldebug.mod"             DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_fldebug_parameters.mod"  DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_linked_lists.mod"        DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_intersection_finder.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_tet_intersection_module.mod"  DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_tri_intersection_module.mod"  DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/include/SuperMesh.pc"                 DESTINATION lib/pkgconfig)
