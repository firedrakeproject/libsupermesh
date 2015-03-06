# CMakeLists files in this project can
# refer to the root source directory of the project as ${SUPERMESH_SOURCE_DIR} and
# to the root binary directory of the project as ${SUPERMESH_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (SUPERMESH)
enable_language (Fortran)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# ToDo TODO todo
# Set installation prefix automatically to custom path if not explicitly specified on the command line
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "/home/ipanourg/fluidity_git/fluidity" CACHE PATH "default install path" FORCE )
endif()

# The version number.
set (SuperMesh_VERSION_MAJOR 0)
set (SuperMesh_VERSION_MINOR 0)
set (SuperMesh_VERSION_MINOR_CHANGE 1)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/include/SuperMeshConfig.h.in"
  "${PROJECT_BINARY_DIR}/SuperMeshConfig.h"
  )

configure_file (
  "${PROJECT_SOURCE_DIR}/config/SuperMesh.pc.in"
  "${PROJECT_BINARY_DIR}/SuperMesh.pc"
  )

#SET(_WARNING_COMPILE_FLAGS "-std=f95 -Wextra -Wall -pedantic -Wextra")
#SET(_WARNING_COMPILE_FLAGS "-Wextra -Wall -pedantic -Wextra")
SET(_WARNING_COMPILE_FLAGS "-Wall")

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${_WARNING_COMPILE_FLAGS}" )
SET( CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} ${_WARNING_COMPILE_FLAGS}" )
SET( CMAKE_Fortran_FLAGS_RELEASE  "${CMAKE_Fortran_FLAGS_RELEASE} ${_WARNING_COMPILE_FLAGS}" )
SET( CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} ${_WARNING_COMPILE_FLAGS}" )


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
include_directories ("${PROJECT_BINARY_DIR}")
include_directories ("${PROJECT_SOURCE_DIR}/spatialindex-1.8.0/include")


# Create a library called "SuperMesh" which includes the source file "Element_Intersection.cpp".
# The extension is already found. Any number of sources could be listed here.
#add_library (SuperMesh STATIC ${PROJECT_SOURCE_DIR}/src/Element_Intersection.cpp ${PROJECT_SOURCE_DIR}/src/Intersection_finder.F90 ${PROJECT_SOURCE_DIR}/src/Linked_Lists.F90)
FILE(GLOB FORTRAN_SRCFILES ${PROJECT_SOURCE_DIR}/src/*.F90)
FILE(GLOB CPP_SRCFILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
#SET( many_srcs "${PROJECT_SOURCE_DIR}/src/Supermesh.F90")

add_library (SuperMesh STATIC ${CPP_SRCFILES} ${FORTRAN_SRCFILES})


install (TARGETS SuperMesh DESTINATION lib)
install (FILES "${PROJECT_BINARY_DIR}/SuperMesh.pc" DESTINATION $ENV{PKG_CONFIG_PATH})
