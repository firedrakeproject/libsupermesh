# CMakeLists files in this project can
# refer to the root source directory of the project as ${SUPERMESH_SOURCE_DIR} and
# to the root binary directory of the project as ${SUPERMESH_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (SUPERMESH C CXX Fortran)
enable_language (Fortran)

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")

FIND_PACKAGE(PETSc REQUIRED)
find_package(MPI REQUIRED)
find_package(Threads)

include_directories(${MPI_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
ADD_DEFINITIONS(${PETSC_DEFINITIONS})

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
#  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -ffree-line-length-0 -pthread -lmpi_f90 -lmpi_f77 -lmpi -ldl -lm -Wl,--export-dynamic -lrt -lnsl -lutil -lm -ldl")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -ffree-line-length-0 -pthread -lmpi_f90 -lmpi_f77 -lmpi -ldl -lm -Wl,--export-dynamic -lrt -lnsl -lutil -lm -ldl")


  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3    -ffree-line-length-0 -fdefault-real-8 -fdefault-double-8 -ldl -lm -Wl,--export-dynamic -lnsl -lutil")
  set (CMAKE_Fortran_FLAGS_DEBUG   "                   -fno-f2c -O0 -g -ffree-line-length-0 -fdefault-real-8 -fdefault-double-8 -ldl -lm -Wl,--export-dynamic -lnsl -lutil")
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
SET( CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE} ${_WARNING_COMPILE_FLAGS}" )
SET( CMAKE_Fortran_FLAGS_RELEASE  "${CMAKE_Fortran_FLAGS_RELEASE} ${_WARNING_COMPILE_FLAGS}" )
SET( CMAKE_Fortran_FLAGS_DEBUG    "${CMAKE_Fortran_FLAGS_DEBUG}   ${_WARNING_COMPILE_FLAGS}" )

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
FILE(GLOB CPP_SRCFILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
FILE(GLOB FORTRAN_SRCFILES ${PROJECT_SOURCE_DIR}/src/*.F90)
#SET( many_srcs "${PROJECT_SOURCE_DIR}/src/Supermesh.F90")

add_library (SuperMesh STATIC ${FORTRAN_SRCFILES} ${CPP_SRCFILES})

if(UNIX AND NOT APPLE)
    target_link_libraries(SuperMesh rt)
endif()

if(MPI_COMPILE_FLAGS)
  set_target_properties(SuperMesh PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(SuperMesh PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

if(CMAKE_THREAD_LIBS_INIT)
  set_target_properties(SuperMesh PROPERTIES
    LINK_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
  target_link_libraries (SuperMesh ${CMAKE_THREAD_LIBS_INIT})
endif()

install (TARGETS SuperMesh DESTINATION lib)
install (FILES "${PROJECT_BINARY_DIR}/SuperMesh.pc" DESTINATION $ENV{PKG_CONFIG_PATH})
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_tet_intersection_module.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_elements.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_vector_tools.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_fields_data_types.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_fields_base.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_fields_allocates.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_fields_manipulation.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_transform_elements.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_linked_lists.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_futils.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_element_numbering.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_construction.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_quadrature.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_shape_functions.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_intersection_finder_module.mod" DESTINATION include)



