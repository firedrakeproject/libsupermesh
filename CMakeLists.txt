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
find_package(LAPACK REQUIRED)
#set(Python_ADDITIONAL_VERSIONS 2.7)
find_package(PythonInterp 2.7 REQUIRED)
find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)

if ("${PETSC_VERSION}" VERSION_LESS 3.4)
  message(FATAL_ERROR "PETSc version 3.4 is required, but not found. (Found version:" ${PETSC_VERSION} ".")
elseif ("${PETSC_VERSION}" VERSION_GREATER 3.5)
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
  set (CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

add_subdirectory(libspud)
add_subdirectory(spatialindex-1.8.0)
add_subdirectory(libjudy)
#add_subdirectory(libjudy)

#include(ExternalProject)
#ExternalProject_Add(
#  libspud
#  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libspud
#  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/libspud/configure --prefix=${CMAKE_CURRENT_SOURCE_DIR}/libspud --enable-shared --enable-static --enable-debugging
#  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/libspud
#  BUILD_COMMAND make -j
#  BUILD_IN_SOURCE 1
#  INSTALL_COMMAND echo "installing"
#)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3    -ffree-line-length-0 -fdefault-real-8 -fdefault-double-8 -ldl -lm -Wl,--export-dynamic -lnsl -lutil")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fcheck=bounds,do,mem,pointer,recursion -frecord-marker=4 -fno-realloc-lhs -finit-integer=-66666 -finit-real=nan -fno-f2c -Og -g -ffree-line-length-0 -fdefault-real-8 -fdefault-double-8 -fimplicit-none -fbacktrace -ldl -lm -Wl,--export-dynamic -lnsl -lutil")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-fbounds-check -frecord-marker=4 -fno-realloc-lhs -finit-integer=-66666 -finit-real=nan -fno-f2c -Og -g -ffree-line-length-0 -fdefault-real-8 -fdefault-double-8 -fimplicit-none -fbacktrace -ldl -lm -Wl,--export-dynamic -lnsl -lutil")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
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
endif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

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
  "${PROJECT_BINARY_DIR}/SuperMeshConfig.h"
  )

configure_file (
  "${PROJECT_SOURCE_DIR}/config/SuperMesh.pc.in"
  "${PROJECT_BINARY_DIR}/SuperMesh.pc"
  )

#SET(_WARNING_COMPILE_FLAGS "-std=f95 -Wextra -Wall -pedantic -Wextra")
SET(_WARNING_COMPILE_FLAGS " ")
SET(_WARNING_COMPILE_FLAGS_FORTRAN_DEBU "-Wextra -Wall -Wimplicit-interface -Wno-surprising")
SET(_WARNING_COMPILE_FLAGS_C "-Wextra -Wall")

SET( CMAKE_CXX_FLAGS             "${CMAKE_CXX_FLAGS} ${_WARNING_COMPILE_FLAGS_C}" )
#SET( CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE} ${_WARNING_COMPILE_FLAGS_FORTRAN}" )
SET( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${_WARNING_COMPILE_FLAGS_RELEASE}" )
SET( CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} ${_WARNING_COMPILE_FLAGS_FORTRAN_DEBUG}" )

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
include_directories ("${PROJECT_SOURCE_DIR}/spatialindex-1.8.0/include/")
include_directories ("${PROJECT_SOURCE_DIR}/spatialindex-1.8.0/")
include_directories ("${PROJECT_SOURCE_DIR}/libwm/")
include_directories ("${PROJECT_SOURCE_DIR}/libspud/")
include_directories ("${PROJECT_SOURCE_DIR}/libspud/include")
include_directories ("${PROJECT_SOURCE_DIR}/libjudy/src/")
include_directories ("${PYTHON_INCLUDE_PATH}")

set(CPP_WM4_FILES "")
set(CPP_WM4_FILES
  ${PROJECT_SOURCE_DIR}/libwm/Wm4Command.cpp
  ${PROJECT_SOURCE_DIR}/libwm/Wm4System.cpp
  ${PROJECT_SOURCE_DIR}/libwm/Wm4Vector2.cpp
  ${PROJECT_SOURCE_DIR}/libwm/Wm4Vector3.cpp
  ${PROJECT_SOURCE_DIR}/libwm/Wm4Math.cpp
)
set (CPP_LIBSPUD_O_FILES "")
set (CPP_LIBSPUD_O_FILES
  ${PROJECT_SOURCE_DIR}/libspud/fspud.o
  ${PROJECT_SOURCE_DIR}/libspud/spud_interfaces.o
  ${PROJECT_SOURCE_DIR}/libspud/spud.o
  ${PROJECT_SOURCE_DIR}/libspud/tinystr.o
  ${PROJECT_SOURCE_DIR}/libspud/tinyxmlerror.o
  ${PROJECT_SOURCE_DIR}/libspud/tinyxml.o
  ${PROJECT_SOURCE_DIR}/libspud/tinyxmlparser.o
)
foreach (test ${CPP_LIBSPUD_O_FILES})
 SET_SOURCE_FILES_PROPERTIES(${test} PROPERTIES GENERATED 1)
endforeach()

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
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/Judy1FreeArray.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/Judy1Set.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/Judy1Count.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/libcount_la-Judy1ByCount.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/Judy1Unset.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/Judy1Test.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/Judy1MallocIF.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/Judy1Tables.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/Judy1Decascade.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/Judy1Cascade.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/Judy1CreateBranch.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/Judy1InsertBranch.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/Judy1First.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/libnext_la-Judy1Next.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/libprev_la-Judy1Prev.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/libnext_la-Judy1NextEmpty.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/Judy1/libprev_la-Judy1PrevEmpty.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/JudyLFreeArray.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/JudyLIns.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/JudyLCount.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/JudyLGet.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/JudyLDel.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/libcount_la-JudyLByCount.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/JudyLTables.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/JudyLMallocIF.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/JudyLDecascade.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/JudyLCascade.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/JudyLInsertBranch.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/JudyLCreateBranch.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/JudyLFirst.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/libnext_la-JudyLNext.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/libprev_la-JudyLPrev.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/libnext_la-JudyLNextEmpty.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyL/libprev_la-JudyLPrevEmpty.o
  ${PROJECT_SOURCE_DIR}/libjudy/src/JudyCommon/JudyMalloc.o
)
foreach (test ${CPP_LIBJUDY_O_FILES})
 SET_SOURCE_FILES_PROPERTIES(${test} PROPERTIES GENERATED 1)
endforeach()

set (CPP_TEST_SOURCE_FILES "")
set (CPP_TEST_SOURCE_FILES
  ${PROJECT_SOURCE_DIR}/src/tests/test_main.cpp
)


# Create a library called "SuperMesh" which includes the source file "Element_Intersection.cpp".
# The extension is already found. Any number of sources could be listed here.
#add_library (SuperMesh STATIC ${PROJECT_SOURCE_DIR}/src/Element_Intersection.cpp ${PROJECT_SOURCE_DIR}/src/Intersection_finder.F90 ${PROJECT_SOURCE_DIR}/src/Linked_Lists.F90)
FILE(GLOB CPP_SRCFILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
FILE(GLOB C_SRCFILES ${PROJECT_SOURCE_DIR}/src/*.c)
FILE(GLOB FORTRAN_SRCFILES ${PROJECT_SOURCE_DIR}/src/*.F90)
FILE(GLOB FORTRAN_TEST_TRI_INTERSECTOR_FILES ${PROJECT_SOURCE_DIR}/src/tests/test_tri_intersector.F90)
FILE(GLOB TEST_DATA_FILES ${PROJECT_SOURCE_DIR}/src/tests/data/*)

add_library (SuperMesh STATIC ${FORTRAN_SRCFILES} ${CPP_SRCFILES} ${C_SRCFILES} ${CPP_WM4_FILES} ${CPP_LIBSPUD_O_FILES} ${CPP_SPATIALINDEX_O_FILES} ${CPP_LIBJUDY_O_FILES})

#add_dependencies(SuperMesh libspud spatialindex spatialindex_c libjudy)
add_dependencies(SuperMesh libjudy libspud spatialindex spatialindex_c)

#add_executable (test_tri_intersector EXCLUDE_FROM_ALL ${FORTRAN_TEST_TRI_INTERSECTOR_FILES} ${CPP_WM4_FILES})
add_executable (test_tri_intersector EXCLUDE_FROM_ALL ${CPP_TEST_SOURCE_FILES} ${FORTRAN_TEST_TRI_INTERSECTOR_FILES})
set_property(SOURCE ${CPP_TEST_SOURCE_FILES} PROPERTY COMPILE_DEFINITIONS "TESTNAME=test_tri_intersector_")

add_custom_command(TARGET test_tri_intersector POST_BUILD
  COMMAND mkdir -p data
  COMMAND cp "${PROJECT_SOURCE_DIR}/src/tests/data/*" data/
)

add_custom_target(test DEPENDS test_tri_intersector)
add_custom_target(distclean 
  COMMAND cd "${PROJECT_SOURCE_DIR}/libspud" && make distclean && rm -rf tmp/
  COMMAND cd "${PROJECT_SOURCE_DIR}/libjudy" && make distclean && rm -rf tmp/ && rm -rf src/libjudy-stamp/libjudy-done && rm -rf src/libjudy-stamp/
  COMMAND cd "${PROJECT_BINARY_DIR}/" && make clean
  COMMAND rm -rf data/
#  COMMAND cd "${PROJECT_SOURCE_DIR}/libjud" && make distcleanlibjudy
)

if(UNIX AND NOT APPLE)
    target_link_libraries(SuperMesh ${MPI_Fortran_LIBRARIES} rt)
    target_link_libraries(test_tri_intersector SuperMesh ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PETSC_LIBRARIES} ${LAPACK_LIBRARIES} ${PYTHON_LIBRARIES} rt)
else ()
    target_link_libraries(SuperMesh ${MPI_Fortran_LIBRARIES})
    target_link_libraries(test_tri_intersector SuperMesh ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PETSC_LIBRARIES} ${LAPACK_LIBRARIES} ${PYTHON_LIBRARIES})
endif()

get_property(supermesh_sources_files TARGET SuperMesh PROPERTY SOURCES)
#message(STATUS "CPP_LIBSPUD_O_FILES:" ${CPP_LIBSPUD_O_FILES} ".")
#message(STATUS "CPP_WM4_FILES:" ${CPP_WM4_FILES} ".")
message(STATUS "SuperMesh SOURCES:" ${supermesh_sources_files} ".")

message(STATUS "PYTHON_LIBRARIES:" ${PYTHON_LIBRARIES} ".")
message(STATUS "PYTHON_INCLUDE_DIRS:" ${PYTHON_INCLUDE_DIRS} ".")

#message(STATUS "MPI_Fortran_FOUND:" ${MPI_Fortran_FOUND} ".")
#message(STATUS "MPI_Fortran_COMPILER:" ${MPI_Fortran_COMPILER} ".")
#message(STATUS "MPI_CXX_FOUND:" ${MPI_CXX_FOUND} ".")
#message(STATUS "MPI_CXX_COMPILER:" ${MPI_CXX_COMPILER} ".")

#message(STATUS "MPI_Fortran_LIBRARIES:" ${MPI_Fortran_LIBRARIES} ".")
#message(STATUS "MPI_CXX_LIBRARIES:" ${MPI_CXX_LIBRARIES} ".")

#message(STATUS "MPI_Fortran_COMPILE_FLAGS:" ${MPI_Fortran_COMPILE_FLAGS} ".")
#message(STATUS "MPI_Fortran_LINK_FLAGS:" ${MPI_Fortran_LINK_FLAGS} ".")
#message(STATUS "MPI_CXX_COMPILE_FLAGS:" ${MPI_CXX_COMPILE_FLAGS} ".")
#message(STATUS "MPI_CXX_LINK_FLAGS:" ${MPI_CXX_LINK_FLAGS} ".")

#message(STATUS "MPI_COMPILE_FLAGS:" ${MPI_COMPILE_FLAGS} ".")

#if (MPI_Fortran_FOUND)
  set_target_properties(SuperMesh            PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
  set_target_properties(SuperMesh            PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${LAPACK_LINKER_FLAGS}")
  target_link_libraries(SuperMesh ${MPI_Fortran_LIBRARIES} rt)
#endif()

if(MPI_Fortran_COMPILE_FLAGS)
  set_target_properties(SuperMesh            PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
  set_target_properties(test_tri_intersector PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
endif()

if(MPI_Fortran_LINK_FLAGS)
  set_target_properties(SuperMesh            PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
  set_target_properties(test_tri_intersector PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
endif()

#if(CMAKE_THREAD_LIBS_INIT)
#  set_target_properties(SuperMesh            PROPERTIES LINK_FLAGS_DEBUG "${CMAKE_THREAD_LIBS_INIT}")
####  set_target_properties(test_tri_intersector PROPERTIES LINK_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
#  target_link_libraries(SuperMesh            ${CMAKE_THREAD_LIBS_INIT})
####  target_link_libraries(test_tri_intersector ${CMAKE_THREAD_LIBS_INIT})
#endif()



install (TARGETS SuperMesh DESTINATION lib)
install (FILES "${PROJECT_BINARY_DIR}/SuperMesh.pc"                         DESTINATION $ENV{PKG_CONFIG_PATH})
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_adjacency_lists.mod"     DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_construction.mod"        DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_elements.mod"            DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_element_numbering.mod"   DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_fields_allocates.mod"    DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_fields_base.mod"         DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_fields_data_types.mod"   DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_fields_manipulation.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_futils.mod"              DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_intersection_finder_module.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_linked_lists.mod"        DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_quadrature.mod"          DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_shape_functions.mod"     DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_sparse_tools.mod"        DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_tet_intersection_module.mod"  DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_transform_elements.mod" DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/libsupermesh_vector_tools.mod"       DESTINATION include)


