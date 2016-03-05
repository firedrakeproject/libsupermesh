cmake_minimum_required(VERSION 2.8)
project(LIBSUPERMESH C CXX Fortran)
enable_language(Fortran)
enable_testing()
set(LIBSUPERMESH_VERSION_MAJOR 0)
set(LIBSUPERMESH_VERSION_MINOR 9)
set(LIBSUPERMESH_VERSION_SUBMINOR 0)
set(LIBSUPERMESH_VERSION_RELEASE 0)
exec_program("git" ${LIBSUPERMESH_SOURCE_DIR} ARGS "rev-parse HEAD" OUTPUT_VARIABLE LIBSUPERMESH_GIT_REVISION RETURN_VALUE git_error)
if(NOT "${git_error}" STREQUAL 0)
  set(LIBSUPERMESH_GIT_REVISION unknown)
endif()
message(STATUS "Git revision: ${LIBSUPERMESH_GIT_REVISION}")
if(LIBSUPERMESH_VERSION_RELEASE STREQUAL 0)
  set(LIBSUPERMESH_VERSION ${LIBSUPERMESH_VERSION_MAJOR}.${LIBSUPERMESH_VERSION_MINOR}.${LIBSUPERMESH_VERSION_SUBMINOR}+)
else()
  set(LIBSUPERMESH_VERSION ${LIBSUPERMESH_VERSION_MAJOR}.${LIBSUPERMESH_VERSION_MINOR}.${LIBSUPERMESH_VERSION_SUBMINOR})
endif()

if(${LIBSUPERMESH_SOURCE_DIR} STREQUAL ${LIBSUPERMESH_BINARY_DIR})
  message(SEND_ERROR "In-source builds are not permitted")
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${LIBSUPERMESH_BINARY_DIR}/include)
include_directories(${LIBSUPERMESH_SOURCE_DIR}/include)
include_directories(${LIBSUPERMESH_BINARY_DIR}/include)
include_directories(${LIBSUPERMESH_BINARY_DIR}/include_local)

string(TOLOWER "${CMAKE_BUILD_TYPE}" lower_cmake_build_type)
if(lower_cmake_build_type STREQUAL debug)
  message(STATUS "Enabling debugging")
  set(LIBSUPERMESH_DEBUG 1)
elseif(lower_cmake_build_type STREQUAL release)
  message(STATUS "Disabling debugging")
  set(LIBSUPERMESH_DEBUG 0)
elseif(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Disabling debugging")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type. Valid values are \"Debug\" and \"Release\"." FORCE)
  set(LIBSUPERMESH_DEBUG 0)
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

option(LIBSUPERMESH_ENABLE_TIMERS "Enable internal timers" OFF)
option(LIBSUPERMESH_OVERLAP_COMPUTE_COMMS "Overlap computation and communication. Efficiency of this depends upon the specific MPI implementation." OFF)

option(LIBSUPERMESH_AUTO_COMPILER_FLAGS "Choose compiler flags automatically. This overrides any other compiler flags defined via CMake variables." ON)
if(LIBSUPERMESH_AUTO_COMPILER_FLAGS)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -pipe")
    if(BUILD_SHARED_LIBS)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
    endif()
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-integer=-66666 -finit-real=nan -fimplicit-none -std=f2008")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -DNDEBUG")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    set(CMAKE_Fortran_FLAGS "-s real64 -fpic")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  endif()

  if(CMAKE_C_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_FLAGS "-pipe -D_GNU_SOURCE -std=c90")
    if(BUILD_SHARED_LIBS)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    endif()
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -DNDEBUG")
  elseif(CMAKE_C_COMPILER_ID MATCHES "Cray")
    set(CMAKE_C_FLAGS "-fpic")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
    set(CMAKE_C_FLAGS_RELEASE "-O2")
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "-pipe -std=c++03")
    if(BUILD_SHARED_LIBS)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -DNDEBUG")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Cray")
    set(CMAKE_CXX_FLAGS "-fpic")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2")
  endif()
else()
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -std=f2008")
    if(BUILD_SHARED_LIBS)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
    endif()
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -s real64 -fpic")
  endif()

  if(CMAKE_C_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_GNU_SOURCE -std=c90")
    if(BUILD_SHARED_LIBS)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    endif()
  elseif(CMAKE_C_COMPILER_ID MATCHES "Cray")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic")
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++03")
    if(BUILD_SHARED_LIBS)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    endif()
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Cray")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")
  endif()
endif()

message(STATUS "CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}")
message(STATUS "CMAKE_Fortran_FLAGS_RELEASE=${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG}")
message(STATUS "CMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}")

include_directories(${LIBSUPERMESH_SOURCE_DIR}/spatialindex-1.8.5/include/)
include_directories(${LIBSUPERMESH_SOURCE_DIR}/spatialindex-1.8.5/src/)
add_subdirectory(spatialindex-1.8.5 EXCLUDE_FROM_ALL)
set(spatialindex_object_files
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tprtree/Index.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tprtree/TPRTree.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tprtree/Leaf.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tprtree/Statistics.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tprtree/Node.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/storagemanager/MemoryStorageManager.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/storagemanager/RandomEvictionsBuffer.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/storagemanager/Buffer.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/storagemanager/DiskStorageManager.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/mvrtree/Index.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/mvrtree/Leaf.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/mvrtree/MVRTree.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/mvrtree/Statistics.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/mvrtree/Node.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/SpatialIndexImpl.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/Point.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/MovingRegion.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/MovingPoint.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/TimeRegion.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/Region.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/LineSegment.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/spatialindex/TimePoint.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tools/Tools.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/rtree/Index.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/rtree/Leaf.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/rtree/RTree.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/rtree/Statistics.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/rtree/BulkLoader.cc${CMAKE_CXX_OUTPUT_EXTENSION}
  ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/rtree/Node.cc${CMAKE_CXX_OUTPUT_EXTENSION}
)
check_function_exists(srand48 HAVE_SRAND48)
if(NOT HAVE_SRAND48)
  list(APPEND spatialindex_object_files
    ${LIBSUPERMESH_BINARY_DIR}/spatialindex-1.8.5/src/CMakeFiles/spatialindex.dir/tools/rand48.cc${CMAKE_CXX_OUTPUT_EXTENSION})
endif()
foreach(spatialindex_object_file ${spatialindex_object_files})
  set_source_files_properties(${spatialindex_object_file} PROPERTIES GENERATED 1)
endforeach()

file(GLOB source_files ${LIBSUPERMESH_SOURCE_DIR}/src/*.F90 ${LIBSUPERMESH_SOURCE_DIR}/src/*.c ${LIBSUPERMESH_SOURCE_DIR}/src/*.cpp)
add_library(supermesh ${source_files} ${spatialindex_object_files})
add_dependencies(supermesh spatialindex)

set(test_main ${LIBSUPERMESH_SOURCE_DIR}/src/tests/test_main.cpp)
file(GLOB unittest_files ${LIBSUPERMESH_SOURCE_DIR}/src/tests/*.F90)
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

file(GLOB benchmark_files ${LIBSUPERMESH_SOURCE_DIR}/tests/*.F90)
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
  "${LIBSUPERMESH_SOURCE_DIR}/config/libsupermesh.pc.in"
  "${LIBSUPERMESH_BINARY_DIR}/config/libsupermesh.pc"
  )
configure_file (
  "${LIBSUPERMESH_SOURCE_DIR}/include/libsupermesh.h.in"
  "${LIBSUPERMESH_BINARY_DIR}/include/libsupermesh.h"
  )
configure_file (
  "${LIBSUPERMESH_SOURCE_DIR}/include/libsupermesh_configuration.h.in"
  "${LIBSUPERMESH_BINARY_DIR}/include_local/libsupermesh_configuration.h"
  )

find_package(LATEX)
if(LATEX_FOUND AND LATEX_PDFLATEX_FOUND AND LATEX_BIBTEX_FOUND)
  include(${LIBSUPERMESH_SOURCE_DIR}/cmake/UseLATEX.cmake)
  set(LATEX_OUTPUT_PATH ${LIBSUPERMESH_BINARY_DIR}/doc)
  add_latex_document(doc/manual.tex INPUTS doc/version.tex CONFIGURE doc/version.tex TARGET_NAME doc FORCE_PDF EXCLUDE_FROM_ALL)# BIBFILES doc/bibliography.bib)
endif()

install(TARGETS supermesh DESTINATION lib)
install(DIRECTORY ${LIBSUPERMESH_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ${LIBSUPERMESH_BINARY_DIR}/config/libsupermesh.pc DESTINATION lib/pkgconfig)

file(GLOB data_files ${LIBSUPERMESH_SOURCE_DIR}/src/tests/data/*.ele
                     ${LIBSUPERMESH_SOURCE_DIR}/src/tests/data/*.node)
file(COPY ${data_files} DESTINATION ${LIBSUPERMESH_BINARY_DIR}/data)
file(GLOB data_files ${LIBSUPERMESH_SOURCE_DIR}/src/tests/data/*.tar.bz2)
foreach(data_file ${data_files})
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvjf ${data_file} WORKING_DIRECTORY ${LIBSUPERMESH_BINARY_DIR}/data)
endforeach()
