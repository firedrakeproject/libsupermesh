cmake_minimum_required(VERSION 3.18)

if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Build type. Valid values are \"Debug\" and \"Release\".")
else()
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type. Valid values are \"Debug\" and \"Release\".")
endif()
string(TOLOWER "${CMAKE_BUILD_TYPE}" lower_cmake_build_type)
if(lower_cmake_build_type STREQUAL debug)
  message(STATUS "Enabling debugging")
  set(LIBSUPERMESH_DEBUG 1)
elseif(lower_cmake_build_type STREQUAL release)
  message(STATUS "Disabling debugging")
  set(LIBSUPERMESH_DEBUG 0)
else()
   message(WARNING "Unrecognised build type")
endif()

project(LIBSUPERMESH C CXX Fortran)
enable_language(Fortran)
enable_testing()
set(LIBSUPERMESH_VERSION_MAJOR 1)
set(LIBSUPERMESH_VERSION_MINOR 0)
set(LIBSUPERMESH_VERSION_PATCH 1)
set(LIBSUPERMESH_VERSION_RELEASE 0)
if(LIBSUPERMESH_VERSION_RELEASE)
  set(LIBSUPERMESH_VERSION ${LIBSUPERMESH_VERSION_MAJOR}.${LIBSUPERMESH_VERSION_MINOR}.${LIBSUPERMESH_VERSION_PATCH})
else()
  set(LIBSUPERMESH_VERSION ${LIBSUPERMESH_VERSION_MAJOR}.${LIBSUPERMESH_VERSION_MINOR}.${LIBSUPERMESH_VERSION_PATCH}+)
endif()
exec_program("git" ${LIBSUPERMESH_SOURCE_DIR} ARGS "rev-parse HEAD" OUTPUT_VARIABLE LIBSUPERMESH_GIT_REVISION RETURN_VALUE git_error)
if(NOT ${git_error} STREQUAL 0)
  set(LIBSUPERMESH_GIT_REVISION unknown)
endif()
message(STATUS "Git revision: ${LIBSUPERMESH_GIT_REVISION}")

if(${LIBSUPERMESH_SOURCE_DIR} STREQUAL ${LIBSUPERMESH_BINARY_DIR})
  message(SEND_ERROR "In-source builds are not permitted")
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${LIBSUPERMESH_BINARY_DIR}/include)
include_directories(${LIBSUPERMESH_SOURCE_DIR}/include)
include_directories(${LIBSUPERMESH_BINARY_DIR}/include)
include_directories(${LIBSUPERMESH_BINARY_DIR}/include_local)
if("${SKBUILD}" STREQUAL 2)
  option(BUILD_SHARED_LIBS "Build shared libraries" ON)
else()
  option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
endif()

find_package(MPI REQUIRED)
set(link_libraries ${MPI_Fortran_LIBRARIES})
include_directories(${MPI_Fortran_INCLUDE_PATH})

option(LIBSUPERMESH_ENABLE_JUDY "Enable use of the Judy library" OFF)
if(LIBSUPERMESH_ENABLE_JUDY)
  # See https://cmake.org/Wiki/CMake:How_To_Find_Libraries
  find_library(JUDY_LIBRARY NAMES Judy)
  if(JUDY_LIBRARY STREQUAL JUDY_LIBRARY-NOTFOUND)
    message(FATAL_ERROR "Judy library not found")
  endif()
  list(APPEND link_libraries ${JUDY_LIBRARY})
  find_path(JUDY_INCLUDE_DIR NAMES Judy.h)
  if(JUDY_INCLUDE_DIR STREQUAL JUDY_INCLUDE_DIR-NOTFOUND)
    message(FATAL_ERROR "Judy header file not found")
  endif()
  include_directories(${JUDY_INCLUDE_DIR})
endif()

find_package(Backtrace)
if(Backtrace_FOUND)
  set(link_libraries ${link_libraries} ${Backtrace_LIBRARIES})
  include_directories(${Backtrace_INCLUDE_DIRS})
endif()

option(LIBSUPERMESH_ENABLE_TIMERS "Enable internal timers" OFF)
option(LIBSUPERMESH_OVERLAP_COMPUTE_COMMS "Overlap computation and communication. Efficiency of this depends upon the specific MPI implementation." OFF)

option(LIBSUPERMESH_DOUBLE_PRECISION "Build with double precision Fortran reals" ON)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -std=f2008")
  if(BUILD_SHARED_LIBS)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -s -fpic")
endif()

if(CMAKE_C_COMPILER_ID MATCHES "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
  if(BUILD_SHARED_LIBS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  endif()
elseif(CMAKE_C_COMPILER_ID MATCHES "Cray")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  if(BUILD_SHARED_LIBS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Cray")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")
endif()

# Use the MPI compilers to compile
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# Print out all the MPI configuration info
message(STATUS "")
message(STATUS "MPI_Fortran_COMPILER=${MPI_Fortran_COMPILER}")
message(STATUS "CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}")
message(STATUS "CMAKE_Fortran_FLAGS_RELEASE=${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "")
message(STATUS "MPI_C_COMPILER=${MPI_C_COMPILER}")
message(STATUS "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG}")
message(STATUS "CMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE}")
message(STATUS "")
message(STATUS "MPI_CXX_COMPILER=${MPI_CXX_COMPILER}")
message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "")
message(STATUS "MPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}")
message(STATUS "")

# Python is not required, but we want to add pip packages to spatialindex search path
set(Python_FIND_VIRTUALENV "FIRST")
find_package(Python COMPONENTS Interpreter Development.Module)
if(Python_FOUND)
    message(STATUS "Python_VERSION_MAJOR " ${Python_VERSION_MAJOR})
    message(STATUS "Python_VERSION_MINOR " ${Python_VERSION_MINOR})
endif()

# Find the spatialindex library
find_library(spatialindexlib
             NAMES spatialindex
             PATHS "/usr/local/lib"
              "/usr/lib"
              "$ENV{VIRTUAL_ENV}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/Rtree.libs"
              "$ENV{VIRTUAL_ENV}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/rtree/lib"
            )

if((NOT spatialindexlib) AND (Python_FOUND))
  message(STATUS "Spatialindex not found, attempting to use the shared library provided by the Rtree Python package")
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import rtree; print(rtree.core.rt._name, end='')"
                  OUTPUT_VARIABLE spatialindexlib
                  RESULT_VARIABLE PythonCommandReturnCode
                  )
  if(PythonCommandReturnCode)
    # Forces build to fail if we still cannot find spatialindex
    unset(spatialindexlib)
    find_library(spatialindexlib NAMES spatialindex REQUIRED)
  else()
    # Get the headers
    execute_process(COMMAND ${Python_EXECUTABLE} -c "import rtree; print(rtree.finder.get_include(), end='')"
                  OUTPUT_VARIABLE spatialindexheaders
                  )
  endif()
endif()

if(spatialindexlib)
  message(STATUS "Found libspatialindex shared library: " ${spatialindexlib})
  set(link_libraries ${link_libraries} ${spatialindexlib})
  if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    # Check whether C++11 strings are in use, I only know how to do this on Linux
    execute_process(COMMAND objdump -T ${spatialindexlib}
                    COMMAND grep -i setproperty
                    OUTPUT_VARIABLE setpropertystring
                    )
    message(STATUS "Set property object signature: " ${setpropertystring})
    if(NOT ${setpropertystring} MATCHES .*__cxx1112basic_string.*)
      message(STATUS "libspatialindex shared library does not use C++11 strings. Adding -D_GLIBCXX_USE_CXX11_ABI=0 to CMAKE_CXX_FLAGS")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
    endif()
  endif()
endif()

# Find the spatialindex headers
find_path(spatialindexheaders
          NAMES spatialindex/SpatialIndex.h
          REQUIRED
          PATHS /usr/local/include /usr/include $ENV{VIRTUAL_ENV}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/rtree/include spatialindexheaders
          )
if(spatialindexheaders)
  message(STATUS "Found libspatialindex headers: " ${spatialindexheaders})
  include_directories(${spatialindexheaders})
endif()

file(GLOB source_files ${LIBSUPERMESH_SOURCE_DIR}/src/*.F90 ${LIBSUPERMESH_SOURCE_DIR}/src/*.c ${LIBSUPERMESH_SOURCE_DIR}/src/*.cpp)
if("${SKBUILD}" STREQUAL 2)
  python_add_library(supermesh MODULE WITH_SOABI ${source_files})
else()
  add_library(supermesh ${source_files})
endif()

target_link_libraries(supermesh PUBLIC ${link_libraries})
#~ set_target_properties(supermesh PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
set_target_properties(supermesh PROPERTIES INSTALL_RPATH "$ORIGIN/../Rtree.libs")
set_target_properties(supermesh PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
if("${SKBUILD}" STREQUAL 2)
  # use, i.e. don't skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

# Only build tests if not built by scikit-build-core
if(NOT "${SKBUILD}" STREQUAL 2)
  set(test_main ${LIBSUPERMESH_SOURCE_DIR}/src/tests/test_main.cpp)
  file(GLOB unittest_files ${LIBSUPERMESH_SOURCE_DIR}/src/tests/*.F90)
  set(unittests "")
  foreach(unittest_file ${unittest_files})
    get_filename_component(unittest ${unittest_file} NAME_WE)
    add_executable(${unittest} ${unittest_file} ${test_main})
    list(APPEND unittests ${unittest})
    set_property(TARGET ${unittest} PROPERTY COMPILE_DEFINITIONS "TESTNAME=${unittest}")
    if(unittest MATCHES parallel)
      add_test(${unittest} ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./${unittest} ${MPIEXEC_POSTFLAGS})
    else()
      add_test(${unittest} ${unittest})
    endif()
    set_tests_properties(${unittest} PROPERTIES FAIL_REGULAR_EXPRESSION "Fail:")
  endforeach()

  set(test_link_libraries supermesh ${link_libraries} ${MPI_CXX_LIBRARIES})
  foreach(test ${unittests})
    target_link_libraries(${test} ${test_link_libraries})
  endforeach()

  foreach(unittest ${unittests})
    set_target_properties(${unittest} PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
  endforeach()
endif()

configure_file (
  "${LIBSUPERMESH_SOURCE_DIR}/config/libsupermesh.pc.in"
  "${LIBSUPERMESH_BINARY_DIR}/config/libsupermesh.pc"
  )
configure_file (
  "${LIBSUPERMESH_SOURCE_DIR}/include/libsupermesh.h.in"
  "${LIBSUPERMESH_BINARY_DIR}/include/libsupermesh.h"
  )
configure_file (
  "${LIBSUPERMESH_SOURCE_DIR}/include/libsupermesh-c.h"
  "${LIBSUPERMESH_BINARY_DIR}/include/libsupermesh-c.h"
  )
configure_file (
  "${LIBSUPERMESH_SOURCE_DIR}/include/libsupermesh_configuration.h.in"
  "${LIBSUPERMESH_BINARY_DIR}/include_local/libsupermesh_configuration.h"
  )

option(ENABLE_DOCS "Enable building of documentation" OFF)
if(ENABLE_DOCS)
  find_package(ImageMagick COMPONENTS convert REQUIRED)
  find_package(LATEX)
  if(NOT LATEX_FOUND OR NOT LATEX_PDFLATEX_FOUND OR NOT LATEX_BIBTEX_FOUND)
    message(SEND_ERROR "Required LaTeX components not found")
  endif()

  include(${LIBSUPERMESH_SOURCE_DIR}/cmake/UseLATEX.cmake)
  set(LATEX_OUTPUT_PATH ${LIBSUPERMESH_BINARY_DIR}/doc)
  add_latex_document(doc/manual.tex INPUTS doc/version.tex
    CONFIGURE doc/version.tex IMAGES doc/hex_ordering.pdf TARGET_NAME doc
    FORCE_PDF EXCLUDE_FROM_ALL BIBFILES doc/bibliography.bib)
endif()

if("${SKBUILD}" STREQUAL 2)
  install(TARGETS supermesh DESTINATION supermesh)
  install(DIRECTORY ${LIBSUPERMESH_BINARY_DIR}/include DESTINATION supermesh)
else()
  install(TARGETS supermesh DESTINATION lib)
  install(DIRECTORY ${LIBSUPERMESH_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(FILES ${LIBSUPERMESH_BINARY_DIR}/config/libsupermesh.pc DESTINATION lib/pkgconfig)
endif()

file(GLOB data_files ${LIBSUPERMESH_SOURCE_DIR}/src/tests/data/*.ele
                     ${LIBSUPERMESH_SOURCE_DIR}/src/tests/data/*.node)
file(COPY ${data_files} DESTINATION ${LIBSUPERMESH_BINARY_DIR}/data)
file(GLOB data_files ${LIBSUPERMESH_SOURCE_DIR}/src/tests/data/*.tar.bz2)
foreach(data_file ${data_files})
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xjf ${data_file} WORKING_DIRECTORY ${LIBSUPERMESH_BINARY_DIR}/data)
endforeach()
