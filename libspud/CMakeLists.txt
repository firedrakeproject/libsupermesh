cmake_minimum_required(VERSION 2.8)

if (SUBPROJECTS_profilling)
  set(local_profile "--enable-profiling")
endif (SUBPROJECTS_profilling)

if (SUBPROJECTS_debugging MATCHES YES)
  set(local_debugging "--enable-debugging")
endif (SUBPROJECTS_debugging MATCHES YES)

string(REPLACE "-fdefault-real-8" "" local_SUBPROJECTS_C_FLAGS ${SUBPROJECTS_C_FLAGS})
string(REPLACE "-fdefault-real-8" "" local_SUBPROJECTS_CXX_FLAGS ${SUBPROJECTS_CXX_FLAGS})
string(REPLACE "-fdefault-real-8" "" local_SUBPROJECTS_Fortran_FLAGS ${SUBPROJECTS_Fortran_FLAGS})
string(REPLACE "-fdefault-double-8" "" local_SUBPROJECTS_C_FLAGS ${local_SUBPROJECTS_C_FLAGS})
string(REPLACE "-fdefault-double-8" "" local_SUBPROJECTS_CXX_FLAGS ${local_SUBPROJECTS_CXX_FLAGS})
string(REPLACE "-fdefault-double-8" "" local_SUBPROJECTS_Fortran_FLAGS ${local_SUBPROJECTS_Fortran_FLAGS})

#message("libspud local_SUBPROJECTS_C_FLAGS:       ${local_SUBPROJECTS_C_FLAGS}")
#message("libspud local_SUBPROJECTS_CXX_FLAGS:     ${local_SUBPROJECTS_CXX_FLAGS}")
#message("libspud local_SUBPROJECTS_Fortran_FLAGS: ${local_SUBPROJECTS_Fortran_FLAGS}")

include(ExternalProject)
ExternalProject_Add(
  libspud
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
#  CONFIGURE_COMMAND FCFLAGS="-fPIC" CXX=mpicxx CC=mpicc FC=mpif90 ${CMAKE_CURRENT_SOURCE_DIR}/configure --prefix=${PROJECT_BINARY_DIR} --enable-shared --enable-static --enable-debugging
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/configure "CFLAGS=${local_SUBPROJECTS_C_FLAGS}" "CXXFLAGS=${local_SUBPROJECTS_CXX_FLAGS}" "FFLAGS=${local_SUBPROJECTS_Fortran_FLAGS}" "FCFLAGS=${local_SUBPROJECTS_Fortran_FLAGS}" --prefix=${PROJECT_BINARY_DIR} --enable-shared --enable-static ${local_profile} ${local_debugging}
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}
  BUILD_COMMAND make 
#  BUILD_DIR ${PROJECT_BINARY_DIR}
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND echo "installing"
)

add_custom_target(distcleanlibspud
#  COMMAND "$(CMAKE)" --build "${CMAKE_BINARY_DIR}" clean
  COMMAND cd "${PROJECT_SOURCE_DIR}" && rm -rf tmp/ && make distclean 
)


