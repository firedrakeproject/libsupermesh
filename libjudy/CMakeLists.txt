cmake_minimum_required(VERSION 2.8)

include(ExternalProject)
ExternalProject_Add(
  libjudy
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
#  CONFIGURE_COMMAND CXX=mpicxx CC=mpicc FC=mpif90 ${CMAKE_CURRENT_SOURCE_DIR}/configure --prefix=${PROJECT_BINARY_DIR} --enable-shared --enable-static --enable-debug
# errors  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/configure "CFLAGS=-O3 -march=native" --prefix=${PROJECT_BINARY_DIR} --enable-shared --enable-static
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/configure --prefix=${PROJECT_BINARY_DIR} --enable-shared --enable-static
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}
  BUILD_COMMAND make
#  BUILD_DIR ${PROJECT_BINARY_DIR}
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND echo "installing"
)

add_custom_target(distcleanlibjudy
#  COMMAND "$(CMAKE)" --build "${CMAKE_BINARY_DIR}" clean
  COMMAND cd "${PROJECT_SOURCE_DIR}" && make distclean 
)


