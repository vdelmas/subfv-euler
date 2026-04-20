cmake_minimum_required(VERSION 3.20)

project(subfv-euler Fortran)

find_package(MPI REQUIRED)

add_subdirectory(external/subfv-core)

add_executable(euler
  src/main.F90
  src/global_data_module.F90
  src/io_module.F90
  src/euler_module.F90
)

target_link_libraries(euler
  subfvcore
  MPI::MPI_Fortran
)

target_include_directories(euler PRIVATE
  ${CMAKE_BINARY_DIR}/mod
)
