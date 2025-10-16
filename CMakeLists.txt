cmake_minimum_required(VERSION 3.10)
project(veluga_ctree LANGUAGES C CXX)

# Standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Options
option(ENABLE_ASAN "Build with AddressSanitizer/UBSan" ON)

# Executable
add_executable(veluga_ctree
  src/main.cxx
  src/global/config.cxx
  src/global/io.cxx
  src/makebr/makebr.cxx
  src/complete_tree/ctree.cxx
  src/utilities/utilities.cxx
)

# ----------------
# MPI (Intel MPI)
# ----------------
find_package(MPI REQUIRED)
if (MPI_FOUND)
  add_compile_definitions(CTREE_USE_MPI)
  target_link_libraries(veluga_ctree PRIVATE MPI::MPI_CXX)
endif()

# ----------------
# OpenMP
# ----------------
find_package(OpenMP REQUIRED COMPONENTS CXX)
if (OpenMP_FOUND)
  add_compile_definitions(CTREE_USE_OMP)
  target_link_libraries(veluga_ctree PRIVATE OpenMP::OpenMP_CXX)
endif()

# ----------------
# HDF5 (Config target or variable)
# ----------------
find_package(HDF5 REQUIRED COMPONENTS C CXX HL)

if (TARGET HDF5::HDF5)
  # For HDF target
  target_link_libraries(veluga_ctree PRIVATE HDF5::HDF5)
elseif (TARGET hdf5::hdf5)
  if (TARGET hdf5::hdf5_cpp)
    target_link_libraries(veluga_ctree PRIVATE hdf5::hdf5_cpp)
  endif()
  if (TARGET hdf5::hdf5_hl)
    target_link_libraries(veluga_ctree PRIVATE hdf5::hdf5_hl)
  endif()
  if (TARGET hdf5::hdf5_hl_cpp)
    target_link_libraries(veluga_ctree PRIVATE hdf5::hdf5_hl_cpp)
  endif()
  target_link_libraries(veluga_ctree PRIVATE hdf5::hdf5)
else()
  # based on find
  target_include_directories(veluga_ctree PRIVATE ${HDF5_INCLUDE_DIRS})
  set(_hdf5_all_libs
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_CXX_HL_LIBRARIES}
  )
  list(REMOVE_DUPLICATES _hdf5_all_libs)
  target_link_libraries(veluga_ctree PRIVATE ${_hdf5_all_libs})
endif()

# ----------------
# GCC version check
# ----------------
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    message(FATAL_ERROR "GCC >= 9 is required. Found ${CMAKE_CXX_COMPILER_VERSION}.")
  endif()
endif()

# ----------------
# optimiaztion & ASAN
# ----------------
if (MSVC)
  target_compile_options(veluga_ctree PRIVATE /W4)
else()
  target_compile_options(veluga_ctree PRIVATE -Wall -Wextra -Wpedantic)
endif()

if (ENABLE_ASAN AND NOT MSVC)
  # optimization crash block
  target_compile_options(veluga_ctree PRIVATE -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined)
  target_link_options(veluga_ctree PRIVATE -fsanitize=address,undefined)
else()
  if (NOT MSVC)
    target_compile_options(veluga_ctree PRIVATE -O3)
  else()
    target_compile_options(veluga_ctree PRIVATE /O2)
  endif()
endif()

# Includes
target_include_directories(veluga_ctree PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Misc
target_compile_definitions(veluga_ctree PRIVATE PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

