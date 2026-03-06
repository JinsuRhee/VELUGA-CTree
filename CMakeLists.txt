# CMakeLists.txt
#
# VELUGA - Ctree
# Written by Jinsu Rhee (jinsu.rhee@gmail.com)
#
# Ver 0.01
#
#-----------------
# Code Structure
#	ctree/
#		├─ CMakeLists.txt
#		└─ src/
#			├─ main.cxx                 # wrapper code
#			├─ global/                  # global headers and utilities
#			│ 	├─ allvar.h
#			│ 	├─ io.h 				# Manage IO
#			│ 	├─ io.cxx 				# IO
#			│ 	└─ config.cxx 			# Read configuration file 
#			├─ makebr/                  # makebr related
#			│ 	├─ makebr.h 			# Makebr header
#			│ 	└─ makebr.cxx 			# Branch making
#			├─ complete_tree/           # complete_tree related
#			│	├─ ctree.h 				# Complete Tree Header
#			│	└─ ctree.cxx 			# Complete Tree part
#			└─ utilities/               # utilities
#				├─ utilities.h
#				└─ utilities.cxx
#-----------------

#-----------------
# User-selectable type options
# Defaults:
#   Catalog ID    = int32
#   Catalog Snap  = int32
#   ParticleID    = int64
#   Merit         = double
#-----------------


cmake_minimum_required(VERSION 3.20)
project(veluga_ctree LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#-----
# Make executable
#-----
add_executable(veluga_ctree
  src/main.cxx
  src/global/config.cxx
  src/global/io.cxx
  src/makebr/makebr.cxx
  src/complete_tree/ctree.cxx
  src/utilities/utilities.cxx
#  src/config.cxx
#  src/utilities.cxx
#  src/makebr.cxx
)


#-----
# Data Type
#-----
set(CTREE_ID_TYPE      "int32"   CACHE STRING "Halo/Node ID type: int32|int64")
set(CTREE_SNAP_TYPE    "int32"   CACHE STRING "Snapshot number type: int32|int64")    
set(CTREE_PARTID_TYPE  "int64"   CACHE STRING "Particle ID type: int32|int64")
set(CTREE_MERIT_TYPE   "float64" CACHE STRING "Merit type: float32|float64")




set_property(CACHE CTREE_ID_TYPE     PROPERTY STRINGS int32 int64)
set_property(CACHE CTREE_SNAP_TYPE   PROPERTY STRINGS int32 int64)
set_property(CACHE CTREE_PARTID_TYPE PROPERTY STRINGS int32 int64)
set_property(CACHE CTREE_MERIT_TYPE  PROPERTY STRINGS float32 float64)

# ---- ID type
if(CTREE_ID_TYPE STREQUAL "int32")
  target_compile_definitions(veluga_ctree PRIVATE CTREE_ID_INT32=1)
elseif(CTREE_ID_TYPE STREQUAL "int64")
  target_compile_definitions(veluga_ctree PRIVATE CTREE_ID_INT64=1)
else()
  message(FATAL_ERROR "Invalid CTREE_ID_TYPE='${CTREE_ID_TYPE}' (use int32|int64)")
endif()

# ---- SNAP type
if(CTREE_SNAP_TYPE STREQUAL "int32")
  target_compile_definitions(veluga_ctree PRIVATE CTREE_SNAP_INT32=1)
elseif(CTREE_SNAP_TYPE STREQUAL "int64")
  target_compile_definitions(veluga_ctree PRIVATE CTREE_SNAP_INT64=1)
else()
  message(FATAL_ERROR "Invalid CTREE_SNAP_TYPE='${CTREE_SNAP_TYPE}' (use int32|int64)")
endif()

# ---- Particle ID type
if(CTREE_PARTID_TYPE STREQUAL "int32")
  target_compile_definitions(veluga_ctree PRIVATE CTREE_PARTID_INT32=1)
elseif(CTREE_PARTID_TYPE STREQUAL "int64")
  target_compile_definitions(veluga_ctree PRIVATE CTREE_PARTID_INT64=1)
else()
  message(FATAL_ERROR "Invalid CTREE_PARTID_TYPE='${CTREE_PARTID_TYPE}' (use int32|int64)")
endif()

# ---- Merit type
if(CTREE_MERIT_TYPE STREQUAL "float32")
  target_compile_definitions(veluga_ctree PRIVATE CTREE_MERIT_FLOAT32=1)
elseif(CTREE_MERIT_TYPE STREQUAL "float64")
  target_compile_definitions(veluga_ctree PRIVATE CTREE_MERIT_FLOAT64=1)
else()
  message(FATAL_ERROR "Invalid CTREE_MERIT_TYPE='${CTREE_MERIT_TYPE}' (use float32|float64)")
endif()

#-----
# MPI
#-----
find_package(MPI REQUIRED)
if (MPI_FOUND)
  add_compile_definitions(CTREE_USE_MPI)
  target_link_libraries(veluga_ctree PRIVATE MPI::MPI_CXX)
endif()

#-----
# OMP
#-----
find_package(OpenMP REQUIRED)
if (OpenMP_FOUND)
  add_compile_definitions(CTREE_USE_OMP)
  target_link_libraries(veluga_ctree PRIVATE OpenMP::OpenMP_CXX)
endif()

#-----
# HDF5
#-----
#find_package(HDF5 REQUIRED COMPONENTS C CXX)
#target_link_libraries(veluga_ctree PRIVATE HDF5::HDF5)
set(HDF5_PREFER_PARALLEL TRUE)

find_package(HDF5 CONFIG QUIET COMPONENTS C CXX)

if (TARGET HDF5::HDF5)
  message(STATUS "Using HDF5 imported target")
  target_link_libraries(veluga_ctree PRIVATE HDF5::HDF5)
  target_include_directories(veluga_ctree PRIVATE ${HDF5_INCLUDE_DIRS})
else()
  find_package(HDF5 REQUIRED COMPONENTS C CXX)
  message(STATUS "Using HDF5 variables (no imported target)")
  target_link_libraries(veluga_ctree PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
  target_include_directories(veluga_ctree PRIVATE ${HDF5_INCLUDE_DIRS})
endif()

#-----
# GCC > 9
#-----
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    message(FATAL_ERROR "GCC >= 9 is required for std::filesystem. Found ${CMAKE_CXX_COMPILER_VERSION}.")
  endif()
endif()



# Add path
target_include_directories(veluga_ctree PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${HDF5_INCLUDE_DIRS}
)

# Link libraries
#target_link_libraries(veluga_ctree
#  PUBLIC MPI::MPI_CXX
#  PUBLIC OpenMP::OpenMP_CXX
#  PUBLIC ${HDF5_CXX_LIBRARIES}
#)

target_compile_definitions(veluga_ctree PRIVATE PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# Optimization & Warning
if (MSVC)
  target_compile_options(veluga_ctree PRIVATE /O2 /W4)
else()
  target_compile_options(veluga_ctree PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()





