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


# Directory structures
#
#	src/
#	|--		allvar.h 			# Farther Header
#	|--		types.h 			# Data Types
#	|-- 	autovec.h 			# Vector Customizing
#	|
#	|
#	|--		main.cxx 			# Main executation file
#	|--		config.cxx 			# Read Configuration
#	|--		utilities.cxx 		# Some utilities ftns
#	|-- 	makebr.cxx 			# make branch main loop
		

cmake_minimum_required(VERSION 3.20)
project(veluga_ctree LANGUAGES CXX)

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
find_package(HDF5 REQUIRED COMPONENTS C CXX)
target_link_libraries(veluga_ctree PRIVATE HDF5::HDF5)

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

# 옵션: 약간의 최적화와 경고
if (MSVC)
  target_compile_options(veluga_ctree PRIVATE /O2 /W4)
else()
  target_compile_options(veluga_ctree PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()





