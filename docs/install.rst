Install Instructions
====================

This document describes how to build **CTree** using CMake.

Prerequisites
--------------

To build **CTree**, the following software and libraries are required.

	- A C++ compiler (**C++17 or later**)
	- GCC (>9)
	- Clang (with C++17 support)
	- CMake (**CMake 3.20 or later**)
	- MPI (OpenMPI, MPICH, Intel MPI)
	- OpenMP
	- HDF5 (**HDF5 with C and C++ interfaces is required**)
	- Operating System (Linux; tested on HPC environments)


Build from Source
-----------------

First, clone the repository:

.. code-block:: bash

   git clone https://github.com/JinsuRhee/CTree.git
   cd CTree

Create a build directory:

.. code-block:: bash

   mkdir build
   cd build

Configure the project with CMake:

.. code-block:: bash

   cmake ../

Build the code:

.. code-block:: bash

   cmake --build .


.. Note ::

	By default, **CTree** is built with the following data types:

	- **Catalog ID**: 32-bit (positive) integer (``int32``)
	- **Snapshot number**: 32-bit integer (``int32``)
	- **Particle ID**: 64-bit integer (``int64``)
	- **Merit value**: 64-bit floating point (``double``)

	These defaults are suitable for most use cases and require no additional
	build options.

	If different data types are required, they can be selected **at compile time**
	via CMake build options. For example:

	.. code-block:: bash

	   cmake \
	     -DCTREE_ID_TYPE=int64 \
	     -DCTREE_SNAP_TYPE=int32 \
	     -DCTREE_PARTID_TYPE=int64 \
	     -DCTREE_MERIT_TYPE=float32 \
	     ../

	   cmake --build build

	Available options are:

	- ``CTREE_ID_TYPE``: ``int32`` | ``int64``
	- ``CTREE_SNAP_TYPE``: ``int32`` | ``int64``
	- ``CTREE_PARTID_TYPE``: ``int32`` | ``int64``
	- ``CTREE_MERIT_TYPE``: ``float32`` | ``float64``

	All type selections are resolved at compile time and applied consistently
	throughout the code base. Changing these options requires a full rebuild
	of the project.


Troubleshooting
---------------
	- If CMake cannot find a dependency, make sure it is installed and visible in your environment.
	- Try removing the build directory and configuring again if you encounter unexpected errors.


