Build Instructions
==================

This document describes how to build the code using CMake.

Prerequisites
--------------

To build **VELUGA-CTree**, the following software and libraries are required.

Compiler
~~~~~~~~

- A C++ compiler with **C++17 support**
- GCC (>9)
- Clang (with C++17 support)

CMake
~~~~~
- **CMake 3.20 or later**

MPI
~~~

- OpenMPI
- MPICH
- Intel MPI

OpenMP
~~~~~~

HDF5
~~~~

- **HDF5 with C and C++ interfaces is required**

Operating System
~~~~~~~~~~~~~~~~

- Linux (tested on HPC environments)


Build from Source
-----------------

First, clone the repository:

.. code-block:: bash

   git clone https://github.com/JinsuRhee/VELUGA-CTree.git
   cd VELUGA-CTree

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


Note
^^^^

By default, **CTree** is built with the following data types:

- **Catalog ID**: 32-bit integer (``int32``)
- **Snapshot number**: 32-bit integer (``int32``)
- **Particle ID**: 64-bit integer (``int64``)
- **Merit value**: 64-bit floating point (``double``)

These defaults are suitable for most use cases and require no additional
build options.

If different data types are required, they can be selected **at compile time**
via CMake build options. For example:

.. code-block:: bash

   cmake -S . -B build \
     -DCTREE_ID_TYPE=int64 \
     -DCTREE_SNAP_TYPE=int32 \
     -DCTREE_PARTID_TYPE=int64 \
     -DCTREE_MERIT_TYPE=float32

   cmake --build build

Available options are:

- ``CTREE_ID_TYPE``: ``int32`` | ``int64``
- ``CTREE_SNAP_TYPE``: ``int32`` | ``int64``
- ``CTREE_PARTID_TYPE``: ``int32`` | ``int64``
- ``CTREE_MERIT_TYPE``: ``float32`` | ``float64``

All type selections are resolved at compile time and applied consistently
throughout the code base. Changing these options requires a full rebuild
of the project.


(Optional) Install
------------------

If installation is supported, you can install the binaries with:

.. code-block:: bash

   cmake --install .

Build Options
-------------

You can pass additional options to CMake if needed. For example:

.. code-block:: bash

   cmake .. -DCMAKE_BUILD_TYPE=Release

Common build types are:

- ``Release``
- ``Debug``
- ``RelWithDebInfo``

Troubleshooting
---------------

- If CMake cannot find a dependency, make sure it is installed and visible in your environment.
- Try removing the build directory and configuring again if you encounter unexpected errors.


