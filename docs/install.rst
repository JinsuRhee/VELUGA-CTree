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


