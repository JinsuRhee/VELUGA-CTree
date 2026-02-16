Output data format
==================

**CTree** writes its main outputs in **binary format** using C++
file streams (``std::ofstream`` with ``std::ios::binary``).  
The binary format is designed for efficient I/O and compact storage of
large tree datasets.

This page documents the exact binary layout of the output files and
provides example routines for reading them.

Overview
--------

Two binary files are written:

- ``ctree_key.dat``: Mapping between ``(Snapshot, Halo ID)`` and branch indices
- ``ctree_tree.dat``: Full tree and branch data

All integer and floating-point values are written in **native binary
representation** without text encoding.

Type Specifiers
---------------

Each binary file begins with one or more **type specifiers**, which encode
the data type used at compile time.

The type specifier is stored as a 32-bit integer with the following meaning:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Value
     - Data type
   * - ``1``
     - 32-bit integer (``int32``)
   * - ``2``
     - 64-bit integer (``int64``)
   * - ``3``
     - 32-bit floating point (``float``)
   * - ``4``
     - 64-bit floating point (``double``)

These specifiers allow the reader to determine the exact binary layout
without prior knowledge of the build configuration.

File: ctree_key.dat
-------------------

Purpose
~~~~~~~

``ctree_key.dat`` provides a mapping between catalog entries and internal
branch indices used in the tree structure.

Binary layout
~~~~~~~~~~~~~

The file is written in the following order:

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Index
     - Type
     - Description
   * - 1
     - ``int32``
     - Type specifier for branch ID (``Tree_BID``)
   * - 2
     - ``int64``
     - Number of key entries ``N``
   * - 3
     - ``Tree_BID``
     - First key value
   * - 4
     - ``Tree_BID[N]``
     - Remaining key values stored as a contiguous array

File: ctree_tree.dat
--------------------

Purpose
~~~~~~~

``ctree_tree.dat`` stores the complete tree structure, including
main branches, progenitors, and merged branches.

Header layout
~~~~~~~~~~~~~

The file header consists of the following elements:

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Index
     - Type
     - Description
   * - 1
     - ``int32``
     - Type specifier for halo ID (GID)
   * - 2
     - ``int32``
     - Type specifier for snapshot number
   * - 3
     - ``int32``
     - Type specifier for branch ID (BID)
   * - 4
     - ``int32``
     - Type specifier for merit value
   * - 5
     - ``int64``
     - Total number of trees ``Ntree``
   * - 6
     - ``int64``
     - Maximum index (``lind``) of the last tree

Per-tree layout
~~~~~~~~~~~~~~~

For each tree, the following data are written sequentially:

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Field
     - Type
     - Description
   * - ``n_branch``
     - ``int32``
     - Number of nodes in the main branch
   * - ``n_numprg``
     - ``int32``
     - Number of merged progenitor branches
   * - ``father_bid``
     - ``int32``
     - Branch ID of the parent tree
   * - ``stat``
     - ``int32``
     - Tree status flag

Main branch arrays
~~~~~~~~~~~~~~~~~

If ``n_branch > 0``, the following arrays of length ``n_branch`` are written:

- ``id``: Halo IDs
- ``snap``: Snapshot numbers
- ``p_id``: Progenitor halo IDs
- ``p_snap``: Progenitor snapshot numbers
- ``p_merit``: Progenitor merit values

Merged progenitor arrays
~~~~~~~~~~~~~~~~~~~~~~~~

If ``n_numprg > 0``, the following arrays of length ``n_numprg`` are written:

- ``m_id``: Merged branch halo IDs
- ``m_snap``: Merged branch snapshot numbers
- ``m_merit``: Merged branch merit values
- ``m_bid``: Merged branch IDs

Reading example
---------------

Python
~~~~~~

.. code-block:: python

   import struct
   import numpy as np

   with open("ctree_tree.dat", "rb") as f:
       tags = struct.unpack("4i", f.read(16))
       n_tree = struct.unpack("q", f.read(8))[0]
       lind   = struct.unpack("q", f.read(8))[0]

       for _ in range(n_tree):
           n_branch, n_numprg, father_bid, stat = struct.unpack("4i", f.read(16))

           if n_branch <= 0:
               continue

           ids   = np.fromfile(f, dtype=np.int64, count=n_branch)
           snaps = np.fromfile(f, dtype=np.int64, count=n_branch)
           p_ids = np.fromfile(f, dtype=np.int64, count=n_branch)
           p_snaps = np.fromfile(f, dtype=np.int64, count=n_branch)
           p_merit = np.fromfile(f, dtype=np.float64, count=n_branch)

           if n_numprg > 0:
               m_ids   = np.fromfile(f, dtype=np.int64, count=n_numprg)
               m_snaps = np.fromfile(f, dtype=np.int64, count=n_numprg)
               m_merit = np.fromfile(f, dtype=np.float64, count=n_numprg)
               m_bid   = np.fromfile(f, dtype=np.int64, count=n_numprg)

IDL
~~~

.. code-block:: IDL
   print, 'hello world'