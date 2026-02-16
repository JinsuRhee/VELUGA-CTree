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

.. note :: Type Specifiers

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

Tree Key File Format
--------------------

``ctree_key.dat`` provides a mapping between catalog entries and internal
branch indices used in the tree structure.

The key array gives a mapping between (ID, Snapshot number) pairs and their indices in the merger tree branches:

.. code-block:: text

   branch_index = key[snapnumber + id * keyvalue]


In **CTree** the default key value is defined as ``the maximum snapshot number + 1``

The file is written in the following order:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Type and Size
     - Description
   * - ``int32 [1]``
     - Type specifier for branch index (``Tree_BID``)
   * - ``int64 [1]``
     - Size of the key array (``N``)
   * - ``Tree_BID [1]``
     - Key value
   * - ``Tree_BID [N]``
     - N elements in the key array

File: ctree_tree.dat
--------------------

``ctree_tree.dat`` stores the whole branch information.
~~~~~~~~~~~~~

The file header consists of the following 6 elements:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Type and Size
     - Description
   * - ``int32 [1]``
     - Type specifier for object ID (``Type_ID``)
   * - ``int32 [1]``
     - Type specifier for snapshot number (``Type_Snap``)
   * - ``int32 [1]``
     - Type specifier for branch index (``Type_BID``)
   * - ``int32 [1]``
     - Type specifier for merit value (``Type_Merit``)
   * - ``int64 [1]``
     - Total number of the branch array (``Nbranch``)

After the header, the file stores **one record per tree**, written
sequentially in a loop over all trees.


.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Type & Size
     - Description
   * - ``int32 [1]``
     - Number of points in the main branch (``N_point``)

       ** Note ** If ``N_point = 0``, the following arrays are not writtien and a new branch data starts

   * - ``int32 [1]``
     - Number of merged progenitor branches (``N_merge``)

   * - ``BID [1]``
     - If this branch is merged into another, the branch index that this branch is merged into

   * - ``int32 [1]``
     - Tree status flag (now it has no information)

   * - ``Type_ID [N_point]``
     - The lisf of IDs of the main branch

   * - ``Type_Snap [N_point]``
     - The list of snapshot numbers of the main branch

   * - ``Type_ID [N_merge]``
     - ** Note ** If ``N_merge = 0``, the following arrays are not writtien and a new branch data starts

       The list of IDs that merged into this branch

   * - ``Type_Snap [N_merge]``
     - The list of Snapshots that merged into this branch

   * - ``Type_Merit [N_merge]``
     - The list of merit scores when merged

   * - ``Type_BID [N_merge]``
     - Branch indices that merged into this branch


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