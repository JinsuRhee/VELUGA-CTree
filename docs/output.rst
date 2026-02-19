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

.. _data_format:
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

Branch File Format
------------------

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

       **Note:** If ``N_point=0``, the following arrays are not writtien and a new branch data starts

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
     - **Note:** If ``N_merge=0``, the following arrays are not writtien and a new branch data starts

       The list of IDs that merged into this branch

   * - ``Type_Snap [N_merge]``
     - The list of Snapshots that merged into this branch

   * - ``Type_Merit [N_merge]``
     - The list of merit scores when merged

   * - ``Type_BID [N_merge]``
     - Branch indices that merged into this branch


.. _reading_example:

Reading example
---------------


Python

``example/rdtree.py`` generates a pickle file containing the tree and key arrays.

The following script shows how to load the tree arrays and extract a branch of an object.
~~~~~~

.. code-block:: python

   import pikcle as pickle
   import numpy as np

   # "ctree.pkl" is what rdtree.py generates

   with open("ctree.pkl", 'rb') as f:
       data = pickle.load(f)

   key = data['key']
   tree = data['tree']

   # Extract a branch of galaxy with (ID=1 & Snap=100)

   id0  = 1
   snap0 = 100

   keyv   = key[snap0 + key[0]*id0]

   if(keyv < 0 ):
       print('no corresponding branch')
   else:
       branch = tree[keyv]

       print(branch['br_len']) # Length of the branch
       print(branch['id']) # List of IDs
       print(branch['snap']) # List of snapshots
       print(branch['merit']) # Merit score of the connections

       print(branch['father_ID']) # If this branch is merged, the father branch index.
                                  # The corresponding branch is tree[branch['father_ID']]
       
       print(branch['n_mergerbr']) # The number of branches that merged into this one
       print(branch['m_id']) # ID list of the merged branches
       print(branch['m_snap']) # Snapshot list of the merged branches
       print(branch['m_merit']) # Merit of mergers
       print(branch['m_bid']) # Indices of the merged branches
                              # , corresponding tree[ branch['m_bid'][:] ]


IDL
~~~

.. code-block:: idl

   ; "ctree.sav" is what rdtree.pro generates
   RESTORE, 'ctree.sav'

   ; Extract a branch of galaxy with (ID=1 & Snap=100)

   id0  = 1L
   snap0 = 100L

   keyv   = tree_key[snap0 + tree_key[0]*id0]

   IF keyv LT 0L THEN BEGIN
      PRINT, 'no corresponding branch'
   ENDIF ELSE BEGIN
      tree0   = tree_data[keyv]   ;; tree_data is a pointer array
      tree0   = *tree0

      PRINT, tree0.br_len     ; Length of the branch
      PRINT, tree0.id         ; List of IDs
      PRINT, tree0.snap       ; List of snapshots
      PRINT, tree0.merit      ; Merit score of the connections
      PRINT, tree0.father_ID  ; If this branch is merged, the father branch index.
                              ; The corresponding branch is tree[branch['father_ID']]
       
      PRINT, tree0.n_mergerbr ; The number of branches that merged into this one
      PRINT, tree0.m_id       ; ID list of the merged branches
      PRINT, tree0.m_snap     ; Snapshot list of the merged branches
      PRINT, tree0.m_merit    ; Merit of mergers
      PRINT, tree0.m_bid      ; Indices of the merged branches
                              ; , corresponding tree[ branch['m_bid'][:] ]
   ENDELSE