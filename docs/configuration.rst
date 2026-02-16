Configuration Options
=====================

This page documents the configuration options used by **CTree**.

General notes:
    - Lines starting with ``#`` are treated as comments.
    - Whitespace around ``=`` is ignored.

I/O related
-----------

.. list-table::
   :header-rows: 0
   :widths: 30 20 20 60

   * - **Option**
     - **Format**  
     - **Default**
     - **Description**

   * - out_dir
     - ``string``
     - *(required)*
     - Path specifier for the output data.

       The output data are "ctree_tree.dat" and "ctree_key.dat"

       - "ctree_tree.dat" has the Lagrangian branch data for the entire branches.
       - "ctree_key.dat" gives the mapping for (#Snap, #ID)->(branch index).
       
       A detailed data format and example routines reading the output are given :output:



Examples
--------

Minimal example:

.. code-block:: ini

   input_file   = catalog.h5
   snapshot_end = 99

Example using a snapshot list:

.. code-block:: ini

   input_file     = catalog.h5
   snapshot_list  = snapshots.txt
   output_dir     = ./ctree_out
   nthreads       = 8
   log_level      = debug

