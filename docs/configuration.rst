Configuration Options
=====================

This page documents the configuration options used by **CTree**.

General notes:
    - Lines starting with ``#`` are treated as comments.
    - Whitespace around ``=`` is ignored.

Option reference
----------------

.. list-table::
   :header-rows: 0
   :widths: 30 20 20 60

   * - **Option**
     - **Format**  
     - **Default**
     - **Description**

   * - **I/O related**
     - 
     - 
     -

   * - out_dir
     - ``string (path)``
     - *(required)*
     - Path specifier for the output data. The output data are "ctree_tree.dat" and "ctree_key.dat"
       "ctree_tree.dat" has the Lagrangian branch 
       data for the entire branches.
       "ctree_key.dat" gives the mapping for (#Snap, #ID)->(branch index).
       A detailed data format and example routines reading the output are given :output:

I/O
~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 60

   * - out_dir
     - ``string (path)``
     - *(required)*
     - Path specifier for the output data. The output data are "ctree_tree.dat" and "ctree_key.dat"
       "ctree_tree.dat" has the Lagrangian branch data for the entire branches.
       "ctree_key.dat" gives the mapping for (#Snap, #ID)->(branch index).
       A detailed data format and example routines reading the output are given :output:

   * - iotype
     - ``string``
     - *(required)*
     - Catalog type. Currently (as of Feb 2026), VELOCIraptor (ascii format), VELUGA (VELOCIraptor post-processed data using https://github.com/JinsuRhee/VELUGA), and HaloMaker catalogs are available.


I/O (VELOCIraptor)
^^^^^^^^^^^^^^^^^^

CTree assumes that a raw VELOCIraptor output (e.g., *.catalog.dat0) is saved in (vr_dir_catalog)/(vr_dir_catalog_prefix)(SSSS)(vr_dir_catalog_suffix).
(SSSS) is the snapshot number with a digit of vr_dir_catalog_snapdigit

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 60
   * - vr_dir_catalog
     - ``string (path)``
     - *(required)*
     - Directory where raw VELOCIraptor outputs are saved. CTree assumes that a raw output (e.g., *.catalog.dat0) is saved in a subdirectory specifying snapshot numbers

  * - vr_dir_catalog_prefix
     - ``string``
     - *(optional)*
     - Prefix for subdirectories

  * - vr_dir_catalog_suffix
     - ``string``
     - *(optional)*
     - Suffix for subdirectories

  * - vr_dir_catalog_snapdigit
     - ``integer``
     - *(optional)*
     - Snapshot number digit


I/O (VELUGA)
^^^^^^^^^^^^

I/O (HaloMaker)
^^^^^^^^^^^^


Snapshots / Time selection
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 60

   * - Option
     - Format
     - Default
     - Description

   * - ``snapshot_start``
     - ``int``
     - ``0``
     - First snapshot number to process (inclusive).

   * - ``snapshot_end``
     - ``int``
     - *(required)*
     - Last snapshot number to process (inclusive).

   * - ``snapshot_step``
     - ``int``
     - ``1``
     - Process every N-th snapshot.

   * - ``snapshot_list``
     - ``string (path)``
     - *(optional)*
     - Alternatively, provide a file containing a list of snapshot numbers (one per line).
       If set, this option overrides ``snapshot_start/end/step``.

Merit / Tree building
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 60

   * - Option
     - Format
     - Default
     - Description

   * - ``core_npart``
     - ``int``
     - ``50``
     - Number of core particles used when computing the merit function.

   * - ``merit_mode``
     - ``string``
     - ``core_count``
     - Merit function mode. For example, ``core_count`` selects the halo containing the
       largest number of core particles.

   * - ``first_progenitor_rule``
     - ``string``
     - ``max_merit``
     - Rule to select the first progenitor (e.g. highest merit score).

Output controls
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 60

   * - Option
     - Format
     - Default
     - Description

   * - ``write_branches``
     - ``bool``
     - ``true``
     - Write intermediate branch products.

   * - ``write_complete_tree``
     - ``bool``
     - ``true``
     - Write the final complete tree output.

   * - ``log_level``
     - ``string``
     - ``info``
     - Logging verbosity. Typical values: ``debug``, ``info``, ``warning``, ``error``.

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

