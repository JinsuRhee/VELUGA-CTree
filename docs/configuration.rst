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

   * - iotype
     - ``string``
     - *(required)*
     - Catalog type.

       Currently (as of Feb 2026), the following options are available.

       - ``VR`` (VELOCIraptor in ascii format)
       - ``VELUGA`` (VELOCIraptor post-processed data using https://github.com/JinsuRhee/VELUGA)
       - ``HM`` (HaloMaker catalogs)

   * - 
     -
     -
     -

   * - vr_dir_catalog
     - ``string``
     - *(required)*
     - Path specifier when ``iotype = VR``

       CTree searches raw VELOCIraptor outputs (*.catalog.dat0) 

       in directories divided by snapshot numbers

       (e.g.,) ``vr_dir_catalog``/``vr_dir_catalog_prefix`` ``SSSS`` ``vr_dir_catalog_suffix``/

   * - vr_dir_catalog_prefix
     - ``string``
     - optional
     - 

   * - vr_dir_catalog_suffix
     - ``string``
     - optional
     -

   * - vr_dir_catalog_snapdigit
     - ``integer``
     - *(required)*
     - The digit number for snapshot number in the directory name

   * - 
     -
     -
     -

   * - veluga_dir_catalog
     - ``string``
     - *(required)*
     - path specifier when ``iotype = VELUGA``

       directory path which includes snap_* directories

       :ref:`horg` option should be given

   * - .. _horg:
     -
     -
     -
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

