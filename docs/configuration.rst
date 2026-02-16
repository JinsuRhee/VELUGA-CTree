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
     - Path specifer for the catalog data ``iotype = VR``

       **CTree** searches raw VELOCIraptor outputs (*.catalog.dat0) in directories divided by snapshot numbers

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
     - Path specifer for the catalog data when ``iotype = VELUGA``

       Path at which ``snap_*`` directories are located

       :ref:`horg <horg-option>` option should be given

   * - 
     -
     -
     -

   * - hm_dir_catalog
     - ``string``
     - *(required)*
     - Path specifer for the catalog data when ``iotype = HM``

       Path at which ``*.bricks`` are located

   * - 
     -
     -
     -

   * - loadtree
     - ``Y or N``
     - *(required)*
     - Flag whether a pre-existing merger tree is used

       If set, **CTree** begins with the given merger tree data and reparis the branches

       If set, ``loadtree_ftree`` and ``loadtree_fkey`` should be given

       The tree data format is described in :doc:`Output data format <output>`

   * - loadtree_ftree
     - ``string``
     - *(required)*
     - Path specifier for the merger tree branch data when ``loadtree=Y``

   * - loadtree_fkey
     - ``string``
     - *(required)*
     - Path specifier for the merger tree key data when ``loadtree=Y``

   * - snapi, snapf
     - ``integer``
     - *(required)*
     - Snapshot range specifiers. **CTree** loops from ``snapf`` to ``snapi``

       A user can also give the snapshot list by a Ascii file with :ref:`snaplist <snaplist-option>`

   * - .. _snaplist-option:
       snaplist
     - ``string``
     - *optional*
     - File specifier including the list of snapshots (line-by-line).

       The list is automatically sorted
   
   * - .. _horg-option:
       horg
     - ``h or g``
     - *optional*
     - A flag for the catalog type (Halo or Galaxy).

       It is only used when reading the catalog data (``iotype==VELUGA``).


Branch related parameters
-------------------------

.. list-table::
   :header-rows: 0
   :widths: 30 20 20 60

   * - **Option**
     - **Format**  
     - **Default**
     - **Description**

   * - ctree_n_search
     - ``integer``
     - 10
     - The number of snapshots to be used in progenitor searches

       It determines the size of list tables of branches where candidate progenitor information is written

   * - ctree_meritfrac
     - ``float``
     - 0.5
     - Fraction to allow a suboptimal later progenitor to be a branch point

       **CTree** first defines the progenitor with the highest merit in the list talbe as the next branch point

       If a later progenitor with a merit > ``ctree_merit_frac`` X (Maximum merit), **CTree** allows this progenitor to be the next one, mainly to avoid frequent and big jumps among the branch points.

   * -
     -
     -
     -

   * - ctree_core_n
     - ``integer``
     - 10
     - The number of snapshots used when extracting core particles.

       ``ctree_core_n`` snapshots with the ``ctree_core_dn`` interval are used

   * - ctree_core_dn
     - ``integer``
     - 5
     - Snapshot interval for the core particle extraction

   * - ctree_minfrac
     - ``float``
     - 0.25
     - Fraction to define core particles.

       A subset of particles with N occurrances are defined as core particles if their fraction is higher than ``ctree_meritfrac``

       Otherwise, **CTree** repeats the core particle searches with a lower occurrance number

   * -
     -
     -
     -

   * - ctree_makecheck
     - ``integer``
     - 10
     - A modulo value for saving checkpoint snapshots.

       If a snapshot number is a modulo of ``ctree_makecheck``, then the data at the corresponding snapshot are saved.

       **CTree** can restart the execution from checkpoint snapshots

       If a negative value is given, then **CTree* does not generate checkpoint snapshots

   * - ctree_loadcheck
     - ``integer``
     - -1
     - A snapshot number for the restart, if exists.

   * -
     -
     -
     -

   * - meritlimit
     - ``float``
     - 0.001
     - The lower limit of the merit for progenitors.

       If the best progenitor has the lower merit than ``meritlimit``, the connection is not allowed

   * - minbranchlength
     - ``integer``
     - 5
     - After all branches are constructed, branches whose length is lower than ``minbranchlength`` are removed

Branch Maker
------------

These options are used when building initial branches based on the connectivity information of objects between snapshots.

**TreeFrog** output is only available at the moment

.. list-table::
   :header-rows: 0
   :widths: 30 20 20 60

   * - **Option**
     - **Format**  
     - **Default**
     - **Description**
     
   * - branchmaker
     - ``Y or N``
     - *(required)*
     - Flag to turn on/off the branchmaker

   * - brtype
     - ``string``
     - *(required)*
     - Tree data catalog type

       - ``TF`` (TreeFrog)

   * - tfdir
     - ``string``
     - *(required)*
     - Path specifier where TreeFrog raw outputs (*.tree) are located

   * - treedir
     - ``des or prg``
     - des
     - Branch making direction

       - ``des`` : descendants searches
       - ``prg`` : progenitors searches 
