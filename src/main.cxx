// src/main.cpp
#ifdef CTREE_USE_MPI
  #include <mpi.h>
#endif

#ifdef CTREE_USE_OMP
  #include <omp.h>
#endif

// Header for global
#include "global/allvar.h"

// Header for Makebr
#include "makebr/makebr.h"

// Header for Ctree
#include "complete_tree/ctree.h"
// Header for utilities
#include "utilities/utilities.h"
//#include <iostream>
//#include "allvar.h"
//#include "logger.h"
//#include "utilities.h"
//#include "config.h"
//#include "autovec.h"
//#include "types.h"


//----- Main code
int main(int argc, char** argv) {
#ifdef CTREE_USE_MPI
  MPI_Init(&argc, &argv);
#endif
  int myrank  = mpi_rank();

#ifdef CTREE_USE_OMP
  int nthreads  = omp_get_max_threads();
  omp_set_num_threads(nthreads);
#endif

  //-----
  // Initial Part
  //-----

  // Logo print
  if(myrank == 0) u_printlogo();

  // Initial Runcheck
  if(myrank == 0){
    if( !u_initialcheck(argc) ){
    	u_stop();
    } else{
    	LOG()  << "Program starts";

    	// TODO 123123
    	// Basic spec info will be printed here
    }
  }

  // Load Configuration
  vctree_set::Settings vh;
  if( !g_load_config(argv[1], vh) ){
  	u_stop();
  } else{
  	if(myrank == 0) LOG() <<"  Configuration Loaded";
  }

  //-----
  // Make Branch Part
  //-----
  Tree::TreeArray tree;
  Tree::TreeKeyArray key;
  if(vh.iotype == "VR"){
    if(myrank == 0) LOG() <<"  Entering to MakeBr";
    //Makebr::mainloop(vh, tree, key);
  }else{
    //TODO 123123
    //For Other OUTPUT
  }




#ifdef CTREE_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  //-----
  // Load Existing tree
  //-----
  if(myrank==0){
    LOG() <<"  Reading to Exsiting Tree";
  }
  //Tree::TreeArray tree2;
  //Tree::TreeKeyArray key2;
  //loadtree(vh, tree, key);

#ifdef CTREE_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  //-----
  // Complete Tree Part
  //-----
  if(myrank==0) LOG() <<"  Entering to Complete Tree";
  Ctree::main(vh, tree, key);

  

LOG()<<" merit comparison check";
LOG()<<" without makebr? save in process, load it. case for single tree";
LOG()<<" commerit2 mpi";

  //-----
  // Save Tree
  //-----
#ifdef CTREE_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if(myrank == 0){
    LOG() <<"  Save the Tree & Key";
    if(myrank == 0) savetree(vh, tree, key);
  }

#ifdef CTREE_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  //-----
  // Program End
  //-----
  if(myrank == 0){
    LOG() <<"Program end";
  }
  MPI_Finalize();
  return 0;
}