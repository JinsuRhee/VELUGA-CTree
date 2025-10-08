#pragma once
#include <algorithm>
#include <cstddef>
#include <chrono>
#include <ctime>
#include <functional>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>


#ifdef CTREE_USE_MPI
	#include <mpi.h>
#endif



// Control variables
namespace vctree_parameters{
	// array size for the elements of tree arrays
	//inline int32_t sarr_size = 1000;
  	//inline int32_t barr_size = 100;
  	//inline int32_t tarr_size = 10000;
  	inline int32_t ctree_nstep 	= 10000;
  	inline int32_t ctree_npid 	= 100000;
  	inline double ctree_minfrac = 0.25;
}


// For Definitions
// [YY:MM:DD:SS] [main.cpp] [123] ' message'
#define LOG()   vctree_log::Line(vctree_log::basename_cstr(__FILE__), __LINE__)
// [YY:MM:DD:SS] [MYCODE] [123] ' message'
#define LOGC(x) vctree_log::Line((x), __LINE__)

//-----
// Name space For Settings definition
//-----
namespace vctree_set{

	// Create Settings structure
	struct Settings {
		std::string iotype 		= "VR";			  // IO type
		std::string	simtype		= "Ramses";		  // Simulation type

		// For VR IO
		std::string vr_dir_catalog 	= "./catalog/";   // VR output directory
		std::string vr_dir_tree    = "";              // TreeFrog data is stored
  		char        horg        = 'g';            	  // 'g' or 'h'
  		
  		// RAMSES related
  		std::string ramses_dir = "";
  		

  		// General settings
  		int32_t snapi   = 70;              // snapshot range
  		int32_t snapf   = 200;
  		double  meritlimit = 0.001;        // melit lower limit to close the branch
  		std::string treedir     = "des";   // Tree direction
  		std::string out_dir     = "";   // Tree direction
  		int32_t minbranchlength = 10;		// Minimum length of brnach

  		// makebr related
  		int64_t treekey     = 1000;  // 10^n; should be larger than the last snapshot number

  		// complete_tree related
  		int32_t ctree_n_search   = 10;     // number of snapshots to be searched simultaneously
  		int32_t ctree_n_step_n   = 10;     // number of branch points (>10) when collecting particles on a existing branch
  		int32_t ctree_n_step_dn  = 5;      // dN between the points (corresponding to 200 MYr seems good)
  		double  ctree_rfact      = 10.;    // Find Candidate galaxies within a sphere of a radius of speed * rfact * dT (10 seems good)
  		int32_t ctree_weighttype = 1; 	   // Weighting type (1 :TF )
  		// TreeFrog I/O related
  		std::string tag_num     = "";
  		std::string tag_off     = "";
  		std::string tag_result  = "";
  		std::string tag_npart   = "";
  		std::string tag_merit   = "";
  		std::string tag_nlink   = "";

  		// Some utils
  		void finalize_paths() {
  			if (iotype == "VR"){
  				vr_dir_catalog = vr_dir_catalog + (horg=='g' ? "Galaxy/VR_Galaxy" : "Halo/VR_Halo");
  				if (vr_dir_tree.empty()) {
      				vr_dir_tree = vr_dir_catalog + "../tree/tfout";
    			}
  			}
  		}

  		// Bring control varibles
  		//int32_t sarr_size = vctree_parameters::sarr_size;
  		//int32_t barr_size = vctree_parameters::barr_size;
  		//int32_t tarr_size = vctree_parameters::tarr_size;
  		int32_t ctree_nstep 	= vctree_parameters::ctree_nstep;
  		int32_t ctree_npid 		= vctree_parameters::ctree_npid;
  		double ctree_minfrac 	= vctree_parameters::ctree_minfrac;

	};
}


//-----
// Name space For Log printing
//-----
namespace vctree_log{
	// Extract Filename
	inline const char* basename_cstr(const char* path) {
  		const char* base = path;
  		for (const char* p = path; *p; ++p) {
    		if (*p == '/' || *p == '\\') base = p + 1;
  		}
  		return base;
	}

	// Time Stamp
	inline std::string timestamp_yy_mm_dd_ss() {
  		using namespace std::chrono;
  		const auto now = system_clock::now();
  		const auto tt  = system_clock::to_time_t(now);
  		std::tm tm{};
#if defined(_WIN32)
  		localtime_s(&tm, &tt);
#else
  		localtime_r(&tt, &tm);
#endif
  		std::ostringstream oss;
  		oss << '[' << std::put_time(&tm, "%y:%m:%d:%S") << ']'; // YY:MM:DD:SS
  		return oss.str();
	}

	// Config
	struct Config {
#ifdef LOGGER_THREAD_SAFE
  		std::mutex mu;
#endif
	};

	inline Config& cfg() { static Config c; return c; }

	// Log by line
	class Line {
		public:
  			Line(const char* code, int line) : code_(code), line_(line) {}
  			~Line() {
#ifdef LOGGER_THREAD_SAFE
    		std::lock_guard<std::mutex> lock(cfg().mu);
#endif
    		std::ostream& out = std::cout;
    		out << timestamp_yy_mm_dd_ss()
        		<< " [" << std::setw(16) << code_ << "] "
        		<< "[" << std::setw(4) << line_ << "] "
        		<< " " << ss_.str() << "\n";
  			}
  			template<typename T>
  			Line& operator<<(const T& v) { ss_ << v; return *this; }

		private:
  			const char* code_;
  			int line_;
  			std::ostringstream ss_;
	};
}



//-----
// Data Type For Tree Array
//-----
namespace Tree{

	using Tree_GID 		= std::int32_t;		// for galaxy ID
	using Tree_Snap 	= std::int32_t;
	using Tree_PID 	 	= std::int64_t;		// for particle ID
	using Tree_merit	= double;
	using Tree_BID 		= std::int64_t;		// Branch Index
	using Tree_I32 		= std::int32_t;
	using Tree_I64 		= std::int64_t;

	inline Tree_I32 tree_stepsize = 100;

	struct TreeSt{
		std::vector<Tree_GID> id = std::vector<Tree_GID>(tree_stepsize);
		std::vector<Tree_Snap> snap = std::vector<Tree_Snap>(tree_stepsize);
		std::vector<Tree_GID> p_id = std::vector<Tree_GID>(tree_stepsize);
		std::vector<Tree_Snap> p_snap = std::vector<Tree_Snap>(tree_stepsize);
		std::vector<Tree_merit> p_merit = std::vector<Tree_merit>(tree_stepsize);
		//std::vector<Tree_GID> d_id = std::vector<Tree_GID>(tree_stepsize);
		//std::vector<Tree_Snap> d_snap = std::vector<Tree_Snap>(tree_stepsize);

		
		std::vector<Tree_GID> m_id; 			// id list merged into this branch
		std::vector<Tree_Snap> m_snap;			// snap list merged into this branch
		std::vector<Tree_merit> m_merit;		// merit when merged
		std::vector<Tree_BID> m_bid;			// tree index for the mergered branch

		Tree_I32 father_bid = -1; 				// if this branch is merged into another, father tree index
		Tree_I32 frag_bid = -1; 				// if this branch is fragmented, originated branch index;

		Tree_I32 numprog = 0;
		Tree_I32 endind = -1; // last index of TreeSt

		Tree_I32 stat = 0; 	// 0 (normal) (-1) natrual end (-2) fragmented

		Tree_BID lind = 0; // last index of vec<TreeST> only stored in the first tree

		Tree_I32 isfree = -1; // only used when MPI used in Ctree
	};

	struct TreeKey{
		Tree_BID ind = -1;
		Tree_I64 key = 0; // only stored in the first tree
	};

	using TreeArray = std::vector<TreeSt>;
	using TreeKeyArray = std::vector<TreeKey>;

	// Vector Input
	template <typename T>
	void treevecinput(std::vector<T>& v, Tree_I32 ind, T val){
		
		if( (double) ind >= (double) v.size()*0.8){
			v.resize(v.size() + tree_stepsize, T{});
		}
		v[ind] = val;

		//LOG()<<ind<<" / "<<val<<" / "<<v[ind];
	}

	// Vector Input
	template <typename T>
	void treevecinput_tofirst(std::vector<T>& v, Tree_I32 ind, T val){
		
		if( (double) ind >= (double) v.size()*0.8){
			v.resize(v.size() + tree_stepsize, T{});
		}
		v.insert(v.begin(), val);
	}

	inline Tree_BID getkey(TreeKeyArray& key, Tree_Snap snap, Tree_GID id){
		return key[ snap + key[0].key * id ].ind;
	}


	inline void treeinit(TreeArray& tree, TreeKeyArray& key, 
		Tree_Snap snap, Tree_GID id){

		tree[0].lind ++;
		Tree_I64 keyval;
		keyval 	= snap + key[0].key * id;

		key[keyval].ind = tree[0].lind;

		TreeSt& treedum 	= tree[ tree[0].lind ];

		treedum.endind ++;
		treevecinput<Tree_GID>(treedum.id, treedum.endind, id); 
		treevecinput<Tree_Snap>(treedum.snap, treedum.endind, snap); 
		treevecinput<Tree_merit>(treedum.p_merit, treedum.endind, 0.);

	}

	inline bool istree(TreeKeyArray& key, Tree_Snap snap, Tree_GID id){
		Tree_I64 keyval = getkey(key, snap, id);
		//keyval 	= snap + key[0].key * id;
		if(keyval > 0){
			return true;
		}
		else{
			return false;	
		} 
	}

	inline TreeSt gettree(TreeArray& tree, TreeKeyArray& key, Tree_Snap snap, Tree_GID id){
		Tree_I64 keyval = getkey(key, snap, id);
		//keyval 	= snap + key[0].key * id;
		TreeSt tree0;

		if(!istree(key, snap, id)){
			LOG()<<"no tree matched";
			int errcode = 1;
    		MPI_Abort(MPI_COMM_WORLD, errcode);
    		std::exit(errcode);
		}else{
			tree0 = tree[keyval];
		}

		return tree0;
	}

	

	inline void treeinput(TreeArray& tree, TreeKeyArray& key, 
		Tree_Snap snap, Tree_GID id, Tree_Snap to_snap, Tree_GID to_id, Tree_merit merit){

		if(!istree(key, snap, id)){
			treeinit(tree, key, snap, id);
		}


		//Tree_Snap p_snap, Tree_GID p_id, Tree_merit p_merit, 
		//Tree_Snap d_snap, Tree_GID d_id){

		Tree_BID keyval_org, keyval_new;
		keyval_org 	= snap + key[0].key * id;
		keyval_new 	= to_snap + key[0].key * to_id;


		key[keyval_new].ind = key[keyval_org].ind;
		TreeSt& treedum 	= tree[ key[keyval_org].ind ];

		
//int myrank  = mpi_rank();
//if(myrank ==0 && snap == 49 && id == 1 && )
//{
//	
//}	
		if(to_snap > treedum.snap[treedum.endind]){
			treedum.endind ++;
			treevecinput<Tree_GID>(treedum.id, treedum.endind, to_id); 
			treevecinput<Tree_Snap>(treedum.snap, treedum.endind, to_snap);
			//<Tree_GID>(treedum.p_id, treedum.endind, p_id); 
			//treevecinput<Tree_Snap>(treedum.p_snap, treedum.endind, p_snap); 
			treevecinput<Tree_merit>(treedum.p_merit, treedum.endind, merit);
			//treevecinput<Tree_GID>(treedum.d_id, treedum.endind, d_id);
			//treevecinput<Tree_Snap>(treedum.d_snap, treedum.endind, d_snap);
	
			//tree[ tree[0].lind ] = treedum;
		}else if(to_snap < treedum.snap[0]){
			treedum.endind ++;
			treevecinput_tofirst<Tree_GID>(treedum.id, treedum.endind, to_id); 
			treevecinput_tofirst<Tree_Snap>(treedum.snap, treedum.endind, to_snap);
			treevecinput_tofirst<Tree_merit>(treedum.p_merit, treedum.endind, merit); 
		}else{
			LOG()<<"Why this happens?";
			int errcode = 1;
    		MPI_Abort(MPI_COMM_WORLD, errcode);
    		std::exit(errcode);
		}

	}

	
	inline void treefree(TreeArray& tree, TreeKeyArray& key, Tree_Snap snap, Tree_GID id){

		
		if(!istree(key, snap, id)){
			LOG()<<"Why?";
			int errcode = 1;
    		MPI_Abort(MPI_COMM_WORLD, errcode);
    		std::exit(errcode);
			return;
		}

		Tree_I64 keyval = getkey(key, snap, id);
		//keyval 	= snap + key[0].key * id;
		TreeSt& tree0 = tree[ keyval ];

		

		// free key
		
		for(Tree_I32 i=0; i<tree0.endind+1; i++){

			keyval 	= tree0.snap[i] + key[0].key * tree0.id[i];
			key[keyval].ind = -1;
		}

		// free this tree
		tree0.id.resize(0);
		tree0.snap.resize(0);
		tree0.p_merit.resize(0);

		tree0.m_id.resize(0);
		tree0.m_snap.resize(0);
		tree0.m_merit.resize(0);
		tree0.m_bid.resize(0);

		tree0.father_bid = -1;
		tree0.numprog = 0;
		tree0.endind = -1;
	}

	// remove connection point before cut_snap
	inline void modifytree(TreeArray& tree, TreeKeyArray& key, Tree_Snap snap, Tree_GID id, Tree_Snap cut_snap){
		Tree_I64 keyval = getkey(key, snap, id);

		TreeSt& tree0	= tree[keyval];

		// modify the tree
		Tree_I32 nnn = 0;
		std::vector<Tree_I32> mod_ind(tree0.endind+1);

		for(Tree_I32 i=0; i<tree0.endind+1; i++){
			if(tree0.snap[i] > cut_snap){
				mod_ind[nnn] = i;
				nnn ++;
			}
		}

		if(nnn == 0){
			if(tree0.endind == 0 && tree0.snap[0] == cut_snap && tree0.snap[tree0.endind] == cut_snap){ // single snapshot case
				treefree(tree, key, snap, id);
			}
			return;
		}

		mod_ind.resize(nnn);

		std::vector<Tree_GID> new_id(tree0.id.size());
		std::vector<Tree_Snap> new_snap(tree0.snap.size());
		std::vector<Tree_merit> new_merit(tree0.p_merit.size());

		for(Tree_I32 i=0; i<nnn; i++){
			new_id[i] 		= tree0.id[mod_ind[i]];
			new_snap[i] 	= tree0.snap[mod_ind[i]];
			new_merit[i] 	= tree0.p_merit[mod_ind[i]];
		}

		for(Tree_I32 i=0; i<tree0.endind+1; i++){
			key[ tree0.snap[i] + key[0].key*tree0.id[i] ].ind = -1;
		}


		tree0.id.clear();
		tree0.snap.clear();
		tree0.p_merit.clear();

		tree0.id 		= std::move(new_id);
		tree0.snap 		= std::move(new_snap);
		tree0.p_merit 	= std::move(new_merit);
		tree0.endind 	= nnn-1;

		for(Tree_I32 i=0; i<nnn; i++){
			key[ tree0.snap[i] + key[0].key*tree0.id[i] ].ind = keyval;
		}

		//

		if(tree0.numprog > 0){
			std::vector<Tree_GID> new_mid(tree0.numprog);
			std::vector<Tree_Snap> new_snap(tree0.numprog);
			std::vector<Tree_merit> new_merit(tree0.numprog);
			std::vector<Tree_BID> new_bid(tree0.numprog);

			Tree_I32 np = 0;

			for(Tree_I32 i=0; i<tree0.numprog; i++){
				if(tree0.m_snap[i] > cut_snap){
					//mod_ind[np] 	= i;
					new_mid[np] 	= tree0.m_id[i];
					new_snap[np] 	= tree0.m_snap[i];
					new_merit[np]	= tree0.m_merit[i];
					new_bid[np]		= tree0.m_bid[i];
					np ++;
				}
			}

			tree0.m_id.clear();
			tree0.m_snap.clear();
			tree0.m_merit.clear();
			tree0.m_bid.clear();
			tree0.numprog = 0;

			if(np > 0){
				//mod_ind.resize(np);

				tree0.m_id 		= std::move(new_mid);
				tree0.m_snap	= std::move(new_snap);
				tree0.m_merit	= std::move(new_merit);
				tree0.m_bid 	= std::move(new_bid);
				tree0.numprog 	= np;
			}
		}

		
	}

}

//-----
// MPI related
//-----
// mpi_types.h
#ifdef CTREE_USE_MPI
namespace mpi_type {
	
	template<typename T>
	MPI_Datatype type();

	// ---- 32-bit signed ----
	template<> inline MPI_Datatype type<std::int32_t>() {
	#ifdef MPI_INT32_T
    	return MPI_INT32_T;
	#else
    	if (sizeof(int)  == 4) return MPI_INT;
    	if (sizeof(long) == 4) return MPI_LONG;
    	if (sizeof(short)== 4) return MPI_SHORT;
    	throw std::logic_error("No 32-bit MPI integer type available");
	#endif
	}

	// ---- 32-bit unsigned ----
	template<> inline MPI_Datatype type<std::uint32_t>() {
	#ifdef MPI_UINT32_T
	    return MPI_UINT32_T;
	#else
	    if (sizeof(unsigned int)  == 4) return MPI_UNSIGNED;
	    if (sizeof(unsigned long) == 4) return MPI_UNSIGNED_LONG;
	    if (sizeof(unsigned short)== 4) return MPI_UNSIGNED_SHORT;
	    throw std::logic_error("No 32-bit MPI unsigned type available");
	#endif
	}

	// ---- 64-bit signed ----
	template<> inline MPI_Datatype type<std::int64_t>() {
	#ifdef MPI_INT64_T
	    return MPI_INT64_T;
	#else
	    if (sizeof(long long) == 8) return MPI_LONG_LONG;
	    if (sizeof(long)      == 8) return MPI_LONG;      // LP64
	    throw std::logic_error("No 64-bit MPI integer type available");
	#endif
	}

	// ---- 64-bit unsigned ----
	template<> inline MPI_Datatype type<std::uint64_t>() {
	#ifdef MPI_UINT64_T
	    return MPI_UINT64_T;
	#else
	    if (sizeof(unsigned long long) == 8) return MPI_UNSIGNED_LONG_LONG;
	    if (sizeof(unsigned long)      == 8) return MPI_UNSIGNED_LONG; // LP64
	    throw std::logic_error("No 64-bit MPI unsigned type available");
	#endif
	}

	// ---- float / double ----
	template<> inline MPI_Datatype type<float>()  { return MPI_FLOAT;  }
	template<> inline MPI_Datatype type<double>() { return MPI_DOUBLE; }



} // namespace mpi_type
#endif // CTREE_USE_MPI

//-----
// Global functions
//-----
bool g_load_config(const std::string& path, vctree_set::Settings& vh);
bool g_set_config(vctree_set::Settings& vh);

