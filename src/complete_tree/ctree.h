#pragma once
#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>   // C++17
#include <fstream>
#include <iomanip>
#include <iostream>
//#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>
#include "global/allvar.h"

#ifdef CTREE_USE_MPI
	#include <mpi.h>
#endif
#ifdef CTREE_USE_OMP
	#include <omp.h>
#endif

namespace Ctree{

	// Data Type
	using CT_I32 	= std::int32_t;
	using CT_I64 	= std::int64_t;
	using CT_float	= float;
	using CT_double = double;

	using CT_Merit 	= CT_float;

	// Specify data type
	using CT_ID 	= IO_dtype::IO_GID;//CT_I32;
	using CT_snap 	= IO_dtype::IO_Snap;//CT_I32;
	using CT_PID 	= IO_dtype::IO_PID;//CT_I64;

	// Branch Controller Array
	struct ListSt{
		CT_Merit merit = 0.;
		CT_ID id = -1;
		CT_snap snap = -1;
	};

	// For log
	struct ctree_num{
		CT_I32 T = 0;
		CT_I32 C = 0;
		CT_I32 B = 0;
	};

	// Particle ID list
	struct PIDSt{
		CT_PID pid;
		CT_ID gid;
		CT_Merit weight;
		CT_I32 n_con;
		CT_I32 maxgid;
	};

	using PIDArray = std::vector<PIDSt>;

	// Snapshot particles
	struct SnapPT{
		PIDArray pid;
		std::vector<CT_I32> gid, n_ptcl, n_ptcl2;
		std::vector<CT_I32> hash;
		std::vector<CT_I32> hash_next;
		CT_I32 maxgid;
		
		//CT_PID maxpid;
	};

	struct ControlSt{
		// ID & Snap at the final snapshot
		CT_ID id0 = 0;
		CT_snap snap0 = 0;

		// Id & Snap at the ending point
		CT_ID id = 0;
		CT_snap snap = 0;

		// Position & Velocity
		CT_double pos[3];
		CT_double vel[3];

		CT_I32 stat = 1; 		// Branch status: 1 (tree exists) 0 (to be connected) -1 (broken)
		CT_I32 detstat = -1; 	// Connection status: 1 (specified) -1 (not yet specified)
		PIDArray p_list; // particle list at the ending point
		CT_I32 n_ptcl = -1; 	// number of particles

		std::vector<ListSt> list; // list of information for the connection candidiates
		CT_I32 list_n = 0;

		CT_I32 last_ind = -1; 	// only stored in the first array

		// Initialize List
    	explicit ControlSt(vctree_set::Settings& vh)
        	: list(vh.ctree_n_search) 
    		{}

    	// Free p_list
    	void Free_plist() {
        	PIDArray().swap(p_list);  // zero capacity
        	n_ptcl = 0;
    	}
	};

	using ControlArray = std::vector<ControlSt>;

	using ControlKey = std::vector<CT_I32>;

	template <class T> struct MeritMatrix {
    	CT_I32 rows{}, cols{};
    	std::vector<T> data;

    	MeritMatrix(CT_I32 r, CT_I32 c, const T& value = T{})
        	: rows(r), cols(c), data(r * c, value) {}

    	// 2D Index
    	T& operator()(CT_I32 r, CT_I32 c){ return data[r * cols + c]; }
    	const T& operator()(CT_I32 r, CT_I32 c) const { return data[r * cols + c]; }
    };


    struct NextSt{
    	CT_Merit merit = 0.;
    	CT_ID id = 0;
    	CT_snap snap = 0;
    };


    using NextArray = std::vector<NextSt>;

    struct CollectSt{
    	CT_snap snap = -1;
    	CT_ID id = -1;
    	CT_I32 lind = 0;
    };

    using CollectArray = std::vector<CollectSt>;

    struct GatherLink{
    	std::vector<CT_I32> ind;
    };

    using GatherLinkArray = std::vector<GatherLink>;

    

    struct CheckSt{
    	CT_Merit merit = 0;
    	CT_ID id = 0;
    	CT_snap snap = 0;
    	CT_ID id0 = 0;
    	CT_snap snap0 = 0;
    	CT_I32 ind = 0;
    };

    using CheckArray = std::vector<CheckSt>;

    struct MeritSt{
    	CT_Merit merit = -1.;
    	CT_ID id = -1;
    };

    struct EndSt{
    	std::vector<Tree::Tree_BID> keyval;
    	std::vector<Tree::Tree_BID> fbid;
    	std::vector<CT_Merit> merit;
    	CT_I32 nn = 0;
    };

    using EndArray = std::vector<EndSt>;

#ifdef CTREE_USE_MPI
    //----- Link MPI
    struct LinkJob{
    	CT_I32 jobnum = -1;

    	CT_I32 ind = -1;
    	CT_I32 snap = -1;
    	CT_I32 id = -1;
    	CT_Merit merit = -1.;
    	CT_I32 dind = -1;
    	Tree::Tree_BID tind = -1;
    	Tree::Tree_BID tind2 = -1;
    	//Tree::Tree_BID bid;
    	//Tree::TreeSt c_tree;

    	//CT_I32 treeinput_ind = -1;
    };
    using JobArray = std::vector<LinkJob>;


    void re_dkey(ControlKey& dkey, CT_I32 ind);
    void in_dkey(ControlKey& dkey, CT_I32 snap, CT_I32 id, CT_I32 ind);
    CT_I32 get_dkey(ControlKey& dkey, CT_I32 snap, CT_I32 id);


    CT_Merit link_commerit(vctree_set::Settings& vh, Tree::TreeSt tree0, CT_I32 snap_to_link, CT_I32 id_to_link, CT_I32 snap_curr);

    void init_job(LinkJob& thisjob);
    LinkJob get_job(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, ControlKey& dkey, CT_I32 ind, CT_I32 snap_to_link, CT_I32 id_to_link, CT_Merit merit_to_link, CT_I32 snap_curr);
    JobArray commque(MPI_Datatype& LINKJOB_T, LinkJob& job, int owner);

    void DoJob1(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, LinkJob& job, CT_I32 snap_curr);
    void DoJob2(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, ControlKey& dkey, LinkJob& job, CT_I32 snap_curr);
    void DoJob3(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, LinkJob& job, CT_I32 snap_curr);
    void DoJob4(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, ControlKey& dkey, LinkJob& job, IO::snapinfo& sinfo, CT_I32 snap_curr);

    void DoJob1a(vctree_set::Settings& vh, ControlArray& data, LinkJob& job, CT_I32 snap_curr);
    void DoJob1b(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, LinkJob& job);

    void DoJob2a(vctree_set::Settings& vh, Tree::TreeArray& tree, ControlArray& data, LinkJob& job, CT_I32 snap_curr);
	void DoJob2b(vctree_set::Settings& vh, ControlArray& data, CT_I32 dind, CT_I32 snap_curr);
	void DoJob2c(Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, LinkJob& job);
	void DoJob2d(Tree::TreeArray& tree, LinkJob& job);

	void DoJob3a(vctree_set::Settings& vh, ControlArray& data, LinkJob& job, CT_I32 snap_curr);
	void DoJob3b(Tree::TreeArray& tree, Tree::Tree_BID org_bid, Tree::Tree_BID com_bid);

	void DoJob4a(vctree_set::Settings& vh, Tree::TreeArray& tree, ControlArray& data, LinkJob& job, CT_I32 snap_curr);
	void DoJob4b(Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, LinkJob& job);
	void DoJob4c(Tree::TreeArray& tree, Tree::TreeKeyArray& key, LinkJob& job);
	void DoJob4d(vctree_set::Settings& vh, ControlArray& data, CT_I32 dind, CT_I32 snap_curr);
	void DoJob4e(vctree_set::Settings& vh, ControlArray& data, CT_I32 dind, CT_I32 snap_curr, CT_I32* jobtype);

	void check_overlap(MPI_Datatype& NEXT_T, LinkJob& thisjob, CT_I32 jobind, CT_I32 ind, std::vector<CT_I32>& cut, NextArray& next_point, std::vector<CT_I32>& islink, CT_I32 rank_index);
	void filter_sort(std::vector<std::pair<CT_I32, CT_I32>>& pairs);
	void syn_data(LinkJob& thisjob, ControlArray& data, ControlKey& dkey);
	void syn_tree(LinkJob& thisjob, Tree::TreeArray& tree, Tree::TreeKeyArray& key);
#endif
	//-----

	void reallocate(vctree_set::Settings& vh, ControlArray& data, CT_I32 nn);

	inline ControlArray allocate(vctree_set::Settings& vh, CT_I32 nn){

		std::vector<ControlSt> controls;
		controls.reserve(nn);

		// OMP HERE?
		for (CT_I32 i = 0; i < nn; i++) controls.emplace_back(vh);


		return controls;
	}

	void inputgal(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, ControlKey& dkey, IO_dtype::GalArray& gal);

	// submoudles
	void classify(vctree_set::Settings& vh, ControlArray& data, IO::snapinfo& sinfo, CT_I32 snap_curr, ctree_num& number);
	
	PIDArray collectpid(vctree_set::Settings& vh, ControlArray& data, Tree::TreeArray& tree, Tree::TreeKeyArray& key);
	PIDArray collectpidalongbranch(vctree_set::Settings& vh, Tree::TreeSt& tree0, CT_I32 snap0, bool opposite);

	SnapPT readsnap(vctree_set::Settings& vh, ControlArray& data, CT_I32 snap_curr);
	SnapPT readsnap2(vctree_set::Settings& vh, CT_I32 snap_curr);

	void commerit(vctree_set::Settings& vh, ControlArray& data, Tree::TreeArray& tree, Tree::TreeKeyArray& key, SnapPT& pid0, CT_I32 snap_curr);

	void link(vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey, Tree::TreeArray& tree, Tree::TreeKeyArray& key, IO::snapinfo& sinfo, CT_I32 snap_curr);

	void addgal(vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey, Tree::TreeArray& tree, Tree::TreeKeyArray& key, CT_I32 snap_curr);
	void delgal(vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey);
	void finalize(vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey, Tree::TreeArray& tree, Tree::TreeKeyArray& key, IO::snapinfo& sinfo, CT_I32 snap_curr, ctree_num& number);
	// ETC
	void makenewbr(vctree_set::Settings& vh, ControlArray& data, CT_I32 ind, CT_snap snap0, CT_ID id0, Tree::TreeArray& tree, Tree::TreeKeyArray& key);
	void expandbr(vctree_set::Settings& vh, ControlArray& data, CT_I32 ind, Tree::TreeArray& tree, Tree::TreeKeyArray& key, CT_I32 id_to_link, CT_I32 snap_to_link, CT_Merit merit_to_link);
	void linkbr(vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey, CT_I32 ind, IO::snapinfo& sinfo, Tree::TreeArray& tree, Tree::TreeKeyArray& key, CT_I32 id_to_link, CT_I32 snap_to_link, CT_Merit merit_to_link, CT_I32 snap_curr);

	CT_Merit get_merit2(std::vector<CT_PID>& pid0, std::vector<CT_PID>& pid, CT_I32 merittype);

	MeritSt get_merit3(SnapPT& pid0, PIDArray& pid, CT_I32 merittype);

	void ctfree(vctree_set::Settings& vh, ControlArray& data, CT_I32 ind, CT_I32 s_end, CT_I32 id_end, CT_I32 snap0);
	CT_I32 wheresnap(IO::snapinfo& sinfo, CT_I32 snap_curr);
	CT_I32 findnextsnap(IO::snapinfo& sinfo, CT_I32 snap_curr);
	std::vector<CT_I32> get_control(ControlArray& data, CT_I32 type);
	void PIDReallocate(vctree_set::Settings& vh, PIDArray& pid, CT_I32 ind);
	std::vector<CT_Merit> get_weight(vctree_set::Settings& vh, std::vector<CT_PID> pid);
	PIDArray get_coreptcl(vctree_set::Settings& vh, PIDArray& pid);

	CT_Merit brcompare(vctree_set::Settings& vh, CT_I32 s0, CT_I32 id0, Tree::TreeSt& tree0, CT_I32 snap0);
	// For main
	void main(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key);
                        //Tree& tree,
                        //Evoldum& evoldum,
                        //Treelog& treelog);

	void findmerge(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key);
	//bool tfcname(vctree_set::Settings& vh);
	//int findnextsnap(vctree_set::Settings& vh, const int snap_curr);
	// For MPI Helper
#ifdef CTREE_USE_MPI
	// ---- raw POD append / read ----
	template<typename T> void append_pod(std::vector<std::uint8_t>& buf, const T& v);

	template<typename T> T read_pod(const std::uint8_t*& p, const std::uint8_t* end);

	template<typename T> void append_vec(std::vector<std::uint8_t>& buf, const std::vector<T>& v);

	template<typename T> std::vector<T> read_vec(const std::uint8_t*& p, const std::uint8_t* end);

	// ---- Serialize / Deserialize for ControlArray ----
	std::vector<std::uint8_t> serialize_d1(const ControlSt& A);
	void deserialize_d1(const std::vector<std::uint8_t>& buf, ControlSt& A);

	// ---- Serialize / Deserialize for ControlArray ----
	std::vector<std::uint8_t> serialize_d2(const ControlSt& A);
	void deserialize_d2(const std::vector<std::uint8_t>& buf, ControlSt& A);


	// ---- Serialize / Deserialize for NextArray ----
	std::vector<std::uint8_t> serialize(const NextSt& x);
	void deserialize(const std::vector<std::uint8_t>& buf, NextSt& out);

	// ---- Serialize / Deserialize for CT_I32 array ----
	std::vector<std::uint8_t> serialize(const CT_I32& x);
	void deserialize(const std::vector<std::uint8_t>& buf, CT_I32& out);

	// ---- Serialize / Deserialize for Job Array ----
	std::vector<std::uint8_t> serialize(const LinkJob& x);
	void deserialize(const std::vector<std::uint8_t>& buf, LinkJob& out);

	// ---- Serialize / Deserialize for Control Array ----
	std::vector<std::uint8_t> serialize(const ControlSt& x);
	void deserialize(const std::vector<std::uint8_t>& buf, ControlSt& out);

	// ---- Serialize / Deserialize for Tree Array ----
	std::vector<std::uint8_t> serialize(const Tree::TreeSt& x);
	void deserialize(const std::vector<std::uint8_t>& buf, Tree::TreeSt& out);

	// ---- Serialize / Deserialize for EndPoints Array ----
	std::vector<std::uint8_t> serialize_ends(const EndSt& x);
	void deserialize_ends(const std::vector<std::uint8_t>& buf, EndSt& out);

	void bcast_blob_from_owner(int owner, std::vector<std::uint8_t>& blob);

    static inline MPI_Datatype make_next_type(){
    	

    	MPI_Datatype CT_MERIT_T;
    	if(sizeof(CT_Merit)==8){
    		CT_MERIT_T = MPI_DOUBLE;
    	}else if(sizeof(CT_Merit)==4){
    		CT_MERIT_T = MPI_FLOAT;
    	}

    	MPI_Datatype CT_ID_T;
    	if(sizeof(CT_ID)==8){
    		CT_ID_T = MPI_INT64_T;
    	}else if(sizeof(CT_ID)==4){
    		CT_ID_T = MPI_INT32_T;
    	}

    	MPI_Datatype CT_Snap_T;
    	if(sizeof(CT_snap)==8){
    		CT_Snap_T = MPI_INT64_T;
    	}else if(sizeof(CT_snap)==4){
    		CT_Snap_T = MPI_INT32_T;
    	}
    	
    	NextSt d;
    	MPI_Aint base, disp[3];
    	int blocklen[3]	= {1, 1, 1};
    	MPI_Datatype types[3] = {CT_MERIT_T, CT_ID_T, CT_Snap_T};

    	MPI_Get_address(&d, 			&base);
    	MPI_Get_address(&d.merit, 		&disp[0]);
    	MPI_Get_address(&d.id, 			&disp[1]);
    	MPI_Get_address(&d.snap, 		&disp[2]);

    	for (int i=0;i<3;++i) disp[i] -= base;

    	MPI_Datatype NEXT_T;
    	MPI_Type_create_struct(3, blocklen, disp, types, &NEXT_T);
    	MPI_Type_commit(&NEXT_T);
    	return NEXT_T;
    	//
    	

    }

	static inline MPI_Datatype make_linkjob_type(){
		      
    	MPI_Datatype CT_I32_T = MPI_INT32_T;

    	// CT_Merit for double-precision
    	MPI_Datatype CT_MERIT_T;
    	if(sizeof(CT_Merit)==8){
    		CT_MERIT_T = MPI_DOUBLE;
    	}else if(sizeof(CT_Merit)==4){
    		CT_MERIT_T = MPI_FLOAT;
    	}

    	MPI_Datatype CT_BID_T;
    	if(sizeof(Tree::Tree_BID)==8){
    		CT_BID_T = MPI_INT64_T;
    	}else if(sizeof(Tree::Tree_BID)==4){
    		CT_BID_T = MPI_INT32_T;
    	}


    	// LinkJob
    	LinkJob dummy;
    	int          blocklen[8] = {1,1,1,1,1,1,1,1};
    	MPI_Datatype types[8]    = {CT_I32_T, CT_I32_T, CT_I32_T, CT_I32_T, CT_I32_T, CT_MERIT_T, CT_BID_T, CT_BID_T};
    	MPI_Aint     disp[8], base;

    	MPI_Get_address(&dummy,           &base);
    	MPI_Get_address(&dummy.jobnum,    &disp[0]);
    	MPI_Get_address(&dummy.ind,       &disp[1]);
    	MPI_Get_address(&dummy.snap,      &disp[2]);
    	MPI_Get_address(&dummy.id,        &disp[3]);
    	MPI_Get_address(&dummy.dind,      &disp[4]);
    	MPI_Get_address(&dummy.merit,     &disp[5]);
    	MPI_Get_address(&dummy.tind,   	  &disp[6]);
    	MPI_Get_address(&dummy.tind2,  	  &disp[7]);
    	

	    for (int i=0;i<8;++i) disp[i] -= base;

    	MPI_Datatype LINKJOB_T;
    	MPI_Type_create_struct(8, blocklen, disp, types, &LINKJOB_T);
    	MPI_Type_commit(&LINKJOB_T);
    	return LINKJOB_T;
	}
#endif

	// For debugging
	void savetree_ctree(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& treekey, CT_I32 snap_curr);
	void loadtree_ctree(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& treekey, CT_I32 snap_curr);
	void savedata(vctree_set::Settings& vh, ControlArray& data, CT_I32 snap_curr);
	ControlArray loaddata(vctree_set::Settings& vh, CT_I32 snap_curr);
	void validate_data(ControlArray& data, ControlArray& data2);
	void validate_tree(Tree::TreeArray& tree, Tree::TreeArray& tree2);
	void validate_treekey(Tree::TreeKeyArray& key, Tree::TreeKeyArray& key2);
}