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

	// Specify data type
	using CT_ID 	= IO_dtype::IO_GID;//CT_I32;
	using CT_snap 	= IO_dtype::IO_Snap;//CT_I32;
	using CT_PID 	= IO_dtype::IO_PID;//CT_I64;

	// Branch Controller Array
	struct ListSt{
		CT_double merit = 0.;
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
		CT_double weight;
		CT_I32 n_con;
		CT_I32 maxgid;
	};

	using PIDArray = std::vector<PIDSt>;

	// Snapshot particles
	struct SnapPT{
		PIDArray pid;
		std::vector<CT_I32> gid, n_ptcl;
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

		// Settings을 받아 list를 '생성 시' 원하는 크기로 초기화
    	explicit ControlSt(const vctree_set::Settings& vh)
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
    	CT_double merit = 0.;
    	CT_ID id = 0;
    	CT_snap snap = 0;
    };


    using NextArray = std::vector<NextSt>;

    struct CheckSt{
    	CT_double merit = 0;
    	CT_ID id = 0;
    	CT_snap snap = 0;
    	CT_ID id0 = 0;
    	CT_snap snap0 = 0;
    	CT_I32 ind = 0;
    };

    using CheckArray = std::vector<CheckSt>;

	//-----

	void reallocate(const vctree_set::Settings& vh, ControlArray& data, CT_I32 nn);

	inline ControlArray allocate(const vctree_set::Settings& vh, CT_I32 nn){

		std::vector<ControlSt> controls;
		controls.reserve(nn);

		// OMP HERE?
		for (CT_I32 i = 0; i < nn; i++) controls.emplace_back(vh);


		return controls;
	}

	void inputgal(const vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key, ControlArray& data, ControlKey& dkey, IO_dtype::GalArray& gal);

	// submoudles
	void classify(const vctree_set::Settings& vh, ControlArray& data, IO::snapinfo& sinfo, CT_I32 snap_curr, ctree_num& number);
	
	PIDArray collectpid(const vctree_set::Settings& vh, ControlArray& data, Tree::TreeArray& tree, Tree::TreeKeyArray& key);
	PIDArray collectpidalongbranch(const vctree_set::Settings& vh, std::vector<CT_snap>& slist, std::vector<CT_ID>& glist);

	SnapPT readsnap(const vctree_set::Settings& vh, ControlArray& data, CT_I32 snap_curr);
	
	void commerit(const vctree_set::Settings& vh, ControlArray& data, PIDArray& pid, SnapPT& pid0, CT_I32 snap_curr);

	void link(const vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey, Tree::TreeArray& tree, Tree::TreeKeyArray& key, IO::snapinfo& sinfo, CT_I32 snap_curr);

	void addgal(const vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey, Tree::TreeArray& tree, Tree::TreeKeyArray& key, CT_I32 snap_curr);

	void finalize(const vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey, Tree::TreeArray& tree, Tree::TreeKeyArray& key, IO::snapinfo& sinfo, CT_I32 snap_curr, ctree_num& number);
	// ETC
	void makenewbr(const vctree_set::Settings& vh, ControlArray& data, CT_I32 ind, CT_snap snap0, CT_ID id0, Tree::TreeArray& tree, Tree::TreeKeyArray& key);
	void expandbr(const vctree_set::Settings& vh, ControlArray& data, CT_I32 ind, Tree::TreeArray& tree, Tree::TreeKeyArray& key, CT_I32 id_to_link, CT_I32 snap_to_link, CT_double merit_to_link);
	void linkbr(const vctree_set::Settings& vh, ControlArray& data, ControlKey& dkey, CT_I32 ind, IO::snapinfo& sinfo, Tree::TreeArray& tree, Tree::TreeKeyArray& key, CT_I32 id_to_link, CT_I32 snap_to_link, CT_double merit_to_link, CT_I32 snap_curr);

	void get_merit(std::vector<CT_PID>& pid_g, std::vector<CT_I32>& gid_g, 
		std::vector<CT_PID>& pid_s, std::vector<CT_I32>& gid_s, 
		std::vector<CT_I32>& hash, std::vector<CT_I32>& hash_next, 
		std::vector<CT_I32>& npart_g, std::vector<CT_I32>& npart_s, 
		std::vector<CT_I32>& match_id, std::vector<CT_double>& match_merit);

	CT_double get_merit2(std::vector<CT_PID>& pid0, std::vector<CT_PID>& pid, std::vector<CT_double>& w0, std::vector<CT_double>& w1, CT_I32 merittype);

	void ctfree(const vctree_set::Settings& vh, ControlArray& data, CT_I32 ind, CT_I32 s_end, CT_I32 id_end, CT_I32 snap0);
	CT_I32 wheresnap(IO::snapinfo& sinfo, CT_I32 snap_curr);
	std::vector<CT_I32> get_control(ControlArray& data, CT_I32 type);
	void PIDReallocate(const vctree_set::Settings& vh, PIDArray& pid, CT_I32 ind);
	std::vector<CT_double> get_weight(const vctree_set::Settings& vh, std::vector<CT_PID> pid);
	PIDArray get_coreptcl(const vctree_set::Settings& vh, PIDArray& pid);

	CT_double brcompare(const vctree_set::Settings& vh, CT_I32 s0, CT_I32 id0, std::vector<CT_I32>& slist, std::vector<CT_I32>& idlist);
	// For main
	void main(const vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key);
                        //Tree& tree,
                        //Evoldum& evoldum,
                        //Treelog& treelog);


	//bool tfcname(const vctree_set::Settings& vh);
	//int findnextsnap(const vctree_set::Settings& vh, const int snap_curr);
	// For MPI Helper
#ifdef CTREE_USE_MPI
	// ---- raw POD append / read ----
	template<typename T> void append_pod(std::vector<std::uint8_t>& buf, const T& v);

	template<typename T> T read_pod(const std::uint8_t*& p, const std::uint8_t* end);

	template<typename T> void append_vec(std::vector<std::uint8_t>& buf, const std::vector<T>& v);

	template<typename T> std::vector<T> read_vec(const std::uint8_t*& p, const std::uint8_t* end);

	// ---- Serialize / Deserialize for ControlArray ----
	std::vector<std::uint8_t> serialize(const ControlSt& A);
	void deserialize(const std::vector<std::uint8_t>& buf, ControlSt& A);

	// ---- Serialize / Deserialize for NextArray ----
	std::vector<std::uint8_t> serialize(const NextSt& x);
	void deserialize(const std::vector<std::uint8_t>& buf, NextSt& out);

	// ---- Serialize / Deserialize for CT_I32 array ----
	std::vector<std::uint8_t> serialize(const CT_I32& x);
	void deserialize(const std::vector<std::uint8_t>& buf, CT_I32& out);

	void bcast_blob_from_owner(int owner, std::vector<std::uint8_t>& blob);


#endif

}