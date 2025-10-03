#pragma once
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>   // C++17
#include <fstream>
#include <hdf5.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <omp.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include "global/allvar.h"
#include "global/io.h"

#ifdef CTREE_USE_MPI
	#include <mpi.h>
#endif


namespace {

    // Extract Integer from tree.snaplist.txt
    inline bool extract_snap_from_line(const std::string& line, int& out) {
        auto p = line.find("snap_");
        if (p == std::string::npos) return false;
        if (p + 5 + 4 > line.size()) return false;            // 최소 4자리
        try {
            out = std::stoi(line.substr(p + 5, 4));
            return true;
        } catch (...) { return false; }
    }

    // Extract Integer "tree.snapshot_0123.tree"
    inline bool extract_number_from_filename(const std::string& filename, int& out) {
        auto us = filename.find_last_of('_');
        auto dot = filename.find_last_of('.');
        if (us == std::string::npos || dot == std::string::npos || us+1 >= dot) return false;
        try {
            out = std::stoi(filename.substr(us + 1, dot - (us + 1)));
            return true;
        } catch (...) { return false; }
    }

    
} // namespace




namespace Makebr{
	// Data Type
	using MBR_GID 		= std::int32_t;		// for galaxy ID
	using MBR_Snap 		= std::int32_t;
	using MBR_PID 	 	= std::int64_t;		// for particle ID
	using MBR_merit		= double;
	using MBR_BID 		= std::int32_t;	


	//-----
	using TF_num		= std::int32_t;
	using TF_off 		= std::int64_t;
	using TF_res 		= std::int64_t;
	using TF_merit		= float;
	using TF_npart 		= std::int64_t;
	using TF_nlink 		= std::int64_t;
	using TF_id 		= std::int64_t;

	struct TFSt{
		std::vector<TF_num> num;
		std::vector<TF_off> off;
		std::vector<TF_res> res;
		std::vector<TF_merit> merit;
		std::vector<TF_npart> npart;
		std::vector<TF_nlink> nlink;
		std::vector<TF_id> id;
	};

	//using TFArr = std::vector<TFSt>;


	// For Makebr Log
	struct Treelog{
		int32_t n_new = 0;
		int32_t n_link = 0;
		int32_t n_link2 = 0;
		int32_t n_link3 = 0;
		int32_t n_broken = 0;
		int32_t n_all = 0;
		int32_t max_id = -1;
		Tree::Tree_BID lind = 0;
	};

	// For Makebr settings
	struct HeaderSt{
		MBR_Snap snapi;
		MBR_Snap snapf;

		std::string tag_id = "ID";
		std::string tag_num;
      	std::string tag_off;
      	std::string tag_result;
      	std::string tag_npart;
      	std::string tag_merit;
      	std::string tag_nlink;
	};

	using MBRHeader = HeaderSt;

	struct EvolSt{
		IO_VR::VRT_GID idc, idn;
		IO_VR::VRT_Snap snapc, snapn;
		IO_VR::VRT_merit merit;
	};

	using EvolArray = std::vector<EvolSt>;

	

	// For Main
	void mainloop(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key);

	

	// For Utilities
	IO_VR::VRT_GID getmaxid(const IO_VR::GalArray& gal);
	MBRHeader setheader(const vctree_set::Settings& vh);
	bool is_snap(const vctree_set::Settings& vh, const IO_VR::VRT_Snap snap_curr);
	MBR_Snap findnextsnap(const vctree_set::Settings& vh, const MBR_Snap& snap_curr);
	bool tfcname(const vctree_set::Settings& vh);

	TFSt readtreefrog(const vctree_set::Settings& vh, const MBRHeader& mbh, const MBR_Snap& snap_curr);

	template <typename T> std::vector<T> slice(const std::vector<T>& v, std::size_t ind1, std::size_t ind2_inclusive);

	// For MPI Helper
#ifdef CTREE_USE_MPI
	// ---- raw POD append / read ----
	template<typename T> void append_pod(std::vector<std::uint8_t>& buf, const T& v);

	template<typename T> T read_pod(const std::uint8_t*& p, const std::uint8_t* end);

	template<typename T> void append_vec(std::vector<std::uint8_t>& buf, const std::vector<T>& v);

	template<typename T> std::vector<T> read_vec(const std::uint8_t*& p, const std::uint8_t* end);

	// ---- Serialize / Deserialize for TFSt ----
	std::vector<std::uint8_t> serialize(const Makebr::TFSt& x);
	void deserialize(const std::vector<std::uint8_t>& buf, Makebr::TFSt& out);

	// ---- Serialize / Deserialize for GalArray ----
	std::vector<std::uint8_t> serialize(const IO_VR::GalArray& A);
	void deserialize(const std::vector<std::uint8_t>& buf, IO_VR::GalArray& A);

	// ---- Serialize / Deserialize for EvolDum ----
	std::vector<std::uint8_t> serialize(const std::vector<std::pair<TF_id, EvolSt>>& upd);
	void deserialize(const std::vector<std::uint8_t>& blob, EvolArray& evoldum);


	void bcast_blob_from_owner(int owner, std::vector<std::uint8_t>& blob);


#endif

}













