#pragma once
#include "global/allvar.h"
//#include "global/io_hdf5.h"
#include "utilities/utilities.h"
#include <algorithm>
#include <filesystem>   // C++17
#include <fstream>
#include <hdf5.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


namespace fs = std::filesystem;
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

    // padding
    inline std::string i4(int x) {
        std::ostringstream oss;
        oss << std::setw(4) << std::setfill('0') << x;
        return oss.str();
    }

    inline std::string i6(int x) {
    	std::ostringstream oss;
        oss << std::setw(6) << std::setfill('0') << x;
        return oss.str();
    }
} // namespace

namespace IO{

	//-----
	// Adjust Configuration File
	//-----
	// Modulate Settings based on the tree direction
	bool set_config_vr(vctree_set::Settings& vh);
	bool set_config(vctree_set::Settings& vh);

	//-----
	// Read Tree Frog Information (if exists) 
	//		If not, null data given and all galaxies are treated as independent tree
	//-----
	auto r_treefrog(const vctree_set::Settings& vh, const int snap_curr);

	//-----
	// read Galaxy Catalog
	//-----
	auto r_gal(const vctree_set::Settings& vh, const int snap_curr, const long long id0, const bool readpart=false);


    //-----
    // Get Snapshot Infomation
    //-----
    auto g_snapinfo(const vctree_set::Settings& vh);


	//-----
	// Find next snapshot (depends on the tree direction)
	//-----
	int findnextsnap(const vctree_set::Settings& vh, const int snap_curr);

	//-----
	// Sink the TreeFrog output name to the snapshot list
	//-----
	bool tfcname(const vctree_set::Settings& vh);
}

namespace rd_vrhdf5{

	struct TreeSetTags {
    	std::string tag_num;
    	std::string tag_off;
    	std::string tag_result;
    	std::string tag_merit;
    	std::string tag_npart;
    	std::string tag_nlink; // attribute on root "."
	};


	// -------------------------------
	// HDF5 native type mapping
	// -------------------------------
	template<typename T> inline hid_t h5_native_type();
	template<> inline hid_t h5_native_type<int8_t>()   { return H5T_NATIVE_INT8;  }
	template<> inline hid_t h5_native_type<uint8_t>()  { return H5T_NATIVE_UINT8; }
	template<> inline hid_t h5_native_type<int16_t>()  { return H5T_NATIVE_INT16; }
	template<> inline hid_t h5_native_type<uint16_t>() { return H5T_NATIVE_UINT16;}
	template<> inline hid_t h5_native_type<int32_t>()  { return H5T_NATIVE_INT32; }
	template<> inline hid_t h5_native_type<uint32_t>() { return H5T_NATIVE_UINT32;}
	template<> inline hid_t h5_native_type<int64_t>()  { return H5T_NATIVE_INT64; }
	template<> inline hid_t h5_native_type<uint64_t>() { return H5T_NATIVE_UINT64;}
	template<> inline hid_t h5_native_type<float>()    { return H5T_NATIVE_FLOAT; }
	template<> inline hid_t h5_native_type<double>()   { return H5T_NATIVE_DOUBLE; }

	//-----
	// Shared Helper
	//-----

	std::string zfill_longlong(long long x, int width);
	bool attribute_exists(hid_t file_or_obj_id, const std::string& object_path, const std::string& attr_name);
	bool attribute_exists(hid_t file_or_obj_id, const std::string& object_path, const std::string& attr_name);

	template<typename T>
	std::vector<T> read_1d_dataset(hid_t file_id, const std::string& path);

	template<typename T>
	T read_scalar_dataset(hid_t file_id, const std::string& path);

	template<typename T>
	std::vector<T> read_attribute_by_name(hid_t file_or_obj_id, const std::string& object_path, const std::string& attr_name);

	template<typename T>
	T read_scalar_dset_or_attr(hid_t file_id, const std::string& object_path, const std::string& name);

	template<typename IdT>
	std::vector<IdT> read_ids_dataset_or_attr(hid_t file_id);

	// -------------------------------
	// Makebr::TreeVec Reader
	// -------------------------------
	//template<typename NumT, typename OffT, typename ResT, typename MerT,
	//         typename NpartT, typename NlinkT, typename IdT>

	template<typename NumT, typename OffT, typename ResT, typename MerT,
         typename NpartT, typename NlinkT, typename IdT>
	struct TreeData {
    	std::vector<NumT>   num;
    	std::vector<OffT>   off;
    	std::vector<ResT>   res;
    	std::vector<MerT>   merit;
    	std::vector<NpartT> npart;
    	std::vector<NlinkT> nlink; // attribute (scalar -> size 1)
    	std::vector<IdT>    id;
	};

	
	template<typename NumT, typename OffT, typename ResT, typename MerT,
	         typename NpartT, typename NlinkT, typename IdT>
	struct TreeArrays_Makebr {
	    Makebr::TreeVec<NumT>   num;
	    Makebr::TreeVec<OffT>   off;
	    Makebr::TreeVec<ResT>   res;
	    Makebr::TreeVec<MerT>   merit;
	    Makebr::TreeVec<NpartT> npart;
	    Makebr::TreeVec<NlinkT> nlink;
	    Makebr::TreeVec<IdT>    id;
	};


	

	template<typename NumT, typename OffT, typename ResT, typename MerT,
	         typename NpartT, typename NlinkT, typename IdT>
	inline TreeArrays_Makebr<NumT,OffT,ResT,MerT,NpartT,NlinkT,IdT> read_treeset_makebr(const std::string& fname, const TreeSetTags& tags, std::size_t growth_step = 1);

	//-----
	// To Read VR Catalog
	//-----
	//==================== r_gal ====================//
	// id0 < 0  → All galaxies
	// id0 >= 0 → galaxy with the corresponding ID
	template<typename IdT, typename NpartT, typename PidT>
	struct GalArrays {
	    Makebr::TreeVec<IdT>    id;
	    Makebr::TreeVec<NpartT> npart;
	    Makebr::TreeVec<PidT>   p_id; // 선택적 로드
	};

	template<typename IdT=int32_t, typename NpartT=int32_t, typename PidT=int64_t>
	inline GalArrays<IdT,NpartT,PidT> r_gal(const std::string& fname, long long id0 = -1, bool read_p_id = false);
}



