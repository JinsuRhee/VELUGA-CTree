#include "global/io.h"
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

//namespace fs = std::filesystem;
//namespace {
//
//    // Extract Integer from tree.snaplist.txt
//    inline bool extract_snap_from_line(const std::string& line, int& out) {
//        auto p = line.find("snap_");
//        if (p == std::string::npos) return false;
//        if (p + 5 + 4 > line.size()) return false;            // 최소 4자리
//        try {
//            out = std::stoi(line.substr(p + 5, 4));
//            return true;
//        } catch (...) { return false; }
//    }
//
//    // Extract Integer "tree.snapshot_0123.tree"
//    inline bool extract_number_from_filename(const std::string& filename, int& out) {
//        auto us = filename.find_last_of('_');
//        auto dot = filename.find_last_of('.');
//        if (us == std::string::npos || dot == std::string::npos || us+1 >= dot) return false;
//        try {
//            out = std::stoi(filename.substr(us + 1, dot - (us + 1)));
//            return true;
//        } catch (...) { return false; }
//    }
//
//    // padding
//    inline std::string i4(int x) {
//        std::ostringstream oss;
//        oss << std::setw(4) << std::setfill('0') << x;
//        return oss.str();
//    }
//
//    inline std::string i6(int x) {
//    	std::ostringstream oss;
//        oss << std::setw(6) << std::setfill('0') << x;
//        return oss.str();
//    }
//} // namespace

namespace IO{

	//-----
	// Adjust Configuration File
	//-----
	// Modulate Settings based on the tree direction
	bool set_config_vr(vctree_set::Settings& vh){
		int32_t snapmin = std::min(vh.snapi, vh.snapf);
  		int32_t snapmax = std::max(vh.snapi, vh.snapf);
		if(vh.treedir == "DES"){
      		vh.snapi  = snapmin;
      		vh.snapf  = snapmax;
  
      		vh.tag_num    = "NumDesc";
      		vh.tag_off    = "DescOffsets";
      		vh.tag_result = "Descendants";
      		vh.tag_npart  = "DescNpart";
      		vh.tag_merit  = "Merits";
      		vh.tag_nlink  = "Nsteps_search_new_links";
    	}
    	else if(vh.treedir == "PRG"){
      		vh.snapi  = snapmax;
      		vh.snapf  = snapmin;
  
      		vh.tag_num    = "NumProgen";
      		vh.tag_off    = "ProgenOffsets";
      		vh.tag_result = "Progenitors";
      		vh.tag_npart  = "ProgenNpart";
      		vh.tag_merit  = "Merits";
      		vh.tag_nlink  = "Nsteps_search_new_links";
  
    	}
    	else{
      		return false;
    	}
    	return true;
	}

	bool set_config(vctree_set::Settings& vh){
		if(vh.iotype == "VR")set_config_vr(vh);
  		else{
    		LOG()<<"Not implemented yet";
    		return false;
  		}
  		return true;
	}

	//-----
	// Read Tree Frog Information (if exists) 
	//		If not, null data given and all galaxies are treated as independent tree
	//-----
	auto r_treefrog(const vctree_set::Settings& vh, const int snap_curr){
		std::string tfname 	= vh.vr_dir_tree + "/tree.snapshot_" + i4(snap_curr) + "VELOCIraptor.tree";

		//using namespace rd_treefrog;

		// 타입은 실제 HDF5 저장 타입에 맞춰 조정하세요.
    	using NumT   = int32_t;
    	using OffT   = int64_t;
    	using ResT   = int64_t; // Descendants가 float/double이면 double로 변경
    	using MerT   = float;
    	using NpartT = int64_t;
    	using NlinkT = int32_t; // attribute
    	using IdT    = int64_t;

  		std::string tag_num     = "";
  		std::string tag_off     = "";
  		std::string tag_result  = "";
  		std::string tag_npart   = "";
  		std::string tag_merit   = "";
  		std::string tag_nlink   = "";


    	rd_vrhdf5::TreeSetTags vt{
        	vh.tag_num,
        	vh.tag_off,
        	vh.tag_result,
        	vh.tag_merit,
        	vh.tag_npart,
        	vh.tag_nlink
    	};

    	
    	auto data = rd_vrhdf5::read_treeset_makebr<NumT,OffT,ResT,MerT,NpartT,NlinkT,IdT>(
        	tfname, vt, /*growth_step=*/1
   		);

   		return data;

	}

	//-----
	// read Galaxy Catalog
	//-----
	auto r_gal(const vctree_set::Settings& vh, const int snap_curr, const long long id0, const bool readpart){
		if(vh.iotype=="VR"){
			std::string fname 	= vh.vr_dir_catalog + "/snap_" + i4(snap_curr) + ".hdf5";
			return rd_vrhdf5::r_gal(fname, id0, readpart);
		}else{
			LOG() << "Not implemented";
			u_stop();
			return rd_vrhdf5::r_gal("", id0, readpart);
		}
	}

	//-----
	// Sink the TreeFrog output name to the snapshot list
	//-----
	bool tfcname(const vctree_set::Settings& vh){

		try {
    		const fs::path tfout = fs::path(vh.vr_dir_tree);

    		// Read snaplist
    		std::vector<int> slist;
    		{
      			std::ifstream in(tfout / "tree.snaplist.txt");
      			if (!in) {
        			std::cerr << "[tfcname] cannot open: " << (tfout / "tree.snaplist.txt") << "\n";
        			return false;
      			}
      			std::string line;
      			while (std::getline(in, line)) {
        			int v;
        			if (extract_snap_from_line(line, v)) slist.push_back(v);
      			}
    		}

    		// list tfile
    		struct Item { int num; fs::path path; };
    		std::vector<Item> items;
    		for (const auto& de : fs::directory_iterator(tfout)) {
      			if (!de.is_regular_file()) continue;
      			const auto fn = de.path().filename().string();
      			if (fn.rfind("tree.snapshot", 0) == 0 && de.path().extension() == ".tree") {
        			int num;
        			if (extract_number_from_filename(fn, num)) {
          				items.push_back({num, de.path()});
        			}
      			}
    		}

    		if (items.size() != slist.size()) {
      			LOG() << "[tfcname] count mismatch: "
                	<< "files=" << items.size()
                	<< ", snaplist=" << slist.size() << "\n";
      			return false;
    		}

    		// Sort
    		std::sort(items.begin(), items.end(),
              [](const Item& a, const Item& b){ return a.num < b.num; });

    		// rename
    		for (std::size_t i = 0; i < items.size(); ++i) {
      			const std::string newbase = "tree.snapshot_" + i4(slist[i]) + "VELOCIraptor";
      			const fs::path target = tfout / newbase;

      			if (items[i].path.filename() == target.filename()) continue; // 이미 같으면 스킵

      			// 충돌 방지: 동일 이름 파일이 있다면 먼저 제거/백업 등 필요 시 처리
      			if (fs::exists(target)) {
        			LOG() << "[tfcname] target already exists, skipping: " << target << "\n";
        			continue;
      			}
      			fs::rename(items[i].path, target);  // mv org -> no-ext name
      			items[i].path = target;             // 경로 갱신
    		}

    		// add .tree
    		for (auto& it : items) {
      			const fs::path target = fs::path(it.path.string() + ".tree");
      			if (fs::exists(target)) {
        			std::cerr << "[tfcname] target already exists, skipping: " << target << "\n";
        			continue;
      			}
      			fs::rename(it.path, target);        // mv name -> name.tree
      			it.path = target;
    		}

    		return true;
  		} catch (const std::exception& e) {
    		std::cerr << "[tfcname] exception: " << e.what() << "\n";
    		return false;
  		}
	}

	//-----
	// Find next snapshot (depends on the tree direction)
	//-----
	int findnextsnap(const vctree_set::Settings& vh, const int snap_curr){
		if(vh.iotype == "VR"){
			int snap_next 	= snap_curr;
			while (true){
				if(vh.treedir == "DES"){
					snap_next ++;
				}else if(vh.treedir == "PRG"){
					snap_next --;
				}else{
					LOG()<<		"Wrong Tree Direction";
					u_stop();
				}

				// 파일명 조립
        		char buf[256];
        		std::snprintf(buf, sizeof(buf), "tree.snapshot_%04dVELOCIraptor.tree", snap_next);
        		std::string dumfname = vh.vr_dir_tree + "/" + buf;
        		
        		// 파일 존재 검사
        		if (fs::exists(dumfname)) {
            		return snap_next;
        		}

        		// 범위 벗어나면 종료
        		if (snap_next < vh.snapi || snap_next > vh.snapf) {
            		return -1;
        		}
			}
		}
		else{
			LOG()<< 		"Not implemented";
			u_stop();
			return -1;
		}
  	}

	//-----
    // Get Snapshot Infomation
    //-----
    auto g_snapinfo(const vctree_set::Settings& vh){

    	// This is the basic structure
    	struct SnapInfo{
    		long 	snap = -1;		// snapshot number
    		double 	aexp = 0.0;		// Scale Factor
    		double 	unit_l = 0.0;	// Length unit
    		double 	unit_t = 0.0; 	// Time unit
    		double 	age = 0.0;		// Age of the Universe at this snapshot
    		double 	kpc	= 0.0;		// 1 kpc as cgs
    	};

    	if(vh.iotype=="VR"){
    		// For the VR (TF type)

    		//----- Get Snapshot list from the TF output
			const fs::path dir = vh.vr_dir_tree;
    		std::vector<long> slist;

    		std::error_code ec;
    		if (!fs::exists(dir, ec) || !fs::is_directory(dir, ec)) {
        		LOG()<<" Directory doesn't exist";
        		u_stop();
    		}

	    	for (const auto& entry : fs::directory_iterator(dir, ec)) {
	        	if (ec) break;
	        	if (!entry.is_regular_file()) continue;

	        	const std::string fname = entry.path().filename().string();

	        	// pattern
	        	constexpr const char* prefix = "tree.snapshot_";
	        	constexpr const char* suffix = "VELOCIraptor.tree";

	        	// pre-check
	        	if (fname.size() < std::char_traits<char>::length(prefix) + std::char_traits<char>::length(suffix)) continue;
	        	if (fname.rfind(prefix, 0) != 0) continue; // prefix로 시작?
	        	if (fname.size() < std::char_traits<char>::length(suffix)) continue;
	        	if (fname.substr(fname.size() - std::char_traits<char>::length(suffix)) != suffix) continue;

	        	const std::size_t i0 = fname.find('_');              // "tree.snapshot_"
	        	const std::size_t i1 = fname.find("VELO");           // "VELOCIraptor.tree"
	        	if (i0 == std::string::npos || i1 == std::string::npos || i1 <= i0 + 1) continue;

	        	const std::string numstr = fname.substr(i0 + 1, i1 - (i0 + 1));
	        	try {
	            	// 숫자(예: "0007" → 7)
	            	long snap = std::stol(numstr);
	            	slist.push_back(snap);
	        	} catch (...) {
	            
	        	}
	    	}

    
    		std::sort(slist.begin(), slist.end());
    		slist.erase(std::unique(slist.begin(), slist.end()), slist.end());

    		for (auto s : slist){
    			LOG()<<s;
    		}




    	}
    	return 1;

    }
}

//-----
// For HDF5 usage
//-----
namespace rd_vrhdf5{
	
	

	

	//-----
	// Shared Helper
	//-----
	inline std::string zfill_longlong(long long x, int width) {
	    std::ostringstream oss;
	    oss << std::setfill('0') << std::setw(width) << x;
	    return oss.str();
	}

	inline bool attribute_exists(hid_t file_or_obj_id,
	                             const std::string& object_path,
	                             const std::string& attr_name) {
	    hid_t obj = H5Oopen(file_or_obj_id, object_path.c_str(), H5P_DEFAULT);
	    if (obj < 0) throw std::runtime_error("Failed to open object: " + object_path);
	    htri_t ex = H5Aexists(obj, attr_name.c_str());
	    H5Oclose(obj);
	    if (ex < 0) throw std::runtime_error("H5Aexists failed on: " + object_path);
	    return ex > 0;
	}

	template<typename T>
	inline std::vector<T> read_1d_dataset(hid_t file_id, const std::string& path) {
	    if (!H5Lexists(file_id, path.c_str(), H5P_DEFAULT))
	        throw std::runtime_error("Dataset not found: " + path);

	    hid_t d = H5Dopen2(file_id, path.c_str(), H5P_DEFAULT);
	    if (d < 0) throw std::runtime_error("Failed to open dataset: " + path);

	    hid_t s = H5Dget_space(d);
	    if (s < 0) { H5Dclose(d); throw std::runtime_error("Failed to get dataspace: " + path); }

	    int nd = H5Sget_simple_extent_ndims(s);
	    if (nd != 1) { H5Sclose(s); H5Dclose(d); throw std::runtime_error("Dataset rank != 1: " + path); }

	    hsize_t dims[1]{};
	    H5Sget_simple_extent_dims(s, dims, nullptr);
	    std::vector<T> buf(static_cast<size_t>(dims[0]));

	    herr_t st = H5Dread(d, h5_native_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
	    H5Sclose(s);
	    H5Dclose(d);
	    if (st < 0) throw std::runtime_error("Failed to read dataset: " + path);
	    return buf;
	}

	// rank 0 스칼라 또는 rank 1 길이 1 허용
	template<typename T>
	inline T read_scalar_dataset(hid_t file_id, const std::string& path) {
	    if (!H5Lexists(file_id, path.c_str(), H5P_DEFAULT))
	        throw std::runtime_error("Dataset not found: " + path);

	    hid_t d = H5Dopen2(file_id, path.c_str(), H5P_DEFAULT);
	    if (d < 0) throw std::runtime_error("Failed to open dataset: " + path);

	    hid_t s = H5Dget_space(d);
	    if (s < 0) { H5Dclose(d); throw std::runtime_error("Failed to get dataspace: " + path); }

	    int nd = H5Sget_simple_extent_ndims(s);
	    if (nd > 1) { H5Sclose(s); H5Dclose(d); throw std::runtime_error("Dataset rank > 1: " + path); }
	    if (nd == 1) {
	        hsize_t dims[1]{};
	        H5Sget_simple_extent_dims(s, dims, nullptr);
	        if (dims[0] != 1) { H5Sclose(s); H5Dclose(d); throw std::runtime_error("Dataset len != 1: " + path); }
	    }

	    T value{};
	    herr_t st = H5Dread(d, h5_native_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
	    H5Sclose(s);
	    H5Dclose(d);
	    if (st < 0) throw std::runtime_error("Failed to read scalar dataset: " + path);
	    return value;
	}

	template<typename T>
	inline std::vector<T> read_attribute_by_name(hid_t file_or_obj_id,
	                                             const std::string& object_path,
	                                             const std::string& attr_name) {
	    hid_t obj = H5Oopen(file_or_obj_id, object_path.c_str(), H5P_DEFAULT);
	    if (obj < 0) throw std::runtime_error("Failed to open object: " + object_path);
	    if (!H5Aexists(obj, attr_name.c_str())) {
	        H5Oclose(obj);
	        throw std::runtime_error("Attribute not found: " + object_path + " / " + attr_name);
	    }
	    hid_t a = H5Aopen(obj, attr_name.c_str(), H5P_DEFAULT);
	    if (a < 0) { H5Oclose(obj); throw std::runtime_error("Failed to open attribute: " + attr_name); }

	    hid_t s = H5Aget_space(a);
	    if (s < 0) { H5Aclose(a); H5Oclose(obj); throw std::runtime_error("Failed to get attr dataspace"); }

	    int nd = H5Sget_simple_extent_ndims(s);
	    hsize_t dims[1]{1};
	    if (nd == 1) H5Sget_simple_extent_dims(s, dims, nullptr);
	    else if (nd > 1) { H5Sclose(s); H5Aclose(a); H5Oclose(obj); throw std::runtime_error("Attr rank > 1 not supported"); }

	    std::vector<T> buf(static_cast<size_t>(dims[0]));
	    herr_t st = H5Aread(a, h5_native_type<T>(), buf.data());
	    H5Sclose(s);
	    H5Aclose(a);
	    H5Oclose(obj);
	    if (st < 0) throw std::runtime_error("Failed to read attribute: " + attr_name);
	    return buf;
	}

	// dataset(객체/이름) 우선, 없으면 동일 위치 attribute 로 대체
	template<typename T>
	inline T read_scalar_dset_or_attr(hid_t file_id,
	                                  const std::string& object_path,
	                                  const std::string& name) {
	    const std::string dset_path = object_path.empty()
	        ? name
	        : (object_path.back() == '/' ? object_path + name : object_path + "/" + name);

	    if (H5Lexists(file_id, dset_path.c_str(), H5P_DEFAULT))
	        return read_scalar_dataset<T>(file_id, dset_path);

	    if (attribute_exists(file_id, object_path, name)) {
	        auto v = read_attribute_by_name<T>(file_id, object_path, name);
	        if (v.size() == 1) return v[0];
	        throw std::runtime_error("Attribute is not scalar: " + object_path + ":" + name);
	    }

	    throw std::runtime_error("Missing dataset and attribute: " + dset_path + " (or attr on " + object_path + ")");
	}

	// 루트의 ID: dataset("ID") 우선, 없으면 attribute("ID")
	template<typename IdT>
	inline std::vector<IdT> read_ids_dataset_or_attr(hid_t file_id) {
	    if (H5Lexists(file_id, "ID", H5P_DEFAULT))
	        return read_1d_dataset<IdT>(file_id, "ID");
	    if (attribute_exists(file_id, ".", "ID"))
	        return read_attribute_by_name<IdT>(file_id, ".", "ID");
	    throw std::runtime_error("Neither dataset 'ID' nor root attribute 'ID' found.");
	}


	// -------------------------------
	// Aggregate return type for TreeFrog
	// -------------------------------
	//template<typename NumT, typename OffT, typename ResT, typename MerT,
    //     typename NpartT, typename NlinkT, typename IdT>
	//struct TreeData {
    //	std::vector<NumT>   num;
    //	std::vector<OffT>   off;
    //	std::vector<ResT>   res;
    //	std::vector<MerT>   merit;
    //	std::vector<NpartT> npart;
    //	std::vector<NlinkT> nlink; // attribute (scalar -> size 1)
    //	std::vector<IdT>    id;
	//};

	// -------------------------------
	// Aggregate return type for VR Catalog
	// -------------------------------
	// id, npart 는 항상 채움. p_id 는 read_p_id==true 일 때만 채움.
	//template<typename IdT, typename NpartT, typename PidT>
	//struct GalArrays {
	//    Makebr::TreeVec<IdT>    id;
	//    Makebr::TreeVec<NpartT> npart;
	//    Makebr::TreeVec<PidT>   p_id; // 선택적 로드
	//};

	// -------------------------------
	// Convert to Makebr::TreeVec
	// -------------------------------
	//template<typename NumT, typename OffT, typename ResT, typename MerT,
	//         typename NpartT, typename NlinkT, typename IdT>
	//struct TreeArrays_Makebr {
	//    Makebr::TreeVec<NumT>   num;
	//    Makebr::TreeVec<OffT>   off;
	//    Makebr::TreeVec<ResT>   res;
	//    Makebr::TreeVec<MerT>   merit;
	//    Makebr::TreeVec<NpartT> npart;
	//    Makebr::TreeVec<NlinkT> nlink;
	//    Makebr::TreeVec<IdT>    id;
	//};

	// -------------------------------
	// Makebr::TreeVec Reader
	// -------------------------------

	template<typename NumT, typename OffT, typename ResT, typename MerT,
	         typename NpartT, typename NlinkT, typename IdT>
	inline TreeArrays_Makebr<NumT,OffT,ResT,MerT,NpartT,NlinkT,IdT>
	read_treeset_makebr(const std::string& fname,
	                    const TreeSetTags& tags,
	                    std::size_t growth_step) {
	
	    hid_t file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	    if (file_id < 0) throw std::runtime_error("Failed to open file: " + fname);
	
	    try {
	        // 1) Load as std::vector
	        auto v_num   = read_1d_dataset<NumT>(file_id, tags.tag_num);
	        auto v_off   = read_1d_dataset<OffT>(file_id, tags.tag_off);
	        auto v_res   = read_1d_dataset<ResT>(file_id, tags.tag_result);
	        auto v_mer   = read_1d_dataset<MerT>(file_id, tags.tag_merit);
	        auto v_np    = read_1d_dataset<NpartT>(file_id, tags.tag_npart);
	        auto v_nl    = read_attribute_by_name<NlinkT>(file_id, ".", tags.tag_nlink);
	        auto v_id    = read_1d_dataset<IdT>(file_id, "ID");
	
	        // 2) Translate into Makebr::TreeVec
	        TreeArrays_Makebr<NumT,OffT,ResT,MerT,NpartT,NlinkT,IdT> out{
	            Makebr::TreeVec<NumT>{NumT{}, growth_step},
	            Makebr::TreeVec<OffT>{OffT{}, growth_step},
	            Makebr::TreeVec<ResT>{ResT{}, growth_step},
	            Makebr::TreeVec<MerT>{MerT{}, growth_step},
	            Makebr::TreeVec<NpartT>{NpartT{}, growth_step},
	            Makebr::TreeVec<NlinkT>{NlinkT{}, growth_step},
	            Makebr::TreeVec<IdT>{IdT{}, growth_step}
	        };
	
	        // Resize & Copy
	        out.num.resize(v_num.size());                  std::copy(v_num.begin(), v_num.end(), out.num.vec().begin());
	        out.off.resize(v_off.size());                  std::copy(v_off.begin(), v_off.end(), out.off.vec().begin());
	        out.res.resize(v_res.size());                  std::copy(v_res.begin(), v_res.end(), out.res.vec().begin());
	        out.merit.resize(v_mer.size());                std::copy(v_mer.begin(), v_mer.end(), out.merit.vec().begin());
	        out.npart.resize(v_np.size());                 std::copy(v_np.begin(), v_np.end(), out.npart.vec().begin());
	        out.nlink.resize(v_nl.size());                 std::copy(v_nl.begin(), v_nl.end(), out.nlink.vec().begin());
	        out.id.resize(v_id.size());                    std::copy(v_id.begin(), v_id.end(), out.id.vec().begin());
	
	        H5Fclose(file_id);
	        return out;
	    } catch (...) {
	        H5Fclose(file_id);
	        throw;
	    }
	}

	//-----
	// To Read VR Catalog
	//-----
	//==================== 메인 함수: r_gal ====================//
	// id0 < 0  → 모든 ID 읽기
	// id0 >= 0 → 해당 ID만 읽기
	// read_p_id 가 true 이면 "ID_xxxxxx/P_prop/P_ID" 도 읽음
	template<typename IdT=int32_t, typename NpartT=int32_t, typename PidT=int64_t>
	inline GalArrays<IdT,NpartT,PidT>
	r_gal(const std::string& fname,
	      long long id0,
	      bool read_p_id) {


		const std::string& id_prefix = "ID_";
		int id_zero_pad = 6;
	    const std::string& g_group = "G_Prop";
	    const std::string& p_group = "P_Prop";
	    const std::string& g_npart_name = "G_npart";
	    const std::string& p_id_name = "P_ID";
	    std::size_t growth_step = 1;

	    hid_t fid = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	    if (fid < 0) throw std::runtime_error("Failed to open file: " + fname);

	    try {
	        // 1) ID 목록 결정
	        std::vector<IdT> ids;
	        if (id0 < 0) {
	            ids = read_ids_dataset_or_attr<IdT>(fid);
	        } else {
	            ids = { static_cast<IdT>(id0) };
	        }

	        // 2) 결과 컨테이너 준비
	        GalArrays<IdT,NpartT,PidT> out;
	        out.id.set_growth_step(growth_step);
	        out.npart.set_growth_step(growth_step);
	        out.p_id.set_growth_step(growth_step);

	        out.id.resize(ids.size());
	        out.npart.resize(ids.size());
	        if (read_p_id) out.p_id.resize(ids.size());  // 읽을 때만 크기 설정

	        // 3) 각 ID 경로에서 읽기
	        for (std::size_t i = 0; i < ids.size(); ++i) {
	            const auto idv = static_cast<long long>(ids[i]);
	            const std::string base = id_prefix + zfill_longlong(idv, id_zero_pad); // "ID_000001"
	            const std::string gobj = base + "/" + g_group; // "ID_000001/G_prop"
	            const std::string pobj = base + "/" + p_group; // "ID_000001/P_prop"

	            // 필수: G_prop/G_npart
	            NpartT np = read_scalar_dset_or_attr<NpartT>(fid, gobj, g_npart_name);

	            // 선택: P_prop/P_ID
	            if (read_p_id) {
	                PidT pv = read_scalar_dset_or_attr<PidT>(fid, pobj, p_id_name);
	                out.p_id[i] = pv;
	            }

	            out.id[i]    = static_cast<IdT>(idv); // G_ID와 동일하므로 별도 읽지 않음
	            out.npart[i] = np;
	        }

	        H5Fclose(fid);
	        return out;
	    } catch (...) {
	        H5Fclose(fid);
	        throw;
	    }
	}

}