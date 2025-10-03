#pragma once
#include "global/allvar.h"
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

//-----
// HDF5 related
//-----

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

	// -------------------------------
	// Aggregate return type for VR Catalog
	// -------------------------------
	// id, npart 는 항상 채움. p_id 는 read_p_id==true 일 때만 채움.
	template<typename IdT, typename NpartT, typename PidT>
	struct GalArrays {
	    Makebr::TreeVec<IdT>    id;
	    Makebr::TreeVec<NpartT> npart;
	    Makebr::TreeVec<PidT>   p_id; // 선택적 로드
	};

	// -------------------------------
	// Convert to Makebr::TreeVec
	// -------------------------------
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

	// -------------------------------
	// Makebr::TreeVec Reader
	// -------------------------------

	template<typename NumT, typename OffT, typename ResT, typename MerT,
	         typename NpartT, typename NlinkT, typename IdT>
	inline TreeArrays_Makebr<NumT,OffT,ResT,MerT,NpartT,NlinkT,IdT>
	read_treeset_makebr(const std::string& fname,
	                    const TreeSetTags& tags,
	                    std::size_t growth_step = 1) {
	
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
	      long long id0 = -1,
	      bool read_p_id = false) {


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
