#pragma once
#include "global/allvar.h"
//#include "global/io_hdf5.h"
#include "utilities/utilities.h"
#include <algorithm>
#include <array>
#include <cstdint>
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



namespace IO_dtype{
	using IO_I32 		= std::int32_t;		
	using IO_I64 		= std::int64_t;
	using IO_float		= float;
	using IO_double		= double;

	using IO_GID 		= std::int32_t; // Catalog galaxy ID type
	using IO_Snap 		= std::int32_t; // Catalog snapshot number type
	using IO_PID 		= std::int64_t; // Catalog particle ID type
	using IO_BID 		= std::int32_t; // Branch ID type

	// For Snapshot info data type
	struct snapSt{
		IO_I32 snum = -1;			// snapshot number
		IO_double aexp;			// scale factor
		IO_double unit_l; 		// sim unit to kpc
		IO_double unit_t;		// sim unit to Gyr
		IO_double h0;
		IO_double om;
		IO_double ol;
		IO_double age;			// age at this zredshift
	};

	using snapinfo = std::vector<snapSt>;

	// Galaxy Catalog Type
	struct GalSt{
		IO_GID id;
		IO_Snap snap;
		IO_I32 npart;
		std::vector<IO_PID> pid;
		IO_double xc, yc, zc;
		IO_double vxc, vyc, vzc;
	};

	using GalArray = std::vector<GalSt>;
}

namespace IO_VR{
	inline std::string zfill_longlong(long long x, int width) {
	    std::ostringstream oss;
	    oss << std::setfill('0') << std::setw(width) << x;
	    return oss.str();
	}

	//-----
	// Data Types
	//-----
	using VRT_GID 		= IO_dtype::IO_GID;//std::int32_t;		// for galaxy ID
	using VRT_Snap 		= IO_dtype::IO_Snap;//std::int32_t;
	using VRT_PID 	 	= IO_dtype::IO_PID;//std::int64_t;		// for particle ID
	using VRT_merit		= IO_dtype::IO_double;//double;
	using VRT_BID 		= IO_dtype::IO_BID;//std::int32_t;		// Branch Index
	using VRT_I32		= IO_dtype::IO_I32;//std::int32_t;
	using VRT_double 	= IO_dtype::IO_double;

	//----- For galaxy array
	using GalSt = IO_dtype::GalSt;
	using GalArray = IO_dtype::GalArray;
	//struct GalSt{
	//	VRT_GID id;
	//	VRT_Snap snap;
	//	VRT_I32 npart;
	//	std::vector<VRT_PID> pid;
	//};
	//using GalArray = std::vector<GalSt>;

	
	//-----
    // HDF5 related
    //-----
    template<typename T> inline hid_t h5_native_type();
	//template<> inline hid_t h5_native_type<int8_t>()   { return H5T_NATIVE_INT8;  }
	//template<> inline hid_t h5_native_type<uint8_t>()  { return H5T_NATIVE_UINT8; }
	//template<> inline hid_t h5_native_type<int16_t>()  { return H5T_NATIVE_INT16; }
	//template<> inline hid_t h5_native_type<uint16_t>() { return H5T_NATIVE_UINT16;}
	template<> inline hid_t h5_native_type<std::int32_t>()  { return H5T_NATIVE_INT32; }
	//template<> inline hid_t h5_native_type<uint32_t>() { return H5T_NATIVE_UINT32;}
	template<> inline hid_t h5_native_type<std::int64_t>()  { return H5T_NATIVE_INT64; }
	//template<> inline hid_t h5_native_type<uint64_t>() { return H5T_NATIVE_UINT64;}
	template<> inline hid_t h5_native_type<float>()    { return H5T_NATIVE_FLOAT; }
	template<> inline hid_t h5_native_type<double>()   { return H5T_NATIVE_DOUBLE; }
	//template<typename T>
	//inline hid_t h5_native_type() {
	//    using std::is_same_v;
	//    if constexpr (is_same_v<T, int8_t>)        return H5T_NATIVE_INT8;
	//    else if constexpr (is_same_v<T, uint8_t>)  return H5T_NATIVE_UINT8;
	//    else if constexpr (is_same_v<T, int16_t>)  return H5T_NATIVE_INT16;
	//    else if constexpr (is_same_v<T, uint16_t>) return H5T_NATIVE_UINT16;
	//    else if constexpr (is_same_v<T, int32_t>)  return H5T_NATIVE_INT32;
	//    else if constexpr (is_same_v<T, uint32_t>) return H5T_NATIVE_UINT32;
	//    else if constexpr (is_same_v<T, int64_t>)  return H5T_NATIVE_INT64;   // LLP64에선 long long
	//    else if constexpr (is_same_v<T, uint64_t>) return H5T_NATIVE_UINT64;
	//    else if constexpr (is_same_v<T, long>)     return H5T_NATIVE_LONG;    // LP64에선 int64_t와 동일
	//    else if constexpr (is_same_v<T, unsigned long>) return H5T_NATIVE_ULONG;
	//    else if constexpr (is_same_v<T, float>)    return H5T_NATIVE_FLOAT;
	//    else if constexpr (is_same_v<T, double>)   return H5T_NATIVE_DOUBLE;
	//    else {
	//        static_assert(sizeof(T) == 0, "Unsupported type for h5_native_type<T>()");
	//    }
	//}

    template<typename T>
    inline std::vector<T> VR_HDF5_rdbyname(const hid_t file_id, const std::string& dname){
    	hid_t d = H5Dopen2(file_id, dname.c_str(), H5P_DEFAULT);
    	if (d < 0) throw std::runtime_error("Failed to open dataset: " + dname);

    	hid_t s = H5Dget_space(d);
    	if (s < 0) { H5Dclose(d); throw std::runtime_error("Failed to get dataspace: " + dname); }

    	int nd = H5Sget_simple_extent_ndims(s);
    	std::vector<T> out;
    	if (nd == 0) {
        	out.resize(1);
    	} else if (nd == 1) {
        	hsize_t dims[1]{};
        	H5Sget_simple_extent_dims(s, dims, nullptr);
        	out.resize(static_cast<size_t>(dims[0]));
    	} else {
        	H5Sclose(s); H5Dclose(d);
        	throw std::runtime_error("Dataset rank > 1 not supported: " + dname);
    	}

    	herr_t st = H5Dread(d, h5_native_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
    	H5Sclose(s);
    	H5Dclose(d);
    	if (st < 0) throw std::runtime_error("Failed to read dataset: " + dname);
    	return out;

    }

    template<typename T>
    inline std::vector<T> VR_HDF5_rdbyattr(const hid_t file_id, const std::string& aname){
    	// open attr 
    	std::string object_path = "/";
    	hid_t a = H5Aopen_by_name(file_id, object_path.c_str(), aname.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    	if (a < 0) throw std::runtime_error("Failed to open attribute: " + object_path + " / " + aname);

    	// dataspace
    	hid_t s = H5Aget_space(a);
    	if (s < 0) { H5Aclose(a); throw std::runtime_error("Failed to get attr dataspace"); }

    	int nd = H5Sget_simple_extent_ndims(s);
    	std::vector<T> out;
    	if (nd == 0) {
        	out.resize(1);
    	} else if (nd == 1) {
 	       hsize_t dims[1]{1};
    	    H5Sget_simple_extent_dims(s, dims, nullptr);
        	out.resize(static_cast<size_t>(dims[0]));
    	} else {
        	H5Sclose(s);
        	H5Aclose(a);
        	throw std::runtime_error("Attribute rank > 1 not supported: " + object_path + " / " + aname);
    	}

    	// read
    	herr_t st = H5Aread(a, h5_native_type<T>(), out.data());
    	H5Sclose(s);
    	H5Aclose(a);
    	if (st < 0) throw std::runtime_error("Failed to read attribute: " + object_path + " / " + aname);

    	return out;

    }
   
	//-----
	// read Galaxy Catalog
	//	
	//-----
	// id0 < 0  → All galaxies
	// id0 >= 0 → A Galaxy corresponding to this ID
	// horg 	→ 'h' or 'g'
	// read_p_id = true → return particle ID

	inline GalArray r_gal(const vctree_set::Settings& vh, const VRT_Snap snap_curr, const VRT_GID id0, const bool readpart=false){
		const std::string& id_prefix = "ID_";
		int id_zero_pad = 6;
	    const std::string& g_group = "G_Prop";
	    const std::string& p_group = "P_Prop";
	    const std::string& g_npart_name = "G_npart";
	    const std::string& p_id_name = "P_ID";
	    //std::size_t growth_step = 1;

	    std::string fname 	= vh.vr_dir_catalog + "/snap_" + i4(snap_curr) + ".hdf5";
		hid_t fid = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	    if (fid < 0) throw std::runtime_error("Failed to open file: " + fname);


	    try {
	        
	        std::vector<VRT_GID> ids;
	        GalArray gal;
	        if (id0 < 0) {
	            ids = VR_HDF5_rdbyname<VRT_GID>(fid, "ID");
	        } else {
	            ids = { static_cast<VRT_GID>(id0) };
	        }

	        gal.resize(ids.size());

	        for(VRT_GID i=0; i < (VRT_GID) ids.size(); i++){
	        	const auto idv = static_cast<VRT_GID>(ids[i]);
	        	const std::string base = id_prefix + zfill_longlong(idv, id_zero_pad); // "ID_000001"
	        	const std::string gobj = base + "/" + g_group; // "ID_000001/G_prop"
	            const std::string pobj = base + "/" + p_group; // "ID_000001/P_prop"

	            gal[i].id 	= static_cast<VRT_GID>(ids[i]);

	            std::vector<VRT_GID> np = VR_HDF5_rdbyname<VRT_GID>(fid, gobj + "/" + g_npart_name);
	            gal[i].npart = static_cast<VRT_I32>(np[0]);

	            // Pos & Vel
	            std::vector<VRT_double> dum;

				dum = VR_HDF5_rdbyname<VRT_double>(fid, gobj + "/G_Xc");
	            gal[i].xc = static_cast<VRT_double>(dum[0]);

	            dum = VR_HDF5_rdbyname<VRT_double>(fid, gobj + "/G_Yc");
	            gal[i].yc = static_cast<VRT_double>(dum[0]);

	            dum = VR_HDF5_rdbyname<VRT_double>(fid, gobj + "/G_Zc");
	            gal[i].zc = static_cast<VRT_double>(dum[0]);

	            dum = VR_HDF5_rdbyname<VRT_double>(fid, gobj + "/G_VXc");
	            gal[i].vxc = static_cast<VRT_double>(dum[0]);

	            dum = VR_HDF5_rdbyname<VRT_double>(fid, gobj + "/G_VYc");
	            gal[i].vyc = static_cast<VRT_double>(dum[0]);

	            dum = VR_HDF5_rdbyname<VRT_double>(fid, gobj + "/G_VZc");
	            gal[i].vzc = static_cast<VRT_double>(dum[0]);


	            if(readpart){
	            	std::vector<VRT_PID> pid = VR_HDF5_rdbyname<VRT_PID>(fid, pobj + "/" + p_id_name);
	            	gal[i].pid 	= pid;
	            }

	            gal[i].snap 	= snap_curr;
	            

	            
	        }

	        H5Fclose(fid);

	        
	        return gal;
	    } catch (...) {
	        H5Fclose(fid);
	        throw;
	    }

	}
}

//-----
// IO for RAMSES
//-----
namespace IO_RAMSES{

	static inline std::string trim(std::string s) {
	    auto not_space = [](unsigned char ch){ return !std::isspace(ch); };
	    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
	    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
	    return s;
	}
	
	template<typename T>
  	bool parse_num(const std::string& s, T& out) {
    	std::istringstream iss(trim(s));
    	iss >> out;
    return (bool)iss && iss.eof();
  	}

  	

	


	using IO_I32 		= IO_dtype::IO_I32;
	using IO_I64 		= IO_dtype::IO_I64;
	using IO_float		= IO_dtype::IO_float;
	using IO_double		= IO_dtype::IO_double;

	IO_dtype::snapinfo get_snapinfo(const vctree_set::Settings& vh);
	bool is_snap(const vctree_set::Settings& vh, const IO_I32 snap_curr);

	IO_dtype::snapSt read_info(const vctree_set::Settings& vh, const IO_I32 snap_curr);

}
//-----
// IO control space
//-----
namespace IO {
	//-----
	// Data Types
	//-----
	using IO_I32 		= IO_dtype::IO_I32;
	using IO_I64 		= IO_dtype::IO_I64;
	using IO_float		= IO_dtype::IO_float;
	using IO_double		= IO_dtype::IO_double;
	using snapinfo 		= IO_dtype::snapinfo;
	//using VRT_Snap 		= std::int32_t;
	//using VRT_PID 	 	= std::int64_t;		// for particle ID
	//using VRT_merit		= double;
	//using VRT_BID 		= std::int32_t;		// Branch Index
	//using VRT_I32		= std::int32_t;

	


	//----- Get SnapInfo
	inline IO_dtype::snapinfo get_snapinfo(const vctree_set::Settings& vh){
		if(vh.simtype == "RAMSES"){
			return IO_RAMSES::get_snapinfo(vh);
		} else{
			LOG()<<"Wrong simulation type: simtype";
			u_stop();
		}
	}

	//----- Read Catalog
	inline IO_dtype::GalArray r_gal(const vctree_set::Settings& vh, const IO_dtype::IO_Snap snap_curr, const IO_dtype::IO_GID id0, const bool readpart=false){
		if(vh.iotype == "VR"){
			return IO_VR::r_gal(vh, snap_curr, id0, readpart);
		} else{
			LOG()<<"Should be implemented for different IO Type";
			u_stop();
		}
	}

	
}




