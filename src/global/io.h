#pragma once
#include "global/allvar.h"
//#include "global/io_hdf5.h"
#include "utilities/utilities.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <filesystem>   // C++17
#include <fstream>
#include <hdf5.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <type_traits>


namespace fs = std::filesystem;



namespace IO_dtype{
	using IO_I32 		= std::int32_t;		
	using IO_I64 		= std::int64_t;
	using IO_float		= float;
	using IO_double		= double;

	using IO_BID 		= std::int32_t; // Branch ID type

//-----
// Compile-time configurable types (set via CMake options)
//
// Defaults (when no CMake option is given):
//   IO_GID     : int32_t
//   IO_Snap    : int32_t
//   IO_PID    	: int64_t
//   IO_merit   : double
//-----

#if defined(CTREE_ID_INT64)
	using IO_GID 		= std::int64_t; // Catalog galaxy ID type
#else
	using IO_GID 		= std::int32_t; // Catalog galaxy ID type
#endif

#if defined(CTREE_SNAP_INT64)
	using IO_Snap 		= std::int64_t; // Catalog snapshot number type
#else
	using IO_Snap 		= std::int32_t; // Catalog snapshot number type
#endif

#if defined(CTREE_PARTID_INT32)
	using IO_PID 		= std::int32_t; // Catalog particle ID type
#else
	using IO_PID 		= std::int64_t; // Catalog particle ID type
#endif

#if defined(CTREE_MERIT_FLOAT32)
	using IO_merit 		= float;
#else
	using IO_merit 		= double;
#endif


// For Snapshot info data type
//	struct snapSt{
//		IO_I32 snum = -1;			// snapshot number
//		IO_double aexp;			// scale factor
//		IO_double unit_l; 		// sim unit to kpc
//		IO_double unit_t;		// sim unit to Gyr
//		IO_double h0;
//		IO_double om;
//		IO_double ol;
//		IO_double age;			// age at this zredshift
//	};
	using snapSt = vctree_set::snapSt;
	using snapinfo = vctree_set::snapinfo;
//	using snapinfo = std::vector<snapSt>;

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

//-----
// FOR VELOCIraptor catalog
//-----
namespace IO_VR{
	//-----
	// Data Types
	//-----
	using VRT_GID 		= IO_dtype::IO_GID;//std::int32_t;		// for galaxy ID
	using VRT_Snap 		= IO_dtype::IO_Snap;//std::int32_t;
	using VRT_PID 	 	= IO_dtype::IO_PID;//std::int64_t;		// for particle ID
	using VRT_merit		= IO_dtype::IO_merit;//double;
	using VRT_BID 		= IO_dtype::IO_BID;//std::int32_t;		// Branch Index
	using VRT_I32		= IO_dtype::IO_I32;//std::int32_t;
	using VRT_I64 		= IO_dtype::IO_I64;
	using VRT_double 	= IO_dtype::IO_double;

	//----- For galaxy array
	using GalSt = IO_dtype::GalSt;
	using GalArray = IO_dtype::GalArray;

	
	inline void in_gpt(std::vector<vctree_set::vr_pointerSt>& gpt, VRT_GID id, VRT_I32 fnum, VRT_I64 offset){
		if((VRT_GID) gpt.size() <= id){
			gpt.resize(id+10000);
		}
		gpt[id].offset	= offset;
		gpt[id].fnum = fnum;
	}

	struct FileSt{
		int N;
		std::string fname;
	};

	inline std::vector<FileSt> r_gal_vrflist(const std::string dir, const std::string prefix){
		std::vector<FileSt> flist;

		if(!fs::exists(dir)){
			throw std::runtime_error("No directory for catalog: " + dir);
			return flist;
		}

		for (const auto& ent : fs::directory_iterator(dir)){

			if( !ent.is_regular_file()) continue;

			std::string fname = ent.path().filename().string();

			if(fname.find(prefix) == std::string::npos) continue;

			std::string mpinum = fname.substr(fname.find_last_of('.')+1);

			if(mpinum == "swp") continue;

			//FileSt dumf;
			//dumf.fname = fname;
			//dumf.N = std::stoi(mpinum);

			flist.push_back(FileSt{std::stoi(mpinum), fname});

		}

		std::sort(flist.begin(), flist.end(), [](const FileSt& a, const FileSt& b){
			return a.N < b.N;
		});

		return flist;
	}

	// assume that 1st and 6th columns are ID and npart, respectively
	static inline bool is_space(char c) { return c == ' ' || c == '\t'; }

	inline void read_prop_fromvr(const std::string& line, IO_dtype::GalSt& gal0){
		const char* s = line.c_str();
		size_t n = line.size(), i=0;

		auto read_token = [&](std::string& tok) -> bool {
			tok.clear();
			while (i < n && is_space(s[i])) ++i;
			if (i >= n) return false;
			size_t start = i;
			while (i < n && !is_space(s[i])) ++i;
			tok.assign(s + start, i - start);
			if (!tok.empty() && tok.back() == '\r') tok.pop_back();
			return true;
		};

		std::string tok;
		int col=0;
		while (read_token(tok)) {
			

			
			if (col == 0){ gal0.id = static_cast<VRT_GID>(std::stol(tok));}
			else if (col == 5){ gal0.npart = static_cast<VRT_I32>(std::stol(tok));}

			col++;
		}

		
	}
	//
	inline GalArray r_gal(vctree_set::Settings& vh, const VRT_Snap snap_curr, const VRT_GID id0, const bool readpart=false){
		const std::string& dir = vh.vr_dir_catalog + "/" + vh.vr_dir_catalog_prefix + iN(snap_curr, vh.vr_dir_catalog_snapdigit) + vh.vr_dir_catalog_suffix;

		GalArray gal;
		VRT_I32 pick_id=0;

		// get filelist
		std::vector<FileSt> flist = r_gal_vrflist(dir, "properties");
		std::vector<FileSt> flist_g = r_gal_vrflist(dir, "catalog_groups");
		std::vector<FileSt> flist_pdum = r_gal_vrflist(dir, "catalog_particles");

		VRT_I32 nfile 	= flist.size();

		std::vector<VRT_I32> nobj;
		nobj.resize(nfile);

		// divide particle files into bound and unbound ones
		std::vector<FileSt> flist_pb, flist_pu;
		flist_pb.resize(flist.size());
		flist_pu.resize(flist.size());

		VRT_I32 ib=0, iu=0;
		for(int i=0; i<nfile*2; i++){
			if(flist_pdum[i].fname.find("unbound") == std::string::npos){
				flist_pb[ib]	= flist_pdum[i];
				ib ++;
			}else{
				flist_pu[iu]	= flist_pdum[i];
				iu ++;
			}
		}

		// allocate
		VRT_I32 ngal = 0;
		for(int i=0; i<nfile; i++){
			std::ifstream fin(dir + "/" + flist[i].fname);
			std::string line;
			while( std::getline(fin, line) ) ngal++;
			ngal -= 3;
		}
		//gal.resize(ngal);


		// get offset vector for this snapshot		
		std::vector<vctree_set::vr_pointerSt>& gpt = vh.vr_gpointer[snap_curr];

		// read all galaxies when id0<0 or never read
		if(id0<0 || gpt.size()==0){
			if(gpt.size()==0){ // if having no prior offset, build it
				VRT_I32 ini_max_id = 100000;
				gpt.resize(ini_max_id);
			}

			//// allocate galaxy array

			gal.resize(ngal);

			VRT_I32 gindex = 0;

			// read file line by line
			for(int i=0; i<nfile; i++){
				std::ifstream f(dir + "/" + flist[i].fname, std::ios::binary);

				// skip first three lines
				std::string line;
				std::getline(f, line);
				std::getline(f, line);
				std::getline(f, line);

				while (true) {
					std::streamoff linestart = f.tellg();
					if( !std::getline(f, line)) break;
					if( line.empty() ) continue;

					read_prop_fromvr(line, gal[gindex]);

					in_gpt(gpt, gal[gindex].id, i, linestart);
					gal[gindex].snap = snap_curr;
					if(id0>0){
						pick_id = gindex;
					}

					gindex++;
				}

				f.close();
				if(i==0) nobj[i] = gindex;
				if(i>0) nobj[i] = gindex - nobj[i-1];
			}

			// read particle files (including bound + unbound)

			gindex = 0;
			for(int i=0; i<nfile; i++){
				std::ifstream fg(dir + "/" + flist_g[i].fname, std::ios::binary);
				std::string line;

				std::getline(fg, line);
				std::getline(fg, line);

				for(int j=0; j<nobj[i]; j++){
					std::getline(fg, line);
					gpt[gal[gindex + j].id].nt 	= std::stol(line);

				}

				std::getline(fg, line);
				VRT_I32 idum = std::stol(line);
				for(int j=0; j<nobj[i]-1; j++){
					std::getline(fg, line);
					gpt[gal[gindex + j].id].nb 	= std::stol(line) - idum;
					idum 	= std::stol(line);
				}

				std::getline(fg, line);
				idum = std::stol(line);
				for(int j=0; j<nobj[i]-1; j++){
					std::getline(fg, line);
					gpt[gal[gindex + j].id].nu 	= std::stol(line) - idum;
					idum 	= std::stol(line);
				}
				// don't forget the last of nb and nu missing here

				fg.close();


				//
				std::ifstream fb(dir + "/" + flist_pb[i].fname, std::ios::binary);
				std::ifstream fu(dir + "/" + flist_pu[i].fname, std::ios::binary);

				std::getline(fb, line); std::getline(fb, line);
				std::getline(fu, line); std::getline(fu, line);


				for(int j=0; j<nobj[i]; j++){
					std::streamoff ls_b = fb.tellg();
					std::streamoff ls_u = fu.tellg();

					gpt[gal[gindex + j].id].offset_pb 	= ls_b;
					gpt[gal[gindex + j].id].offset_pu 	= ls_u;

					
					std::vector<VRT_PID> pid;
					VRT_PID pindex = 0;
					if(readpart){
						pid.resize(gal[gindex+j].npart);
					}

					// update nb and nu here
					if(j<nobj[i]-1){
						for(int k=0; k<gpt[gal[gindex + j].id].nb; k++){
							std::getline(fb, line);
							if(readpart){pid[pindex] = std::stol(line); pindex++;}
						}

						for(int k=0; k<gpt[gal[gindex + j].id].nu; k++){
							std::getline(fu, line);
							if(readpart){pid[pindex] = std::stol(line); pindex++;}
						}
					}else{
						gpt[gal[gindex + j].id].nb = 0;
						gpt[gal[gindex + j].id].nu = 0;

						while(std::getline(fb, line)){
							gpt[gal[gindex + j].id].nb ++;
							if(readpart){pid[pindex] = std::stol(line); pindex++;}
						}
						while(std::getline(fu, line)){
							gpt[gal[gindex + j].id].nu ++;
							if(readpart){pid[pindex] = std::stol(line); pindex++;}
						}
					}

					if(readpart){gal[gindex+j].pid 	= pid; pid.clear();}

				}
				
				fb.close();
				fu.close();



				gindex += nobj[i];
			}
		}
		// read one galaxy with given ID
		else{
			ngal 	= 1;
			gal.resize(ngal);
			pick_id = 0;

			VRT_I32 mpinum = gpt[id0].fnum;
			VRT_I64 offset = gpt[id0].offset;

			std::ifstream f(dir + "/" + flist[mpinum].fname, std::ios::binary);
			f.seekg(offset);
			std::string line;
			std::getline(f, line);

			read_prop_fromvr(line, gal[0]);
			gal[0].snap = snap_curr;

			if(readpart){
				std::ifstream fb(dir + "/" + flist_pb[mpinum].fname, std::ios::binary);
				std::ifstream fu(dir + "/" + flist_pu[mpinum].fname, std::ios::binary);

				fb.seekg(gpt[id0].offset_pb);
				fu.seekg(gpt[id0].offset_pu);

				std::string line;

				std::vector<VRT_PID> pid;
				pid.resize(gal[0].npart);

				VRT_PID pindex = 0;
				for(VRT_PID k=0; k<gpt[id0].nb; k++){

					std::getline(fb, line);
					pid[pindex] = std::stol(line);

					pindex ++;
				}

				for(VRT_PID k=0; k<gpt[id0].nu; k++){
					std::getline(fu, line);
					pid[pindex] = std::stol(line);

					pindex ++;
				}


				gal[0].pid 	= pid;
				fb.close();
				fu.close();
			}

			f.close();
		}

		if(id0>0){
			return {gal[pick_id]};
		}
		else{
			return gal;
		}


		return gal;

//		const std::string& id_prefix = "ID_";
//		int id_zero_pad = 6;
//	    const std::string& g_group = "G_Prop";
//	    const std::string& p_group = "P_Prop";
//	    const std::string& g_npart_name = "G_npart";
//	    const std::string& p_id_name = "P_ID";
//	    //std::size_t growth_step = 1;
//
//	    std::string fname 	= vh.veluga_dir_catalog + "/snap_" + i4(snap_curr) + ".hdf5";
//		hid_t fid = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
//
//	    if (fid < 0) throw std::runtime_error("Failed to open file: " + fname);
//
//
//	    try {
//	        
//	        std::vector<VRT_GID> ids;
//	        GalArray gal;
//	        if (id0 < 0) {
//	            ids = VR_HDF5_rdbyname<VRT_GID>(fid, "ID");
//	        } else {
//	            ids = { static_cast<VRT_GID>(id0) };
//	        }
//
//	        gal.resize(ids.size());
//
//	        for(VRT_GID i=0; i < (VRT_GID) ids.size(); i++){
//	        	const auto idv = static_cast<VRT_GID>(ids[i]);
//	        	const std::string base = id_prefix + zfill_longlong(idv, id_zero_pad); // "ID_000001"
//	        	const std::string gobj = base + "/" + g_group; // "ID_000001/G_prop"
//	            const std::string pobj = base + "/" + p_group; // "ID_000001/P_prop"
//
//	            gal[i].id 	= static_cast<VRT_GID>(ids[i]);
//
//	            std::vector<VRT_GID> np = VR_HDF5_rdbyname<VRT_GID>(fid, gobj + "/" + g_npart_name);
//	            gal[i].npart = static_cast<VRT_I32>(np[0]);
//
//	            // Pos & Vel
//	            std::vector<VRT_double> dum;
//
//				dum = VR_HDF5_rdbyname<VRT_double>(fid, gobj + "/G_Xc");
//	            gal[i].xc = static_cast<VRT_double>(dum[0]);
//
//	            dum = VR_HDF5_rdbyname<VRT_double>(fid, gobj + "/G_Yc");
//	            gal[i].yc = static_cast<VRT_double>(dum[0]);
//
//	            dum = VR_HDF5_rdbyname<VRT_double>(fid, gobj + "/G_Zc");
//	            gal[i].zc = static_cast<VRT_double>(dum[0]);
//
//	            dum = VR_HDF5_rdbyname<VRT_double>(fid, gobj + "/G_VXc");
//	            gal[i].vxc = static_cast<VRT_double>(dum[0]);
//
//	            dum = VR_HDF5_rdbyname<VRT_double>(fid, gobj + "/G_VYc");
//	            gal[i].vyc = static_cast<VRT_double>(dum[0]);
//
//	            dum = VR_HDF5_rdbyname<VRT_double>(fid, gobj + "/G_VZc");
//	            gal[i].vzc = static_cast<VRT_double>(dum[0]);
//
//
//	            if(readpart){
//	            	std::vector<VRT_PID> pid = VR_HDF5_rdbyname<VRT_PID>(fid, pobj + "/" + p_id_name);
//	            	gal[i].pid 	= pid;
//	            }
//
//	            gal[i].snap 	= snap_curr;
//	            
//
//	            
//	        }
//
//	        H5Fclose(fid);
//
//	        
//	        return gal;
//	    } catch (...) {
//	        H5Fclose(fid);
//	        throw;
//	    }
//
	}
}


//-----
// FOR VELUGA catalog
//-----
namespace IO_VELUGA{
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
	using VRT_merit		= IO_dtype::IO_merit;//double;
	using VRT_BID 		= IO_dtype::IO_BID;//std::int32_t;		// Branch Index
	using VRT_I32		= IO_dtype::IO_I32;//std::int32_t;
	using VRT_I64 		= IO_dtype::IO_I64;
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
	//    else if constexpr (is_same_v<T, int64_t>)  return H5T_NATIVE_INT64;   
	//    else if constexpr (is_same_v<T, uint64_t>) return H5T_NATIVE_UINT64;
	//    else if constexpr (is_same_v<T, long>)     return H5T_NATIVE_LONG;    
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
	// id0 < 0  All galaxies
	// id0 >= 0 A Galaxy corresponding to this ID
	// horg 	'h' or 'g'
	// read_p_id = true return particle ID

	inline GalArray r_gal(const vctree_set::Settings& vh, const VRT_Snap snap_curr, const VRT_GID id0, const bool readpart=false){
		const std::string& id_prefix = "ID_";
		int id_zero_pad = 6;
		const std::string& g_group = "G_Prop";
		const std::string& p_group = "P_Prop";
		const std::string& g_npart_name = "G_npart";
		const std::string& p_id_name = "P_ID";
		//std::size_t growth_step = 1;

		std::string fname 	= vh.veluga_dir_catalog + "/snap_" + i4(snap_curr) + ".hdf5";
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
// For HaloMaker catalog
//-----
namespace IO_HM{
	inline std::string zfill_longlong(long long x, int width) {
		std::ostringstream oss;
		oss << std::setfill('0') << std::setw(width) << x;
		return oss.str();
	}

	// Endian swap
	inline uint16_t bswap16(uint16_t x) { return __builtin_bswap16(x); }
	inline uint32_t bswap32(uint32_t x) { return __builtin_bswap32(x); }
	inline uint64_t bswap64(uint64_t x) { return __builtin_bswap64(x); }
	template <typename T>
	T bswap(T v) {
		if constexpr (sizeof(T) == 2) {
			uint16_t x; std::memcpy(&x, &v, 2);
			x = bswap16(x);
			std::memcpy(&v, &x, 2);
		} else if constexpr (sizeof(T) == 4) {
			uint32_t x; std::memcpy(&x, &v, 4);
			x = bswap32(x);
			std::memcpy(&v, &x, 4);
		} else if constexpr (sizeof(T) == 8) {
			uint64_t x; std::memcpy(&x, &v, 8);
			x = bswap64(x);
			std::memcpy(&v, &x, 8);
		}
		return v;
	}

	inline std::vector<char> read_f77_record(std::ifstream& ifs) {
		int32_t len1 = 0, len2 = 0;

		if (!ifs.read(reinterpret_cast<char*>(&len1), 4))
			throw std::runtime_error("Failed to read record header");

		std::vector<char> buf(len1);
		if (!ifs.read(buf.data(), len1))
			throw std::runtime_error("Failed to read record payload");

		if (!ifs.read(reinterpret_cast<char*>(&len2), 4))
			throw std::runtime_error("Failed to read record footer");

		if (len1 != len2)
			throw std::runtime_error("Record length mismatch");

		return buf;
	}

	// ---- SKIP record  ----
	inline void skip_f77_record(std::ifstream& ifs) {
		int32_t len1 = 0, len2 = 0;

		// header
		if (!ifs.read(reinterpret_cast<char*>(&len1), 4))
			throw std::runtime_error("Failed to read record header");

		// skip payload
		if (!ifs.seekg(len1, std::ios::cur))
			throw std::runtime_error("Failed to skip record payload");

		// footer
		if (!ifs.read(reinterpret_cast<char*>(&len2), 4))
			throw std::runtime_error("Failed to read record footer");

		if (len1 != len2)
			throw std::runtime_error("Record length mismatch on skip");
	}

	//-----
	// Data Types
	//-----
	using HM_GID 		= IO_dtype::IO_GID;//std::int32_t;		// for galaxy ID
	using HM_Snap 		= IO_dtype::IO_Snap;//std::int32_t;
	using HM_PID 	 	= IO_dtype::IO_PID;//std::int64_t;		// for particle ID
	using HM_merit		= IO_dtype::IO_merit;//double;
	using HM_BID 		= IO_dtype::IO_BID;//std::int32_t;		// Branch Index
	using HM_I32		= IO_dtype::IO_I32;//std::int32_t;
	using HM_I64 		= IO_dtype::IO_I64;
	using HM_double 	= IO_dtype::IO_double;

	//----- For galaxy array
	using GalSt = IO_dtype::GalSt;
	using GalArray = IO_dtype::GalArray;

	inline void in_gpt(std::vector<int64_t>& gpt, HM_I32 id, HM_I64 pointer){
		if((HM_I32) gpt.size() <= id){
			gpt.resize(id*2);
		}
		gpt[id]	= pointer;
	}
	//-----
	// read Galaxy Catalog
	//	
	//-----
	// id0 < 0  All galaxies
	// id0 >= 0 A Galaxy corresponding to this ID
	// horg 	'h' or 'g'
	// read_p_id = true return particle ID

	inline GalArray r_gal(vctree_set::Settings& vh, const HM_Snap snap_curr, const HM_GID id0, const bool readpart=false){

		std::string fname 	= vh.hm_dir_catalog + "/tree_bricks" + i5(snap_curr);

		std::ifstream ifs(fname, std::ios::binary);

		if(!ifs){
			throw std::runtime_error("Failed to open file: " + fname);
		}

		//check the presence of pointer
		std::vector<int64_t>& gpt = vh.hm_gpointer[snap_curr];

		if(gpt.size()==0){ // have no information about file pointer

			HM_I32 ini_max_id 	= 100000;
			gpt.resize(ini_max_id);
			GalArray gal;

			HM_I32 nbodies, nmain, nsub, nall, nbin;
			HM_I64 curr_pt = 0;
			nbin = 0;

			// READ nbodies
			auto rec 	= read_f77_record(ifs);
			std::memcpy(&nbodies, rec.data(), 4);
			curr_pt += (4+4+4);

			// SKIP FOR DOUBLES (massp, aexp, omega_t, age_univ)
			for(int k=0; k<4; k++) skip_f77_record(ifs);
			curr_pt += (4+8+4)*4;

			// Read n_main & n_sub
			rec 	= read_f77_record(ifs);
			std::memcpy(&nmain, rec.data(), 4);
			std::memcpy(&nsub, rec.data() + 4, 4);
			curr_pt += (4+4*2+4);


			nall 	= nmain + nsub;

			if(id0>0){
				gal.resize(1);
			}else{
				gal.resize(nall);
			}

			// Read with loop
			HM_I32 numpart, gid;
			std::vector<HM_I32> pid;
			std::vector<char> rec2;

			HM_I64 curr_pt_old;
			for(HM_I32 i=0; i<nall; i++){

				curr_pt_old = curr_pt;

				// numpart
				rec 	= read_f77_record(ifs);
				std::memcpy(&numpart, rec.data(), 4);
				curr_pt 	+= (4+4+4);

				// read Particle ID
				if(readpart){
					rec2 	= read_f77_record(ifs);
					//pid.resize(numpart);
				}else{
					skip_f77_record(ifs);
				}


				curr_pt 	+= (4+4*numpart+4);
				
				// read ID
				rec 	= read_f77_record(ifs);
				std::memcpy(&gid, rec.data(), 4);
				curr_pt 	+= (4+4+4);

				// skip timestep
				skip_f77_record(ifs);
				curr_pt 	+= (4+4+4);

				// skip level ...
				skip_f77_record(ifs);
				curr_pt 	+= (4+4*5+4);

				// skip mass
				skip_f77_record(ifs);
				curr_pt 	+= (4+8+4);

				// skip x, y, z
				skip_f77_record(ifs);
				curr_pt 	+= (4+8*3+4);

				// skip vx, vy, vz
				skip_f77_record(ifs);
				curr_pt 	+= (4+8*3+4);

				// skip lx, ly, lz
				skip_f77_record(ifs);
				curr_pt 	+= (4+8*3+4);

				// skip shapes
				skip_f77_record(ifs);
				curr_pt 	+= (4+8*4+4);

				// skip energies
				skip_f77_record(ifs);
				curr_pt 	+= (4+8*3+4);

				// skip spin
				skip_f77_record(ifs);
				curr_pt 	+= (4+8+4);

				// skip sigma
				skip_f77_record(ifs);
				if(vh.horg == 'h'){
					curr_pt 	+= (4+8+4);
				}else if(vh.horg == 'g'){
					curr_pt 	+= (4+8*3+4);
				}
				

				// skip virial
				skip_f77_record(ifs);
				curr_pt 	+= (4+8*4+4);

				// skip profile
				skip_f77_record(ifs);
				curr_pt 	+= (4+8*2+4);

				// additional skip for GalaxyMaker Data
				if(vh.horg == 'g'){
					// # of bin
					rec 	= read_f77_record(ifs);
					std::memcpy(&nbin, rec.data(), 4);
					curr_pt 	+= (4+4*1+4);

					// radial bin
					skip_f77_record(ifs);
					curr_pt 	+= (4+8*nbin+4);

					// density
					skip_f77_record(ifs);
					curr_pt 	+= (4+8*nbin+4);
				}
				// input
				in_gpt(gpt, gid, curr_pt_old);
				
				if(id0>0){
					if(id0 == gid){
						gal[0].snap 	= snap_curr;
						gal[0].id 		= gid;
						gal[0].npart 	= numpart;

						if(readpart){
							gal[0].pid.resize(numpart);
							std::memcpy(gal[0].pid.data(), rec2.data(), numpart * sizeof(HM_PID));
						}
					}
				}else{
					gal[i].snap 	= snap_curr;
					gal[i].id 		= gid;
					gal[i].npart 	= numpart;

					if(readpart){
						gal[i].pid.resize(numpart);
						std::memcpy(gal[i].pid.data(), rec2.data(), numpart * sizeof(HM_PID));
					}
				}
			}


			return gal;
		}else{
			GalArray gal;

			HM_I32 nbodies, nmain, nsub, nall;
			HM_I32 numpart, gid;
			std::vector<char> rec2;

			if(id0>0){
				gal.resize(1);
				ifs.seekg(gpt[id0], std::ios::cur);

				// numpart
				auto rec 	= read_f77_record(ifs);
				std::memcpy(&numpart, rec.data(), 4);
				
				// read Particle ID
				if(readpart){
					rec2 	= read_f77_record(ifs);
					//pid.resize(numpart);
				}else{
					skip_f77_record(ifs);
				}
								
				// read ID
				rec 	= read_f77_record(ifs);
				std::memcpy(&gid, rec.data(), 4);
			
				gal[0].snap 	= snap_curr;
				gal[0].id 		= gid;
				gal[0].npart 	= numpart;

				if(readpart){
					gal[0].pid.resize(numpart);
					std::memcpy(gal[0].pid.data(), rec2.data(), numpart * sizeof(HM_PID));
				}

			}else{
				// READ nbodies
				auto rec 	= read_f77_record(ifs);
				std::memcpy(&nbodies, rec.data(), 4);


				// SKIP FOR DOUBLES (massp, aexp, omega_t, age_univ)
				for(int k=0; k<4; k++) skip_f77_record(ifs);
				
				// Read n_main & n_sub
				rec 	= read_f77_record(ifs);
				std::memcpy(&nmain, rec.data(), 4);
				std::memcpy(&nsub, rec.data() + 4, 4);


				nall 	= nmain + nsub;

				gal.resize(nall);


				// Read with loop
				for(HM_I32 i=0; i<nall; i++){

					// numpart
					rec 	= read_f77_record(ifs);
					std::memcpy(&numpart, rec.data(), 4);

					// read Particle ID
					if(readpart){
						rec2 	= read_f77_record(ifs);
					}else{
						skip_f77_record(ifs);
						}
					// read ID
					rec 	= read_f77_record(ifs);
					std::memcpy(&gid, rec.data(), 4);

					// skip timestep
					skip_f77_record(ifs);

					// skip level ...
					skip_f77_record(ifs);

					// skip mass
					skip_f77_record(ifs);

					// skip x, y, z
					skip_f77_record(ifs);

					// skip vx, vy, vz
					skip_f77_record(ifs);

					// skip lx, ly, lz
					skip_f77_record(ifs);

					// skip shapes
					skip_f77_record(ifs);

					// skip energies
					skip_f77_record(ifs);

					// skip spin
					skip_f77_record(ifs);

					// skip sigma
					skip_f77_record(ifs);

					// skip virial
					skip_f77_record(ifs);

					// skip profile
					skip_f77_record(ifs);

					// additional skip for GalaxyMaker Data
					if(vh.horg == 'g'){
						// # of bin
						skip_f77_record(ifs);

						// radial bin
						skip_f77_record(ifs);
					
						// density
						skip_f77_record(ifs);
					}

				
					gal[i].snap 	= snap_curr;
					gal[i].id 		= gid;
					gal[i].npart 	= numpart;

					if(readpart){
						gal[i].pid.resize(numpart);
						std::memcpy(gal[i].pid.data(), rec2.data(), numpart * sizeof(HM_PID));
					}
				}
							
			}
			return gal;
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

	//IO_dtype::snapinfo get_snapinfo(const vctree_set::Settings& vh); // unflagged
	bool is_snap(const vctree_set::Settings& vh, const IO_I32 snap_curr);

	//IO_dtype::snapSt read_info(const vctree_set::Settings& vh, const IO_I32 snap_curr);

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
//	inline IO_dtype::snapinfo get_snapinfo(vctree_set::Settings& vh){
//
//		std::vector<IO::IO_I32> snaplist;
//		IO::IO_I32 maxsnap = -1;
//
//		if(is_file(vh)){
//			std::ifstream file(vh.snaplist);
//		    std::string line;
//    		IO_I32 value;
//
//    		while (std::getline(file, line)) {
//        		if (line.empty()) continue;
//       			value = std::stoi(line);
//       			snaplist.push_back(value);
//
//       			if(value > maxsnap) maxsnap = value;
//        	}
//
//        	std::sort(snaplist.begin(), snaplist.end());
//        	vh.snapi 	= snaplist[0];
//        	vh.snapf 	= maxsnap;
//
//        	file.close();
//		}else{
//			for(IO::IO_I32 i=vh.snapi; i<vh.snapf+1; i++){
//				snaplist.push_back(i);
//				if(i > maxsnap) maxsnap = i;
//			}
//		}
//
//		//----- Allocate Snap
//		IO_dtype::snapinfo sinfo;
//		sinfo.resize(maxsnap+1);
//
//		//----- Read Info
//		for (IO::IO_I32 s : snaplist) {
//			sinfo[s].snum 	= s;
//		}
//
//		//----- Allocate key variables
//		if(vh.iotype == "HM"){
//   			vh.hm_gpointer.resize(vh.snapf+1);
//  		}
//
//		if(vh.iotype == "VR"){
//			vh.vr_gpointer.resize(vh.snapf+1);
//		}
//
//		return sinfo;
//	}

	//----- Read Catalog
	inline IO_dtype::GalArray r_gal(vctree_set::Settings& vh, const IO_dtype::IO_Snap snap_curr, const IO_dtype::IO_GID id0, const      bool readpart=false){
			IO_dtype::GalArray gal;
			if(vh.iotype == "VELUGA"){
					gal = IO_VELUGA::r_gal(vh, snap_curr, id0, readpart);
			}else if(vh.iotype == "HM"){
					gal = IO_HM::r_gal(vh, snap_curr, id0, readpart);
			}else if(vh.iotype == "VR"){
					gal = IO_VR::r_gal(vh, snap_curr, id0, readpart);
			//}else if(vh.iotype == "ANY"){
					//return IO_ANY::r_gal(vh, snap_curr, id0, readpart);
			}else{
					LOG()<<"Should be implemented for different IO Type";
					u_stop();
			}
			// Wrong ID check
			if(id0 <0){
					IO_dtype::IO_GID maxid = -1;
					for(IO_dtype::IO_I32 i=0; i< (IO_dtype::IO_I32) gal.size(); i++){
							if(gal[i].id >= maxid) maxid = gal[i].id;
					}
					std::vector<IO_dtype::IO_GID> idcheck(maxid+1, static_cast<IO_dtype::IO_GID>(-1));
					for(IO_dtype::IO_I32 i=0; i< (IO_dtype::IO_I32) gal.size(); i++){
							if(idcheck[gal[i].id]>0){
									LOG()<<"Duplicate ID exists (ID = "<<gal[i].id<<" )";
									u_stop();
							}
			
							idcheck[gal[i].id]=1;
					}
			
			}
	 
			return gal;
	}


}




