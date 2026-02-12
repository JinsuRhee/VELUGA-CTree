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


namespace IO_RAMSES{

//	IO_dtype::snapinfo get_snapinfo(const vctree_set::Settings& vh){
//
//
//		//----- Gather snapshot list
//		std::vector<IO::IO_I32> snaplist;
//		IO::IO_I32 maxsnap = -1;
//		for(IO::IO_I32 i=vh.snapi; i<vh.snapf+1; i++){
//
//			if(!is_snap(vh, i)) continue;
//			snaplist.push_back(i);
//			if(i > maxsnap) maxsnap = i;		
//		}
//
//		//----- Allocate Snap
//		IO_dtype::snapinfo sinfo;
//		sinfo.resize(maxsnap+1);
//
//
//		//----- Read Info
//		for (IO::IO_I32 s : snaplist) {
//			sinfo[s]	= read_info(vh, s);
//		}
//
//		//----- Compute the age of the universe
//
//		return sinfo;
//	}

//	IO_dtype::snapSt read_info(const vctree_set::Settings& vh, const IO_I32 snap_curr){
//
//		std::string infoname 	= vh.ramses_dir + "/output_" + i5(snap_curr) + "/info_" + i5(snap_curr) + ".txt";
//
//		std::ifstream fin(infoname);
//    	if (!fin){
//    		throw std::runtime_error("cannot open file: " + infoname);	
//    	} 
//
//    	// initialize pointer
//    	fin.clear();
//		fin.seekg(0);
//    	//if (auto* f = dynamic_cast<std::ifstream*>(&fin)) { f->clear(); f->seekg(0); }
//    	//if (auto* s = dynamic_cast<std::istringstream*>(&fin)) { s->clear(); s->seekg(0); }
//
//    	std::string line;
//    	for(int i=1; i<10; i++) std::getline(fin, line);
//
//    	// aexp
//		std::getline(fin, line);
//    	std::string t1 = trim(line);
//
//    	// H0
//    	std::getline(fin, line);
//    	std::string t_h0 = trim(line);
//
//    	// om
//    	std::getline(fin, line);
//    	std::string t_om = trim(line);
//
//    	// oL
//    	std::getline(fin, line);
//    	std::string t_ol = trim(line);
//
//    	for(int i=14; i<16; i++) std::getline(fin, line);
//    	
//    	// unit_l
//    	std::getline(fin, line);
//    	std::string t2 = trim(line);
//
//    	for(int i=17; i<18; i++) std::getline(fin, line);
//
//    	// unit_t
//    	std::getline(fin, line);
//    	std::string t3 = trim(line);
//
//    	// extract
//    	auto eq1 	= t1.find('=');
//    	auto eq2 	= t2.find('=');
//    	auto eq3 	= t3.find('=');
//
// 		auto eq_h0	= t_h0.find('=');
// 		auto eq_om 	= t_om.find('=');
// 		auto eq_ol 	= t_ol.find('=');
//
//    	// numbers
//    	std::string rhs1, rhs2, rhs3, rhs_h0, rhs_om, rhs_ol;
//    	rhs1 	= trim(t1.substr(eq1+1));
//    	rhs2 	= trim(t2.substr(eq2+1));
//    	rhs3 	= trim(t3.substr(eq3+1));
//
//    	rhs_h0	= trim(t_h0.substr(eq_h0+1));
//    	rhs_om	= trim(t_om.substr(eq_om+1));
//    	rhs_ol	= trim(t_ol.substr(eq_ol+1));
//
//    	// input
//    	IO_dtype::snapSt thissnap;
//
//    	parse_num<IO_dtype::IO_double>(rhs1, thissnap.aexp);
//    	parse_num<IO_dtype::IO_double>(rhs2, thissnap.unit_l);
//    	parse_num<IO_dtype::IO_double>(rhs3, thissnap.unit_t);
//    	parse_num<IO_dtype::IO_double>(rhs_h0, thissnap.h0);
//    	parse_num<IO_dtype::IO_double>(rhs_om, thissnap.om);
//    	parse_num<IO_dtype::IO_double>(rhs_ol, thissnap.ol);
//
//    	thissnap.snum 	= snap_curr;
//
//    	return thissnap;
//	}


//	bool is_snap(const vctree_set::Settings& vh, const IO::IO_I32 snap_curr){
//		char buf[256], buf2[256];
//    	std::snprintf(buf, sizeof(buf), "output_%05d", snap_curr);
//    	std::snprintf(buf2, sizeof(buf2), "info_%05d.txt", snap_curr);
//    	std::string dumfname = vh.ramses_dir + "/" + buf + "/" + buf2;
//    		
//    	// Is File?
//    	if (fs::exists(dumfname)) {
//        	return true;
//    	} else{
//    		return false;
//    	}
//	}
	
	
}
