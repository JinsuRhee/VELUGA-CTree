#pragma once
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>   // C++17
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include "global/allvar.h"

#ifdef CTREE_USE_MPI
    #include <mpi.h>
#endif


// Logo Printing
bool u_printlogo();

// Run Stat check
bool u_initialcheck(int argc);

// String
bool is_file(std::string file);
std::string i4(int x);
std::string i5(int x);
std::string i6(int x);
std::string iN(int x, int N);

// Get Memory usage
template <class T> inline double how_big(const std::vector<T>& v) {
    return static_cast<double>(v.capacity()) * sizeof(T) / (1024.0L * 1024.0L * 1024.0L);
}



// Get Mpi rank
int mpi_rank();

// Stop the program
[[noreturn]] void u_stop();

// Save the tree
template <typename T>
constexpr std::int32_t savetree_gettype() {
    if constexpr (std::is_same_v<T, std::int32_t>)  return 1;
    else if constexpr (std::is_same_v<T, std::int64_t>) return 2;
    else if constexpr (std::is_same_v<T, float>) return 3;
    else if constexpr (std::is_same_v<T, double>) return 4;
    else return -1;
}

void savetree_base(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& key, bool usename);



//----- Load related (should be modulated in the future)
std::int32_t loadtree_getsize(std::int32_t tag);
std::string loadtree_gettypename(std::int32_t tag);
template <typename T> std::string loadtree_gettypenamebysize(T var);

void loadtree_base(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& treekey, bool usename);


template <typename T>
static inline void loadtree_read(std::istream& in, T& v) {
    in.read(reinterpret_cast<char*>(&v), sizeof(T));
    if (!in) throw std::runtime_error("binary read failed");
}
//
template <typename T>
inline void loadtree_vecread(std::istream& in, std::vector<T>& vec, std::size_t n) {
    vec.resize(n);
    if (n == 0) return;
    in.read(reinterpret_cast<char*>(vec.data()), sizeof(T)*n);
    if (!in) throw std::runtime_error("binary vector read failed");
}



//template <class SrcT, class BidT>
//void read_bid_scatter(std::istream& in, std::vector<Tree::TreeKey>& treekey, std::size_t n) {
//    constexpr std::size_t BLOCK = 1 << 15; // 32K elements
//    std::vector<SrcT> buf(std::min(n, BLOCK));
//
//    std::size_t i = 0;
//    while (i < n) {
//        std::size_t m = std::min(BLOCK, n - i);
//        in.read(reinterpret_cast<char*>(buf.data()), sizeof(SrcT) * m);
//        if (!in) throw std::runtime_error("binary read failed");
//
//        auto* tk = treekey.data() + i;
//        for (std::size_t j = 0; j < m; ++j)
//            tk[j].bid = static_cast<BidT>(buf[j]);
//        i += m;
//    }
//}
//
//
//inline void load_treekey_from_stream(std::istream& in, Tree::TreeKeyArray& treekey) {
//    
//    // Check BID data type
//    std::int32_t tag = -1;
//    loadtree_read<std::int32_t>(in, tag);
//
//
//    const std::int32_t TAG_BID = savetree_gettype<Tree::Tree_BID>();
//    if (tag != TAG_BID){
//    	LOG()<<"should be considered here";
//    	u_stop();
//    }
//
//    // Get the number of key elements
//    std::int64_t nkey = -1;
//    loadtree_read<std::int64_t>(in, nkey);
//
//
//    // Read The tree Key
//    std::int64_t keyval = -1;
//    loadtree_read<std::int64_t>(in, keyval);
//    
//
//
//    // Read Key
//    //auto keyind = loadtree_mkvector(tag, nkey);
//    if(tag == 1){ std::vector<std::int32_t> keyind; }
//    else if(tag == 2){ std::vector<std::int64_t> keyind; }
//    else{
//    	std::vector<std::int32_t> keyind;
//    	LOG()<<"Wrong Tag for tree type : tag is beyond the allowed range";
//    	u_stop();
//    }
//    
//	keyind.resize(nkey);    
//	loadtree_vecread(in, keyind, nkey);
//
//    // Allocate and copy
//    treekey.resize(nkey);
//
//    treekey[0].key 	= keyval;
//    for (std::int64_t i = 0; i < nkey; i++) {
//        treekey[i].bid = keyind[i];
//	}
//}
//
////-----
//// MPI related
////-----
//#if CTREE_USE_MPI
//inline std::vector<char> read_file_bytes(const std::string& path) {
//    std::ifstream in(path, std::ios::binary);
//    if (!in) throw std::runtime_error("cannot open: " + path);
//    in.seekg(0, std::ios::end);
//    const std::streamsize sz = in.tellg();
//    in.seekg(0, std::ios::beg);
//    std::vector<char> buf(static_cast<size_t>(sz));
//    if (sz > 0) in.read(buf.data(), sz);
//    if (!in) throw std::runtime_error("file read failed: " + path);
//    return buf;
//}
//
//inline void mpi_bcast_bytes(std::vector<char>& buf, int root, MPI_Comm comm) {
//    // Size broad cast
//    std::int64_t n = static_cast<long long>(buf.size());
//    MPI_Bcast(&n, 1, MPI_LONG_LONG, root, comm);
//    // Resize in Sender
//    if (int rank = 0; (MPI_Comm_rank(comm, &rank), true), rank != root) {
//        buf.resize(static_cast<size_t>(n));
//    }
//    // Send the byte
//    if (n > 0) {
//        MPI_Bcast(buf.data(), static_cast<int>(n), MPI_BYTE, root, comm);
//    }
//}
//#endif

