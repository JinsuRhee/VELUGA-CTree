
#include "utilities/utilities.h"
#include "global/io.h"

#include <variant>







// For logo print
bool u_printlogo() {

    const std::string logopath = std::string(PROJECT_SOURCE_DIR) + "/src/veluga_logo.txt";


    std::ifstream in(logopath);
    if (!in) {
        std::cerr << "failed to open: " << logopath << '\n';
        return false;
    }
    std::string line;
    while (std::getline(in, line)) {
        std::cout << line << '\n';
    }
    return true;
}

// For Initial Run stat Check
bool u_initialcheck(int argc) {
    if (argc < 2) {
        LOG() << "config file is required";
        return false;
    } else{
        return true;
    }
}

// padding
std::string i4(int x) {
    std::ostringstream oss;
    oss << std::setw(4) << std::setfill('0') << x;
    return oss.str();
}

std::string i5(int x) {
    std::ostringstream oss;
    oss << std::setw(5) << std::setfill('0') << x;
    return oss.str();
}

std::string i6(int x) {
    std::ostringstream oss;
    oss << std::setw(6) << std::setfill('0') << x;
    return oss.str();
}

// padding
std::string iN(int x, int N) {
    std::ostringstream oss;
    oss << std::setw(N) << std::setfill('0') << x;
    return oss.str();
}

// Is file
bool is_file(std::string file){
    if (std::filesystem::exists(file)) {
        return true;
    } else{
        return false;
    }
}

// Get MPI rank
int mpi_rank() {
#ifdef CTREE_USE_MPI
  static int r = []{ int rr=0; MPI_Comm_rank(MPI_COMM_WORLD,&rr); return rr; }();
  return r;
#else
  return 0;
#endif
}

// Stop the Program
void u_stop() {
    int errcode = 1;
    MPI_Abort(MPI_COMM_WORLD, errcode);
    std::exit(errcode);
}

// Get Memory usage


// Save Tree
void savetree_base(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& treekey, bool usename){

    

    namespace fs = std::filesystem;

    constexpr std::int32_t TAG_GID = savetree_gettype<Tree::Tree_GID>();
    constexpr std::int32_t TAG_BID = savetree_gettype<Tree::Tree_BID>();
    constexpr std::int32_t TAG_SNAP = savetree_gettype<Tree::Tree_Snap>();
    constexpr std::int32_t TAG_MERIT = savetree_gettype<Tree::Tree_merit>();



    //-----
    // Save TreeKey
    //-----
    std::string path;
    if(usename == true) path = vh.loadtree_fkey;
    else path = vh.out_dir + "/ctree_key.dat";
    LOG()<<"    Writing TreeKey in "<<path;

    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("save_treekey_bin: cannot open " + path);

    // First element as type specifer
    //      I32 type
    //      1 I32
    //      2 I64
    //      3 float
    //      4 double
    out.write(reinterpret_cast<const char*>(&TAG_BID), sizeof(TAG_BID));

    // Second elements Key size in I64
    std::int64_t n = static_cast<std::int64_t>(treekey.size());
    out.write(reinterpret_cast<const char*>(&n), sizeof(n));

    // Third elements as Tree Key
    std::int64_t tk = treekey[0];
    out.write(reinterpret_cast<const char*>(&tk), sizeof(tk));

    // Left for elements
    std::vector<Tree::Tree_BID> inds;
    inds.reserve(treekey.size());
    for (const auto& k : treekey) inds.push_back(k);

    if (!inds.empty()) {
        out.write(reinterpret_cast<const char*>(inds.data()),
                  inds.size() * sizeof(Tree::Tree_BID));
    }


    //-----
    // Save Tree
    //-----
    std::string path2;
    if(usename == true) path2 = vh.loadtree_ftree;
    else path2 = vh.out_dir + "/ctree_tree.dat";
    LOG()<<"    Writing Tree in "<<path2;

    std::ofstream out2(path2, std::ios::binary);
    if (!out2) throw std::runtime_error("save_tree_bin: cannot open " + path2);

    // 1 GID specifer
    out2.write(reinterpret_cast<const char*>(&TAG_GID), sizeof(TAG_GID));
    // 2 Snap specifer
    out2.write(reinterpret_cast<const char*>(&TAG_SNAP), sizeof(TAG_SNAP));
    // 3 BID specifer
    out2.write(reinterpret_cast<const char*>(&TAG_BID), sizeof(TAG_BID));
    // 4 merit specifer
    out2.write(reinterpret_cast<const char*>(&TAG_MERIT), sizeof(TAG_MERIT));


    // 5 Tree Size in I64
    //std::int64_t n2 = static_cast<std::int64_t>(tree.size());
    //out2.write(reinterpret_cast<const char*>(&n2), sizeof(n2));

    // 6 End ind in I64
    std::int64_t n3 = static_cast<std::int64_t>(tree[0].lind);
    out2.write(reinterpret_cast<const char*>(&n3), sizeof(n3));

    // loop for tree
    
    for(auto t : tree){
        
        Tree::Tree_I32 n_branch = t.endind + 1;
        Tree::Tree_I32 n_numprg = t.numprog;
        Tree::Tree_BID father_bid = t.father_bid;
        Tree::Tree_I32 t_stat   = t.stat;
      
        
        // Size for the main branch
        out2.write(reinterpret_cast<const char*>(&n_branch), sizeof(n_branch));

        // skip empty tree
        if(n_branch <= 0) continue;

        
        // Size for the merged branch
        out2.write(reinterpret_cast<const char*>(&n_numprg), sizeof(n_numprg));

        // Branch ID of father
        out2.write(reinterpret_cast<const char*>(&father_bid), sizeof(father_bid));

        // Tree Stat
        out2.write(reinterpret_cast<const char*>(&t_stat), sizeof(t_stat));

        // Branch_ID
        for(std::int32_t i=0; i<n_branch; i++){
            out2.write(reinterpret_cast<const char*>(&t.id[i]), sizeof(t.id[i]));
        }

        // Branch_snap
        for(std::int32_t i=0; i<n_branch; i++){
            out2.write(reinterpret_cast<const char*>(&t.snap[i]), sizeof(t.snap[i]));
        }

        // Branch_ID (progenitor)
//        for(std::int32_t i=0; i<n_branch; i++){
//            out2.write(reinterpret_cast<const char*>(&t.p_id[i]), sizeof(t.p_id[i]));
//        }
//
//        // Branch_Snap (progenitor)
//        for(std::int32_t i=0; i<n_branch; i++){
//            out2.write(reinterpret_cast<const char*>(&t.p_snap[i]), sizeof(t.p_snap[i]));
//        }
//
        // Branch_Merit (progenitor)
        for(std::int32_t i=0; i<n_branch; i++){
            out2.write(reinterpret_cast<const char*>(&t.p_merit[i]), sizeof(t.p_merit[i]));
        }

        if(n_numprg>0){

            // Merged Branch ID
            for(std::int32_t i=0; i<n_numprg; i++){
                out2.write(reinterpret_cast<const char*>(&t.m_id[i]), sizeof(t.m_id[i]));
            }

            // Merged Branch Snap
            for(std::int32_t i=0; i<n_numprg; i++){
                out2.write(reinterpret_cast<const char*>(&t.m_snap[i]), sizeof(t.m_snap[i]));
            }

            // Merged Branch Merit
            for(std::int32_t i=0; i<n_numprg; i++){
                out2.write(reinterpret_cast<const char*>(&t.m_merit[i]), sizeof(t.m_merit[i]));
            }

            // Merged Branch BID
            for(std::int32_t i=0; i<n_numprg; i++){
                out2.write(reinterpret_cast<const char*>(&t.m_bid[i]), sizeof(t.m_bid[i]));
            }
        }
    }  
}

//-----
// Load Part
//-----
std::int32_t loadtree_getsize(std::int32_t tag){
    switch(tag){
        case 1:
            return (std::int32_t) sizeof(std::int32_t);
        case 2:
            return (std::int32_t) sizeof(std::int64_t);
        case 3:
            return (std::int32_t) sizeof(float);
        case 4:
            return (std::int32_t) sizeof(double);
    }

    return -1;
}

std::string loadtree_gettypename(std::int32_t tag){
    switch(tag){
        case 1:
            return "int32";
        case 2:
            return "int64";
        case 3:
            return "float";
        case 4:
            return "double";
    }

    return "no matched case";
}

template <typename T> std::string loadtree_gettypenamebysize(T var){
    std::int32_t varsize = sizeof(var);

    if(varsize == sizeof(std::int32_t) || varsize == sizeof(float)){
        return "int32 or float";
    }
    if(varsize == sizeof(std::int64_t) || varsize == sizeof(double)){
        return "int64 or double";
    }

    return "no matched case";
}

void loadtree_base(vctree_set::Settings& vh, Tree::TreeArray& tree, Tree::TreeKeyArray& treekey, bool usename){

    int myrank  = mpi_rank();
    // Simple version used
    {
        // file check
        

        std::string path; 
        
        if(usename == true) path = vh.loadtree_fkey;
        else path = vh.out_dir + "/ctree_key.dat";

        if(myrank==0){
            LOG() << "    Reading TreeKey from " << path;
        }

        std::ifstream in(path, std::ios::binary);
        if (!in) {
            LOG()<<"load_treekey_bin: cannot open " + path;
            u_stop();
        }

        // Read Data Type
        std::int32_t bidtag;
        loadtree_read(in, bidtag);
    
        // Read Nkey
        std::int64_t nbid;
        loadtree_read(in, nbid);

        // TreeKey
        std::int64_t keyval;
        loadtree_read(in, keyval);

        // Validate the type
        if(sizeof(IO_dtype::IO_BID) != loadtree_getsize(bidtag)){
            LOG()<<"Wrong type match for branch index";
            LOG()<<"    given : "<<loadtree_gettypename(bidtag);
            LOG()<<" compiled : "<<loadtree_gettypenamebysize((Tree::Tree_BID) 1);
            u_stop();
        }

        // Allocate & Read
        //Tree::TreeKeyArray treekey;
        treekey.resize(nbid);

        loadtree_vecread<Tree::Tree_BID>(in, treekey, nbid);


        treekey[0]  = keyval;
//        if(bidtag == 1){
//            std::vecotr<std::int32_t> tree_keyload;
//            treekey_load.resize(nbid);
//            loadtree_vecread<std::int32_t>(in, treekey_load, nbid);
//        }else if(bidtag == 2){
//            std::vecotr<std::int64_t> tree_keyload;
//            treekey_load.resize(nbid);
//            loadtree_vecread<std::int64_t>(in, treekey_load, nbid);
//        }
//        else{
//            LOG()<<"Wrong Tag for branch index given (current): "<<bidtag;
//            LOG()<<"    (allowed): 1 or 2";
//            u_stop();
//        }
//
//        if(sizeof(IO_dtype::IO_BID) != sizeof(tree_keyload[0])){
//            LOG()<"Incorrect matching for branch index type :";
//            LOG()<<"        Given type = "<<sizeof(tree_keyload[0]);
//            LOG()<<"        Compiled type = "<<Sizeof(IO_dtype::IO_BID);
//
//            u_stop();
//        }
        //auto treekey_load = loadtree_mkvector(bidtag, nbid);
        //std::vector<std::int32_t> treekey_load;
        
//        LoadVector x;
//
//        if(bidtag == 1){
//            auto& v = std::get<std::vector<int32_t>>(treekey_load);
//            x = v[0];
//        }else if(bidtag == 2){
//            auto& v = std::get<std::vector<int64_t>>(treekey_load);
//            x = v[0];
//        }else{
//            LOG()<<"Wrong Tag for branch index given (current): "<<bidtag;
//            LOG()<<"    (allowed): 1 or 2";
//            u_stop();
//        }
//            
//        std::visit([](auto val) {}, x);

//        // Copy
//        Tree::TreeKeyArray treekey;
//        treekey.resize(nbid);
//
//        for(std::int64_t i=0; i<nbid; i++){
//            treekey[i]  = x[i];
//        }
//
//        treekey[0] = keyval;


    } 


    {
        
        std::string path;
        if(usename == true) path = vh.loadtree_ftree;
        else path = vh.out_dir + "/ctree_tree.dat";

        if(myrank==0){
            LOG() << "    Reading Tree from " << path;
        }

        std::ifstream in(path, std::ios::binary);
        if (!in) {
            LOG()<<"load_treekey_bin: cannot open " + path;
            u_stop();
        }

        // Read Data Type
        std::int32_t gidtag, snaptag, bidtag, mertag;
        std::int32_t nbranch, nmerge, fid, t_stat;



        loadtree_read(in, gidtag);
        loadtree_read(in, snaptag);
        loadtree_read(in, bidtag);
        loadtree_read(in, mertag);

        //----- Validate
        if(sizeof(IO_dtype::IO_GID) != loadtree_getsize(gidtag)){
            LOG()<<"Wrong type match for ID";
            LOG()<<"    given : "<<loadtree_gettypename(gidtag);
            LOG()<<" compiled : "<<loadtree_gettypenamebysize((IO_dtype::IO_GID) 1);
            u_stop();
        }

        if(sizeof(IO_dtype::IO_Snap) != loadtree_getsize(snaptag)){
            LOG()<<"Wrong type match for snapshot number";
            LOG()<<"    given : "<<loadtree_gettypename(snaptag);
            LOG()<<" compiled : "<<loadtree_gettypenamebysize((IO_dtype::IO_Snap) 1);
            u_stop();
        }

        if(sizeof(IO_dtype::IO_BID) != loadtree_getsize(bidtag)){
            LOG()<<"Wrong type match for branch index";
            LOG()<<"    given : "<<loadtree_gettypename(bidtag);
            LOG()<<" compiled : "<<loadtree_gettypenamebysize((IO_dtype::IO_BID) 1);
            u_stop();
        }

        if(sizeof(IO_dtype::IO_merit) != loadtree_getsize(mertag)){
            LOG()<<"Wrong type match for merit";
            LOG()<<"    given : "<<loadtree_gettypename(mertag);
            LOG()<<" compiled : "<<loadtree_gettypenamebysize((IO_dtype::IO_merit) 1);
            u_stop();
        }

        //-----
    
        // Read Ntree
        //std::int64_t ntree;
        //loadtree_read(in, ntree);

        std::int64_t lastind;
        loadtree_read(in, lastind);

        //tree.resize(ntree);
        
        tree.resize(lastind+vh.ctree_nstep);
        tree[0].lind = lastind;

        //for(std::int64_t i=0; i<ntree; i++){
        for(std::int64_t i=0; i<tree[0].lind+1; i++){

            loadtree_read(in, nbranch);

            tree[i].endind = nbranch-1;

            if(nbranch <= 0) continue;

            loadtree_read(in, nmerge);
            loadtree_read(in, fid);
            loadtree_read(in, t_stat);

            tree[i].father_bid = fid;
            tree[i].numprog = nmerge;
            tree[i].stat    = t_stat;

            tree[i].id.resize(nbranch);
            tree[i].snap.resize(nbranch);
//            tree[i].p_id.resize(nbranch);
//            tree[i].p_snap.resize(nbranch);
            tree[i].p_merit.resize(nbranch);

            loadtree_vecread<IO_dtype::IO_GID>(in, tree[i].id, nbranch);
            loadtree_vecread<IO_dtype::IO_Snap>(in, tree[i].snap, nbranch);
            //loadtree_vecread<std::int32_t>(in, tree[i].p_id, nbranch);
            //loadtree_vecread<std::int32_t>(in, tree[i].p_snap, nbranch);
            //loadtree_vecread<double>(in, tree[i].p_merit, nbranch);

            if(nmerge >= 1){
                tree[i].m_id.resize(nmerge);
                tree[i].m_snap.resize(nmerge);
                tree[i].m_merit.resize(nmerge);
                tree[i].m_bid.resize(nmerge);

                loadtree_vecread<IO_dtype::IO_GID>(in, tree[i].m_id, nmerge);
                loadtree_vecread<IO_dtype::IO_Snap>(in, tree[i].m_snap, nmerge);
                loadtree_vecread<IO_dtype::IO_merit>(in, tree[i].m_merit, nmerge);
                loadtree_vecread<IO_dtype::IO_BID>(in, tree[i].m_bid, nmerge);
            }


        }




    }
    // -----
    // Load TreeKey
    // -----
//    std::vector<char> key_buf, tree_buf;
//    Tree::TreeArray tree;
//    Tree::TreeKeyArray treekey;
//    int myrank  = mpi_rank();
//
//#ifdef CTREE_USE_MPI
//    MPI_Barrier(MPI_COMM_WORLD);
//
//
//    if (myrank == 0) {
//        key_buf  = read_file_bytes(vh.out_dir + "/ctree_key.dat");
//        tree_buf = read_file_bytes(vh.out_dir + "/ctree_tree.dat");
//        LOG() << "    Reading Tree and TreeKey from ";
//        LOG() << "      tree : "<<tree_buf;
//        LOG() << "      key  : "<<key_buf;
//    }
//
//    mpi_bcast_bytes(key_buf,  0, MPI_COMM_WORLD);
//    mpi_bcast_bytes(tree_buf, 0, MPI_COMM_WORLD);
//
//    {
//        std::string key_str(key_buf.begin(), key_buf.end());
//        std::istringstream key_in(key_str, std::ios::binary);
//        load_treekey_from_stream(key_in, treekey);
//    }
//    {
//        std::string tree_str(tree_buf.begin(), tree_buf.end());
//        std::istringstream tree_in(tree_str, std::ios::binary);
//        //load_tree_from_stream(tree_in, tree);
//    }
//
//    MPI_Barrier(MPI_COMM_WORLD);
//
//#else
//    {
//        std::ifstream in(vh.out_dir + "/ctree_key.dat", std::ios::binary);
//        load_treekey_from_stream(in, treekey);
//    }
//    {
//        std::ifstream in(vh.out_dir + "/ctree_tree.dat", std::ios::binary);
//        //load_tree_from_stream(in, tree);
//    }
//
//#endif
}


