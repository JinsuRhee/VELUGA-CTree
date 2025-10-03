#include <string>
//-----
// Run Control Variables
//    - not communicated by the configuration file
//-----
namespace Settings_int{
  inline int32_t sarr_size = 1000;
  inline int32_t barr_size = 100;
  inline int32_t tarr_size = 10000;
}

//-----
// Configuration related
//  - communicated by the configuration file
//-----
// Settings for the entire run
struct Settings {
  
  std::string dir_catalog = "./catalog/";   // VR output directory
  std::string dir_tree    = "";             // TreeFrog data is stored
  char        horg        = 'g';            // 'g' or 'h'
  std::string treedir     = "des";

  // General settings
  int32_t snapi   = 70;              // snapshot range
  int32_t snapf   = 200;
  double  meritlimit = 0.001;        // melit lower limit to close the branch

  // Branch making settings
  //int32_t makebr_nprog      = 5000;
  
  int64_t makebr_bidkey     = 1000;         // 10^n; should be larger than the last snapshot number


  // Ctree
  int32_t ctree_n_search   = 10;             // number of snapshots to be searched simultaneously
  int32_t ctree_n_step_N   = 10;             // number of branch points (>10) when collecting particles on a existing branch
  int32_t ctree_n_step_dn  = 5;              // dN between the points (corresponding to 200 MYr seems good)
  double  ctree_rfact      = 10.;            // Find Candidate galaxies within a sphere of a radius of speed * rfact * dT (10 seems good)
  
  //int32_t ctree_rerunmod   = 5;
  //int32_t ctree_rerun      = -1; // -1: off

  // TreeFrog I/O related
  std::string tag_num     = "";
  std::string tag_off     = "";
  std::string tag_result  = "";
  std::string tag_npart   = "";
  std::string tag_merit   = "";
  std::string tag_nlink   = "";

  // 유틸
  void finalize_paths() {
    if (dir_tree.empty()) {
      dir_tree = dir_catalog + (horg=='g' ? "Galaxy/tree" : "Halo/tree");
    }
  }

  // Bring control varibles
  int32_t sarr_size = Settings_int::sarr_size;
  int32_t barr_size = Settings_int::barr_size;
  int32_t tarr_size = Settings_int::tarr_size;
};

// Load Config File
bool load_config(const std::string& path, Settings& vh);
bool set_config(Settings& vh);