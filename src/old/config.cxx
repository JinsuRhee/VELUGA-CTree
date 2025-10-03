#include "config.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>

namespace {

  inline std::string trim(std::string s) {
    auto not_space = [](int ch){ return !std::isspace(ch); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
  }
  inline std::string tolower_copy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return (char)std::tolower(c); });
    return s;
  }
//inline bool parse_bool(const std::string& v, bool& out) {
//  auto s = tolower_copy(trim(v));
//  if (s=="1"||s=="true"||s=="yes"||s=="on")  { out=true;  return true; }
//  if (s=="0"||s=="false"||s=="no" ||s=="off"){ out=false; return true; }
//  return false;
//}
//inline bool parse_treedir(const std::string& v, TreeDir& out) {
//  auto s = tolower_copy(trim(v));
//  if (s=="des" || s=="desc" || s=="descendants") { out = TreeDir::Descendants; return true; }
//  if (s=="prg" || s=="prog" || s=="progenitors") { out = TreeDir::Progenitors; return true; }
//  return false;
//}

  template<typename T>
  bool parse_num(const std::string& s, T& out) {
    std::istringstream iss(trim(s));
    iss >> out;
    return (bool)iss && iss.eof();
  }


} // namespace

bool load_config(const std::string& path, Settings& vh) {
  std::ifstream in(path);
  if (!in) {
    std::cerr << "config open failed: " << path << "\n";
    return false;
  }

  std::string line;
  int lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    // 주석/빈줄 제거 (#, ;, // 지원)
    auto pos_hash = line.find('#');
    auto pos_semi = line.find(';');
    auto pos_sl   = line.find("//");
    size_t cut = std::min({pos_hash, pos_semi, pos_sl,
                           std::string::npos});
    if (cut != std::string::npos) line = line.substr(0, cut);
    line = trim(line);
    if (line.empty()) continue;

    auto eq = line.find('=');
    if (eq == std::string::npos) {
      std::cerr << "config: ignore line " << lineno << " (no '=')\n";
      continue;
    }

    std::string key = trim(line.substr(0, eq));
    std::string val = trim(line.substr(eq+1));
    auto lkey = tolower_copy(key);

    // General (String)
    if (lkey=="dir_catalog") vh.dir_catalog = val;
    else if (lkey=="dir_tree") vh.dir_tree = val;
    else if (lkey=="horg") {
      if (!val.empty()) vh.horg = (char)std::tolower((unsigned char)val[0]);
    }
    else if(lkey=="treedir") vh.treedir = val;
      
    // General (number)
    else if (lkey=="snapi")  { parse_num(val, vh.snapi); }
    else if (lkey=="snapf")    { parse_num(val, vh.snapf); }

    // Make Branch
    else if (lkey=="meritlimit"){ parse_num(val, vh.meritlimit); }
    //else if (lkey=="snap_step")   { parse_num(val, vh.snap_step); }
    //else if (lkey=="makebr_nprog"){ parse_num(val, vh.makebr_nprog); }
    else if (lkey=="makebr_bidkey"){ parse_num(val, vh.makebr_bidkey); }

    // Ctree
    else if (lkey=="ctree_n_search")   { parse_num(val, vh.ctree_n_search); }
    else if (lkey=="ctree_n_step_n")   { parse_num(val, vh.ctree_n_step_N); }
    else if (lkey=="ctree_n_step_dn")  { parse_num(val, vh.ctree_n_step_dn); }
    else if (lkey=="ctree_rfact")      { parse_num(val, vh.ctree_rfact); }
    //else if (lkey=="ctree_rerunmod")   { parse_num(val, vh.ctree_rerunmod); }
    //else if (lkey=="ctree_rerun")      { parse_num(val, vh.ctree_rerun); }
    else {
      std::cerr << "config: unknown key '" << key << "' at line " << lineno << "\n";
    }
  }



  vh.finalize_paths();
  return true;
}

// Modulate Settings based on the tree direction
bool set_config(Settings& vh){

  int32_t snapmin = std::min(vh.snapi, vh.snapf);
  int32_t snapmax = std::max(vh.snapi, vh.snapf);
  if(vh.treedir == "des"){

    vh.snapi  = snapmin;
    vh.snapf  = snapmax;

    vh.tag_num    = "NumDesc";
    vh.tag_off    = "DescOffsets";
    vh.tag_result = "Descendants";
    vh.tag_npart  = "DescNpart";
    vh.tag_merit  = "Merits";
    vh.tag_nlink  = "Nsteps_search_new_links";

  }
  else if(vh.treedir == "prg"){
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