#include "global/allvar.h"
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

std::string normalize_iotype(std::string s, bool nocap) {
    auto not_space = [](unsigned char ch){ return !std::isspace(ch); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());

    if (s.size() >= 2) {
        char a = s.front(), b = s.back();
        if ((a=='"' && b=='"') || (a=='\'' && b=='\'')) {
            s = s.substr(1, s.size()-2);
        }
    }

    // To Captial letters
    if(!nocap)std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::toupper(c); });

    return s;
}

bool g_load_config(const std::string& path, vctree_set::Settings& vh) {
  std::ifstream in(path);
  if (!in) {
    LOG() << "config open failed: " << path << "\n";

    return false;
  }

  std::string line;
  int lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
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

    // IO settings
    if (lkey=="iotype") vh.iotype = normalize_iotype(val, false);
    //else if (lkey=="simtype") vh.simtype = normalize_iotype(val, false);
    else if (lkey=="out_dir") vh.out_dir = normalize_iotype(val, true);

    // VR related
    else if (lkey=="vr_dir_catalog") vh.vr_dir_catalog = normalize_iotype(val, true);
    //else if (lkey=="vr_dir_tree") vh.vr_dir_tree = normalize_iotype(val, true);

    // HMaker related
    else if (lkey=="hm_dir_catalog") vh.hm_dir_catalog = normalize_iotype(val, true);

    // Ramses related
    //else if (lkey=="ramses_dir") vh.ramses_dir = normalize_iotype(val, true);


    // Branch related
    else if (lkey=="branchmaker") vh.branchmaker = normalize_iotype(val, false);
    else if (lkey=="brtype") vh.brtype = normalize_iotype(val, false);
    else if (lkey=="treedir") vh.treedir = normalize_iotype(val, false);

    //// -- TreeFrog
    else if (lkey=="tf_dir") vh.tf_dir = normalize_iotype(val, true);

    // Load tree
    else if (lkey=="loadtree") vh.loadtree = normalize_iotype(val, false);
    else if (lkey=="loadtree_ftree") vh.loadtree_ftree = normalize_iotype(val, true);
    else if (lkey=="loadtree_fkey") vh.loadtree_fkey = normalize_iotype(val, true);
    
    // General (String)
    else if (lkey=="horg") {
      if (!val.empty()) vh.horg = (char)std::tolower((unsigned char)val[0]);
    }
    
      
    // General (number)
    else if (lkey=="snapi")  { parse_num(val, vh.snapi); }
    else if (lkey=="snapf")    { parse_num(val, vh.snapf); }
    else if (lkey=="snaplist")    { vh.snaplist = normalize_iotype(val, true); }
    else if (lkey=="minbranchlength") { parse_num(val, vh.minbranchlength); }

    // Make Branch
    else if (lkey=="meritlimit"){ parse_num(val, vh.meritlimit); }
    //else if (lkey=="snap_step")   { parse_num(val, vh.snap_step); }
    //else if (lkey=="makebr_nprog"){ parse_num(val, vh.makebr_nprog); }
   

    // Ctree
    else if (lkey=="ctree_minfrac")    { parse_num(val, vh.ctree_minfrac); }
    else if (lkey=="ctree_meritfrac")    { parse_num(val, vh.ctree_meritfrac); }
    else if (lkey=="ctree_n_search")   { parse_num(val, vh.ctree_n_search); }
    else if (lkey=="ctree_core_n")   { parse_num(val, vh.ctree_core_n); }
    else if (lkey=="ctree_core_dn")  { parse_num(val, vh.ctree_core_dn); }
    else if (lkey=="ctree_rfact")      { parse_num(val, vh.ctree_rfact); }
    else if (lkey=="ctree_weighttype") { parse_num(val, vh.ctree_weighttype); }
    else if (lkey=="ctree_makecheck") { parse_num(val, vh.ctree_makecheck); }
    else if (lkey=="ctree_loadcheck") { parse_num(val, vh.ctree_loadcheck); }
    //else if (lkey=="ctree_rerunmod")   { parse_num(val, vh.ctree_rerunmod); }
    //else if (lkey=="ctree_rerun")      { parse_num(val, vh.ctree_rerun); }
    else {
      std::cerr << "config: unknown key '" << key << "' at line " << lineno << "\n";
    }
  }

  vh.finalize_paths();

  if(vh.iotype == "HM"){
    vh.hm_gpointer.resize(vh.snapf+1);
  }
  return true;
}

