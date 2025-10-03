#include "utilities.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <filesystem>   // C++17
#include "utilities.h"





// For logo print
bool util_printlogo() {

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

// For modifying TreeFrog Result file name
// src/tfcname.cxx
namespace fs = std::filesystem;

namespace {

    // Extract Integer from tree.snaplist.txt
    inline bool extract_snap_from_line(const std::string& line, int& out) {
        auto p = line.find("snap_");
        if (p == std::string::npos) return false;
        if (p + 5 + 4 > line.size()) return false;            // 최소 4자리
        try {
            out = std::stoi(line.substr(p + 5, 4));
            return true;
        } catch (...) { return false; }
    }

    // Extract Integer "tree.snapshot_0123.tree"
    inline bool extract_number_from_filename(const std::string& filename, int& out) {
        auto us = filename.find_last_of('_');
        auto dot = filename.find_last_of('.');
        if (us == std::string::npos || dot == std::string::npos || us+1 >= dot) return false;
        try {
            out = std::stoi(filename.substr(us + 1, dot - (us + 1)));
            return true;
        } catch (...) { return false; }
    }

    // padding
    inline std::string i4(int x) {
        std::ostringstream oss;
        oss << std::setw(4) << std::setfill('0') << x;
        return oss.str();
    }
} // namespace


bool util_tfcname(const std::string& dir_tree) {
  try {
    const fs::path tfout = fs::path(dir_tree) / "tfout";

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
      std::cerr << "[tfcname] count mismatch: "
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
        std::cerr << "[tfcname] target already exists, skipping: " << target << "\n";
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


