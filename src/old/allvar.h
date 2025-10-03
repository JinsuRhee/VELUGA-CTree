#pragma once
#include <string>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <cstddef>
#include <functional>
#include "autovec.h"

//-----
// General control variables
//-----
//namespace veluga_parameters{
//  inline int32_t tree_snapsize    = 1000;           // Initial snapshot size of branch arrays
//  inline int32_t tree_brsize  = 10;                 // Initial branch size of branch arrays
//}





//-----
// Make Branch Tree Data
//-----
//namespace Makebr{
//
//  
//  // Adjust Policy
//  struct BranchGrow {
//    virtual ~BranchGrow() = default;
//  };
//
//
//  // Branch Struct
//  struct Branch {
//    // snap array
//    AutoVec<int>    ID{-1}, snap{-1}, p_snap{-1}, p_id{-1}, d_snap{-1}, d_id{-1};
//    AutoVec<double> p_merit{-1.0};
//
//    // progenitor arrays
//    AutoVec<int>    m_id{-1}, m_snap{-1}, m_BID{-1};
//    AutoVec<double> m_merit{-1.0};
//
//    int f_BID   = -1;   // Father branch Index that this branch is merged
//
//    int endind  = -1;
//    int numprog = 1;
//
//    // 초기 추정치로 내부 길이 잡기 (필요 시)
//    Branch() = default;
//    Branch(int nsnap_guess, int nprog_guess) {
//      // 초기 길이를 잡되, 값은 -1/-1.0으로 채움
//      ID.resize(nsnap_guess);    snap.resize(nsnap_guess);
//      p_snap.resize(nsnap_guess); p_id.resize(nsnap_guess);
//      d_snap.resize(nsnap_guess); d_id.resize(nsnap_guess);
//      p_merit.resize(nsnap_guess);
//
//      m_id.resize(nprog_guess);  m_snap.resize(nprog_guess);
//      m_BID.resize(nprog_guess); m_merit.resize(nprog_guess);
//
//      endind = -1;
//      numprog = 1;
//    }
//  };
//
//  // ---------------------------------------------
//  // Branch fractory: Extend tree
//  // ---------------------------------------------
//  struct Branchfactory {
//    int nsnap_guess = 0;
//    int nprog_guess = 0;
//    Branch operator()(std::size_t) const {
//      return Branch(nsnap_guess, nprog_guess);
//    }
//  };
//
//  // ---------------------------------------------
//  // tree container: Make as AutoVec<Branch> , autumatic extendsion accessing to tree[i]
//  // ---------------------------------------------
//  using Tree = AutoVec<Branch>;
//
//  // 헬퍼: 초기 guess 크기의 tree 만들기 (초기 개수도 정하고 싶으면 resize 호출)
//  inline Tree make_tree(int nsnap_guess, int nprog_guess, std::size_t initial_tree_guess = 0) {
//    Tree t{ BranchFactory{nsnap_guess, nprog_guess} }; // 팩토리 기반
//    if (initial_tree_guess) t.resize(initial_tree_guess);
//    return t;
//  }
//
//};
