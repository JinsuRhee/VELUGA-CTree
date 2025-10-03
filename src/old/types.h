#pragma once
#include "autovec.h"


namespace Makebr {
  // Tree Log
  struct Treelog{
    std::int32_t n_new = 0;
    std::int32_t n_link = 0;
    std::int32_t n_link2 = 0;
    std::int32_t n_link3 = 0;
    std::int32_t n_broken = 0;
    std::int32_t n_all = 0;
    std::int32_t max_id = -1;
  };

  // Evolution check dummy array
  struct EvoItem {
    std::int32_t ID   = -1;   // LONARR → 32-bit 정수
    std::int32_t snap = -1;   // LONARR
    double       merit= -1.0; // DBLARR
  };

  using Evoldum = AutoVec<EvoItem>;

  // Generating helper
  inline Evoldum make_evoldum(std::size_t initial_guess = 0,
                            std::size_t growth_step   = 1000) {
    // Extension
    Evoldum e{EvoItem{}, growth_step};
    if (initial_guess) e.resize(initial_guess);
    return e;
  }


  // Tree Array
  struct Branch {
    // snaparr 계열 (기본 -1/-1.0)
    AutoVec<int>    ID{-1}, snap{-1}, p_snap{-1}, p_id{-1}, d_snap{-1}, d_id{-1};
    AutoVec<double> p_merit{-1.0};

    // galarr 계열
    AutoVec<int>    m_id{-1}, m_snap{-1}, m_BID{-1};
    AutoVec<double> m_merit{-1.0};

    int f_BID = -1;
    int endind  = -1;
    int numprog = 1;

    // tree Stat
    //  1   : processing
    //  <0  : broken
    //  0   : completed

    int stat = 1;

    Branch() = default;

    // 초기 추정치 + 증가 단위(step) 지정
    Branch(int nsnap_guess, int nprog_guess,
           std::size_t snap_step = 1, std::size_t prog_step = 1)
    {
      set_snap_growth_step(snap_step);
      set_prog_growth_step(prog_step);

      // 초기 길이는 step의 배수로 맞춰서 확보됨
      ID.resize(nsnap_guess);    snap.resize(nsnap_guess);
      p_snap.resize(nsnap_guess); p_id.resize(nsnap_guess);
      d_snap.resize(nsnap_guess); d_id.resize(nsnap_guess);
      p_merit.resize(nsnap_guess);

      m_id.resize(nprog_guess);  m_snap.resize(nprog_guess);
      m_BID.resize(nprog_guess); m_merit.resize(nprog_guess);
    }

    // 나중에 동적으로 step 바꾸고 싶을 때
    void set_snap_growth_step(std::size_t step) {
      ID.set_growth_step(step); snap.set_growth_step(step);
      p_snap.set_growth_step(step); p_id.set_growth_step(step);
      d_snap.set_growth_step(step); d_id.set_growth_step(step);
      p_merit.set_growth_step(step);
    }
    void set_prog_growth_step(std::size_t step) {
      m_id.set_growth_step(step); m_snap.set_growth_step(step);
      m_BID.set_growth_step(step); m_merit.set_growth_step(step);
    }
  };

  // 새 Branch 생성용 팩토리(outer tree가 커질 때 사용)
  struct BranchFactory {
    int nsnap_guess = 0;
    int nprog_guess = 0;
    std::size_t snap_step = 1;
    std::size_t prog_step = 1;

    Branch operator()(std::size_t) const {
      return Branch(nsnap_guess, nprog_guess, snap_step, prog_step);
    }
  };

  // tree: Branch의 자동 확장 배열
  using Tree = AutoVec<Branch>;

  // tree 생성: 각 step 지정 가능 (tree 자체 step도 설정)
  inline Tree make_tree(int nsnap_guess, int nprog_guess,
                      std::size_t initial_tree_guess = 0,
                      std::size_t snap_step = 1,
                      std::size_t prog_step = 1,
                      std::size_t tree_step = 1)
  {
    Tree t{ BranchFactory{nsnap_guess, nprog_guess, snap_step, prog_step}, tree_step };
    if (initial_tree_guess) t.resize(initial_tree_guess);
    return t;
  }

  // Main Loop
  struct Settings;
  int main_loop(Tree& tree, Evoldum& evoldum, ::Settings& vh, Treelog& treelog);

} // namespace Makebr
