// src/main.cpp
#include <mpi.h>
#include <omp.h>
#include <iostream>
//#include "allvar.h"
#include "logger.h"
#include "utilities.h"
#include "config.h"
#include "autovec.h"
#include "types.h"


//#include "types.hpp"
//#include "veluga_makebr.hpp"
//#include "veluga_utilites.hpp"
//#include "veluga_ctree.hpp"

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  //-----
  // Program Start
  //-----

  // Logo print
  util_printlogo();

  //-----
  // Initial check
  //-----
  if (argc < 2) {
    LOG() << "config file is required";
    return 1;
  }

  // TODO 123123
  // Basic spec info will be printed here
  LOG()  << "Program start";
  

  //------
  // Read Configuration and Make a settings
  //-----
  Settings vh;
  if (!load_config(argv[1], vh)) {
    LOG() << "Failed to load config. Using defaults where possible.\n";
    return 1;
  }
  if (!set_config(vh)){
    LOG() << "Wrong Tree direction given (treedir = des or prg)";
    return 1;
  }

  //-----
  // TreeFrog output modification
  //-----
  util_tfcname(vh.dir_tree);
  

  //-----
  // Make Branch
  //-----

  // Allocate Tree Array
  auto tree = Makebr::make_tree(vh.sarr_size, vh.barr_size, vh.tarr_size);
  auto evoldum = Makebr::make_evoldum(vh.tarr_size, vh.tarr_size);
  Makebr::Treelog treelog;


  // Main Loop
  int rcheck = Makebr::main_loop(tree, evoldum, vh, treelog);



  //std::size_t i = 150000;        // 크기를 초과해도 접근 시 자동으로 확장
  //evoldum[i].ID    = 123456;
  //evoldum[i].snap  = 117;
  //evoldum[i].merit = 0.83;
//  // TODO: 인자/설정 파일에서 세팅 로드
//  s.dir_catalog = "/path/to/catalog/";
//  s.dir_tree    = s.dir_catalog + (s.horg=='g' ? "Galaxy/tree" : "Halo/tree");
//  s.makebr_treedir = TreeDir::Descendants; // or Progenitors
//  s.snap_start = 70; s.snap_end = 200; s.snap_step = 1; // 예시
//  s.num_thread = omp_get_max_threads();
//
//  omp_set_num_threads(s.num_thread);
//
//  // 1) TreeFrog 파일명 정규화 (IDL: veluga_makebr_tfcname)
//  io::normalize_treefrog_filenames(s);
//
//  // 2) TreeSet 초기화
//  auto ts = makebr::init_treeset(s);
//
//  // 3) 메모리 할당
//  const i32 max_ngal = 100000;
//  auto tree = makebr::allocate_branches(s, ts, max_ngal);
//  std::vector<std::shared_ptr<Branch>> complete(max_ngal);
//  BranchWork wk; wk.id.assign(max_ngal, -1); wk.snap.assign(max_ngal, -1); wk.merit.assign(max_ngal, -1.0);
//
//  i32 gind = -1, n_comp = 0;
//  i64 max_id_seen = -1;
//
//  // 4) 스냅샷 루프
//  for (i32 snap = ts.N0; (ts.DN>0 ? snap<=ts.N1 : snap>=ts.N1); snap += ts.DN) {
//
//    // 빈 스냅샷 건너뛰기 (IDL의 파일 체크 자리)
//    // if (!io::file_exists(...)) continue;
//
//    const i32 snap_curr = snap;
//    i32 snap_next = makebr::find_next_snap(s, snap_curr);
//
//    // 마지막 스냅샷 처리
//    if (snap == ts.N1 || snap_next < 0) {
//      makebr::process_last_snapshot(s, tree, complete, n_comp, wk, gind, snap_curr);
//      // 남은 브랜치들 마무리(자리)
//      break;
//    }
//
//    // 데이터 로드 (TreeFrog HDF5 & 갤/헤일로)
//    // const auto tf = io::read_treefrog_h5(...);  // TODO
//    // const auto g_curr = io::read_galaxies(snap_curr, s.horg, s);
//    // const auto g_next = io::read_galaxies(snap_next, s.horg, s);
//
//    // max_id_seen = std::max<i64>(max_id_seen, ...);
//
//    // 매칭 단계 (OpenMP 병렬화 지점 예시)
//    // makebr::match_step(s, tree, complete, n_comp, ...);
//
//    // 가지 정리
//    // makebr::prune_stale(s, wk, tree, gind, complete, n_comp, snap_curr, tf.nlink);
//
//    // 로그/저장 체크포인트 (자리)
//    if (snap % 5 == 0) {
//      // TODO: 중간 저장
//    }
//
//    // 재할당 필요시
//    if ((i32)tree.size() - gind < (i32)(0.2 * max_ngal)) {
//      makebr::ensure_capacity(tree, wk, gind, max_ngal);
//    }
//    if ((i32)complete.size() - n_comp < (i32)(0.2 * max_ngal)) {
//      makebr::ensure_capacity_complete(complete, n_comp, max_ngal);
//    }
//
//    if (MPI::COMM_WORLD.Get_rank() == 0) {
//      std::cout << snap << " / " << ts.N1 << " processed\n";
//    }
//  }
//
//  // 5) 키 생성 + 검증 + 최종 저장(자리)
//  // auto key = makebr::build_key(s, complete, max_id_seen);
//  // TODO: 검증 & 저장
//
  MPI_Finalize();
  return 0;
}
