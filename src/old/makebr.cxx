// makebr.cxx
#include "types.h"
#include "config.h"     // Settings의 "정의"가 필요하니까 여기서 include
#include "logger.h"
// #include "utilities.h"      // util_tfcname 사용 시
// #include "run_make_branch.h"// run_make_branch_loop 사용 시

#include <algorithm>
#include <string>

namespace Makebr {

// Settings가 어떤 필드를 갖는지에 맞춰 수정하세요.
// 예: snap_i/snap_f/snap_step을 쓰는 경우
static inline int snap_start(const Settings& vh) { return vh.snap_i; }
static inline int snap_end  (const Settings& vh) { return vh.snap_f; }
static inline int snap_step (const Settings& vh) { return vh.snap_step; }

// 만약 Settings가 std::array<int,3> makebr_snap라면 위 3개를 이렇게 바꾸세요:
// static inline int snap_start(const Settings& vh){ return vh.makebr_snap[0]; }
// static inline int snap_end  (const Settings& vh){ return vh.makebr_snap[1]; }
// static inline int snap_step (const Settings& vh){ return vh.makebr_snap[2]; }

int main_loop(Tree& tree, Evoldum& evoldum, Settings& vh, Treelog& treelog)
{
  if (vh.treedir != "des" && vh.treedir != "prg") {
    LOG() << "Wrong Tree direction given (treedir = des or prg)";
    return 1;
  }
  if (snap_step(vh) == 0) {
    LOG() << "snap_step cannot be 0";
    return 1;
  }

  LOG() << "Program start";
  LOG() << "dir_tree=" << vh.dir_tree
        << " treedir=" << vh.treedir
        << " snap=[" << snap_start(vh) << "," << snap_end(vh)
        << "] step=" << snap_step(vh);

  // 필요 시:
  // util_tfcname(vh.dir_tree);

  int rc = 0;
  // rc = run_make_branch_loop(tree, evoldum, vh, treelog);

  LOG() << "Program end (rc=" << rc << ")";
  return rc;
}

} // namespace Makebr
