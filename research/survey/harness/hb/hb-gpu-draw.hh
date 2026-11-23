// SHIM, not HarfBuzz code.  Stands in for src/hb-gpu-draw.hh so that an
// unmodified copy of src/hb-gpu-cu2qu.hh can be compiled alone: the few
// things that header uses, with the accumulator calls recorded.
#ifndef SHIM_HB_GPU_DRAW_HH
#define SHIM_HB_GPU_DRAW_HH

#include <cstdint>
#include <vector>

#define unlikely(x) (__builtin_expect (!!(x), 0))
#define likely(x) (__builtin_expect (!!(x), 1))

template <typename T> static inline T hb_max (T a, T b) { return a > b ? a : b; }

#define HB_BUDGET_1 1u

static inline bool
hb_budget_spend (int64_t &budget, unsigned int cost, unsigned int mult = 1)
{
  budget -= (int64_t) cost * mult;
  return budget >= 0;
}

struct shim_seg_t
{
  uint32_t kind; /* 1 line, 2 quad */
  double x0, y0, cx, cy, x1, y1;
};

struct hb_gpu_draw_t
{
  bool success = true;
  int64_t budget = INT64_MAX / 2;
  double cur_x = 0, cur_y = 0;
  std::vector<shim_seg_t> segs;

  int64_t *get_budget_remaining () { return &budget; }

  void acc_line_to (double x, double y)
  {
    segs.push_back ({1u, cur_x, cur_y, (cur_x + x) * .5, (cur_y + y) * .5, x, y});
    cur_x = x; cur_y = y;
  }
  void acc_conic_to (double cx, double cy, double x, double y)
  {
    segs.push_back ({2u, cur_x, cur_y, cx, cy, x, y});
    cur_x = x; cur_y = y;
  }
};

#endif
