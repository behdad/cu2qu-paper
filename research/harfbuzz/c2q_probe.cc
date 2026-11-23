/* Budget units and blob bytes for a few single cubics through the public
 * draw funcs, to compare two builds. */

#include <hb.h>
#include <hb-gpu.h>

#include <stdint.h>
#include <stdio.h>

struct cubic_t { const char *name; float p[8]; };

int
main ()
{
  static const cubic_t cubics[] = {
    {"quarter circle r=100",          {100, 0,  100, 55.23f,  55.23f, 100,  0, 100}},
    {"quarter circle r=1000",         {1000, 0, 1000, 552.3f, 552.3f, 1000, 0, 1000}},
    {"last control on end",           {0, 0,  0, 100,  100, 100,  100, 100}},
    {"first control on start",        {0, 0,  0, 0,    0, 100,    100, 100}},
    {"both controls on their ends",   {0, 0,  0, 0,    100, 100,  100, 100}},
    {"s-curve",                       {0, 0,  100, 0,  0, 100,    100, 100}},
    {"half circle (parallel tangents)", {0, 0, 0, 133, 200, 133,  200, 0}},
    {"loop back to start",            {0, 0,  200, 100, -200, 100, 0, 0}},
    {"exact quadratic",               {0, 0,  0, 66.6666667f, 33.3333333f, 100, 100, 100}},
  };
  hb_gpu_draw_t *draw = hb_gpu_draw_create_or_fail ();
  hb_draw_funcs_t *funcs = hb_gpu_draw_get_funcs (draw);
  for (const cubic_t &c : cubics)
  {
    hb_gpu_draw_reset (draw);
    hb_draw_set_budget (funcs, draw, 1 << 20);
    hb_draw_state_t st = HB_DRAW_STATE_DEFAULT;
    hb_draw_move_to (funcs, draw, &st, c.p[0], c.p[1]);
    int64_t before = hb_draw_get_budget_remaining (funcs, draw);
    hb_draw_cubic_to (funcs, draw, &st, c.p[2], c.p[3], c.p[4], c.p[5], c.p[6], c.p[7]);
    int64_t spent = before - hb_draw_get_budget_remaining (funcs, draw);
    hb_draw_close_path (funcs, draw, &st);
    hb_blob_t *blob = hb_gpu_draw_encode (draw, nullptr);
    printf ("%-34s units %4lld  blob %5u bytes\n", c.name, (long long) spent,
	    blob ? hb_blob_get_length (blob) : 0);
    hb_blob_destroy (blob);
  }
  hb_gpu_draw_destroy (draw);
  return 0;
}
