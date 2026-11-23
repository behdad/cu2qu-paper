// Runs HarfBuzz's hb_gpu_cubic_to_quadratics () from an unmodified copy of
// src/hb-gpu-cu2qu.hh (copied by build_hb.sh), against the shim in hb/.
//
//   hb_c2q <corpus.bin> <tol> <out.bin> [time-reps]

#include "hb-gpu-cu2qu.hh"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>

int main (int argc, char **argv)
{
  if (argc < 4) return 2;
  FILE *f = fopen (argv[1], "rb");
  if (!f) { perror ("corpus"); return 1; }
  char magic[8];
  uint32_t count;
  if (fread (magic, 1, 8, f) != 8 || memcmp (magic, "C2QCOR01", 8)) return 1;
  if (fread (&count, 4, 1, f) != 1) return 1;
  std::vector<double> pts ((size_t) count * 8);
  if (fread (pts.data (), 8, pts.size (), f) != pts.size ()) return 1;
  fclose (f);

  double tol = atof (argv[2]);
  int reps = argc > 4 ? atoi (argv[4]) : 0;

  FILE *o = fopen (argv[3], "wb");
  if (!o) { perror ("out"); return 1; }
  fwrite ("C2QRES01", 1, 8, o);
  fwrite (&count, 4, 1, o);
  hb_gpu_draw_t g;
  for (uint32_t i = 0; i < count; i++)
  {
    const double *c = &pts[(size_t) i * 8];
    g.segs.clear ();
    g.cur_x = c[0]; g.cur_y = c[1];
    hb_gpu_cubic_to_quadratics (&g, {c[0], c[1]}, {c[2], c[3]}, {c[4], c[5]}, {c[6], c[7]}, tol);
    uint32_t n = g.segs.size ();
    fwrite (&n, 4, 1, o);
    for (auto &s : g.segs)
    {
      fwrite (&s.kind, 4, 1, o);
      double v[6] = {s.x0, s.y0, s.cx, s.cy, s.x1, s.y1};
      fwrite (v, 8, 6, o);
    }
  }
  fclose (o);

  if (reps > 0)
  {
    double best = 1e30;
    size_t nq = 0;
    double check = 0;
    for (int r = 0; r < reps; r++)
    {
      auto t0 = std::chrono::steady_clock::now ();
      nq = 0;
      for (uint32_t i = 0; i < count; i++)
      {
	const double *c = &pts[(size_t) i * 8];
	g.segs.clear ();
	g.cur_x = c[0]; g.cur_y = c[1];
	hb_gpu_cubic_to_quadratics (&g, {c[0], c[1]}, {c[2], c[3]}, {c[4], c[5]}, {c[6], c[7]}, tol);
	nq += g.segs.size ();
	check += g.segs.back ().cx;
      }
      auto t1 = std::chrono::steady_clock::now ();
      double dt = std::chrono::duration<double> (t1 - t0).count ();
      if (dt < best) best = dt;
    }
    printf ("{\"method\":\"harfbuzz_hb_gpu_cu2qu\",\"tol\":%g,\"ns_per_cubic\":%.1f,\"quads\":%zu,\"check\":%g}\n",
	    tol, best * 1e9 / count, nq, check);
  }
  return 0;
}
