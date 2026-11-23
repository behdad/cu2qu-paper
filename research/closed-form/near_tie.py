#!/usr/bin/env python3
# Usage: near_tie.py <path to a COPY of hb-gpu-cu2qu.hh with the trial loop>
# Instruments the trial loop: when the n that the trials settle on differs
# from the n of the closed form, print how close the closed form's n was to
# the tolerance.  For a throwaway export, not the worktree.
import sys
p = sys.argv[1]
s = open(p).read()
a = '''    if (!hb_gpu_cubic_approx_spline (cubic, n, tolerance, off, on))
      continue;
'''
b = '''    if (!hb_gpu_cubic_approx_spline (cubic, n, tolerance, off, on))
      continue;
    {
      auto Phi = [] (unsigned m) {
	if (m == 2) return 4. / 27.;
	double gg = (m - 2) / (4. * (m - 1));
	double first = 4. / 27. * (1 + gg) * (1 + gg) * (1 + gg);
	double w2 = (.25 - .5 / (m - 1)) / 3;
	double second = (m - 3) / (4. * (m - 1)) + 2 * w2 * sqrt (w2);
	return first > second ? first : second;
      };
      unsigned m = 2;
      while (m <= HB_GPU_CU2QU_MAX_N && !(d3 * Phi (m) <= tolerance * m * m * m)) m++;
      static unsigned long long total, differ;
      total++;
      if (m != n)
      {
	differ++;
	fprintf (stderr, "trial n %u, formula n %u: error at formula n is %.9f of the tolerance (%llu of %llu differ)\\n",
		 n, m, d3 * Phi (m) / (tolerance * m * m * m), differ, total);
      }
    }
'''
assert s.count(a) == 1
s = s.replace(a, b)
if '#include <stdio.h>' not in s:
    s = s.replace('#include <cmath>\n', '#include <cmath>\n#include <stdio.h>\n', 1)
open(p, 'w').write(s)
