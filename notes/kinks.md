# Kinks from the single quadratic: two measurements

Behdad's requirement (2026-11-23): "The equal parameters for fonts is more
than folklore. It's what avoids kinks. Condition for kink-free interpolation
is either equal angle or equal ratio of handles." And exact end tangents are
"nonnegotiable for fonts".

For n >= 2 cu2qu is linear in the control points and keeps both. For n = 1
(the intersection of the end tangents) it is not: the quadratic's handles are
the cubic's divided by fractions that differ from master to master.

## The data

The masters of 24 variable-font sources, extracted to arrays by
`research/kinks/` (stage 1; the arrays are 272 MB and are not in this
repository; on the machine they were made on they are in
`~/.cache/cu2qu-work/agent-kinks/data`): Cantarell, Noto Sans CJK and Noto
Serif CJK as CFF2 masters, and 21 Noto design spaces, of 2 to 16 masters.
3,589,652 cubic segments; converted together at 0.001 em, 447,956 become
single quadratics (12.5%), 428,843 of them in Noto Serif CJK.

## Second measurement (in the note): `scripts/kinks.py`, `data/kinks.txt`

An independent implementation, written from a description of what to measure
and not from the first one's code.

- Smooth joint: a side is a cubic, and the two directions at the joint are
  exactly parallel and point the same way in every master.
- Instances: for every pair of masters, a quarter, a half and three quarters
  of the way.
- Change by conversion: the largest difference, over the instances, between
  the kink of the quadratic instance and that of the cubic instance. Change by
  rounding: the same between rounded and unrounded quadratics.

| | |
|---|---|
| smooth joints | 1,056,932 |
| next to at least one single quadratic | 194,219 |
| of those, not horizontal or vertical in every master | 4,939 |
| changed by conversion by more than 0.1 / 0.5 / 1 degree | 616 / 70 / 26, largest 4.43 |
| changed by rounding by more than 0.1 / 0.5 / 1 degree | 2,175 / 1,932 / 1,850, largest 16.3 |
| smooth joints not next to a single quadratic | 862,713 |
| largest change by conversion there | 3.3e-12 degrees, leaving out 69 joints |

Things to know about these numbers:

- All 616 are among the 4,939 joints that are not horizontal or vertical in
  every master. At a joint that is horizontal (or vertical) in every master
  the ratio of the handles cannot make a kink.
- The two CJK sources hold 189,599 of the 194,219 joints and 613 of the 616.
- The 69 joints: in an instance, the handles of two masters cancel and leave
  a handle of no length. The direction then comes from a farther control
  point, which the conversion does not keep for any n. They come from pairs
  of masters that are not neighbours.
- One joint (Noto Serif CJK `glyph36167`) has its tangent turned around by
  the single quadratic in a master itself: a hook one unit long, smaller than
  the tolerance, whose tangents meet behind a control point.
- With signed angles the counts are 619 / 73 / 26 and 2,175 / 1,949 / 1,869.
- Taking joints that are smooth to within one degree instead of exactly:
  224,065 next to a single quadratic; 4,951 / 547 / 33 changed by conversion
  (largest still 4.43), 28,885 / 22,525 / 18,960 by rounding.
- Worst joint next to a single quadratic: Noto Sans CJK `uni5FDB`, between
  the two weight masters at three quarters: 0 degrees as cubics, 4.43 as
  quadratics (3.59 rounded). Both sides became single quadratics.

## First measurement: `research/kinks/`

645 / 77 / 27 joints changed by conversion, of 194,218; rounding 3,982 /
3,387 / 3,062 over all smooth joints. The differences from the second are
all in the definitions, and the second reproduces the first exactly when
given them: it counted collinear line-to-line joints as smooth (24,852 more),
sampled tenths between neighbouring masters (and midpoints on several axes)
instead of quarters between all pairs, and used signed angles. Its other
findings, not measured again:

- Forbidding single quadratics for masters would cost 2.4% more points (8%
  for Noto Naskh Arabic, about 25% for masters that came from TrueType).
- tan(kink) = (R - 1) sin(theta) / (2 sqrt(R) + (1 + R) cos(theta)) for two
  masters halfway, R the ratio of their ratios, theta the angle between their
  tangents.
- A single quadratic passes the equal-parameter test only with both of its
  first differences e1, e2 within 2 sqrt(3) tolerance (9/4 for a right angle
  between the end tangents), which is what keeps the fractions near 2/3 and
  the kinks small unless the handles are a few units long.

## What was decided

Leave it. The exception is real and small beside rounding. The note says so
in section 6 and in "Loose ends".
