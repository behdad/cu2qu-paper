#!/bin/sh
# Convergence check of the instrument on three adversarial sets: outer maximum on 256 points per quadratic instead of 64.
W=/home/behdad/hb/cu2qu-work/research-error
for name in cusp loop straight; do
  $W/gt $W/sets/$name.bin $W/sets/$name.gt_ramp.J256 8 0 0 256 512
done
