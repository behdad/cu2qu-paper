#!/bin/bash
# Regenerate every table of the report (outputs next to the scripts, *.out).
# The long steps are not rerun here:
#   python3 corpus_run.py          closed-form schemes on the corpus (about a minute)
#   python3 yardstick_run.py       numerical optimum on the corpus (about 20 minutes on 14 cores)
#   python3 chd.py                 instrument self-test (about 4 minutes; chd_selftest.out)
cd "$(dirname "$0")"
VENV=/home/behdad/hb/cu2qu-work/research-weights/venv/bin/python
PY="python3 -W ignore -u"
$PY common.py                > common.out 2>&1
$PY q1_lp.py                 > q1_lp.out 2>&1
$PY q1_proofs.py             > q1_proofs.out 2>&1
$VENV -W ignore -u q1_odd_exact.py > q1_odd_exact.out 2>&1
$PY q1_relax.py              > q1_relax.out 2>&1
$PY q1_blocks.py             > q1_blocks.out 2>&1
$PY q1_angle_budget.py       > q1_angle_budget.out 2>&1
$PY q2_third_order_check.py  > q2_third_order_check.out 2>&1
$PY q2_fourth_order.py       > q2_fourth_order.out 2>&1
$PY q2_model_check.py        > q2_model_check.out 2>&1
$PY q2_balance_weights.py    > q2_balance_weights.out 2>&1
$PY q2_structure.py          > q2_structure.out 2>&1
$PY q2_density_check.py      > q2_density_check.out 2>&1
$PY q2_qm_affine.py          > q2_qm_affine.out 2>&1
$PY q3_slide_ep.py           > q3_slide_ep.out 2>&1
$PY q3_t_equal_param.py      > q3_t_equal_param.out 2>&1
$PY q3_slide_corpus.py       > q3_slide_corpus.out 2>&1     # builds a table of LPs on first run (a few minutes, 12 processes)
$PY canon.py                 > canon.out 2>&1
$PY tquad.py                 > tquad.out 2>&1
$PY gest.py                  > gest.out 2>&1
$PY balanced.py              > balanced.out 2>&1
$PY corpus_stats.py          > corpus_stats.out 2>&1
$PY yardstick_stats.py       > yardstick_stats.out 2>&1
$PY bt_vs_yard.py            > bt_vs_yard.out 2>&1
$PY opt_check.py             > opt_check.out 2>&1
$PY geo2.py                  > geo2.out 2>&1
SAFETY=1.0 $PY geo2.py       > geo2_safety1.0.out 2>&1
$PY q5_metrics.py            > q5_metrics.out 2>&1
echo done
