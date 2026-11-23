#!/usr/bin/env python3
"""Characterize the differences between two ft_equiv.py dumps.

Usage: ft_diff.py old.txt new.txt
For each differing line: its class, and how the two results differ (number
of points; exception; or same count with different coordinates).
"""
import ast, collections, sys

def parse(s):
    if s.startswith('['):
        return ast.literal_eval(s.replace('nan', 'None').replace('inf', '1e999'))
    return s

old = open(sys.argv[1]).read().splitlines()
new = open(sys.argv[2]).read().splitlines()
assert len(old) == len(new)
by_tag = collections.Counter()
kinds = collections.Counter()
total = collections.Counter()
examples = collections.defaultdict(list)
for a, b in zip(old, new):
    tag, ra = a.split(' ', 1)
    tag2, rb = b.split(' ', 1)
    assert tag == tag2
    total[tag] += 1
    if ra == rb:
        continue
    by_tag[tag] += 1
    try:
        pa, pb = parse(ra), parse(rb)
    except Exception:
        pa, pb = ra[:20], rb[:20]
    if isinstance(pa, str) or isinstance(pb, str):
        kind = f'{pa if isinstance(pa, str) else "spline"} -> {pb if isinstance(pb, str) else "spline"}'
    elif pa and isinstance(pa[0], list):
        # several curves: compare the point count of the first
        na, nb = len(pa[0]), len(pb[0])
        kind = 'multi: same count, other points' if na == nb else f'multi: {na - 2} -> {nb - 2} quadratics' if nb < na else f'multi: MORE {na - 2} -> {nb - 2}'
    else:
        na, nb = len(pa), len(pb)
        kind = 'same count, other points' if na == nb else ('one quadratic fewer' if nb == na - 1 else f'{na - 2} -> {nb - 2} quadratics')
    kinds[(tag, kind)] += 1
    if len(examples[(tag, kind)]) < 1:
        examples[(tag, kind)].append((ra[:110], rb[:110]))
print('class: differing / total')
for tag in sorted(total):
    print(f'  {tag:10s} {by_tag[tag]:6d} / {total[tag]}')
print('how they differ:')
for (tag, kind), n in sorted(kinds.items()):
    print(f'  {tag:10s} {kind:42s} {n}')
for key in sorted(examples):
    for ra, rb in examples[key]:
        print(key, '\n   old', ra, '\n   new', rb)
