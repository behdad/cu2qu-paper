#!/usr/bin/env python3
# Usage: gdiff_larger.py <scratch> <tag>
# Where did blobs grow?  By scale, by how much, and a few examples.
import sys, re, collections
S, T = sys.argv[1], sys.argv[2]
base = open(f'{S}/out.base.{T}.txt').read().splitlines()
fix = open(f'{S}/out.fix.{T}.txt').read().splitlines()
pat = re.compile(r' (draw|paint)=(NULL|(\d+):[0-9a-f]+)')
by_scale = collections.defaultdict(lambda: [0, 0, 0, 0, 0])  # blobs, smaller, larger, bytes_b, bytes_f
grow = collections.Counter()
examples = []
for b, f in zip(base, fix):
    head = b.split(':', 1)[0]
    if ' cubic s' not in head:
        continue
    scale = head.split(' cubic s')[1].split(' ')[0]
    upem_rel = None
    for (tb, vb, lb), (tf, vf, lf) in zip(pat.findall(b), pat.findall(f)):
        if not (lb and lf):
            continue
        lb, lf = int(lb), int(lf)
        st = by_scale[scale]
        st[0] += 1; st[3] += lb; st[4] += lf
        if lf < lb: st[1] += 1
        if lf > lb:
            st[2] += 1
            d = lf - lb
            grow[d] += 1
            if len(examples) < 400:
                examples.append((d, lb, lf, head, tb))
print('growth in bytes -> count:', sorted(grow.items())[:12], '... max', max(grow) if grow else 0)
tot = collections.defaultdict(lambda: [0, 0, 0, 0, 0])
for scale, st in by_scale.items():
    k = 's30' if scale == '30' else 'other'
    for i in range(5): tot[k][i] += st[i]
for k, st in tot.items():
    print(f'{k}: blobs {st[0]}, smaller {st[1]}, larger {st[2]}, bytes {st[3]} -> {st[4]} ({(st[4]/st[3]-1)*100:+.1f}%)')
examples.sort(reverse=True)
for d, lb, lf, head, tb in examples[:8]:
    print(f'  +{d} bytes ({lb} -> {lf}) {tb} {head[-90:]}')
