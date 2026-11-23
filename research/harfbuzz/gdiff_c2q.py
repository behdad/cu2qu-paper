#!/usr/bin/env python3
# Usage: gdiff_c2q.py <scratch> <tag>
# Compares two gpu_dump outputs where one side converts cubics differently:
# faces without CFF/CFF2 must be identical; for the others, report how many
# glyphs change, whether any changes between encoded and not, and the bytes.
import sys, re
S, T = sys.argv[1], sys.argv[2]
base = open(f'{S}/out.base.{T}.txt').read().splitlines()
fix = open(f'{S}/out.fix.{T}.txt').read().splitlines()
assert len(base) == len(fix), (len(base), len(fix))
pat = re.compile(r' (draw|paint)=(NULL|(\d+):[0-9a-f]+)')
st = {'quad': [0, 0], 'cubic': [0, 0]}
to_null = from_null = 0
bytes_b = bytes_f = 0
larger = smaller = 0
fonts_cubic = set(); fonts_changed = set(); fonts_quad = set()
bad = []
synth = [0, 0]
for b, f in zip(base, fix):
    if b.startswith('synth'):
        synth[0] += 1
        if b != f: synth[1] += 1
        continue
    kind = 'cubic' if ' cubic s' in b.split(':', 1)[0] else 'quad'
    st[kind][0] += 1
    font = b.split('#', 1)[0]
    if kind == 'cubic': fonts_cubic.add(font)
    else: fonts_quad.add(font)
    if b != f:
        st[kind][1] += 1
        if kind == 'quad': bad.append((b, f))
        else: fonts_changed.add(font)
    if kind == 'cubic':
        mb, mf = pat.findall(b), pat.findall(f)
        for (tb, vb, lb), (tf, vf, lf) in zip(mb, mf):
            if vb != 'NULL' and vf == 'NULL': to_null += 1
            if vb == 'NULL' and vf != 'NULL': from_null += 1
            if lb and lf:
                bytes_b += int(lb); bytes_f += int(lf)
                if int(lf) > int(lb): larger += 1
                if int(lf) < int(lb): smaller += 1
print(f"faces without cubics: {len(fonts_quad)} fonts, {st['quad'][0]} glyph lines, {st['quad'][1]} differ")
print(f"faces with CFF/CFF2:  {len(fonts_cubic)} fonts, {st['cubic'][0]} glyph lines, {st['cubic'][1]} differ, in {len(fonts_changed)} fonts")
print(f"  blobs that became NULL: {to_null}; NULL that became blobs: {from_null}")
if bytes_b:
    print(f"  bytes of blobs encoded on both sides: {bytes_b} -> {bytes_f} ({(bytes_f / bytes_b - 1) * 100:+.1f}%); smaller {smaller}, larger {larger}")
if synth[0]:
    print(f"synthetic outlines: {synth[0]}, {synth[1]} differ")
for b, f in bad[:4]:
    print('<', b[:200]); print('>', f[:200])
