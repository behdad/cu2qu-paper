#!/usr/bin/env python3
"""Find the inputs on which a kdiff binary dies (stack overflow), by feeding
it lines in chunks and bisecting the chunks that kill it.

Usage: find_crash.py <binary> <inputs> <first line> <last line> [max crashes]
Prints each crashing input decoded.
"""
import struct, subprocess, sys

binary, inputs = sys.argv[1], sys.argv[2]
first, last = int(sys.argv[3]), int(sys.argv[4])
max_crashes = int(sys.argv[5]) if len(sys.argv) > 5 else 5
lines = open(inputs).read().splitlines()[first - 1 : last]


def dies(chunk):
    r = subprocess.run([binary], input="\n".join(chunk) + "\n", capture_output=True, text=True)
    return r.returncode != 0


def unbits(s):
    return struct.unpack(">d", bytes.fromhex(s))[0]


found = []


def search(lo, hi):
    if len(found) >= max_crashes or not dies(lines[lo:hi]):
        return
    if hi - lo == 1:
        found.append(lo)
        return
    mid = (lo + hi) // 2
    search(lo, mid)
    search(mid, hi)


CHUNK = 2000
total = 0
for start in range(0, len(lines), CHUNK):
    search(start, min(start + CHUNK, len(lines)))
    if len(found) >= max_crashes:
        break
for i in found:
    f = lines[i].split()
    v = [unbits(s) for s in f[1:]]
    print(f"line {first + i}: accuracy {v[0]!r}  cubic {[(v[1], v[2]), (v[3], v[4]), (v[5], v[6]), (v[7], v[8])]!r}")
    dx = v[7] - 3 * v[5] + 3 * v[3] - v[1]
    dy = v[8] - 3 * v[6] + 3 * v[4] - v[2]
    print(f"   third difference ({dx!r}, {dy!r}); over accuracy: {(dx * dx + dy * dy) ** 0.5 / v[0]!r}")
print("crashing inputs found:", len(found), "(stopped at", max_crashes, ")" if len(found) >= max_crashes else ")")
