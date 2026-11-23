#!/usr/bin/env python3
"""Compare two result files of ft_run.py / the kurbo harness, line by line,
and say how they differ.

Usage: cmp_out.py <a> <b> [--classes general|adversarial|ties] [--show N]
       [--sample-of <inputs_adversarial.txt>]

Kinds of difference:
  no end -> N quadratics   one side did not end (R: RecursionError, CRASH)
  count                    different numbers of quadratics (a decision)
  no result                one side found no approximation
  value                    same numbers of points, some coordinate differs:
                           the largest relative difference is reported
  zero sign                same values, a zero of the other sign
"""
import collections, struct, sys

GENERAL = [(120000, "random float"), (240000, "random int"), (280000, "third difference at a limit"),
           (340000, "nearly quadratic"), (380000, "groups"), (535066, "font"), (555066, "font pairs")]
ADVERSARIAL = [(42000, "three equal points"), (84000, "three equal points, integers"), (126000, "one handle exact"),
               (168000, "farthest on an axis"), (189000, "flat in the middle"), (210000, "S"), (252000, "nearly quadratic")]


def unbits(s):
    return struct.unpack(">d", bytes.fromhex(s))[0]


def parse(line):
    line = line.strip()
    if line in ("-", "R", "CRASH"):
        return line
    return [[unbits(s) for s in part.split()] for part in line.split(" | ")]


def main():
    args = sys.argv[1:]
    a, b = args[0], args[1]
    classes = None
    show = 0
    stride = 1
    if "--classes" in args:
        classes = {"general": GENERAL, "adversarial": ADVERSARIAL, "ties": [(10**9, "single quadratic at its limit")]}[args[args.index("--classes") + 1]]
    if "--show" in args:
        show = int(args[args.index("--show") + 1])
    if "--sampled" in args:
        stride = 10  # every tenth cubic (7 lines each) of the adversarial set

    def cls(i):
        if classes is None:
            return "all"
        if stride != 1:
            i = (i // 7) * 7 * stride + i % 7
        for end, name in classes:
            if i < end:
                return name
        return "?"

    kinds = collections.defaultdict(collections.Counter)
    totals = collections.Counter()
    worst = collections.defaultdict(float)
    shown = collections.Counter()
    n = 0
    for i, (la, lb) in enumerate(zip(open(a), open(b))):
        n += 1
        c = cls(i)
        totals[c] += 1
        if la == lb:
            continue
        pa, pb = parse(la), parse(lb)
        if pa in ("R", "CRASH") or pb in ("R", "CRASH"):
            other = pb if pa in ("R", "CRASH") else pa
            if other in ("R", "CRASH"):
                continue  # both did not end
            kind = "no end -> %s" % ("no result" if other == "-" else "%d quadratic(s)" % (len(other[0]) // 2 - 2))
        elif pa == "-" or pb == "-":
            kind = "no result on one side"
        elif [len(s) for s in pa] != [len(s) for s in pb]:
            kind = "count %d -> %d" % (len(pa[0]) // 2 - 2, len(pb[0]) // 2 - 2)
        else:
            d = max(abs(x - y) for sa, sb in zip(pa, pb) for x, y in zip(sa, sb))
            m = max(max(abs(x) for s in pa for x in s), 1e-300)
            rel = d / m
            worst[c] = max(worst[c], rel)
            kind = "zero sign" if d == 0 else "value"
        kinds[c][kind] += 1
        if shown[(c, kind)] < show:
            shown[(c, kind)] += 1
            print(f"   line {i + 1} [{c}] {kind}")
    print(f"{a.split('/')[-1]} vs {b.split('/')[-1]}: {n} lines")
    any_diff = False
    for c in totals:
        k = kinds[c]
        if k:
            any_diff = True
        rest = ", ".join(f"{name}: {v}" for name, v in sorted(k.items()))
        extra = f"  (largest relative difference of a coordinate {worst[c]:.2e})" if worst[c] else ""
        print(f"  {c:32s} {totals[c]:7d}  identical {totals[c] - sum(k.values()):7d}  {rest}{extra}")
    if not any_diff:
        print("  identical")


main()
