"""Which of the real master groups come out differently from fontTools
upstream main and from the branch, and how (see oldnew.sh)."""

import sys

d = sys.argv[1]
fewer = more = same = 0
single_changed = 0
for i, (a, b) in enumerate(zip(open(d + "/py-main.txt"), open(d + "/py.txt"))):
    if a == b:
        continue
    na = [len(s.split()) // 2 - 2 for s in a.strip().split(" | ")]
    nb = [len(s.split()) // 2 - 2 for s in b.strip().split(" | ")]
    if nb[0] < na[0]:
        fewer += 1
    elif nb[0] > na[0]:
        more += 1
    else:
        same += 1
    if na[0] == 1 or nb[0] == 1:
        single_changed += 1
    print("group %d: quadratics per curve %s -> %s" % (i, na, nb))
print(
    "branch has fewer quadratics in %d groups, more in %d, the same number in %d; a single quadratic is involved in %d"
    % (fewer, more, same, single_changed)
)
