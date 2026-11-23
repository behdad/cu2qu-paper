#!/usr/bin/env python3
"""Check the page of the note against the scripts' output.

Usage: check_page.py page.html DATA_DIR [--links [--repository-pending]]

Checks that the page is well formed (tags closed, ids unique, links within
the page have a target, formula delimiters paired), that its tables hold the
numbers the formulas and the measurements give, and that the numbers quoted
in the text are the ones in DATA_DIR (the output of measure.py,
compare_revs.py, perf_fonttools.py and render_figures.js). With --links it
also fetches every external link; --repository-pending lets the links to this
repository on GitHub fail, for before it is there.
"""

import html.parser
import hashlib
import json
import math
import os
import re
import sys
import urllib.request

import cu2qu_math as M

VOID = {"meta", "link", "br", "hr", "img", "input", "path", "circle", "rect", "source", "col", "wbr"}


class Page(html.parser.HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack = []
        self.ids = []
        self.internal = []
        self.external = []
        self.errors = []
        self.text = []
        self.tables = {}
        self.table = None
        self.row = None
        self.cell = None
        self.skip = 0

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            self.ids.append(attrs["id"])
        href = attrs.get("href")
        if tag == "a" and href:
            (self.internal if href.startswith("#") else self.external).append(href)
        if tag in ("script", "style", "title", "desc"):
            self.skip += 1
        if tag == "table":
            self.table = self.tables.setdefault(attrs.get("id", "table%d" % len(self.tables)), [])
        elif tag == "tr" and self.table is not None:
            self.row = []
            self.table.append(self.row)
        elif tag in ("td", "th") and self.row is not None:
            self.cell = []
        if tag not in VOID:
            self.stack.append(tag)

    def handle_startendtag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            self.ids.append(attrs["id"])

    def handle_endtag(self, tag):
        if tag in VOID:
            return
        if tag in ("script", "style", "title", "desc"):
            self.skip -= 1
        if tag in ("td", "th") and self.cell is not None:
            self.row.append(" ".join("".join(self.cell).split()))
            self.cell = None
        elif tag == "table":
            self.table = self.row = None
        if not self.stack or self.stack[-1] != tag:
            self.errors.append("closing </%s> at line %d, open: %s" % (tag, self.getpos()[0], self.stack[-3:]))
            if tag in self.stack:
                while self.stack and self.stack.pop() != tag:
                    pass
        else:
            self.stack.pop()

    def handle_data(self, data):
        if self.skip:
            return
        self.text.append(data)
        if self.cell is not None:
            self.cell.append(data)


def fmt(x, digits):
    """A number as the page writes it: thousands separators, no trailing zeros
    for whole numbers."""
    if abs(x - round(x)) < 1e-9:
        return "{:,}".format(int(round(x)))
    return "{:,.{d}f}".format(x, d=digits)


def main():
    path, data = sys.argv[1], sys.argv[2]
    source = open(path, encoding="utf-8").read()
    page = Page()
    page.feed(source)
    page.close()
    problems = list(page.errors)
    if page.stack:
        problems.append("not closed: %s" % page.stack)
    seen = set()
    for i in page.ids:
        if i in seen:
            problems.append("id used twice: %s" % i)
        seen.add(i)
    for href in page.internal:
        if href[1:] not in seen:
            problems.append("no target for %s" % href)
    body = source[source.index("<main>"):source.index("</main>")]
    for a, b in ((r"\(", r"\)"), (r"\[", r"\]")):
        if body.count(a) != body.count(b):
            problems.append("%s and %s are not paired: %d, %d" % (a, b, body.count(a), body.count(b)))
    for name in ("construction", "error", "count", "kink", "pace"):
        if not re.search(r"<!-- figure:%s -->\n<svg" % name, source):
            problems.append("figure %s is missing" % name)
    text = " ".join("".join(page.text).split())

    def need(what, *strings):
        for s in strings:
            if s not in text:
                problems.append("%s: the page does not say %r" % (what, s))

    # --- Table 1, from the formulas
    rows = page.tables["table-limits"][1:]
    for row in rows:
        n = int(row[0])
        want = [str(n), "%.5f" % M.error_factor(n), fmt(M.limit(n), 2), fmt(n**3 / M.best_error_factor(n), 2)]
        if row != want:
            problems.append("Table 1, n = %d: page %s, computed %s" % (n, row, want))
    need("limits", "\\(36/\\sqrt3 \\approx 20.8\\)", "20.78")
    assert abs(36 / math.sqrt(3) - 20.78) < 0.005
    assert abs(0.25 + math.sqrt(3) / 36 - 0.2981) < 5e-5
    assert abs((4 / 27) / (math.sqrt(3) / 36) - 3.08) < 0.005 and abs(((4 / 27) / (math.sqrt(3) / 36)) ** (1 / 3) - 1.45) < 0.005
    assert abs(((0.25 + math.sqrt(3) / 36) / (4 / 27)) ** (1 / 3) - 1.26) < 0.005
    assert ["%.4f" % M.best_error_factor(n) for n in (5, 7, 9)] == ["0.1519", "0.1494", "0.1482"]
    need("constants", "0.2981", "3.08", "1.45", "1.26", "giving 0.1519", "0.1494 for seven and 0.1482 for nine")

    # --- the measurements on the three fonts
    measure = open(os.path.join(data, "measure.txt")).read()

    def line(prefix):
        for ln in measure.split("\n"):
            if ln.strip().startswith(prefix):
                return ln[len(ln) - len(ln.lstrip()) + len(prefix):]
        raise KeyError(prefix)

    def floats(prefix):
        return [float(x) for x in re.findall(r"-?\d+\.\d+", line(prefix))]

    def counts(prefix):
        return re.findall(r"(\d+) \((\d+\.\d)%\)", line(prefix))

    names = (
        ("cu2qu", "cu2qu: one if it fits at equal parameters, else the ramp"),
        ("the same with the best weights (section 4)", "the same with the best weights (sec. 4)"),
        ("HarfBuzz’s converter: one quadratic by distance (section 7)", "HarfBuzz: one by the bound of sec. 7, else the ramp"),
        ("mid-point spline, end tangents not kept", "mid-point (kurbo to_quads, lyon)"),
        ("cu2qu’s construction, accepted by measured distance", "cu2qu's constructions, least n within tolerance"),
    )
    table = {row[0]: row[1:] for row in page.tables["table-corpus"][1:]}
    values = {}
    for label, prefix in names:
        values[label] = floats(prefix)
        want = ["%.3f" % v for v in values[label]]
        if table.get(label) != want:
            problems.append("Table 2, %s: page %s, measured %s" % (label, table.get(label), want))
    cu2qu = values["cu2qu"]
    best = values["the same with the best weights (section 4)"]
    least = values["cu2qu’s construction, accepted by measured distance"]
    hb = values["HarfBuzz’s converter: one quadratic by distance (section 7)"]
    got = ["%.1f" % (100 * (1 - b / c)) for b, c in zip(best, cu2qu)]
    need("best weights save", "%s%% of the quadratics at a tolerance of one unit, %s%% at half a unit, and %s%% at a tenth" % tuple(got))
    by_distance = floats("HarfBuzz, its splines accepted by distance")
    saved = ["%.0f" % (100 * (1 - a / b)) for a, b in zip(by_distance, hb)]
    need("a renderer could save", "would take %s%%, %s%% and %s%% off the quadratics that HarfBuzz’s converter makes now" % tuple(saved))
    need("by distance", "1.70, 2.11 and 3.69 quadratics per cubic where cu2qu takes 2.34, 2.84 and 4.80")
    if ["%.2f" % v for v in least] != ["1.70", "2.11", "3.69"] or ["%.2f" % v for v in cu2qu] != ["2.34", "2.84", "4.80"]:
        problems.append("by distance: measured %s %s" % (least, cu2qu))
    quarter = [1 - a / b for a, b in zip(least, cu2qu)]
    if not all(0.22 < q < 0.28 for q in quarter):
        problems.append("'about a quarter fewer': %s" % quarter)
    fill = floats("cu2qu's output: distance / tolerance, mean")
    if ["%.2f" % v for v in fill] != ["0.27", "0.34", "0.45"]:
        problems.append("distance over tolerance: %s" % fill)
    need("fill", "0.27, 0.34 and 0.45 of the tolerance")
    single = counts("fits at equal parameters (cu2qu's test)")
    ruled = counts("ruled out by |d3| > 36/sqrt(3) x tolerance")
    bound = counts("within tolerance by the bound on the distance of sec. 7")
    one = counts("one quadratic within tolerance")
    if [s[1] for s in single] != ["6.3", "2.4", "0.3"]:
        problems.append("singles: %s" % single)
    if [str(round(float(r[1]))) for r in ruled] != ["61", "80", "97"]:
        problems.append("ruled out: %s" % ruled)
    need("singles", "6.3%, 2.4% and 0.3% of all cubics", "61% of the cubics at a tolerance of one unit, 80% at half a unit and 97% at a tenth")
    hb_one = counts("HarfBuzz's quadratic within tolerance")
    if bound[1][1] != "24.4" or hb_one[1][1] != "33.3" or one[0][1] != "45.8":
        problems.append("one quadratic by distance: %s %s %s" % (bound, one, hb_one))
    need("HarfBuzz", "24.4% become one quadratic, where 2.4% did at equal parameters; 33.3% are in fact that close to this quadratic", "go from 2.84 to 2.59", "45.8% are within a unit")
    share = floats("share of that saving in cubics that take one, percent")
    if ["%.0f" % v for v in share] != ["69", "52", "22"]:
        problems.append("share of the saving in single quadratics: %s" % share)
    need("share of the saving", "69% of that saving is in cubics for which one quadratic would do", "At half a unit it is 52% of the saving, and at a tenth 22%")
    lines = counts("cubics  ")
    with_lines = floats("segments per cubic with those as one line")
    if lines[1][1] != "2.1" or "%.3f" % with_lines[1] != "2.526" or "%.3f" % hb[1] != "2.585":
        problems.append("HarfBuzz's lines: %s %s" % (lines, with_lines))
    need("HarfBuzz's lines", "That is 2.1% of these cubics at half a unit, and takes its 2.585 down to 2.526 segments per cubic")
    zero = re.search(r"a control point on its end point: \d+ \((\d+\.\d)\d%\)", measure)
    if zero.group(1) != "1.1":
        problems.append("handles of no length: %s" % zero.group(0))
    need("handles of no length", "as 1.1% of the cubics above have")
    if "%.2f" % hb[1] != "2.59":
        problems.append("HarfBuzz at 0.5: %s" % hb)
    m = re.search(r"bound / measured median (\d+\.\d+)", measure)
    if m.group(1) != "1.10":
        problems.append("bound over measured: %s" % m.group(1))
    need("bound", "its median is 1.10 times it")
    if "below it by at most" not in measure or float(re.search(r"below it by at most (\S+) units", measure).group(1)) > 1e-9:
        problems.append("the bound is below the measured distance somewhere")
    share = floats("tolerance 1.0:")[:5]
    if ["%.1f" % v for v in share] != ["6.3", "65.9", "17.6", "8.7", "1.6"] or "(most: 8)" not in line("tolerance 1.0:"):
        problems.append("shares: %s" % line("tolerance 1.0:"))
    need("shares", "6.3% of the cubics become one quadratic, 65.9% two, 17.6% three, 8.7% four, and 1.6% more, up to eight")
    pct = floats("|d3| percentiles 50 / 90 / 99 / max:")
    if pct != [27.7, 131.6, 336.0, 1574.0]:
        problems.append("percentiles: %s" % pct)
    need("percentiles", "shorter than 27.7 units, and the longest is 1,574", "half of them are below 27.7, 90% below 131.6 and 99% below 336")
    angle = line("mid-point spline: angle")  # heading only; the numbers are on the next lines
    m = re.search(r"tolerance 1\.0: more than 0\.1 / 1 / 5 degrees on \d+ \([\d.]+%\) / \d+ \((\d+)\.\d%\)", measure)
    if m.group(1) != "57":
        problems.append("mid-point tangents: %s" % m.group(0))
    need("mid-point", "more than a degree for 57% of the cubics")
    need("corpus", "77,533 cubics")
    if "corpus: 77533 cubics" not in measure:
        problems.append("corpus size")

    # --- fontTools before and after
    compare = open(os.path.join(data, "compare-fonttools.txt")).read()
    differ = re.findall(r"results that differ: (\d+) of 77533", compare)
    if differ != ["3", "7", "0"] or compare.count("|d3| / tolerance = 128: 4 -> 3 quadratics") != 10:
        problems.append("fontTools before and after: %s" % differ)
    if re.findall(r"not the closed form's: before (\d+), after (\d+)", compare) != [("3", "0"), ("7", "0"), ("0", "0")]:
        problems.append("fontTools against the closed form")
    need("ties", "it rejects 3 cubics of 77,533 at a tolerance of one unit and 7 at half a unit", "all but 3 cubics at a tolerance of one unit, 7 at half a unit and none at a tenth")
    # the closed form of the single quadratic's error decides as the test does
    if re.findall(r"single quadratics after: (\d+); cubics that the closed form of its error decides otherwise: (\d+)", compare) != [(s_[0], "0") for s_ in single]:
        problems.append("the single quadratic's closed form against the test")
    need("single quadratic, closed form", "On the fonts above it decides every cubic as the recursive test does")

    # what the old test did at a limit
    ties = {}
    for ln in open(os.path.join(data, "ties.txt")):
        m = re.match(r"(\S+) x tolerance\s+all\s+(\d+) \|\s+(\d+)\s+(\d+)\s+(\d+) \|", ln)
        if m:
            ties[m.group(1)] = tuple(int(x) for x in m.group(2, 3, 4, 5))
    if ties != {"54": (3000, 3000, 0, 0), "128": (3000, 855, 2145, 0), "843.75": (4200, 1849, 2290, 61)}:
        problems.append("ties: %s" % ties)
    need("ties at a limit", "it accepted all 3,000 at 54 times the tolerance, rejected 2,145 of 3,000 at 128 times, and at 843.75 times rejected 2,290 of 4,200 and did not return for 61")

    # the same comparison on groups of masters
    masters = open(os.path.join(data, "compare-fonttools-masters.txt")).read()
    m = re.search(r"^total\s+(\d+)\s+(\d+)\s+4 -> 3: (\d+)$", masters, re.M)
    differing = [ln.split()[0] for ln in masters.split("\n")[3:] if len(ln.split()) > 3 and ln.split()[0] != "total"]
    if not m or m.group(1, 2, 3) != ("3589652", "142", "142") or differing != ["notosanscjk", "notoserifcjk"]:
        problems.append("masters before and after: %s %s" % (m and m.groups(), differing))
    need("masters before and after", "the 3,589,652 groups of compatible masters", "for all but 142, all in the two CJK fonts and each from four quadratics to three")

    # the inputs that do not return
    endless = open(os.path.join(data, "endless.txt")).read()
    blocks = endless.split("curve_to_quadratic(")[1:]
    if len(blocks) != 6 or not all("(compiled)" in b and "killed by SIGSEGV" in b and b.count("RecursionError") == 1 for b in blocks):
        problems.append("endless.txt: fontTools before")
    if not all(re.search(r"after \S+ \(pure Python\)\s+(returns \d+ points|ApproxNotFoundError)", b) for b in blocks):
        problems.append("endless.txt: fontTools after")
    kurbo = open(os.path.join(data, "endless-kurbo.txt")).read()
    outcome = {}
    for b in kurbo.split("CubicBez(")[1:]:
        head, before, after = b.split("\n")[:3]
        outcome[head.split(")).approx_spline(")[0] + ")@" + head.split("approx_spline(")[1].split(")")[0]] = ("stack overflow" in before, "returns" in after)
    want = {
        "(0, 0), (0, 0), (0, 0), (25, 0)@100/27": (False, True),
        "(-45, -3), (-30, -2), (0, 5), (0, 9)@20/3": (False, True),
        "(-258, -394), (-148, -194), (-155, -13), (105, 661)@5": (False, True),
        "(-22, 289), (31, 324), (-137, 476), (1499, 3445)@4": (False, True),
        "(386, -158), (527, -28), (480, -99), (2270, 2329)@4": (True, True),
        "(-20, 7), (-20, -27.5), (-12, -32.5), (4, -32)@3.5555555555555554": (True, True),
        "(-31, 97), (-51, 93), (-91, 76), (-139, 52)@1.987615979999813": (True, True),
        "(0, 0), (nan, 0), (1, 1), (2, 2)@1": (True, True),
        "(0, 0), (0, 100), (100, 200), (354, 300)@nan": (True, True),
    }
    if outcome != want:
        problems.append("endless-kurbo.txt: %s" % outcome)
    assert 32 / 9 == 3.5555555555555554
    need(
        "inputs that do not return",
        "curve_to_quadratic([(0, 0), (0, 0), (0, 0), (25, 0)], 100 / 27)",
        "curve_to_quadratic([(-45, -3), (-30, -2), (0, 5), (0, 9)], 20 / 3)",
        "curve_to_quadratic([(-258, -394), (-148, -194), (-155, -13), (105, 661)], 5)",
        "curve_to_quadratic([(-22, 289), (31, 324), (-137, 476), (1499, 3445)], 4)",
        "whose third differences are 128 × 5 and 843.75 × 4",
        "returns from these four and overflows its stack on others",
        "(386, −158), (527, −28), (480, −99), (2270, 2329) at an accuracy of 4, again 843.75 × 4",
        "(−20, 7), (−20, −27.5), (−12, −32.5), (4, −32) at 32/9",
    )
    for c, r in (([(-258, -394), (-148, -194), (-155, -13), (105, 661)], 128 * 5), ([(-22, 289), (31, 324), (-137, 476), (1499, 3445)], 843.75 * 4), ([(386, -158), (527, -28), (480, -99), (2270, 2329)], 843.75 * 4)):
        if abs(M.third_difference([complex(*q) for q in c])) != r:
            problems.append("third difference of %s" % c)

    # kurbo: the same bits as fontTools, and its instructions
    kurbo = open(os.path.join(data, "compare-kurbo.txt")).read()
    if kurbo.count("splines with a different number of points: 0; with the same number and a different bit: 0") != 3 or "kurbo 3 quadratics, fontTools 4" not in kurbo:
        problems.append("compare-kurbo.txt")
    need("kurbo and fontTools", "gives fontTools’ coordinates bit for bit on the fonts above", "Within rounding of a limit the two can still differ by a quadratic")
    perf_kurbo = dict(((m.group(1), m.group(2)), int(m.group(3))) for m in re.finditer(r"^(before|after)\s+(\S+)\s+(\d+)\s", open(os.path.join(data, "perf-kurbo.txt")).read(), re.M))
    need(
        "kurbo's instructions",
        "go from {:,} to {:,} at a tolerance of one unit, and from {:,} to {:,} at half a unit".format(
            perf_kurbo["before", "1"], perf_kurbo["after", "1"], perf_kurbo["before", "0.5"], perf_kurbo["after", "0.5"]
        ),
    )

    # the fonts
    fonts = open(os.path.join(data, "cubics.fonts.txt")).read()
    for name, version, glyphs in (("SourceCodePro-Regular", "2.042", "1568 of 1568"), ("STIXTwoMath-Regular", "2.13 ", "3000 of"), ("NotoSansSC-Regular", "2.004", "3000 of")):
        block = fonts.split(name)[1].split("/usr/share")[0]
        if "Version " + version not in block or glyphs not in block or "units per em 1000" not in block:
            problems.append("fonts: %s" % name)
    need("fonts", "all of Source Code Pro Regular 2.042, and the first 3,000 glyphs of STIX Two Math 2.13 and of Noto Sans SC Regular 2.004")

    # --- the joints of section 6
    kinks = open(os.path.join(data, "kinks.txt")).read()

    def total(heading, text=kinks):
        part = text.split(heading)[1]
        return re.search(r"^  total\s+(.*)$", part, re.M).group(1).split()

    if total("1. Sources and conversion")[-4:] != ["9877525", "3589652", "447956", "12.48%"] or kinks.count("\n  noto") + kinks.count("\n  cantarell") < 24:
        problems.append("kinks: sources")
    sources = re.findall(r"^  (\S+)\s+(\d+)\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+[\d.]+%$", kinks, re.M)
    if len(sources) != 24 or min(int(m_) for _, m_ in sources) != 2 or max(int(m_) for _, m_ in sources) != 16:
        problems.append("kinks: %d sources" % len(sources))
    joints = total("2. Joints")
    if joints != ["9877525", "5491830", "1056932", "7599", "194219", "18.38%", "4939", "1391"]:
        problems.append("kinks: joints %s" % joints)
    cjk = [int(re.search(r"^  %s\s.*\s(\d+)\s+\d+$" % name, kinks.split("2. Joints")[1].split("3. Smooth")[0], re.M).group(1)) for name in ("notosanscjk", "notoserifcjk")]
    if sum(cjk) < 0.99 * 4939:
        problems.append("kinks: joints off the axes in the CJK fonts: %s" % cjk)
    beside = total("3. Smooth joints next to a single quadratic")
    if beside != ["194219", "616", "70", "26", "4.4304", "2175", "1932", "1850", "16.3139"]:
        problems.append("kinks: beside a single quadratic %s" % beside)
    other = total("4. Smooth joints not next to a single quadratic")
    if other[:4] != ["862713", "180.0000", "3.3e-12", "69"]:
        problems.append("kinks: other joints %s" % other)
    if "At 50 of them no master's cubic has a control point at the joint" not in kinks or "At 19 every master's cubic has one there" not in kinks or "At 0 some masters have one" not in kinks:
        problems.append("kinks: the 69 joints")
    if "joints with a change by conversion above 0.1, 0.5, 1 degree: 29, 29, 28 over all instances;\n    0, 0, 0 over regular instances." not in kinks:
        problems.append("kinks: changes at the 69 joints")
    in_master = re.findall(r"^  \S+\s+\S+\s+\d+\s+2/2\s+\d+\s+0\s+([\d.]+)\s+0\s+no$", kinks, re.M)
    if in_master != ["26.5651", "11.3099"]:
        problems.append("kinks: changes in a master itself: %s" % in_master)
    within = open(os.path.join(data, "kinks-within-one-degree.txt")).read()
    if "of the 34785 that are not H/V (see 2): conversion 4951, 547, 33 " not in within or "rounding 28885, 22525, 18960 " not in within:
        problems.append("kinks: joints smooth to within a degree")
    need(
        "kinks",
        "the masters of 24 variable fonts",
        "the sources of 21 design spaces of other Noto families, of 2 to 16 masters each",
        "They have 3.59 million cubics, of which 12.5% become single quadratics",
        "Of 1,056,932 such joints, 194,219 have a single quadratic on at least one side. All but 4,939",
        "a tenth of a degree at 616, by more than half a degree at 70 and by more than a degree at 26; the largest change is 4.4 degrees",
        "at 2,175, 1,932 and 1,850 of the same joints",
        "34,785 that are not horizontal or vertical have a single quadratic beside them; the conversion changes 4,951, 547 and 33 of them by more than those amounts, and rounding 28,885, 22,525 and 18,960",
        "At the 862,713 exactly smooth joints",
        "except at 69 where an instance has a handle of no length: at 50 because the handles of two masters cancel there, at 19 because no master has one",
        "at 29 of the 69 the angle changes by more than a tenth of a degree",
        "by 27 and 11 degrees",
        "One joint of the 194,219 of section 6 has it",
    )
    if "above 0.1, 0.5, 1 degree at 1, 1, 1 of these joints, largest 180.0000" not in kinks:
        problems.append("kinks: the tangent turned around")
    assert abs(2 * math.sqrt(3) - 3.5) < 0.05
    need("turned around", "\\(2\\sqrt3 \\approx 3.5\\) times the tolerance")

    # --- constants of the related work
    assert abs(M.error_factor(2) / 8 - 1 / 54) < 1e-15
    need("related work", "it is 1/54 of the length of the third difference", "\\(\\lvert\\mathbf{d}_3\\rvert/32\\)", "432")

    perf = open(os.path.join(data, "perf-fonttools.txt")).read()
    got = {}
    for ln in perf.split("\n"):
        m = re.match(r"(before|after), (compiled|pure Python)\s+(instr\.|microsec\.)\s+(\S+)\s+(\S+)\s+(\S+)", ln)
        if m:
            got[m.group(1), m.group(2), m.group(3)] = m.group(4, 5, 6)
    rows = page.tables["table-speed"][1:]
    for i, row in enumerate(rows):
        want = []
        for kind in ("compiled", "pure Python"):
            for when in ("before", "after"):
                want.append("%s (%s)" % ("{:,}".format(int(got[when, kind, "instr."][i])), got[when, kind, "microsec."][i]))
        if row[1:] != want:
            problems.append("Table 3, row %d: page %s, measured %s" % (i + 1, row[1:], want))

    # --- the figures' numbers, as cu2qu-figures.js computes them
    figures = json.load(open(os.path.join(data, "figures.json")))
    pace = figures["pace"]
    if (round(pace["equal"]), "%.1f" % pace["apart"], "%.1f" % pace["bound"]) != (71, "0.7", "2.4"):
        problems.append("figure 5: %s" % pace)
    need("figure 5", "up to 71 apart, on a chord of 520; as curves they are 0.7 apart, and the bound says 2.4", "0.27 and 0.97")
    if abs(figures["kink"]["angles"][0]) > 1e-9 or round(abs(figures["kink"]["angles"][1])) != 16:
        problems.append("figure 4: %s" % figures["kink"])
    if [v for _, v in figures["count"]["percentiles"]] != pct[:3]:
        problems.append("figure 3: the marks are at %s, the percentiles are %s" % (figures["count"]["percentiles"], pct))
    if figures["count"]["ends"] != [8, 12, 15] or figures["construction"]["n"] != 3:
        problems.append("figures 1 and 3: %s" % figures)
    if abs(figures["construction"]["formula"] - figures["construction"]["measured"]) > 1e-9:
        problems.append("figure 1: formula and measurement differ")

    # the figure script's address changes when the script does
    script = open(os.path.join(os.path.dirname(os.path.abspath(path)), "cu2qu-figures.js"), "rb").read()
    blob = hashlib.sha1(b"blob %d\0" % len(script) + script).hexdigest()[:10]
    if 'src="cu2qu-figures.js?v=%s"' % blob not in source:
        problems.append("the figure script is linked with another version than %s" % blob)

    if "--links" in sys.argv:
        for href in sorted(set(page.external)):
            if not href.startswith("http"):
                target = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(path)), href.split("?")[0].split("#")[0]))
                if not os.path.exists(target):
                    problems.append("link to a file that is not there: %s" % href)
                continue
            try:
                request = urllib.request.Request(href, headers={"User-Agent": "Mozilla/5.0 (link check)"})
                with urllib.request.urlopen(request, timeout=30) as response:
                    status = response.status
                    answer = response.read(4000)
                # a page can answer 200 and say that it is not there
                if status == 200 and (len(answer) < 400 or b"cannot be displayed" in answer):
                    status = "200, with %d bytes that are not the page" % len(answer)
            except Exception as e:  # noqa: BLE001
                status = str(e)
            print("  %s  %s" % (status, href))
            # the DOI resolver's target turns scripts away; a browser gets through
            if status != 200 and not (href.startswith("https://doi.org/") and "403" in str(status)):
                if "--repository-pending" in sys.argv and href.startswith("https://github.com/behdad/cu2qu-paper"):
                    continue  # this repository, before it is on GitHub
                problems.append("link: %s: %s" % (href, status))

    print("%d formulas, %d internal links, %d external links, %d tables" % (body.count(r"\(") + body.count(r"\["), len(page.internal), len(set(page.external)), len(page.tables)))
    if problems:
        print("\n".join("PROBLEM: " + p for p in problems))
        sys.exit(1)
    print("the page agrees with the data")


if __name__ == "__main__":
    main()
