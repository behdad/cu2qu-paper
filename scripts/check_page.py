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
    for a, b in ((r"\(", r"\)"), (r"\[", r"\]")):
        if source.count(a) != source.count(b):
            problems.append("%s and %s are not paired: %d, %d" % (a, b, source.count(a), source.count(b)))
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
    assert abs(M.FIVE_BLOCK_FACTOR / (4 / 27) - 1.025) < 0.001
    need("constants", "0.2981", "3.08", "1.45", "1.26", "0.1519", "within 2.5%")

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
        ("HarfBuzz: one quadratic by distance (section 7)", "HarfBuzz: one by the bound of sec. 7, else the ramp"),
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
    hb = values["HarfBuzz: one quadratic by distance (section 7)"]
    need("best weights save", "0.6% of the quadratics at a tolerance of one unit, 2.6% at half a unit, and 10.3% at a tenth")
    got = ["%.1f" % (100 * (1 - b / c)) for b, c in zip(best, cu2qu)]
    if got != ["0.6", "2.6", "10.3"]:
        problems.append("best weights save %s" % got)
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
    if bound[1][1] != "24.4" or one[1][1] != "32.6" or one[0][1] != "45.8":
        problems.append("one quadratic by distance: %s %s" % (bound, one))
    need("HarfBuzz", "24.4% become one quadratic, where 2.4% did at equal parameters and 32.6% are in fact that close", "from 2.84 to 2.59", "45.8% are within a unit")
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
    need("ties", "3 cubics of 77,533 at a tolerance of one unit and to 7 at half a unit", "all but 3 cubics at a tolerance of one unit, 7 at half a unit and none at a tenth")

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

    # --- the figures' numbers, as cu2qu.js computes them
    figures = json.load(open(os.path.join(data, "figures.json")))
    pace = figures["pace"]
    if (round(pace["equal"]), "%.1f" % pace["apart"], "%.1f" % pace["bound"]) != (71, "0.7", "2.4"):
        problems.append("figure 5: %s" % pace)
    need("figure 5", "up to 71 apart, on a chord of 520; as curves they are 0.7 apart, and the bound says 2.4", "0.27 and 0.97")
    if abs(figures["kink"]["angles"][0]) > 1e-9 or round(abs(figures["kink"]["angles"][1])) != 16:
        problems.append("figure 4: %s" % figures["kink"])
    if figures["count"]["ends"] != [8, 12, 15] or figures["construction"]["n"] != 3:
        problems.append("figures 1 and 3: %s" % figures)
    if abs(figures["construction"]["formula"] - figures["construction"]["measured"]) > 1e-9:
        problems.append("figure 1: formula and measurement differ")

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
            except Exception as e:  # noqa: BLE001
                status = str(e)
            print("  %s  %s" % (status, href))
            # the DOI resolver's target turns scripts away; a browser gets through
            if status != 200 and not (href.startswith("https://doi.org/") and "403" in str(status)):
                if "--repository-pending" in sys.argv and href.startswith("https://github.com/behdad/cu2qu-paper"):
                    continue  # this repository, before it is on GitHub
                problems.append("link: %s: %s" % (href, status))

    print("%d formulas, %d internal links, %d external links, %d tables" % (source.count(r"\(") + source.count(r"\["), len(page.internal), len(set(page.external)), len(page.tables)))
    if problems:
        print("\n".join("PROBLEM: " + p for p in problems))
        sys.exit(1)
    print("the page agrees with the data")


if __name__ == "__main__":
    main()
