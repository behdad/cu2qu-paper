#!/usr/bin/env python3
"""Lay the page out at several widths in a headless browser and look for
things that do not fit.

Usage: check_layout.py page.html [--verbose]

At each width (three phones, a printed page's width, and a desktop window),
after MathJax has typeset the formulas:

- the page must not be wider than the window: anything that sticks out must
  be inside a box that scrolls sideways on its own;
- at the printed page's width and on the desktop no display formula, table,
  figure or listing may be wider than its box (paper does not scroll; the
  listings are let off here, since the page prints them smaller);
- no formula may have failed to typeset.

Needs google-chrome. The page is served from its directory's parent, as
test_interaction.py does. A headless window cannot be narrower than 500
pixels, so the page is laid out in frames of the widths wanted, side by side
in one wide window.
"""

import functools
import http.server
import json
import os
import re
import subprocess
import sys
import threading

# name, window width, and whether boxes that scroll sideways are an error
WIDTHS = [("phone, 320", 320, False), ("phone, 390", 390, False), ("phone, 480", 480, False), ("printed page, 700", 700, True), ("desktop, 1280", 1280, True)]

DRIVER = r"""
<script>
window.addEventListener("load", async () => {
  const out = { errors: [] };
  try {
    if (window.MathJax && MathJax.startup) await MathJax.startup.promise;
    await new Promise((r) => setTimeout(r, 300));
    const width = document.documentElement.clientWidth;
    out.window = width;
    out.page = document.documentElement.scrollWidth;
    const scrolls = (el) => {
      for (let p = el.parentElement; p && p !== document.body; p = p.parentElement) {
        const o = getComputedStyle(p).overflowX;
        if (o === "auto" || o === "scroll" || o === "hidden") return true;
      }
      return false;
    };
    const name = (el) => {
      const text = (el.textContent || "").replace(/\s+/g, " ").trim().slice(0, 60);
      return el.tagName.toLowerCase() + (el.className && el.className.baseVal === undefined ? "." + String(el.className).replace(/\s+/g, ".") : "") + ": " + text;
    };
    out.sticking = [];
    for (const el of document.querySelectorAll("main *")) {
      const r = el.getBoundingClientRect();
      if (r.width && r.right > width + 1 && !scrolls(el) && !el.closest("mjx-container > *")) out.sticking.push(name(el) + " [" + Math.round(r.width) + " wide, right edge " + Math.round(r.right) + "]");
    }
    out.sticking = out.sticking.slice(0, 20);
    out.scrolling = [];
    out.listings = [];
    for (const el of document.querySelectorAll("main div.eq, main .table-scroll, main pre, main figure")) {
      if (el.scrollWidth > el.clientWidth + 1) (el.tagName === "PRE" ? out.listings : out.scrolling).push(name(el) + " [" + el.scrollWidth + " in " + el.clientWidth + "]");
    }
    out.formulas = document.querySelectorAll("mjx-container").length;
    out.failed = document.querySelectorAll("mjx-merror").length;
    out.raw = (document.querySelector("main").textContent.match(/\\\(|\\\[/g) || []).length;
  } catch (e) {
    out.errors.push(String(e));
  }
  const pre = document.createElement("pre");
  pre.id = "layoutlog";
  pre.textContent = JSON.stringify(out);
  document.body.append(pre);
});
</script>
"""

WRAPPER = r"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>layout</title></head>
<body style="margin:0">
%(frames)s
<script>
const frames = [...document.querySelectorAll("iframe")];
const timer = setInterval(() => {
  const logs = frames.map((f) => { try { return f.contentDocument.getElementById("layoutlog"); } catch (e) { return null; } });
  if (logs.some((l) => !l)) return;
  clearInterval(timer);
  const pre = document.createElement("pre");
  pre.id = "layoutlogs";
  pre.textContent = JSON.stringify(logs.map((l) => JSON.parse(l.textContent)));
  document.body.append(pre);
}, 200);
</script>
</body></html>
"""


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def main():
    page = os.path.realpath(sys.argv[1])
    directory = os.path.dirname(page)
    root = os.path.dirname(directory)
    test = os.path.join(directory, ".cu2qu-layout-test.html")
    wrapper = os.path.join(directory, ".cu2qu-layout-frames.html")
    source = open(page, encoding="utf-8").read()
    open(test, "w", encoding="utf-8").write(source.replace("</body>", DRIVER + "</body>"))
    frames = "\n".join('<iframe src="%s" style="width:%dpx;height:900px;border:0"></iframe>' % (os.path.basename(test), width) for _, width, _ in WIDTHS)
    open(wrapper, "w", encoding="utf-8").write(WRAPPER % {"frames": frames})
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=root))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    url = "http://127.0.0.1:%d/%s/%s" % (server.server_address[1], os.path.basename(directory), os.path.basename(wrapper))
    try:
        dom = subprocess.run(
            ["google-chrome", "--headless=new", "--no-sandbox", "--disable-gpu", "--hide-scrollbars", "--window-size=4000,1000", "--virtual-time-budget=40000", "--dump-dom", url],
            capture_output=True, text=True, timeout=300,
        ).stdout
    finally:
        os.unlink(test)
        os.unlink(wrapper)
        server.shutdown()
    m = re.search(r'<pre id="layoutlogs">(.*?)</pre>', dom, re.S)
    if not m:
        print("the driver did not run")
        sys.exit(1)
    logs = json.loads(m.group(1).replace("&gt;", ">").replace("&lt;", "<").replace("&amp;", "&"))
    bad = False
    for (label, width, strict), out in zip(WIDTHS, logs):
        problems = list(out["errors"])
        if out["window"] != width:
            problems.append("the window is %d wide, not %d" % (out["window"], width))
        if out["page"] > out["window"] + 1:
            problems.append("the page is %d wide in a window of %d" % (out["page"], out["window"]))
        problems += ["sticks out: " + t for t in out["sticking"]]
        if strict:
            problems += ["does not fit its box: " + t for t in out["scrolling"]]
        if out["failed"] or out["raw"]:
            problems.append("%d formulas failed to typeset, %d left as TeX" % (out["failed"], out["raw"]))
        boxes = len(out["scrolling"]) + len(out["listings"])
        print("%-18s %d formulas; %d boxes scroll sideways; %s" % (label + ":", out["formulas"], boxes, "ok" if not problems else "%d problems" % len(problems)))
        for t in problems:
            print("    " + t)
        if "--verbose" in sys.argv:
            for t in (out["scrolling"] if not strict else []) + out["listings"]:
                print("    scrolls: " + t)
        bad = bad or bool(problems)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
