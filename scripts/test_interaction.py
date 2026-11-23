#!/usr/bin/env python3
"""Drive the live figures in a headless browser and check what they say.

Usage: test_interaction.py page.html

The page is served from its own directory's parent (so that ../style.css and
cu2qu.js resolve), with a script added that moves every slider, drags
control points and moves the pointer over the plot, and writes down what the
figures report. Needs google-chrome.
"""

import functools
import http.server
import os
import re
import subprocess
import sys
import threading

DRIVER = r"""
<script>
window.addEventListener("load", () => setTimeout(() => {
  const log = [];
  const fig = (n) => document.querySelector(`figure[data-fig="${n}"]`);
  const readout = (n) => fig(n).querySelector(".fig-readout").textContent;
  const slide = (input, v) => { input.value = v; input.dispatchEvent(new Event("input")); };
  const point = (type, x, y, target) => target.dispatchEvent(new PointerEvent(type, { clientX: x, clientY: y, pointerId: 1, bubbles: true }));
  const centre = (el) => { const r = el.getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; };
  try {
    log.push("live figures: " + document.querySelectorAll("figure.live").length);

    log.push("construction 0: " + readout("construction"));
    const tol = fig("construction").querySelector("input[type=range]");
    slide(tol, 0);
    log.push("construction tolerance least: " + readout("construction"));
    slide(tol, 100);
    log.push("construction tolerance most: " + readout("construction"));
    slide(tol, 47);
    let svg = fig("construction").querySelector("svg");
    let [x, y] = centre(svg.querySelector('[data-handle="1"]'));
    point("pointerdown", x, y, svg.querySelector('[data-handle="1"]'));
    point("pointermove", x + 120, y + 90, svg);
    point("pointerup", x + 120, y + 90, svg);
    log.push("construction after dragging c1: " + readout("construction"));
    point("pointermove", x + 300, y + 300, svg);
    log.push("construction after a move without a button: " + readout("construction"));

    slide(fig("error").querySelector("input[type=range]"), 12);
    log.push("error: " + [...fig("error").querySelectorAll("text")].map((t) => t.textContent).filter((t) => t.includes("Φ")).join("; "));

    svg = fig("count").querySelector("svg");
    const r = svg.getBoundingClientRect();
    const at = (vx, vy) => [r.left + (vx / 760) * r.width, r.top + (vy / 400) * r.height];
    // 44 + (log10(128.5) - 1) / 3 * 558 in the figure's own coordinates
    [x, y] = at(44 + ((Math.log10(128.5) - 1) / 3) * 558, 150);
    point("pointermove", x, y, svg);
    log.push("count at 128.5: " + readout("count"));
    [x, y] = at(44 + ((Math.log10(127.5) - 1) / 3) * 558, 150);
    point("pointermove", x, y, svg);
    log.push("count at 127.5: " + readout("count"));
    svg.dispatchEvent(new PointerEvent("pointerleave", {}));
    log.push("count after leaving: [" + readout("count") + "] hidden " + fig("count").querySelector(".fig-readout").hidden);

    const mu = fig("kink").querySelector("input[type=range]");
    const corners = () => [...fig("kink").querySelectorAll("text")].map((t) => t.textContent).filter((t) => /smooth|corner|way|master|halfway/.test(t)).join("; ");
    log.push("kink 0.5: " + corners());
    slide(mu, 0);
    log.push("kink 0: " + corners());
    slide(mu, 0.25);
    log.push("kink 0.25: " + corners());

    log.push("pace 0: " + readout("pace"));
    const [first, second] = fig("pace").querySelectorAll("input[type=range]");
    slide(first, 2 / 3);
    slide(second, 2 / 3);
    log.push("pace 2/3, 2/3: " + readout("pace"));
    svg = fig("pace").querySelector("svg");
    [x, y] = centre(svg.querySelector('[data-handle="1"]'));
    point("pointerdown", x, y, svg.querySelector('[data-handle="1"]'));
    point("pointermove", x - 60, y + 60, svg);
    point("pointerup", x - 60, y + 60, svg);
    log.push("pace after dragging c1 back: " + readout("pace") + " sliders " + first.value + " " + second.value);
  } catch (e) {
    log.push("ERROR " + e + "\n" + e.stack);
  }
  const pre = document.createElement("pre");
  pre.id = "testlog";
  pre.textContent = log.join("\n");
  document.body.append(pre);
}, 300));
</script>
"""


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def main():
    page = os.path.realpath(sys.argv[1])
    directory = os.path.dirname(page)
    root = os.path.dirname(directory)
    test = os.path.join(directory, ".cu2qu-interaction-test.html")
    source = open(page, encoding="utf-8").read()
    open(test, "w", encoding="utf-8").write(source.replace("</body>", DRIVER + "</body>"))
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=root))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    url = "http://127.0.0.1:%d/%s/%s" % (server.server_address[1], os.path.basename(directory), os.path.basename(test))
    try:
        dom = subprocess.run(
            ["google-chrome", "--headless=new", "--no-sandbox", "--disable-gpu", "--virtual-time-budget=12000", "--dump-dom", url],
            capture_output=True, text=True, timeout=180,
        ).stdout
    finally:
        os.unlink(test)
        server.shutdown()
    m = re.search(r'<pre id="testlog">(.*?)</pre>', dom, re.S)
    if not m:
        print("no log: the driver did not run")
        sys.exit(1)
    log = m.group(1).replace("&gt;", ">").replace("&lt;", "<").replace("&amp;", "&")
    print(log)
    if "ERROR" in log:
        sys.exit(1)


if __name__ == "__main__":
    main()
