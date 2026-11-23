#!/usr/bin/env python3
"""Drive the live figures in a headless browser and check what they say.

Usage: test_interaction.py page.html

The page is served from its own directory's parent (so that ../style.css and
cu2qu-figures.js resolve), with a script added that moves every slider, drags
control points, presses the reset buttons and moves the pointer over the plot,
and writes down what the figures report. It fails if a reset does not bring a
figure back, if a drag starts or goes on when it should not, if a control
point leaves its drawing, or if Figure 3 shows a number that does not go with
its counts. The events are made up by the script; test_touch.js sends real
touches. Needs google-chrome.
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
  const point = (type, x, y, target, more) => target.dispatchEvent(new PointerEvent(type, Object.assign({ clientX: x, clientY: y, pointerId: 1, bubbles: true }, more || {})));
  const centre = (el) => { const r = el.getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; };
  const reset = (n) => fig(n).querySelector(".fig-controls button").click();
  const same = (now, before) => (now === before ? "as at first" : "DIFFERS: " + now);
  // the number a readout of figure 3 shows must give the count it shows
  const consistent = (text) => {
    const m = text.match(/= ([\d.]+): cu2qu (\d+), best weights (\d+), mid-point (\d+)/);
    if (!m) return " INCONSISTENT (no numbers)";
    const F = Cu2quFigures.math;
    const r = Number(m[1]);
    return F.count(r, 1, 1000) === Number(m[2]) && F.countBest(r) === Number(m[3]) && F.countMidpoint(r) === Number(m[4]) ? "" : " INCONSISTENT";
  };
  try {
    log.push("live figures: " + document.querySelectorAll("figure.live").length);

    const first = readout("construction");
    log.push("construction 0: " + first);
    const tol = fig("construction").querySelector("input[type=range]");
    log.push("construction slider starts at " + tol.value);
    slide(tol, 0);
    log.push("construction tolerance least: " + readout("construction"));
    slide(tol, 100);
    log.push("construction tolerance most: " + readout("construction"));
    slide(tol, 50);
    log.push("construction with the slider back: " + same(readout("construction"), first));
    let svg = fig("construction").querySelector("svg");
    let [x, y] = centre(svg.querySelector('[data-handle="1"]'));
    point("pointerdown", x, y, svg.querySelector('[data-handle="1"]'), { button: 2, buttons: 2 });
    point("pointermove", x + 50, y + 50, svg, { buttons: 2 });
    point("pointerup", x + 50, y + 50, svg);
    log.push("construction after the other button: " + same(readout("construction"), first));
    point("pointerdown", x, y, svg.querySelector('[data-handle="1"]'));
    point("pointermove", x + 120, y + 90, svg);
    point("pointerup", x + 120, y + 90, svg);
    const dragged = readout("construction");
    log.push("construction after dragging c1: " + dragged);
    point("pointermove", x + 300, y + 300, svg);
    log.push("construction after a move without a button: " + (readout("construction") === dragged ? "unchanged" : "DIFFERS: " + readout("construction")));
    // a release that never arrives: a mouse that moves with no button down lets go
    [x, y] = centre(svg.querySelector('[data-handle="2"]'));
    point("pointerdown", x, y, svg.querySelector('[data-handle="2"]'), { pointerType: "mouse", buttons: 1 });
    point("pointermove", x - 40, y + 40, svg, { pointerType: "mouse", buttons: 1 });
    const held = readout("construction");
    point("pointermove", x - 200, y + 200, svg, { pointerType: "mouse", buttons: 0 });
    point("pointermove", x - 220, y + 220, svg, { pointerType: "mouse", buttons: 0 });
    log.push("construction after a lost release: " + (held !== dragged && readout("construction") === held ? "let go" : "STILL DRAGGING"));
    reset("construction");
    log.push("construction after reset: " + same(readout("construction"), first) + ", slider " + tol.value);

    const phis = () => [...fig("error").querySelectorAll("text")].map((t) => t.textContent).filter((t) => t.includes("Φ")).join("; ");
    const phis0 = phis();
    slide(fig("error").querySelector("input[type=range]"), 12);
    log.push("error: " + phis());
    reset("error");
    log.push("error after reset: " + same(phis(), phis0));

    log.push("count at first: " + readout("count"));
    const cs = fig("count").querySelector("input[type=range]");
    slide(cs, 369);
    log.push("count, slider at 369: " + readout("count") + consistent(readout("count")));
    slide(cs, 370);
    log.push("count, slider at 370: " + readout("count") + consistent(readout("count")));
    svg = fig("count").querySelector("svg");
    const r = svg.getBoundingClientRect();
    const at = (vx, vy) => [r.left + (vx / 760) * r.width, r.top + (vy / 400) * r.height];
    // 44 + (log10(128.5) - 1) / 3 * 558 in the figure's own coordinates
    [x, y] = at(44 + ((Math.log10(128.5) - 1) / 3) * 558, 150);
    point("pointermove", x, y, svg);
    log.push("count at 128.5: " + readout("count") + consistent(readout("count")) + ", slider " + cs.value);
    [x, y] = at(44 + ((Math.log10(127.5) - 1) / 3) * 558, 150);
    point("pointermove", x, y, svg);
    log.push("count at 127.5: " + readout("count") + consistent(readout("count")));
    // a hair above a limit, by a press as a finger makes: the number must be above it too
    [x, y] = at(44 + ((Math.log10(128.02) - 1) / 3) * 558, 150);
    point("pointerdown", x, y, svg);
    log.push("count a hair above 128: " + readout("count") + consistent(readout("count")));
    const before = readout("count");
    svg.dispatchEvent(new PointerEvent("pointerleave", {}));
    log.push("count after leaving: " + (readout("count") === before ? "stays" : "DIFFERS: " + readout("count")));

    const mu = fig("kink").querySelector("input[type=range]");
    const corners = () => [...fig("kink").querySelectorAll("text")].map((t) => t.textContent).filter((t) => /smooth|corner|way|master|halfway/.test(t)).join("; ");
    const corners0 = corners();
    log.push("kink 0.5: " + corners0);
    slide(mu, 0);
    log.push("kink 0: " + corners());
    slide(mu, 0.25);
    log.push("kink 0.25: " + corners());
    reset("kink");
    log.push("kink after reset: " + same(corners(), corners0));

    const pace0 = readout("pace");
    log.push("pace 0: " + pace0);
    const [l1, l2] = fig("pace").querySelectorAll("input[type=range]");
    log.push("pace sliders go to " + l1.max + " and " + l2.max);
    slide(l1, 2 / 3);
    slide(l2, 2 / 3);
    log.push("pace 2/3, 2/3: " + readout("pace"));
    slide(l1, 0);
    slide(l2, 1.33);
    log.push("pace 0, 1.33: " + readout("pace"));
    slide(l1, 1.33);
    svg = fig("pace").querySelector("svg");
    const box = svg.viewBox.baseVal;
    const inside = [...svg.querySelectorAll("[data-handle]")].every((h) => {
      const cx = Number(h.getAttribute("cx")), cy = Number(h.getAttribute("cy"));
      return cx > 10 && cx < box.width - 10 && cy > 10 && cy < box.height - 10;
    });
    log.push("pace 1.33, 1.33: " + readout("pace") + " Control points " + (inside ? "inside" : "OUTSIDE") + " the drawing.");
    slide(l1, 2 / 3);
    slide(l2, 2 / 3);
    [x, y] = centre(svg.querySelector('[data-handle="1"]'));
    point("pointerdown", x, y, svg.querySelector('[data-handle="1"]'));
    point("pointermove", x - 60, y + 60, svg);
    point("pointerup", x - 60, y + 60, svg);
    log.push("pace after dragging c1 back: " + readout("pace") + " sliders " + l1.value + " " + l2.value);
    reset("pace");
    log.push("pace after reset: " + same(readout("pace"), pace0) + ", sliders " + l1.value + " " + l2.value);
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
    bad = [word for word in ("ERROR", "DIFFERS", "STILL DRAGGING", "INCONSISTENT", "OUTSIDE") if word in log]
    if bad or "live figures: 5" not in log:
        print("FAILED: %s" % (", ".join(bad) or "not every figure is live"))
        sys.exit(1)


if __name__ == "__main__":
    main()
