#!/usr/bin/env node
// Drag the figures' control points with real touches and a real mouse, in a
// headless Chrome, and see that the figure follows and the page stays put.
//
// Usage: node test_touch.js SITE [PAGE]      (PAGE defaults to doc/cu2qu.html)
//
// test_interaction.py makes its events up in the page, which says nothing
// about what a browser does with a finger: scroll the page, and cancel the
// pointer. Here Chrome emulates a phone 390 pixels wide and is sent touches
// through its debugging port:
//
// - a finger on a control point drags it: the curve changes, the page and
//   the figure's own sideways scroll do not move, no pointer is cancelled;
// - a finger anywhere else on a figure still scrolls the page;
// - a tap on the plot of Figure 3 reads it.
//
// Then, as a desktop, the same drag with the mouse, and a press of the other
// button, which must not drag.
//
// Needs node 22 or later (for its WebSocket) and google-chrome. Temporary
// files go to the system's temporary directory (TMPDIR).

const { spawn } = require("child_process");
const fs = require("fs");
const http = require("http");
const os = require("os");
const path = require("path");

const [site, page = "doc/cu2qu.html"] = process.argv.slice(2);
if (!site) {
  console.error("usage: test_touch.js SITE [PAGE]");
  process.exit(2);
}

const TYPES = { ".html": "text/html; charset=utf-8", ".js": "text/javascript", ".css": "text/css", ".ico": "image/x-icon", ".png": "image/png", ".svg": "image/svg+xml" };

function serve(root) {
  const server = http.createServer((request, response) => {
    const file = path.join(root, decodeURIComponent(request.url.split("?")[0]));
    if (!file.startsWith(path.resolve(root)) || !fs.existsSync(file) || fs.statSync(file).isDirectory()) {
      response.writeHead(404);
      response.end();
      return;
    }
    response.writeHead(200, { "Content-Type": TYPES[path.extname(file)] || "application/octet-stream" });
    fs.createReadStream(file).pipe(response);
  });
  return new Promise((ok) => server.listen(0, "127.0.0.1", () => ok(server)));
}

function startChrome() {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), "cu2qu-touch-"));
  const child = spawn("google-chrome", ["--headless=new", "--no-sandbox", "--disable-gpu", "--remote-debugging-port=0", `--user-data-dir=${dir}`, "about:blank"], { stdio: ["ignore", "ignore", "pipe"] });
  return new Promise((ok, bad) => {
    let text = "";
    child.stderr.on("data", (d) => {
      text += d;
      const m = text.match(/DevTools listening on ws:\/\/127\.0\.0\.1:(\d+)\//);
      if (m) ok({ child, dir, port: Number(m[1]) });
    });
    child.on("exit", () => bad(new Error("chrome exited: " + text)));
    setTimeout(() => bad(new Error("chrome did not start: " + text)), 30000);
  });
}

async function open(port) {
  const target = await (await fetch(`http://127.0.0.1:${port}/json/new?about:blank`, { method: "PUT" })).json();
  const ws = new WebSocket(target.webSocketDebuggerUrl);
  await new Promise((ok, bad) => { ws.onopen = ok; ws.onerror = bad; });
  let id = 0;
  const pending = new Map();
  const waiting = [];
  ws.onmessage = (message) => {
    const msg = JSON.parse(message.data);
    if (msg.id) {
      const p = pending.get(msg.id);
      pending.delete(msg.id);
      if (msg.error) p.bad(new Error(JSON.stringify(msg.error)));
      else p.ok(msg.result);
    } else {
      for (const w of waiting.slice()) {
        if (w.method === msg.method) {
          waiting.splice(waiting.indexOf(w), 1);
          w.ok(msg.params);
        }
      }
    }
  };
  const send = (method, params) => new Promise((ok, bad) => {
    pending.set(++id, { ok, bad });
    ws.send(JSON.stringify({ id, method, params: params || {} }));
  });
  const event = (method) => new Promise((ok) => waiting.push({ method, ok }));
  const evaluate = async (expression) => {
    const r = await send("Runtime.evaluate", { expression, returnByValue: true, awaitPromise: true });
    if (r.exceptionDetails) throw new Error(JSON.stringify(r.exceptionDetails.exception || r.exceptionDetails));
    return r.result.value;
  };
  const goto = async (url) => {
    const loaded = event("Page.loadEventFired");
    await send("Page.navigate", { url });
    await loaded;
    // until the figures are live and the formulas are typeset, so that nothing moves afterwards
    await evaluate(`new Promise((ok) => { const t = setInterval(() => {
      if (document.querySelectorAll("figure.live").length === 5 && window.MathJax && MathJax.startup && MathJax.startup.promise) {
        clearInterval(t); MathJax.startup.promise.then(() => setTimeout(ok, 500));
      } }, 100); })`);
    await evaluate(`(() => { window.__cancelled = 0; window.addEventListener("pointercancel", () => { window.__cancelled++; }, true); })()`);
  };
  await send("Page.enable");
  await send("Runtime.enable");
  return { send, evaluate, goto, close: () => ws.close() };
}

const sleep = (ms) => new Promise((ok) => setTimeout(ok, ms));

// Where a handle is, after bringing it into view, and what the figure shows.
const LOOK = (fig, handle) => `(() => {
  const f = document.querySelector('figure[data-fig="${fig}"]');
  const svg = f.querySelector("svg");
  f.scrollIntoView({ block: "center", behavior: "instant" });
  let el = svg.querySelector('[data-handle="${handle}"]');
  let b = el.getBoundingClientRect();
  f.scrollLeft += b.left + b.width / 2 - (f.getBoundingClientRect().left + f.clientWidth / 2);
  b = el.getBoundingClientRect();
  const x = b.left + b.width / 2, y = b.top + b.height / 2;
  const hit = document.elementFromPoint(x, y);
  return { x, y, onHandle: !!(hit && hit.closest("[data-handle]")), scrollY, scrollLeft: f.scrollLeft,
    curve: svg.querySelector("path.f-cubic").getAttribute("d"), cancelled: window.__cancelled, dragging: svg.classList.contains("dragging") };
})()`;
const STATE = (fig) => `(() => {
  const f = document.querySelector('figure[data-fig="${fig}"]');
  const svg = f.querySelector("svg");
  const cubic = svg.querySelector("path.f-cubic");
  return { scrollY, scrollLeft: f.scrollLeft, curve: cubic ? cubic.getAttribute("d") : "", cancelled: window.__cancelled,
    dragging: svg.classList.contains("dragging"), readout: f.querySelector(".fig-readout").textContent };
})()`;

let failures = 0;
function report(name, ok, detail) {
  console.log(`${ok ? "ok    " : "FAILED"} ${name}: ${detail}`);
  if (!ok) failures++;
}

async function touches(tab, url) {
  await tab.send("Emulation.setDeviceMetricsOverride", { width: 390, height: 844, deviceScaleFactor: 2, mobile: true });
  await tab.send("Emulation.setTouchEmulationEnabled", { enabled: true, maxTouchPoints: 5 });
  await tab.goto(url);
  const width = await tab.evaluate("[innerWidth, document.documentElement.scrollWidth]");
  report("a phone 390 wide", width[0] === 390 && width[1] <= 391, `the window is ${width[0]} wide and the page ${width[1]}`);
  const touch = (type, points) => tab.send("Input.dispatchTouchEvent", { type, touchPoints: points.map(([x, y]) => ({ x, y, id: 1, radiusX: 8, radiusY: 8, force: 1 })) });
  const path = [[2, -12], [4, -24], [8, -40], [14, -56], [22, -70], [34, -80], [50, -86], [70, -90]];

  for (const [fig, handle] of [["construction", 1], ["construction", 3], ["pace", 1], ["pace", 2]]) {
    const before = await tab.evaluate(LOOK(fig, handle));
    await sleep(200);
    await touch("touchStart", [[before.x, before.y]]);
    await sleep(60);
    for (const [dx, dy] of path) {
      await touch("touchMove", [[before.x + dx, before.y + dy]]);
      await sleep(30);
    }
    await touch("touchEnd", []);
    await sleep(300);
    const after = await tab.evaluate(STATE(fig));
    const ok = before.onHandle && after.curve !== before.curve && Math.abs(after.scrollY - before.scrollY) < 1 &&
      Math.abs(after.scrollLeft - before.scrollLeft) < 1 && after.cancelled === before.cancelled && !after.dragging;
    report(`a finger on control point ${handle} of figure "${fig}"`, ok,
      `the curve ${after.curve !== before.curve ? "follows" : "does not move"}; the page scrolls ${Math.round(after.scrollY - before.scrollY)}, the figure ${Math.round(after.scrollLeft - before.scrollLeft)}; ` +
      `pointers cancelled: ${after.cancelled - before.cancelled}`);
  }

  // a finger on a figure, off its control points, scrolls the page as anywhere
  for (const fig of ["construction", "pace"]) {
    const at = await tab.evaluate(`(() => {
      const f = document.querySelector('figure[data-fig="${fig}"]');
      f.scrollIntoView({ block: "center", behavior: "instant" });
      f.scrollLeft = 0;
      const b = f.querySelector("svg").getBoundingClientRect();
      const x = b.left + 150, y = b.top + b.height * 0.55;
      const hit = document.elementFromPoint(x, y);
      return { x, y, onHandle: !!(hit && hit.closest("[data-handle]")), inFigure: !!(hit && hit.closest("figure")), scrollY, curve: f.querySelector("path.f-cubic").getAttribute("d") };
    })()`);
    await sleep(200);
    await touch("touchStart", [[at.x, at.y]]);
    await sleep(60);
    for (let i = 1; i <= 10; i++) {
      await touch("touchMove", [[at.x, at.y - 20 * i]]);
      await sleep(30);
    }
    await touch("touchEnd", []);
    await sleep(500);
    const after = await tab.evaluate(STATE(fig));
    report(`a finger on figure "${fig}", off its control points`, at.inFigure && !at.onHandle && after.scrollY - at.scrollY > 100 && after.curve === at.curve,
      `the page scrolls ${Math.round(after.scrollY - at.scrollY)}; the curve ${after.curve === at.curve ? "stays" : "MOVES"}`);
  }

  // a tap on the plot of figure 3 reads it
  const plot = await tab.evaluate(`(() => {
    const f = document.querySelector('figure[data-fig="count"]');
    f.scrollIntoView({ block: "center", behavior: "instant" });
    f.scrollLeft = 0;
    const b = f.querySelector("svg").getBoundingClientRect();
    return { x: b.left + b.width * 0.35, y: b.top + b.height * 0.4, readout: f.querySelector(".fig-readout").textContent, top: f.querySelector(".fig-controls").getBoundingClientRect().top };
  })()`);
  await sleep(200);
  await touch("touchStart", [[plot.x, plot.y]]);
  await sleep(60);
  await touch("touchEnd", []);
  await sleep(300);
  const read = await tab.evaluate(`(() => { const f = document.querySelector('figure[data-fig="count"]'); return { readout: f.querySelector(".fig-readout").textContent, top: f.querySelector(".fig-controls").getBoundingClientRect().top }; })()`);
  report("a tap on the plot of figure 3", /third difference \/ tolerance = [\d.]+: cu2qu \d+/.test(read.readout) && !/third difference/.test(plot.readout),
    `it says "${read.readout}"; what is under it moves ${Math.round(read.top - plot.top)}`);
}

async function mouse(tab, url) {
  await tab.send("Emulation.setDeviceMetricsOverride", { width: 1200, height: 900, deviceScaleFactor: 1, mobile: false });
  await tab.goto(url);
  const move = (type, x, y, button, buttons) => tab.send("Input.dispatchMouseEvent", { type, x, y, button, buttons, clickCount: type === "mouseMoved" ? 0 : 1, pointerType: "mouse" });
  for (const [button, buttons, drags] of [["right", 2, false], ["left", 1, true]]) {
    const before = await tab.evaluate(LOOK("construction", 2));
    await move("mouseMoved", before.x, before.y, "none", 0);
    await move("mousePressed", before.x, before.y, button, buttons);
    for (let i = 1; i <= 6; i++) await move("mouseMoved", before.x - 12 * i, before.y + 10 * i, button, buttons);
    await move("mouseReleased", before.x - 72, before.y + 60, button, 0);
    await sleep(200);
    const after = await tab.evaluate(STATE("construction"));
    report(`the ${button} mouse button on a control point`, before.onHandle && (after.curve !== before.curve) === drags && !after.dragging,
      `the curve ${after.curve !== before.curve ? "follows" : "does not move"}`);
    // afterwards the pointer moves freely
    await move("mouseMoved", before.x - 200, before.y + 150, "none", 0);
    await sleep(100);
    const later = await tab.evaluate(STATE("construction"));
    report("the pointer moving on with no button", later.curve === after.curve, `the curve ${later.curve === after.curve ? "stays" : "MOVES"}`);
  }
}

(async () => {
  const server = await serve(path.resolve(site));
  const chrome = await startChrome();
  const url = `http://127.0.0.1:${server.address().port}/${page}`;
  try {
    let tab = await open(chrome.port);
    await touches(tab, url);
    tab.close();
    tab = await open(chrome.port);
    await mouse(tab, url);
    tab.close();
  } finally {
    chrome.child.removeAllListeners("exit");
    chrome.child.kill();
    server.close();
    await sleep(300);
    fs.rmSync(chrome.dir, { recursive: true, force: true });
  }
  console.log(failures ? `${failures} failed` : "touch and mouse dragging work");
  process.exit(failures ? 1 : 0);
})().catch((e) => {
  console.error("FAILED", e);
  process.exit(1);
});
