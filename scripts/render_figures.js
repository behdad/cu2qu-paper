#!/usr/bin/env node
// Put the figures into the page, as they are before anything is moved.
//
// Usage: node render_figures.js page.html cu2qu-figures.js [numbers.json] [--check]
//
// The page holds each figure between two comments,
//
//     <!-- figure:NAME -->
//     ...
//     <!-- /figure:NAME -->
//
// and this script replaces what is between them with the SVG that cu2qu-figures.js
// draws for its default state, so that the page shows the figures without
// scripts too. With scripts, cu2qu-figures.js redraws them in place and makes them
// live. The numbers the captions quote are written to numbers.json, where
// check_page.py compares them with the text.

const fs = require("fs");
const path = require("path");

const args = process.argv.slice(2);
const check = args.includes("--check");
const [page, lib, numbers] = args.filter((a) => !a.startsWith("--"));
if (!page || !lib) {
  console.error("usage: render_figures.js page.html cu2qu-figures.js [numbers.json] [--check]");
  process.exit(2);
}
const F = require(path.resolve(lib));

const html = fs.readFileSync(page, "utf8");
let out = html;
const report = {};
for (const name of F.names) {
  const { svg, readout, drawn } = F.figure(name);
  const block = svg + (readout ? `\n<p class="fig-readout">${readout}</p>` : "");
  const re = new RegExp(`(<!-- figure:${name} -->\\n)[\\s\\S]*?(<!-- /figure:${name} -->)`);
  if (!re.test(out)) {
    console.error(`no place for figure ${name} in ${page}`);
    process.exitCode = 1;
    continue;
  }
  out = out.replace(re, (m, a, b) => a + block + "\n" + b);
  if (name === "construction") report.construction = { n: drawn.n, d3: drawn.d3, formula: drawn.formula, measured: drawn.measured };
  if (name === "count") report.count = { ends: drawn.ends };
  if (name === "kink") report.kink = { angles: drawn.angles };
  if (name === "pace") report.pace = drawn.numbers;
}
if (numbers) fs.writeFileSync(numbers, JSON.stringify(report, null, 1) + "\n");
if (check) {
  if (out !== html) {
    console.error(`${page}: figures are not up to date`);
    process.exit(1);
  }
  console.log(`${page}: figures up to date`);
} else if (out !== html) {
  fs.writeFileSync(page, out);
  console.log(`${page}: figures updated`);
} else {
  console.log(`${page}: figures unchanged`);
}
