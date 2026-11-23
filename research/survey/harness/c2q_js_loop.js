// How many segments does cubic2quad's loop
//   for (var t = 0; t < 1; t += (1 / segmentsCount)) approximation.push(...)
// produce for each segmentsCount, and what does the library then return for
// inflection-free cubics over the corpus (counts by number of segments)?
//
//   node c2q_js_loop.js <path/to/cubic2quad/index.js> <corpus.json> <tol>
'use strict'
const fs = require('fs')
const path = require('path')
for (let segmentsCount = 1; segmentsCount <= 8; segmentsCount++) {
  let k = 0
  let last = 0
  for (let t = 0; t < 1; t += (1 / segmentsCount)) { k++; last = t }
  console.log(`segmentsCount ${segmentsCount}: loop body runs ${k} times; last segment starts at t = ${last}`)
}
const [lib, corpusPath, tolStr] = process.argv.slice(2)
if (lib) {
  const cubic2quad = require(path.resolve(lib))
  const tol = parseFloat(tolStr)
  const cubics = JSON.parse(fs.readFileSync(corpusPath, 'utf8'))
  const hist = {}
  let noInfl = 0
  for (const c of cubics) {
    // inflection test as in the library's solveInflections (p t^2 + q t + r = 0)
    const [x1, y1, x2, y2, x3, y3, x4, y4] = [c[0][0], c[0][1], c[1][0], c[1][1], c[2][0], c[2][1], c[3][0], c[3][1]]
    const p = -(x4 * (y1 - 2 * y2 + y3)) + x3 * (2 * y1 - 3 * y2 + y4) + x1 * (y2 - 2 * y3 + y4) - x2 * (y1 - 3 * y3 + 2 * y4)
    const q = x4 * (y1 - y2) + 3 * x3 * (-y1 + y2) + x2 * (2 * y1 - 3 * y3 + y4) - x1 * (2 * y2 - 3 * y3 + y4)
    const r = x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1 + y3)
    const roots = cubic2quad.quadSolve(p, q, r).filter(t => t > 1e-8 && t < 1 - 1e-8)
    if (roots.length) continue
    noInfl++
    const res = cubic2quad(x1, y1, x2, y2, x3, y3, x4, y4, tol)
    const n = (res.length - 2) / 4
    hist[n] = (hist[n] || 0) + 1
  }
  console.log(`tolerance ${tol}: ${noInfl} cubics without inflection; segments returned:`, JSON.stringify(hist))
}
