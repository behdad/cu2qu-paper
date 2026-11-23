// Run the cubic2quad npm library (fontello/cubic2quad, used by svg2ttf) over the corpus.
//
//   node c2q_js.js <path/to/cubic2quad/index.js> <corpus.json> <tol> <out.bin>
//
// Output: "C2QRES01", u32 count, per cubic u32 nseg, nseg * {u32 kind=2, 6 f64}.
'use strict'
const fs = require('fs')
const path = require('path')

const [lib, corpusPath, tolStr, out] = process.argv.slice(2)
const cubic2quad = require(path.resolve(lib))
const tol = parseFloat(tolStr)
const cubics = JSON.parse(fs.readFileSync(corpusPath, 'utf8'))

function convert (c) {
  return cubic2quad(c[0][0], c[0][1], c[1][0], c[1][1], c[2][0], c[2][1], c[3][0], c[3][1], tol)
}

// warm up, then time
for (let i = 0; i < 2000; i++) convert(cubics[i % cubics.length])
let best = Infinity
let results
for (let rep = 0; rep < 3; rep++) {
  const t0 = process.hrtime.bigint()
  results = cubics.map(convert)
  const dt = Number(process.hrtime.bigint() - t0)
  if (dt < best) best = dt
}

let nq = 0
for (const r of results) nq += (r.length - 2) / 4
const buf = Buffer.alloc(12 + results.length * 4 + nq * 52)
let pos = 0
buf.write('C2QRES01', pos, 'latin1'); pos += 8
buf.writeUInt32LE(results.length, pos); pos += 4
for (const r of results) {
  const n = (r.length - 2) / 4
  buf.writeUInt32LE(n, pos); pos += 4
  for (let k = 0; k < n; k++) {
    buf.writeUInt32LE(2, pos); pos += 4
    for (let j = 0; j < 6; j++) {
      buf.writeDoubleLE(r[4 * k + j], pos); pos += 8
    }
  }
}
fs.writeFileSync(out, buf)
console.log(JSON.stringify({
  method: 'cubic2quad_js',
  tol,
  ns_per_cubic: Math.round(best / cubics.length * 10) / 10,
  quads: nq
}))
