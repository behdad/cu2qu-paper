#!/usr/bin/env python3
"""Write the cubics of some CFF fonts to a JSON file, in font units.

Usage: corpus.py out.json max_glyphs font...

The corpus of the note is the first 3,000 glyphs of each of three fonts as
Fedora 43 ships them:

    corpus.py data/cubics.json 3000 \\
        /usr/share/fonts/adobe-source-code-pro-fonts/SourceCodePro-Regular.otf \\
        /usr/share/fonts/stix-fonts/STIXTwoMath-Regular.otf \\
        /usr/share/fonts/google-noto-sans-sc-fonts/NotoSansSC-Regular.otf

A second file, out.json with ".fonts.txt" in place of ".json", records what
was read: file, version string, units per em, glyphs and cubics taken.
"""
import hashlib
import json
import sys

from fontTools.pens.recordingPen import RecordingPen
from fontTools.ttLib import TTFont

out, max_glyphs = sys.argv[1], int(sys.argv[2])
cubics = []
notes = []
for path in sys.argv[3:]:
    font = TTFont(path, fontNumber=0)
    glyph_set = font.getGlyphSet()
    before = len(cubics)
    names = font.getGlyphOrder()[:max_glyphs]
    for name in names:
        pen = RecordingPen()
        glyph_set[name].draw(pen)
        current = None
        for op, args in pen.value:
            if op in ("moveTo", "lineTo"):
                current = args[0]
            elif op == "curveTo":
                if len(args) == 3 and current is not None:
                    cubics.append([current, args[0], args[1], args[2]])
                current = args[-1]
            elif op == "qCurveTo":
                current = args[-1]
    version = font["name"].getDebugName(5)
    sha = hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]
    notes.append(
        "%s\n  %s\n  units per em %d, %d of %d glyphs, %d cubics, sha256 %s..."
        % (
            path,
            version,
            font["head"].unitsPerEm,
            len(names),
            font["maxp"].numGlyphs,
            len(cubics) - before,
            sha,
        )
    )
with open(out, "w") as f:
    json.dump(cubics, f)
notes.append("total %d cubics" % len(cubics))
with open(out.replace(".json", ".fonts.txt"), "w") as f:
    f.write("\n".join(notes) + "\n")
print("\n".join(notes))
