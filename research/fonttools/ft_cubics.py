#!/usr/bin/env python3
"""Write the cubics of some CFF fonts to a JSON file (font units).
Usage: ft_cubics.py out.json max_glyphs font..."""
import json, sys
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import RecordingPen

out, max_glyphs = sys.argv[1], int(sys.argv[2])
cubics = []
for path in sys.argv[3:]:
    font = TTFont(path, fontNumber=0)
    gs = font.getGlyphSet()
    n0 = len(cubics)
    for name in font.getGlyphOrder()[:max_glyphs]:
        pen = RecordingPen()
        try:
            gs[name].draw(pen)
        except Exception:
            continue
        cur = None
        for op, args in pen.value:
            if op in ('moveTo', 'lineTo'):
                cur = args[0]
            elif op == 'curveTo':
                # may be a super-bezier with more than 3 points; keep plain cubics
                if len(args) == 3 and cur is not None:
                    cubics.append([cur, args[0], args[1], args[2]])
                cur = args[-1]
            elif op == 'qCurveTo':
                cur = args[-1]
    print(path, len(cubics) - n0, 'cubics', file=sys.stderr)
json.dump(cubics, open(out, 'w'))
print('total', len(cubics), file=sys.stderr)
