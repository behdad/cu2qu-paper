"""List designspaces under a directory with their axes and sources."""

import glob
import os
import sys

from fontTools.designspaceLib import DesignSpaceDocument

for path in sorted(glob.glob(os.path.join(sys.argv[1], "*.designspace"))):
    try:
        doc = DesignSpaceDocument.fromfile(path)
    except Exception as e:
        print(os.path.basename(path), "ERROR", e)
        continue
    full = [s for s in doc.sources if not s.layerName]
    missing = [s.filename for s in full if not os.path.isdir(s.path)]
    size = 0
    if full and not missing:
        size = len(os.listdir(os.path.join(full[0].path, "glyphs")))
    print(
        "%-45s axes=%s sources=%d full=%d missing=%d glyphs~%d"
        % (
            os.path.basename(path),
            ",".join(a.tag for a in doc.axes),
            len(doc.sources),
            len(full),
            len(missing),
            size,
        )
    )
