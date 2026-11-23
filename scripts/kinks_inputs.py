#!/usr/bin/env python3
"""List the masters that kinks.py and compare_revs_masters.py read.

Usage: kinks_inputs.py MASTERS_DIR

MASTERS_DIR holds NAME.npz with NAME.npz.meta.json beside it, one for each
font, written by research/kinks/stage1.py: the point-compatible cubic masters
of a variable font (a CFF2 font as shipped, or a design space from its
sources) as arrays. This prints, for each, where it came from and what it
holds, and the version of each origin: the package of a font file, the
commit of a checkout. The arrays themselves (272 MB) are not in this
repository.
"""

import glob
import json
import os
import subprocess
import sys


def run(*command):
    try:
        return subprocess.run(command, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main():
    directory = sys.argv[1]
    home = os.path.expanduser("~")
    origins = {}
    print("%-24s %7s %7s %9s %5s  %s" % ("name", "masters", "glyphs", "segments", "upem", "from"))
    for path in sorted(glob.glob(os.path.join(directory, "*.npz.meta.json"))):
        m = json.load(open(path))
        name = os.path.basename(path)[: -len(".npz.meta.json")]
        source = m["path"].replace(home, "~")
        if m.get("font_number") not in (None, 0, -1):
            source += " (font %s)" % m["font_number"]
        skipped = m.get("skipped") or {}
        note = "" if not skipped else "; skipped: %s" % json.dumps(skipped)
        print("%-24s %7d %7s %9s %5s  %s%s" % (name, len(m["locations"]), m.get("glyphs", ""), m.get("segments", ""), m["upem"], source, note))
        top = run("git", "-C", os.path.dirname(m["path"]), "rev-parse", "--show-toplevel")
        if top:
            origins[top.replace(home, "~")] = "%s at %s (%s)" % (
                run("git", "-C", top, "remote", "get-url", "origin"),
                run("git", "-C", top, "rev-parse", "HEAD"),
                run("git", "-C", top, "log", "-1", "--format=%cs"),
            )
        else:
            origins[m["path"]] = "package %s" % run("rpm", "-qf", "--qf", "%{NAME}-%{VERSION}-%{RELEASE}.%{ARCH}", m["path"])
    print()
    for origin in sorted(origins):
        print("%s: %s" % (origin, origins[origin]))


if __name__ == "__main__":
    main()
