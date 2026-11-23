#!/usr/bin/env python3
"""Print sha/date/author/subject from a saved GitHub 'commits' API page."""
import json
import sys

d = json.load(open(sys.argv[1]))
last = int(sys.argv[2]) if len(sys.argv) > 2 else 60
print(len(d), "commits on page")
for c in d[-last:]:
    print(
        c["sha"][:10],
        c["commit"]["author"]["date"],
        c["commit"]["author"]["name"],
        "|",
        c["commit"]["message"].split("\n")[0][:100],
    )
