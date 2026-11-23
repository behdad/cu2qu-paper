#!/usr/bin/env python3
"""json2bin.py in.json out.bin -- write a JSON list of cubics as N x 8 doubles."""
import json, sys
import numpy as np
c = np.array(json.load(open(sys.argv[1])), dtype=np.float64).reshape(-1, 8)
c.tofile(sys.argv[2])
print(len(c), 'cubics ->', sys.argv[2])
