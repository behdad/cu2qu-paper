#!/bin/bash
# Run the shared virtualenv's python on the fontTools tree of the agent clone
# (branch as checked out there), pure Python unless FT_TREE says otherwise.
W=/home/behdad/.cache/cu2qu-work
export PYTHONPATH=${FT_TREE:-$W/agent-kinks/fonttools}/Lib
export PYTHONDONTWRITEBYTECODE=1
exec $W/venv/bin/python "$@"
