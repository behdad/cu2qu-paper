#!/bin/bash
# Everything, on the current working tree of ~/fonttools: black, rebuild the
# compiled export, cu2qu tests pure and compiled, odd inputs, differential
# dumps, whole test suite pure and compiled.
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=/home/behdad/hb/cu2qu-work
echo "=== black"
(cd /home/behdad/fonttools && "$W/venv/bin/black" --check Lib/fontTools/cu2qu/cu2qu.py Tests/cu2qu/cu2qu_test.py 2>&1 | tail -2)
echo "=== build"
bash "$S/ft_build_cy.sh" new 2>&1 | tail -1
echo "=== cu2qu tests"
bash "$S/ft_check_cy.sh" 2>&1 | grep -E "^==|passed|failed|error" | head -12
echo "=== odd inputs"
bash "$S/ft_odd.sh" 2>&1 | cut -c1-200 | head -24
echo "=== dumps"
cp "$W/out.new-cy.txt" "$W/out.new-cy.prev.txt"
bash "$S/ft_dumps.sh" 2>&1 | tail -5
cmp -s "$W/out.new-cy.txt" "$W/out.new-cy.prev.txt" && echo "new compiled dump: same as the previous one" || echo "new compiled dump: CHANGED from the previous one"
rm -f "$W/out.new-cy.prev.txt"
echo "=== whole suite"
bash "$S/ft_fulltests.sh" all 2>&1 | tail -6
