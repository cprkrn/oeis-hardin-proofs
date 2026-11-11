# Proofs of 887 empirical recurrences of R. H. Hardin in the OEIS

**Draft note, Nov 2026.** Author: _Your Name_. The proofs were computer-assisted, using an AI agent.

## Summary

R. H. Hardin contributed thousands of OEIS entries counting arrays under local constraints. Each has an "Empirical:" linear recurrence that was fitted to computed terms but never proved. We prove 887 of these recurrences, in two families:

| Family (OEIS name pattern) | Proved |
|---|---|
| "Number of nXk 0..m arrays with every element equal/unequal to S king-move (or horizontally, vertically or antidiagonally) adjacent elements, with upper left element zero" | 752 |
| "Number of (n+1) X (k+1) 0..m arrays with every 2 X 2 subblock having its diagonal sum differing from its antidiagonal sum by c [, with no adjacent elements equal]" (constant-stress tilings) | 135 |

None of the 887 entries mentioned a proof at the time of writing. No empirical recurrence we tested was false. 137 further entries in these families were too large for this run.

## Method

Fix the width. Generate arrays row by row, with state = the last two rows (first family) or the last row (second family). Then a(n) = uᵀ Tⁿ⁻ᶜ v for a transfer matrix T on D states. By Cayley–Hamilton, (a(n)) satisfies a linear recurrence of order at most D.

Let Q be Hardin's recurrence. Then b = Q(E)a also satisfies that order-D recurrence. If b(n) = 0 for D consecutive n, b vanishes for all later n. We check Q exactly, in integer arithmetic, on Hardin's whole stated range up to D consecutive terms beyond it. That proves the formula as stated.

**Interpretation check.** Each definition is encoded from the OEIS name. Every term listed in the OEIS (%S/%T/%U) is reproduced exactly, so the encoding matches the intended problem.

## Reproducibility

Files:
- `hardin_prove.py` and `hardin_stress.py`: the programs.
- `proved.csv`: per entry, the A-number, the number of states D, the recurrence order and the formula.
- `oeis_edits.txt`: proposed OEIS edits.

Data: the OEIS mirror at github.com/oeis/oeisdata, Nov 2026. Each program re-derives every result from that data.

**Independent verification.** A second AI model wrote fresh code, sharing none with the original. For 677 of the neighbour-count entries (all of widths 2–4) it ran three independent checks:
- direct brute-force enumeration of the small arrays;
- its own transfer matrix;
- its own Cayley–Hamilton proof.

All 677 were proved, none refuted, and all 17,185 OEIS terms matched. Each recurrence holds from exactly the start Hardin stated. The "antidiagonally adjacent" family means 6 neighbours, not 8, and the brute-force counts confirm this reading. The 30 "Empirical g.f." lines by C. Barker on these entries follow from the proved recurrences.

The other 210 proofs (75 of width 5, plus the 135 constant-stress tilings) use the same programs and the same argument, and reproduce every OEIS term. They have not been re-implemented independently.
