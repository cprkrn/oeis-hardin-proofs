# Proofs of R. H. Hardin's empirical OEIS recurrences

- `hardin_prove.py`: neighbour-count families ("every element equal/unequal to S king-move or horizontally/vertically/antidiagonally adjacent elements").
- `hardin_stress.py`: constant-stress 2 X 2 tilings.
- `proved.csv`: every proved entry, with its A-number, number of states, recurrence order and formula.
- `note.md`: the method and the independent verification.

Rerun (needs Python 3 and the OEIS data mirror):

    git clone --depth 1 https://github.com/oeis/oeisdata
    python3 hardin_prove.py oeisdata/seq results1.json
    python3 hardin_stress.py oeisdata/seq results2.json

## Batch 2 (`batch2/`)

A generic prover (`engine.py`, `families.py`, `prove2.py`) covers 30 further Hardin families. It uses exact state merging to bound the recurrence order and works in either orientation.

- **Proved:** 2,067 empirical formulas. 22 of them already had proofs by Robert Israel, so **2,045 are new** (`batch2/proved_new.csv`).
- **False:** 2 (see `ERRATA.md`).
- **Not attempted:** 561 were too large, and 291 have formulas that are not linear recurrences or polynomials.
- **Checks:** every listed OEIS term is reproduced. For every proved entry, adding 1 to one coefficient makes the check fail, as it should.
- **Independent audit:** an independent re-implementation (separate code, no shared files) confirmed all 2,045 proofs. It also confirmed both refutations, which it found to be truncated recurrences (see `ERRATA.md`).

Total new proofs in this repository: 887 + 2,045 = 2,932.

## Batch 3 (`batch3_gf/`): empirical generating functions

- **Checked:** 3,601 "Empirical g.f." formulas, 3,445 of them by Colin Barker.
- **Proved:** 705. Of these, 603 follow from the proved recurrences, and 102 were computed directly. 61 of the 102 come from 5 new array families.
- **Incorrect:** 26 (listed in `ERRATA.md`).
- **Term check:** all listed terms of 2,968 computable entries were rechecked, and no wrong terms were found.

## Independent verification (all batches)

An independent auditor rewrote everything from scratch and had no disagreements:
- 677 + 75 + 135 recurrences in batch 1;
- all 2,045 recurrences in batch 2;
- 704 of the 705 proved g.f.s (A234139 was not checked because of a compute limit);
- all 28 errata.
