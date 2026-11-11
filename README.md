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
- **Not yet done:** batch 2 has not been independently re-implemented. The two refutations have been confirmed independently.

Total new proofs in this repository: 887 + 2,045 = 2,932.
