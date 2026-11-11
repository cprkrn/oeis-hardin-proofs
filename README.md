# Proofs of R. H. Hardin's empirical OEIS recurrences

- `hardin_prove.py`: neighbour-count families ("every element equal/unequal to S king-move or horizontally/vertically/antidiagonally adjacent elements").
- `hardin_stress.py`: constant-stress 2 X 2 tilings.
- `proved.csv`: every proved entry, with its A-number, number of states, recurrence order and formula.
- `note.md`: the method and the independent verification.

Rerun (needs Python 3 and the OEIS data mirror):

    git clone --depth 1 https://github.com/oeis/oeisdata
    python3 hardin_prove.py oeisdata/seq results1.json
    python3 hardin_stress.py oeisdata/seq results2.json
