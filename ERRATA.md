# Errata: two empirical recurrences that are false

Both entries are counts of 0..4 arrays in which each element x equals the number of its horizontal and vertical neighbours equal to f(x).

For each entry:
- every listed OEIS term is reproduced exactly;
- the stated recurrence fails from the first n it claims;
- the failure was found by `batch2/prove2.py`;
- it was confirmed by an independent program, `batch2/refuted_independent_check.py`.

| A-number | Definition | Claim | First failure | Actual a(n) | Recurrence gives |
|---|---|---|---|---|---|
| A197230 | n X 3, f = (4,3,0,1,2) | order-22 recurrence from n=23 | n=23 | 1327965802062332 | 1327965802062198 |
| A196074 | n X 4, f = (0,3,2,1,4) | order-37 recurrence from n=38 | n=38 | 299231374661523 | 299231374661916 |

The true minimal recurrence orders are about 25 (A197230) and 43 (A196074). These estimates come from Berlekamp–Massey mod a prime and are unproved. Suggested OEIS action: correct or remove the "Empirical:" formula and cite these values.
