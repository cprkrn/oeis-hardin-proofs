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

## Ready-to-paste OEIS edits

A197230: add
  %F A197230 The empirical recurrence above fails at n=23: a(23) = 1327965802062332, while the recurrence gives 1327965802062198. See link. - _Cooper Kernan_, Nov 12 2026
  %H A197230 Cooper Kernan, <a href="https://github.com/cprkrn/oeis-hardin-proofs/blob/main/ERRATA.md">Counterexample to the empirical recurrence</a>

A196074: add
  %F A196074 The empirical recurrence above fails at n=38: a(38) = 299231374661523, while the recurrence gives 299231374661916. See link. - _Cooper Kernan_, Nov 12 2026
  %H A196074 Cooper Kernan, <a href="https://github.com/cprkrn/oeis-hardin-proofs/blob/main/ERRATA.md">Counterexample to the empirical recurrence</a>


# Errata: 26 incorrect "Empirical g.f." formulas

Found by expanding each g.f. and comparing it with the listed terms; see `batch3_gf/gf_results.csv` for details. A sample of six was independently re-expanded with sympy and confirmed.

- **Wrong by a constant factor:** A194476, A194481 (1/3 of true terms); A232591 (1/2); A280227 (a(n+1)/2).
- **Belong to a different sequence or are wrong outright:** A201081 (gives 3,5,15,21,... vs 3,7,13,25,...); A208410 (gives 2,15,182,... vs 5,187,10682,...).
- **Missing constant term:** A183777 (a(0)=16, but the g.f. gives 0; correct for n>=1).
- **Offset (off by a factor of x):** A197274, A297953 (proved); A068778, A187452, A200456, A201875, A202883, A209603, A209962, A220172, A221241, A230470, A231081, A239532, A268133, A274099, A274233, A284076 (x*G(x) matches all listed terms); A144135 (G(x)/x matches).
