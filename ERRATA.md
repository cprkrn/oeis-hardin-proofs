# Errata: incorrect "Empirical" formulas in the OEIS

How each formula below was checked:
- The sequence was recomputed exactly from its OEIS definition.
- Every listed OEIS term is reproduced, so the terms are correct and the error is in the formula.
- Every item was confirmed by two independent programs written separately (`batch2/` and `batch3_gf/`, plus an independent audit).

## 1. Recurrences that are truncated (2)

The printed coefficients are exactly the leading coefficients of the true recurrence, with the final terms missing. The minimal orders are proved (exact computation).

| A-number | Printed order | True minimal order | First failure | Missing terms to append |
|---|---|---|---|---|
| A197230 | 22 | 25 | n=23: a(23) = 1327965802062332, but the printed recurrence gives 1327965802062198 | `+28*a(n-23) -82*a(n-24) +40*a(n-25)` |
| A196074 | 37 | 43 | n=38: a(38) = 299231374661523, but the printed recurrence gives 299231374661916 | `+274*a(n-38) +112*a(n-39) +84*a(n-40) -56*a(n-41) +8*a(n-42) -8*a(n-43)` |

With the missing terms appended, each recurrence was checked to hold over a further 59 consecutive terms beyond its order.

**Ready-to-paste OEIS edits**

    %F A197230 The empirical recurrence above is incomplete: the full recurrence (order 25, minimal) adds +28*a(n-23) -82*a(n-24) +40*a(n-25). Without these terms it fails at n=23 (a(23) = 1327965802062332). See link. - _Cooper Kernan_, Nov 12 2026
    %H A197230 Cooper Kernan, <a href="https://github.com/cprkrn/oeis-hardin-proofs/blob/main/ERRATA.md">Errata and full recurrence</a>

    %F A196074 The empirical recurrence above is incomplete: the full recurrence (order 43, minimal) adds +274*a(n-38) +112*a(n-39) +84*a(n-40) -56*a(n-41) +8*a(n-42) -8*a(n-43). Without these terms it fails at n=38 (a(38) = 299231374661523). See link. - _Cooper Kernan_, Nov 12 2026
    %H A196074 Cooper Kernan, <a href="https://github.com/cprkrn/oeis-hardin-proofs/blob/main/ERRATA.md">Errata and full recurrence</a>

## 2. Generating functions that are substantively wrong (6), plus one wrong constant term

| A-number | Problem |
|---|---|
| A194476, A194481 | The g.f.'s coefficients are exactly 1/3 of the true terms. |
| A232591 | The g.f.'s coefficients are exactly 1/2 of the true terms. |
| A280227 | The g.f. gives exactly a(n)/2 for every n. |
| A201081 | The g.f. expands to 3,5,15,21,... but the true terms are 3,7,13,25,47,... |
| A208410 | The g.f. is character-for-character the g.f. of A208409, which it generates correctly; it was evidently copied onto the wrong entry. |
| A183777 | The constant term is wrong: a(0)=16 but the g.f. gives 0. The g.f. is proved correct for n>=1. |

## 3. Offset corrections only (20)

These rational functions are mathematically right, but one power of x is off for the entry's offset. Two of these entries (A144135, A187452) carry later comments calling the g.f. verified, so editors may treat this group as cosmetic.

- **x·G(x) matches all terms:** A197274 and A297953 (proved), and A068778, A187452, A200456, A201875, A202883, A209603, A209962, A220172, A221241, A230470, A231081, A239532, A268133, A274099, A274233, A284076.
- **G(x)/x matches:** A144135.

## Suggested order of submission

Section 1 first, then section 2, then section 3 as a single batch note.
