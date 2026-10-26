# Certified oblivious strategies

`eps_nNN_sS.json`: the strategy achieving the certified lower bound on eps*(n, s) in
data/oblivious_v1.csv, written by `spr obl grid`. Cells are 0-based; each entry is
`[tuple, "p/q"]` with an exact rational probability. Re-verify with
`spr obl verify FILE --claim p/q`. Instances whose certified lower bound is 0 have no file.

`golden_sS.json`: the golden coin-driven strategy on F(s+3) - 2 cells
(scripts/golden.py), with masses in Q(sqrt 5) written as `[a, b]` for `a + b*sqrt5`
(a, b rational strings). Its eps is exactly (3 - sqrt 5)/2; the check is exact in
Q(sqrt 5) in scripts/golden.py.
