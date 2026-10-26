# Reproducing the results

This file lists, for every table and every certified number in the data, the command that
regenerates or verifies it, its measured cost, and the committed file it is checked against.
All commands run from the repository root.

## Quick start

```
scripts/verify_all.sh
```

This is the default target. It builds the workspace, runs the tests, checks the solver
against the Python reference, regenerates the value tables up to moderate `n`, re-runs the
gap-lemma and erasure checks on reduced ranges, re-verifies every stored oblivious strategy
exactly (twice, with independent verifiers), re-certifies 52 of the 59 oblivious instances,
re-derives ten closed forms, checks the golden constructions in Q(sqrt 5), and checks several
derived statements against the data. It compares every result with the committed files and prints
`verify_all: PASS` or the list of failures.

Measured on an Apple silicon laptop with 8 GB RAM and other applications open: about 300 s
wall time with the workspace already built (add about 35 s for a cold build), 123 checks, and
at most about 330 MB for any single step. It writes only to a temporary directory.

The long runs are opt-in, one at a time: `scripts/verify_long.sh TARGET` (table below).

## Toolchain

| tool | version used |
|---|---|
| rustc, cargo | 1.99.1 (stable; `rust-toolchain.toml`) |
| HiGHS bindings | `highs` crate 2.4, `highs-sys` 1.15 (HiGHS is built by cargo) |
| Python | 3.9.6 |
| sympy | 1.14.0 (closed forms) |
| scipy | 1.13.1 (only the search scripts `golden.py` and `cube_compress.py`) |

`scripts/verify_golden.py`, `compare_values.py`, `recertify_obl.py` and `make_tables.py` use
only the standard library. The reference oracle
`reference/spr_ref.py` is run, never modified.

## Every result, its command and its data

Results are named by what they state; a write-up with numbered statements is in preparation.

Times are wall-clock on the machine above; memory is peak resident size of the solver
process. "default" means the step is part of `scripts/verify_all.sh`.

### Exact values of the adaptive game

| result | what | command | cost | checked against |
|---|---|---|---|---|
| value tables | value tables (default: forced n <= 18, full n <= 16, QV n <= 19) | `spr tabulate --variant V --max-n N --max-states 40000000 --out F.csv`, then `python3 scripts/compare_values.py F.csv data/values_V_v1.csv` | default: 5 s, 3 s, 3 s; 330, 165, 165 MB | `data/values_{forced,full,qv}_v1.csv` |
| same, full ranges | forced n <= 20, full n <= 18, QV n <= 22 | `scripts/verify_long.sh forced20`, `full18`, `qv22` | see the long-run table | same |
| value tables, cross-checks | reference implementation and v0 tables | `spr check-oracle` | default: 17 s, 85 MB | `reference/spr_ref.py` (run live), `data/spr_values_v0.md` |
| same | brute force, invariants, reductions as identities | `cargo test --release` | default: about 25 s | (self-contained tests in `crates/*/tests`) |
| Qiao-Valiant threshold values | `opt_QV(2^t - 1, t) = t`, `opt_QV(2^t - 2, t) = t - 1`, t <= 5 | `spr known-values --max-t 5` | default: < 0.1 s | Qiao and Valiant, Lemma 5(1) |
| reuse advantage | reuse advantage (217 of 342; forced > QV on 279 of 420) | `python3 scripts/make_tables.py OUTDIR` (macros in OUTDIR/numbers.tex: `ReuseAdv`, `ForcedGtQV`) | default: < 1 s | the three value CSVs |
| binary-search bound | tightness of B at step-up cells (42 of 64 forced n <= 15; 62 of 106, 54 of 88, 68 of 83) | same script, which also reads the out-of-table values | default: < 1 s | value CSVs and `data/values_small_long_v1.csv` |
| same | out-of-table values opt(g, r), r > 2g | `spr tabulate --variant {full,forced} --max-n 9 --s-hi 40 --max-states 40000000 --out F.csv`, then `compare_values.py F.csv data/values_small_long_v1.csv --where variant=V` | default: < 0.1 s, 2 MB | `data/values_small_long_v1.csv` |
| maximal erasure | maximal erasure, n <= 14, `|x| + r <= 2n` | `scripts/verify_long.sh erasure14` (`spr find-erasure-counterexample --max-n 14 --max-states 120000000`) | 851 s, 1.28 GB | prints `no counterexample`; per-n counts in the output |
| same, reduced | n <= 12 | `spr find-erasure-counterexample --max-n 12 --max-states 40000000` | default: 14 s, 210 MB | 12,223,144 pairs, no counterexample |

### Oblivious no-reuse Pointers

| result | what | command | cost | checked against |
|---|---|---|---|---|
| oblivious optima (lower bounds) | every stored strategy re-verified exactly, by both verifiers, against its certified eps | `spr obl verify FILE --claim P/Q --naive` for each `data/strategies/eps_*.json` | default: about 1.5 s for all 41 files | `info.eps_certified` in each file |
| same (upper bounds and ov*) | Farkas certificates recomputed and checked exactly; ov* exact | `python3 scripts/recertify_obl.py --max-seconds 10` (52 instances) | default: 50 s | `data/oblivious_v1.csv` |
| same, all 59 instances | including (13,5) | `scripts/verify_long.sh obl-grid` | about 25 min (the original grid took 1236 s of solve time; (13,5) alone 984 s) | `data/oblivious_v1.csv` |
| same, classes; the thresholds 2s - 1 and 2^s - 1 against the data | 18 instances with n < 2s - 1 have bracket [0, 1.0002e-9]; 26 with n >= 2^s - 1 are exactly 1/2; 15 in between, widest bracket 5.6e-8 | `make_tables.py` (macros `NObl*`, `OblMaxWidth`; fails if the data disagree with either lemma) | default | `data/oblivious_v1.csv` |
| closed forms | closed forms, ten instances | `spr obl active-set N S --out B.json`, then `python3 scripts/closed_form.py B.json` | default: about 20 s in all | prints `lower <= closed <= upper: True` |
| same, (10,5), (12,5), (13,5) | higher-degree and slow instances | `scripts/verify_long.sh closed-forms` | (10,5): 13 s; (12,5): 198 s; (13,5): about 17 min (solve) | same |
| product construction | products (3,2)x(3,2) and (3,2)x(7,3) give exactly 1/4 | `spr obl product A.json B.json` | default: < 0.1 s | prints `lemma holds` |
| binomial cube of Dagan et al. | binomial cube exactly 2^-k for (3,1), (4,1), (4,2), (5,2) | `spr obl cube D K` | default: < 1 s | |
| oblivious binary search | oblivious binary search on 15 cells has eps = 1/2 | `spr obl binary-search 4` | default | |
| golden family | golden family, s = 3..7, exact in Q(sqrt 5) | `python3 scripts/verify_golden.py` | default: 150 s together with the next row | `data/strategies/golden_s*.json` |
| golden cube | golden and fair cube, eight (d, k) | same script (enumerates every coin outcome) | (included above; (5,2) dominates) | |
| cube compression | fewest cells for the cube's survival events | `scripts/verify_long.sh cube-compress` (`python3 scripts/cube_compress.py D K {fair,golden}`) | about 2 s each, but up to 2.2 GB for (5,3) | sizes quoted in the paper; the MILP is a search device, the expanded strategy is checked exactly |
| small instances | 21 pairs pass the caps, the threshold 2s - 1 removes 8, 13 remain (first (15,8)), needing eps between 0.437 and 0.500 | `make_tables.py` (macros `Range*`) | default | `data/values_qv_v1.csv` |

### Adaptive structure and further measurements

| result | what | command | cost | checked against |
|---|---|---|---|---|
| gap lemma check and slack | gap lemma on every forced-reachable board, n <= 14; slack pattern | `spr gap-check --max-n 14 --max-states 40000000` | default: 0.3 s, 32 MB | prints `no violation: 245759 boards, 4734974 checks`; the default target also checks the tight total 2,572,341 and the slack pattern |
| tensorization exponents | tensorization exponents; `opt_QV(20, 20) = 6` gives beta = 0.4182 | `make_tables.py` (`OUTDIR/tensor.tex`, macros `BetaAtOne`, `TensorBestPhi`) | default | `data/values_qv_v1.csv` |
| root values | full = forced on all 342 cells | `make_tables.py` (macro `FullEqForced`) | default | value CSVs |
| refills | refills along principal variations, n <= 16 | `spr gap-pv --max-n 16 --out F --summary G --max-states 40000000`, then `cmp` | default: 1 s, 85 MB | `data/gap_decomposition_v1.csv`, `data/gap_pv_summary_v1.csv` (byte-identical) |
| same | counts and survival rates (843 refills on 192 lines, 160 of 163, 498 kept; 0.591 against 0.436 reweighted) | `make_tables.py` (macros `Refill*`, `OtherRateReweighted`) | default: < 1 s | same files |
| long lines | principal variations at s = n and s = 2n, n <= 18 (forced) | `spr gap-pv --diag-only --max-n 18 --out F --summary G --max-states 40000000`, then `cmp` | default: 5 s, 330 MB | `data/long_game_pv_v1.csv`, `data/long_game_pv_summary_v1.csv` (byte-identical) |
| policy values | fixed-policy values | `spr policy-eval --policy P [--full] --max-n N --s-mode {both,all} --out F --max-states 40000000`, then `compare_values.py F data/policy_values_v1.csv --where policy=P --where labeler=L` (the exact list is in `verify_all.sh`) | default: < 1 s each | `data/policy_values_v1.csv` |
| staircase recursions | staircase recursions consistent on all 420 forced cells (t = 3, 4, 5) | `make_tables.py` (macro `StairConsistent`; the script fails if the count changes) | default | `data/values_forced_v1.csv` |
| first moves | first moves, forced, n <= 15 | `spr first-moves --variant forced --max-n 15 --out F.csv`, then `cmp` | default: 0.3 s, 43 MB | `data/first_moves_forced_v1.csv` (byte-identical) |
| same, claims | every first move worth at least opt(n, s - 1); centred interval at step-up columns | `make_tables.py` (fails if either claim fails on the data) | default | `data/first_moves_forced_v1.csv` |
| same, QV remark | a corner first move at QV (3, 3) is worth 1 | `spr principal-variation 3 3 qv` | default | prints `0:1 1:2 2:1` |

The golden family placements were found by `python3 scripts/golden.py S --save F` (a MILP search
device using scipy); the default target verifies the stored results without scipy.

### Long runs

Run one at a time with `scripts/verify_long.sh TARGET`, on an otherwise idle machine.

| target | what | measured | memory |
|---|---|---|---|
| `forced20` | forced table n <= 20 | 43 s in the original run (row n = 20 swap-affected) | 53.3M states; table 0.87 GB, about 1.3 GB while the table grows (estimate from the table layout, not measured) |
| `full18` | full table n <= 18 | 1002 s in the original run, of which 991 s for row n = 18 under heavy swap; without swap it should be far less (not measured) | 45.1M states; table about 0.87 GB, about 1.3 GB while growing (estimate) |
| `qv22` | QV table n <= 22 | 359 s in the original run (rows n = 21, 22 swap-affected) | 70.2M states; table about 1.75 GB, up to about 2.6 GB while growing (estimate). Above the 2 GB guideline: run only with that much memory free |
| `erasure14` | maximal erasure on n <= 14 | 851 s | 1.28 GB |
| `obl-grid` | all 59 oblivious instances | about 25 min | under 0.2 GB per instance |
| `closed-forms` | (10,5), (12,5), (13,5) | about 20 min, dominated by (13,5) | under 0.2 GB |
| `cube-compress` | compression sizes | under 10 s | up to 2.2 GB (scipy) |

The solver's values are independent of timing and of swap.

## Not reproduced by a script

* The Farkas vectors behind the certified upper bounds are not stored. They are recomputed and
  re-checked exactly by `spr obl solve`, which `recertify_obl.py` runs.
* Exponents of asymptotic families (for example 0.5438958 for the binomial cube of Dagan et
  al.) are evaluations of closed-form expressions, not certified computations.
