# spr-calibration-code

Code and data for exact computations on the sign preservation with reuse game.
Author: Tave La Rocque.

## The game and why it matters

In online binary prediction a forecaster announces a probability before each outcome, and its
calibration error measures how far its predictions are from the frequencies they claim. The
best known bounds on the optimal calibration error both pass through a two-player combinatorial
game introduced by Dagan, Daskalakis, Fishelson, Golowich, Kleinberg and Okoroafor ("Breaking
the T^{2/3} Barrier for Sequential Calibration", STOC 2025, arXiv:2406.13668), building on the
sign preservation game of Qiao and Valiant. In sign preservation with reuse (SPR), one player,
the Pointer, repeatedly names an empty cell of a row of n cells; the other, the Labeler, writes a
plus or a minus sign there and may first erase minus signs to the left of the cell and plus
signs to its right. The Pointer wants many signs to remain after s rounds, the Labeler wants
few. How the optimal number of surviving signs grows with n and s is tied to the optimal
calibration rate in both directions, and the strongest known lower bound comes from an
oblivious Pointer strategy whose signs each survive with a guaranteed probability.

## What the code computes

An exact memoised minimax solver in Rust computes the value of the game in three variants: the
full game (up to 18 cells), the game in which the Labeler must erase every erasable sign (up to
20 cells), and the Qiao-Valiant game without reuse (up to 22 cells), each for every number of
rounds up to twice the number of cells. The values are exact integers, and the solver's state
reductions are cross-checked against an unreduced brute force and an independent Python
reference implementation (`reference/spr_ref.py`).

For small boards, the best worst-case survival probability of an oblivious no-reuse Pointer is
found with a linear program (HiGHS) and then certified exactly. Each lower bound is the exact
value of an explicit stored strategy, computed in rational arithmetic by two independently
written verifiers, and each upper bound is an exact rational Farkas certificate. Closed forms
are derived symbolically where possible.

The code also computes the exact value of a fixed Pointer policy against a best-responding
Labeler, checks a gap lemma and the optimality of maximal erasure on every state of small
boards, and checks golden-ratio constructions exactly.

The results are in `data/`: value tables (`values_*_v1.csv`, `tables_v1.md`), certified
oblivious brackets and closed forms (`oblivious_v1.csv`), the certified strategies themselves
(`data/strategies/`), and the outputs of the other computations.

## Checking everything

You need a stable Rust toolchain (`rust-toolchain.toml`), CMake and a C++ compiler (the
`highs-sys` crate builds HiGHS from source), and Python 3 with sympy. Then run

```
scripts/verify_all.sh
```

from the repository root. It builds the workspace and runs the test suite (including the
brute-force and reference comparisons), checks the solver against the Python reference,
recomputes the value tables up to moderate size and compares them cell by cell with the
committed files, reruns the gap-lemma and maximal-erasure checks on reduced ranges, re-verifies
every stored oblivious strategy exactly with both verifiers, re-certifies the faster oblivious
instances and ten closed forms, checks the golden constructions exactly in Q(sqrt 5), and checks
several derived counts and statements against the data. It takes about five minutes and under
0.5 GB of memory on an 8 GB laptop, writes only to a temporary directory, and ends with
`verify_all: PASS` or a list of the checks that failed. The longer runs that produced the full
tables are opt-in through `scripts/verify_long.sh`; `REPRODUCE.md` lists every command, its
measured cost and the file it checks.

## Status

A write-up of the results is in preparation.

Parts of this code were developed with AI assistance (Claude, Anthropic).

## License

MIT; see `LICENSE`.
