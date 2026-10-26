# Exact SPR values, first pass (25 Oct 2026)

Produced by spr_ref.py (memoised integer minimax, reflection symmetry). Rows n, columns s = 1..16.

## opt(n, s), full SPR (Labeler may erase any subset), n <= 13, s <= 10

```
n= 1   1  1  1  1  1  1  1  1  1  1
n= 2   1  1  1  1  1  1  1  1  1  1
n= 3   1  2  2  2  2  2  2  2  2  2
n= 4   1  2  2  2  2  2  2  2  2  2
n= 5   1  2  2  3  3  3  3  3  3  3
n= 6   1  2  2  3  3  3  3  3  3  3
n= 7   1  2  3  3  3  3  4  4  4  4
n= 8   1  2  3  3  3  3  4  4  4  4
n= 9   1  2  3  3  3  4  4  4  4  4
n=10   1  2  3  3  3  4  4  4  4  4
n=11   1  2  3  3  4  4  4  4  5  5
n=12   1  2  3  3  4  4  4  5  5  5
n=13   1  2  3  3  4  4  4  5  5  5
```

## SPR with forced maximal erasure, n <= 18, s <= 16
Identical to full SPR on every instance where both were computed (n <= 13, s <= 10).

```
n= 1   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
n= 2   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
n= 3   1  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
n= 4   1  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
n= 5   1  2  2  3  3  3  3  3  3  3  3  3  3  3  3  3
n= 6   1  2  2  3  3  3  3  3  3  3  3  3  3  3  3  3
n= 7   1  2  3  3  3  3  4  4  4  4  4  4  4  4  4  4
n= 8   1  2  3  3  3  3  4  4  4  4  4  4  4  4  4  4
n= 9   1  2  3  3  3  4  4  4  4  4  5  5  5  5  5  5
n=10   1  2  3  3  3  4  4  4  4  4  5  5  5  5  5  5
n=11   1  2  3  3  4  4  4  4  5  5  5  5  5  5  5  6
n=12   1  2  3  3  4  4  4  5  5  5  5  5  5  5  5  6
n=13   1  2  3  3  4  4  4  5  5  5  5  5  6  6  6  6
n=14   1  2  3  3  4  4  4  5  5  5  5  5  6  6  6  6
n=15   1  2  3  4  4  4  5  5  5  5  6  6  6  6  6  6
n=16   1  2  3  4  4  4  5  5  5  5  6  6  6  6  6  6
n=17   1  2  3  4  4  5  5  5  5  6  6  6  6  6  6  7
n=18   1  2  3  4  4  5  5  5  5  6  6  6  6  6  7  7
```

## opt_QV(n, s), Qiao-Valiant game (no reuse, automatic erasure), n <= 13, s <= 10

```
n= 1   1  1  1  1  1  1  1  1  1  1
n= 2   1  1  1  1  1  1  1  1  1  1
n= 3   1  2  2  2  2  2  2  2  2  2
n= 4   1  2  2  2  2  2  2  2  2  2
n= 5   1  2  2  2  2  2  2  2  2  2
n= 6   1  2  2  3  3  3  3  3  3  3
n= 7   1  2  3  3  3  3  3  3  3  3
n= 8   1  2  3  3  3  3  3  3  3  3
n= 9   1  2  3  3  3  3  3  3  3  3
n=10   1  2  3  3  3  4  4  4  4  4
n=11   1  2  3  3  4  4  4  4  4  4
n=12   1  2  3  3  4  4  4  4  4  4
n=13   1  2  3  3  4  4  4  4  4  4
```

## Observations (empirical, unproven)
1. Reproduces Qiao-Valiant Lemma 5: opt(3,2)=2, opt(7,3)=3, opt(15,4)=4; and opt(14,4)=3, so 2^t - 1 cells remains the threshold for keeping all t signs even with reuse, for t <= 4.
2. Partial erasure never helps the Labeler (full SPR = forced maximal erasure, n <= 13, s <= 10). Under forced erasure every board satisfies: every '+' lies left of every '-'.
3. Reuse strictly helps the Pointer: e.g. opt(5,4)=3 vs opt_QV(5,4)=2, opt(7,7)=4 vs 3, opt(11,9)=5 vs 4. opt >= opt_QV on every computed instance.
