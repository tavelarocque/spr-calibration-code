"""Exact reference solver for Sign Preservation with Reuse (SPR) and the
Qiao-Valiant restriction (QV). Brute force, memoised, small n only.
This file is the correctness oracle for the Rust solver. Do not optimise it.
Cell codes: 0 empty, 1 '+', 2 '-', 3 empty-and-burned (QV only, cannot be reused).
"""
import sys
from functools import lru_cache
from itertools import combinations

sys.setrecursionlimit(100000)

def count(x):
    return sum(1 for c in x if c in (1, 2))

def canon(x):
    # reflection symmetry: reverse strip and swap + / -
    sw = {0: 0, 1: 2, 2: 1, 3: 3}
    y = tuple(sw[c] for c in reversed(x))
    return min(x, y)

def erasable(x, j):
    return [i for i in range(len(x)) if (i < j and x[i] == 2) or (i > j and x[i] == 1)]

def subsets(lst, mode):
    if mode == "any":
        for k in range(len(lst) + 1):
            for c in combinations(lst, k):
                yield c
    elif mode == "allornone":
        yield ()
        if lst:
            yield tuple(lst)
    elif mode == "all":
        yield tuple(lst)

def make_solver(game, lmode):
    """game: 'spr' (reuse allowed) or 'qv' (no reuse, erased cells burned).
    lmode: Labeler erasure options: 'any', 'allornone', 'all'."""
    @lru_cache(maxsize=None)
    def V(x, r):
        best = count(x)  # Pointer may stop
        if r == 0:
            return best
        for j in range(len(x)):
            if x[j] != 0:
                continue
            worst = None
            er = erasable(x, j)
            for E in subsets(er, lmode):
                y = list(x)
                for i in E:
                    y[i] = 3 if game == "qv" else 0
                for sgn in (1, 2):
                    y[j] = sgn
                    v = V(canon(tuple(y)), r - 1)
                    if worst is None or v < worst:
                        worst = v
                        if worst <= best:
                            break
                if worst is not None and worst <= best:
                    break
            if worst > best:
                best = worst
        return best
    return lambda n, s: V(tuple([0] * n), s)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    S = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    variants = {
        "SPR": make_solver("spr", "any"),
        "SPR_allornone": make_solver("spr", "allornone"),
        "SPR_forced": make_solver("spr", "all"),
        "QV": make_solver("qv", "all"),
    }
    for name, f in variants.items():
        print(f"== {name}: rows n, columns s = 1..{S}")
        for n in range(1, N + 1):
            print(f"n={n:2d}  " + " ".join(f"{f(n, s):2d}" for s in range(1, S + 1)), flush=True)
