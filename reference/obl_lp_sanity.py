"""Sanity check: optimal worst-case preservation probability eps*(n, s) of an
oblivious no-reuse Pointer (DDFGKO Def 5.1), by LP feasibility + bisection.
Floating point, exploratory only. Also the average-case value ov(n, s).
"""
import sys, math, itertools
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix

def tuples(n, s):
    return list(itertools.permutations(range(n), s))

def build(n, s):
    T = tuples(n, s)
    idx = {t: i for i, t in enumerate(T)}
    # prefixes of length 1..s-1
    rows = []  # (ext indices, up indices, down indices)
    pref = {}
    for i, t in enumerate(T):
        for L in range(1, s):
            h = t[:L]
            d = pref.setdefault(h, ([], [], []))
            d[0].append(i)
            last = h[-1]
            rest = t[L:]
            if all(x > last for x in rest):
                d[1].append(i)
            if all(x < last for x in rest):
                d[2].append(i)
    return T, pref

def feasible(T, pref, eps):
    m = len(pref) * 2
    A = lil_matrix((m, len(T)))
    r = 0
    for h, (ext, up, dn) in pref.items():
        for j in ext:
            A[r, j] += eps
            A[r + 1, j] += eps
        for j in up:
            A[r, j] -= 1
        for j in dn:
            A[r + 1, j] -= 1
        r += 2
    res = linprog(np.zeros(len(T)), A_ub=A.tocsr(), b_ub=np.zeros(m),
                  A_eq=np.ones((1, len(T))), b_eq=[1], bounds=(0, None), method="highs")
    return res.status == 0

def eps_star(n, s, iters=30):
    T, pref = build(n, s)
    lo, hi = 0.0, 1.0
    for _ in range(iters):
        mid = (lo + hi) / 2
        if feasible(T, pref, mid):
            lo = mid
        else:
            hi = mid
    return lo

def E(n, s, eps):
    if eps <= 0 or n * n * s <= 1:
        return float("nan")
    return math.log(eps * n * s) / math.log(n * n * s)

if __name__ == "__main__":
    N = int(sys.argv[1]); S = int(sys.argv[2])
    for n in range(2, N + 1):
        out = []
        for s in range(2, min(S, n) + 1):
            e = eps_star(n, s)
            out.append(f"s={s}: eps={e:.5f} E={E(n, s, e):.4f}")
        print(f"n={n:2d}  " + " | ".join(out), flush=True)
