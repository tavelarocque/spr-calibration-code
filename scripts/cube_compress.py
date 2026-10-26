"""Fewest cells for the cube's survival events, fair or golden (a search device).

Usage: python3 scripts/cube_compress.py d k {fair|golden}

The cube (scripts/golden_cube.py) fixes, for every sign t, survival events in terms of
coins: fair, "all later above" iff xi_{u_j} = + for every zero j of t's label, "below" iff
all are -; golden, "below" iff g_{u_j} = J for all j, "above" iff g_{u_j} = C and
g_{u_j.1} = C (second coin where it exists) for all j. Here the coin-to-cell map is not
fixed by coordinates. Instead:

* cell(t) is an integer that may depend only on coins that belong to an event of some
  earlier sign and to no event of sign t or any later sign (so, conditioned on the
  cells up to time t'', the coins of sign t'''s events are still fresh);
* for every t'' < t and every joint coin assignment under which sign t'''s "above"
  ("below") event holds, cell(t) > cell(t'') (cell(t) < cell(t''));
* cells of different times are always distinct (they co-occur in some tuple).

A MILP (scipy.optimize.milp, HiGHS) minimises the number of cells. The result is then
expanded over all coin outcomes and verified EXACTLY (Fractions, or Q(sqrt 5)). The MILP
is only a search device; only the exact check counts.
"""

import itertools
import math
import sys
from fractions import Fraction

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

sys.path.insert(0, __import__("os").path.dirname(__file__))
from golden import PHI_INV, PHI_INV2, QS5, exact_eps  # noqa: E402
from golden_cube import eps_fair, labels  # noqa: E402


def events(W, mode):
    """For each time t: (coins in its events, above(assignment), below(assignment))."""
    d = len(W[0])
    trie1 = {w[:m] for w in W for m in range(d) if w[m] == 1}  # nodes with a 1-child
    out = []
    for t, w in enumerate(W):
        later = W[t + 1:]
        levels = []
        for l in range(d):
            if w[l] != 0:
                continue
            u = w[:l]
            if not any(x[:l] == u and x[l] == 1 for x in later):
                continue  # vacuous level: no later label branches here
            second = u + (1,) if (u + (1,)) in trie1 else None
            levels.append((u, second))
        coins = set()
        for u, sec in levels:
            coins.add(u)
            if mode == "golden" and sec is not None:
                coins.add(sec)
        out.append((levels, coins))
    return out


def solve(d, k, mode):
    W = labels(d, k)
    s = len(W)
    ev = events(W, mode)
    ev_coins = [c for _, c in ev]
    D = []
    for t in range(s):
        earlier = set().union(*ev_coins[:t]) if t else set()
        later = set().union(*ev_coins[t:])
        D.append(sorted(earlier - later))
    # variables: one integer per (t, assignment of D[t])
    var = {}
    for t in range(s):
        for a in itertools.product((0, 1), repeat=len(D[t])):
            var[(t, a)] = len(var)
    V = len(var)

    def above(levels, asg):
        return all(asg[u] == 1 and (mode == "fair" or sec is None or asg[sec] == 1) for u, sec in levels)

    def below(levels, asg):
        return all(asg[u] == 0 for u, _ in levels)

    rows, lo, hi = [], [], []
    ordered = set()
    pairs = []
    for t2 in range(s):
        for t in range(t2 + 1, s):
            coins = sorted(set(D[t2]) | set(D[t]) | ev_coins[t2])
            for bits in itertools.product((0, 1), repeat=len(coins)):
                asg = dict(zip(coins, bits))
                a2 = tuple(asg[c] for c in D[t2])
                a = tuple(asg[c] for c in D[t])
                i2, i = var[(t2, a2)], var[(t, a)]
                if above(ev[t2][0], asg):
                    rows.append((i, i2, 1, np.inf))  # x_i - x_i2 >= 1
                    ordered.add((i, i2))
                elif below(ev[t2][0], asg):
                    rows.append((i, i2, -np.inf, -1))
                    ordered.add((i, i2))
                pairs.append((i, i2))
    free = sorted({p for p in pairs if p not in ordered})
    big = 4 * V + 8
    nvar = V + 1 + len(free)
    A, L, H = [], [], []

    def add(coefs, l, h):
        r = np.zeros(nvar)
        for j, c in coefs:
            r[j] += c
        A.append(r)
        L.append(l)
        H.append(h)

    for i, i2, l, h in set(rows):
        add([(i, 1), (i2, -1)], l, h)
    for z, (i, i2) in enumerate(free):
        zi = V + 1 + z
        add([(i, 1), (i2, -1), (zi, -big)], 1 - big, np.inf)
        add([(i, 1), (i2, -1), (zi, -big)], -np.inf, -1)
    for j in range(V):
        add([(V, 1), (j, -1)], 0, np.inf)
    c = np.zeros(nvar)
    c[V] = 1
    res = milp(c, constraints=LinearConstraint(np.array(A), L, H), integrality=np.ones(nvar),
               bounds=Bounds(np.zeros(nvar), np.concatenate([np.full(V + 1, big), np.ones(len(free))])))
    if not res.success:
        raise SystemExit(f"MILP failed: {res.message}")
    x = np.round(res.x).astype(int)
    n = int(x[V]) + 1
    # expand over all coins used anywhere and verify exactly
    allc = sorted(set().union(*ev_coins))
    tuples = {}
    for bits in itertools.product((0, 1), repeat=len(allc)):
        asg = dict(zip(allc, bits))
        tup = tuple(int(x[var[(t, tuple(asg[cc] for cc in D[t]))]]) for t in range(s))
        if mode == "fair":
            m = Fraction(1, 2 ** len(allc))
        else:
            m = QS5(1)
            for cc in allc:
                m = m * (PHI_INV if asg[cc] else PHI_INV2)
        tuples[tup] = tuples[tup] + m if tup in tuples else m
    tuples = list(tuples.items())
    assert all(len(set(t)) == s for t, _ in tuples), "cells repeat within a tuple"
    eps = eps_fair(s, tuples) if mode == "fair" else exact_eps(s, tuples)
    return n, s, eps, len(tuples), len(allc)


def main():
    d, k, mode = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    n, s, eps, sup, coins = solve(d, k, mode)
    ef = float(eps)
    cube_n = math.comb(d, k) * 2 ** (d - k)
    E = math.log(ef * n * s) / math.log(n * n * s)
    print(f"{mode:6s} d={d} k={k}: compressed n={n} (coordinate cube {cube_n}), s={s}, coins={coins}, "
          f"support={sup}, exact eps={eps} (~{ef:.6f}), E={E:.4f}")


if __name__ == "__main__":
    main()
