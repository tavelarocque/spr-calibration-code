"""The binomial cube with fair coins (DDFGKO Lemma 6.1) and its golden variant.

Usage: python3 scripts/golden_cube.py d k [--save-fair path] [--save-golden path]

Both use the cube's time labels (binary strings w of length d with exactly k zeros, in
lexicographic order) and the cube's cells (strings in {-1,0,1}^d with exactly k zeros, in
lexicographic order, -1 < 0 < 1); a cell copies its label's zero pattern. They differ only
in the value written at a 1-position m of a label w. Let v = w_1 .. w_{m-1} be the prefix
before m (a node of the label trie).

  fair   (the cube):  value = xi_v, one fair sign per trie node.
  golden:             value = -1 if g_v = J;
                      if g_v = C: +1, except when w_{m+1} = 1, in which case the value is
                      +1 if g_{v.1} = C and -1 if g_{v.1} = J.
  The golden coins are i.i.d. with P(C) = 1/phi, P(J) = 1/phi^2.

Survival (golden). A sign at time t whose label has zeros at l_1 < ... < l_k (prefixes
u_j before them) ends up with all later cells below it if g_{u_j} = J for every j, and
with all later cells above it if g_{u_j} = C and g_{u_j.1} = C for every j (the second
coin only where the node u_j.1 has a 1-child). So eps >= min(phi^-2k, phi^-2k) = phi^-2k.
This script does not rely on that argument: it enumerates every coin outcome and
computes eps exactly (Fractions for fair, Q(sqrt 5) for golden).
"""

import json
import math
import sys
from fractions import Fraction
from itertools import combinations

sys.path.insert(0, __import__("os").path.dirname(__file__))
from golden import PHI_INV, PHI_INV2, QS5, exact_eps  # noqa: E402


def labels(d, k):
    out = []
    for zeros in combinations(range(d), k):
        out.append(tuple(0 if i in zeros else 1 for i in range(d)))
    return sorted(out)


def cells(d, k):
    out = []

    def rec(pre, z):
        if len(pre) == d:
            if z == k:
                out.append(tuple(pre))
            return
        for a in (-1, 0, 1):
            if a == 0 and z == k:
                continue
            rec(pre + [a], z + (a == 0))

    rec([], 0)
    return {c: i for i, c in enumerate(sorted(out))}


def build(d, k, mode):
    W = labels(d, k)
    index = cells(d, k)
    # coins: one per trie node with a 1-child
    nodes = sorted({w[:m] for w in W for m in range(d) if w[m] == 1})
    pos = {v: i for i, v in enumerate(nodes)}
    c = len(nodes)
    if c > 22:
        raise SystemExit(f"{c} coins: 2^{c} outcomes is too many to enumerate")
    out = {}
    for outcome in range(1 << c):
        bit = lambda v: (outcome >> pos[v]) & 1  # noqa: E731  (1 = '+' / 'C')
        tup = []
        for w in W:
            q = []
            for m in range(d):
                if w[m] == 0:
                    q.append(0)
                    continue
                v = w[:m]
                if mode == "fair":
                    q.append(1 if bit(v) else -1)
                else:
                    if not bit(v):  # J
                        q.append(-1)
                    elif m + 1 < d and w[m + 1] == 1:
                        q.append(1 if bit(v + (1,)) else -1)
                    else:
                        q.append(1)
            tup.append(index[tuple(q)])
        tup = tuple(tup)
        if mode == "fair":
            mass = Fraction(1, 1 << c)
        else:
            mass = QS5(1)
            for v in nodes:
                mass = mass * (PHI_INV if bit(v) else PHI_INV2)
        out[tup] = out[tup] + mass if tup in out else mass
    return len(index), len(W), c, list(out.items())


def eps_fair(s, tuples):
    pref = {}
    for t, m in tuples:
        for i in range(1, s):
            h, last = t[:i], t[i - 1]
            mu, u, dn = pref.get(h, (0, 0, 0))
            pref[h] = (mu + m, u + (m if all(x > last for x in t[i:]) else 0),
                       dn + (m if all(x < last for x in t[i:]) else 0))
    return min(min(u, dn) / mu for mu, u, dn in pref.values())


def E(n, s, eps):
    return math.log(eps * n * s) / math.log(n * n * s)


def main():
    d, k = int(sys.argv[1]), int(sys.argv[2])
    for mode in ("fair", "golden"):
        n, s, c, tuples = build(d, k, mode)
        if mode == "fair":
            eps = eps_fair(s, tuples)
            ef = float(eps)
            want = Fraction(1, 2 ** k)
        else:
            eps = exact_eps(s, tuples)
            ef = float(eps)
            want = QS5(1)
            for _ in range(k):
                want = want * PHI_INV2
        print(f"{mode:6s} d={d} k={k}: n={n} s={s} coins={c} support={len(tuples)} "
              f"exact eps={eps} (~{ef:.6f}); expected {'2^-k' if mode == 'fair' else 'phi^-2k'}"
              f" ~{float(want):.6f}, eps >= expected: {not (eps < want)}; E={E(n, s, ef):.4f}")
        key = f"--save-{mode}"
        if key in sys.argv:
            path = sys.argv[sys.argv.index(key) + 1]
            if mode == "fair":
                data = [[list(t), f"{m.numerator}/{m.denominator}"] for t, m in tuples]
            else:
                data = [[list(t), [str(m.a), str(m.b)]] for t, m in tuples]
            json.dump({"n": n, "s": s, "mode": mode, "tuples": data}, open(path, "w"))


if __name__ == "__main__":
    main()
