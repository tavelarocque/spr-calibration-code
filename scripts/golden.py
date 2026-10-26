"""The golden coin-driven oblivious strategy, generalising the optima at (6,3), (11,4).

Usage: python3 scripts/golden.py s [n_max] [--save path.json]

Structure read off the exact optima:
  * A fair coin picks the mirror image; then steps 2..s are driven by i.i.d. coins,
    C with probability 1/phi and J with probability 1/phi^2 (phi = golden ratio).
  * Every node (coin prefix) has an outward direction. The root's outward direction is
    "down" (the mirror coin flips everything). A C-move goes to the outward side of the
    current cell and keeps the direction; a J-move goes to the inward side and reverses
    it (the new outward direction points away from the cell it came from).
  * Placement constraints, for every earlier cell x_i:
      - if the step after i is J, all later cells lie on the inward side of x_i;
      - if the two steps after i are C, C, all later cells lie on the outward side;
      - if i = s - 1 and the last step is C, the last cell lies on the outward side.
    Each sign then survives on its J side with probability 1/phi^2 and on its C side with
    probability 1/phi^2 (1/phi for the second-to-last cell), so eps = 1/phi^2 = (3 - sqrt 5)/2.

This script finds the smallest n on which such a placement exists, as a small integer
program (scipy.optimize.milp, which calls HiGHS): one integer cell per node, the strict
side constraints above, and for every other ancestor pair on a common path a binary
choice of side (cells on one path must be distinct). The optimum is then verified
EXACTLY in Q(sqrt 5); the MILP itself is only a search device.
"""

import json
import sys
from fractions import Fraction


class QS5:
    """a + b sqrt(5) with rational a, b; exact arithmetic and comparison."""

    __slots__ = ("a", "b")

    def __init__(self, a, b=0):
        self.a, self.b = Fraction(a), Fraction(b)

    def __add__(self, o):
        o = o if isinstance(o, QS5) else QS5(o)
        return QS5(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __sub__(self, o):
        o = o if isinstance(o, QS5) else QS5(o)
        return QS5(self.a - o.a, self.b - o.b)

    def __mul__(self, o):
        o = o if isinstance(o, QS5) else QS5(o)
        return QS5(self.a * o.a + 5 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def inv(self):
        d = self.a * self.a - 5 * self.b * self.b
        return QS5(self.a / d, -self.b / d)

    def __truediv__(self, o):
        o = o if isinstance(o, QS5) else QS5(o)
        return self * o.inv()

    def sign(self):
        # sign of a + b sqrt5, exactly
        a, b = self.a, self.b
        if b == 0:
            return (a > 0) - (a < 0)
        if a == 0:
            return (b > 0) - (b < 0)
        if (a > 0) == (b > 0):
            return 1 if a > 0 else -1
        # opposite signs: compare a^2 with 5 b^2
        c = a * a - 5 * b * b
        s = (c > 0) - (c < 0)
        return s if a > 0 else -s

    def __lt__(self, o):
        return (self - o).sign() < 0

    def __le__(self, o):
        return (self - o).sign() <= 0

    def __eq__(self, o):
        o = o if isinstance(o, QS5) else QS5(o)
        return self.a == o.a and self.b == o.b

    def __repr__(self):
        return f"{self.a} + {self.b}*sqrt5"

    def __float__(self):
        return float(self.a) + float(self.b) * 5 ** 0.5


PHI_INV = QS5(Fraction(-1, 2), Fraction(1, 2))  # (sqrt5 - 1)/2
PHI_INV2 = QS5(Fraction(3, 2), Fraction(-1, 2))  # (3 - sqrt5)/2


def nodes_for(s):
    """Coin prefixes (strings over 'CJ') for steps 2..s, root = ''."""
    out = [""]
    for depth in range(1, s):
        out += [p + c for p in out if len(p) == depth - 1 for c in "CJ"]
    return out


def outward(prefix):
    """Outward direction (-1 down, +1 up) of the node reached by `prefix`."""
    d = -1
    for c in prefix:
        if c == "J":
            d = -d
    return d


def constraints(prefix):
    """For the node at `prefix` (its own cell), the list of (ancestor_prefix, side)
    meaning: this cell must lie strictly on `side` (+1 above / -1 below) of the
    ancestor's cell."""
    out = []
    k = len(prefix)
    for i in range(k):  # ancestor at prefix[:i]; the step after it is prefix[i]
        a = prefix[:i]
        oa = outward(a)
        rest = prefix[i:]
        if rest[0] == "J":
            out.append((a, -oa))
        elif len(rest) >= 2 and rest[1] == "C":
            out.append((a, oa))
        elif len(rest) == 1:  # the C-child itself lies outward
            out.append((a, oa))
    return out


def search(s, n):
    """Backtracking placement of all nodes on cells 0..n-1 (one mirror half)."""
    order = nodes_for(s)
    cons = {p: constraints(p) for p in order}
    cell = {}

    def rec(idx):
        if idx == len(order):
            return True
        p = order[idx]
        lo, hi = 0, n - 1
        for a, side in cons[p]:
            if side > 0:
                lo = max(lo, cell[a] + 1)
            else:
                hi = min(hi, cell[a] - 1)
        used = {cell[p[:i]] for i in range(len(p))}
        # try cells nearest the middle of the allowed interval first
        mid = (lo + hi) / 2
        for x in sorted(range(lo, hi + 1), key=lambda v: (abs(v - mid), v)):
            if x in used:
                continue
            cell[p] = x
            if rec(idx + 1):
                return True
            del cell[p]
        return False

    return dict(cell) if rec(0) else None


def search_milp(s):
    """Minimal n and a placement, by MILP. Returns (n, cell)."""
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp

    order = nodes_for(s)
    idx = {p: i for i, p in enumerate(order)}
    V = len(order)
    pairs = []  # (v, a) ancestor pairs with no side constraint
    sided = []  # (v, a, side)
    for p in order:
        cons = dict(constraints(p))
        for i in range(len(p)):
            a = p[:i]
            if a in cons:
                sided.append((idx[p], idx[a], cons[a]))
            else:
                pairs.append((idx[p], idx[a]))
    nvar = V + 1 + len(pairs)  # cells, N (max cell), one binary per free pair
    big = 4 * 2 ** s
    rows, lo, hi = [], [], []

    def row(coefs, l, h):
        r = np.zeros(nvar)
        for j, c in coefs:
            r[j] += c
        rows.append(r)
        lo.append(l)
        hi.append(h)

    for v, a, side in sided:  # side * (x_v - x_a) >= 1
        row([(v, side), (a, -side)], 1, np.inf)
    for k, (v, a) in enumerate(pairs):
        z = V + 1 + k  # z = 1: x_v > x_a, z = 0: x_v < x_a
        row([(v, 1), (a, -1), (z, -big)], 1 - big, np.inf)  # x_v - x_a >= 1 - big(1 - z)
        row([(v, 1), (a, -1), (z, -big)], -np.inf, -1)      # x_v - x_a <= -1 + big z
    for v in range(V):
        row([(V, 1), (v, -1)], 0, np.inf)  # N >= x_v
    c = np.zeros(nvar)
    c[V] = 1
    integrality = np.ones(nvar)
    lb = np.zeros(nvar)
    ub = np.concatenate([np.full(V + 1, big), np.ones(len(pairs))])
    res = milp(c, constraints=LinearConstraint(np.array(rows), lo, hi), integrality=integrality,
               bounds=Bounds(lb, ub))
    if not res.success:
        raise RuntimeError(res.message)
    x = np.round(res.x).astype(int)
    cell = {p: int(x[idx[p]]) for p in order}
    return int(x[V]) + 1, cell


def strategy(s, n, cell):
    """All tuples with exact masses in Q(sqrt5), both mirror halves."""
    leaves = [p for p in nodes_for(s) if len(p) == s - 1]
    tuples = []
    for p in leaves:
        t = tuple(cell[p[:i]] for i in range(s))
        m = QS5(Fraction(1, 2))
        for c in p:
            m = m * (PHI_INV if c == "C" else PHI_INV2)
        tuples.append((t, m))
        tuples.append((tuple(n - 1 - x for x in t), m))
    return tuples


def exact_eps(s, tuples):
    total = QS5(0)
    for _, m in tuples:
        total = total + m
    assert total == QS5(1), total
    assert all(len(set(t)) == s for t, _ in tuples)
    pref = {}
    for t, m in tuples:
        for i in range(1, s):
            h, last = t[:i], t[i - 1]
            mu, u, d = pref.get(h, (QS5(0), QS5(0), QS5(0)))
            up = all(c > last for c in t[i:])
            down = all(c < last for c in t[i:])
            pref[h] = (mu + m, u + (m if up else QS5(0)), d + (m if down else QS5(0)))
    eps = None
    for mu, u, d in pref.values():
        for v in (u / mu, d / mu):
            if eps is None or v < eps:
                eps = v
    return eps


def main():
    s = int(sys.argv[1])
    n, cell = search_milp(s)
    tuples = strategy(s, n, cell)
    eps = exact_eps(s, tuples)
    print(f"s={s}: smallest n with a golden placement (MILP optimum): n={n}; exact eps = {eps} "
          f"(= (3 - sqrt5)/2: {eps == PHI_INV2}); support {len(tuples)} tuples")
    for p in sorted(cell, key=lambda q: (len(q), q)):
        print(f"  {p or '(root)':8s} -> cell {cell[p]}")
    if "--save" in sys.argv:
        path = sys.argv[sys.argv.index("--save") + 1]
        json.dump({"n": n, "s": s, "eps": "(3 - sqrt5)/2",
                   "tuples": [[list(t), [str(m.a), str(m.b)]] for t, m in tuples],
                   "masses": "a + b*sqrt5 as [a, b]"}, open(path, "w"))
    return 0 if eps == PHI_INV2 else 1


if __name__ == "__main__":
    sys.exit(main())
