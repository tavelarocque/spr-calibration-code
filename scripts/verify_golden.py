"""Exact checks of the golden-ratio constructions (golden family and golden cube).

Usage: python3 scripts/verify_golden.py

1. For every data/strategies/golden_sS.json: the masses (in Q(sqrt 5)) are nonnegative and
   sum to 1, the tuples have distinct cells, the cell count is F(s+3) - 2, and eps is
   exactly (3 - sqrt 5)/2.
2. For eight small (d, k): the binomial cube has eps exactly 2^-k and the golden
   cube has eps exactly phi^-2k, by enumeration of every coin outcome.

Exits non-zero on any failure. Uses only the standard library.
"""

import glob
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from golden import PHI_INV2, QS5, exact_eps  # noqa: E402
import golden_cube  # noqa: E402

ROOT = os.path.dirname(HERE)


def fib(m):
    a, b = 0, 1
    for _ in range(m):
        a, b = b, a + b
    return a


def check_family():
    ok = True
    files = sorted(glob.glob(os.path.join(ROOT, "data", "strategies", "golden_s*.json")))
    if not files:
        print("no golden strategy files found")
        return False
    for path in files:
        d = json.load(open(path))
        n, s = d["n"], d["s"]
        tuples = [(tuple(t), QS5(Fraction(m[0]), Fraction(m[1]))) for t, m in d["tuples"]]
        good = all(m.sign() >= 0 for _, m in tuples)
        good &= all(len(set(t)) == s and all(0 <= c < n for c in t) for t, _ in tuples)
        good &= n == fib(s + 3) - 2
        eps = exact_eps(s, tuples)  # also asserts that the masses sum to 1
        good &= eps == PHI_INV2
        print(f"golden family s={s}: n={n} (F(s+3)-2={fib(s + 3) - 2}), eps={eps}: "
              f"{'OK' if good else 'FAIL'}")
        ok &= good
    return ok


def check_cubes():
    ok = True
    for d, k in [(2, 1), (3, 1), (4, 1), (3, 2), (4, 2), (4, 3), (5, 2), (5, 3)]:
        n, s, _, tuples = golden_cube.build(d, k, "fair")
        fair = golden_cube.eps_fair(s, tuples)
        n2, s2, _, tuples2 = golden_cube.build(d, k, "golden")
        gold = exact_eps(s2, tuples2)
        want = QS5(1)
        for _ in range(k):
            want = want * PHI_INV2
        good = fair == Fraction(1, 2 ** k) and gold == want and (n, s) == (n2, s2)
        print(f"cube d={d} k={k}: n={n} s={s} fair eps={fair}, golden eps={gold}: "
              f"{'OK' if good else 'FAIL'}")
        ok &= good
    return ok


if __name__ == "__main__":
    a = check_family()
    b = check_cubes()
    print("verify_golden:", "PASS" if a and b else "FAIL")
    sys.exit(0 if a and b else 1)
