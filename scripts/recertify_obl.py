"""Re-run the oblivious certification and compare with data/oblivious_v1.csv.

Usage: python3 scripts/recertify_obl.py [--max-seconds X] [--spr PATH]

For every instance of data/oblivious_v1.csv whose recorded solve time is at most X seconds
(default 10; use a large X for the whole grid, about 25 minutes), runs `spr obl solve n s`
and `spr obl ov n s`. The solver checks every certificate exactly itself; this script then
checks that the new certified bracket for eps* overlaps the stored one (new lower <= stored
upper and stored lower <= new upper) and that ov* is certified exact and equal to the
stored value. The Farkas vectors are recomputed, not read from disk: they are not stored.
"""

import csv
import os
import re
import subprocess
import sys
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def arg(name, default):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default


def main():
    max_s = float(arg("--max-seconds", "10"))
    spr = arg("--spr", os.path.join(ROOT, "target", "release", "spr"))
    rows = list(csv.DictReader(open(os.path.join(ROOT, "data", "oblivious_v1.csv"))))
    ok, count = True, 0
    for r in rows:
        if float(r["seconds"]) > max_s:
            continue
        n, s = r["n"], r["s"]
        out = subprocess.run([spr, "obl", "solve", n, s], capture_output=True, text=True).stdout
        lo = re.search(r"certified lower bound: (\S+)", out)
        up = re.search(r"upper bound: (\S+) \(", out)
        if not lo or not up or "checked exactly" not in out:
            print(f"({n},{s}): could not parse solver output\n{out}")
            ok = False
            continue
        new_lo, new_up = Fraction(lo.group(1)), Fraction(up.group(1))
        old_lo, old_up = Fraction(r["eps_lower"]), Fraction(r["eps_upper"])
        good = new_lo <= new_up and new_lo <= old_up and old_lo <= new_up
        ov = subprocess.run([spr, "obl", "ov", n, s], capture_output=True, text=True).stdout
        m = re.search(r"certified \[(\S+), (\S+)\] EXACT", ov)
        good_ov = bool(m) and Fraction(m.group(1)) == Fraction(m.group(2)) == Fraction(r["ov_star"])
        print(f"({n},{s}): eps bracket [{float(new_lo):.10f}, {float(new_up):.10f}] "
              f"vs stored [{float(old_lo):.10f}, {float(old_up):.10f}] {'OK' if good else 'FAIL'}; "
              f"ov* {m.group(1) if m else '?'} {'OK' if good_ov else 'FAIL'}")
        ok &= good and good_ov
        count += 1
    print(f"recertify_obl: {count} instances, {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
