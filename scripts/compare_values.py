"""Compare a freshly computed value CSV against a committed one.

Usage: python3 scripts/compare_values.py NEW.csv COMMITTED.csv [--where COL=VAL ...]

Both files have columns n, s, value (extra columns are ignored). Each --where COL=VAL keeps
only the rows whose column COL, where the file has it, equals VAL (for example variant=full in
values_small_long_v1.csv, or policy=binary-search and labeler=full in policy_values_v1.csv).
Every (n, s) in NEW must appear in COMMITTED with the same value. Exits non-zero otherwise.
"""

import csv
import sys


def load(path, where):
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            if any(c in r and r[c] != v for c, v in where):
                continue
            out[(int(r["n"]), int(r["s"]))] = int(r["value"])
    return out


def main():
    argv = sys.argv[1:]
    where, args = [], []
    i = 0
    while i < len(argv):
        if argv[i] == "--where":
            c, v = argv[i + 1].split("=", 1)
            where.append((c, v))
            i += 2
        else:
            args.append(argv[i])
            i += 1
    new, old = load(args[0], where), load(args[1], where)
    missing = [k for k in new if k not in old]
    diff = [(k, new[k], old[k]) for k in new if k in old and new[k] != old[k]]
    if not new:
        print(f"{args[0]}: no rows")
        return 1
    if missing or diff:
        print(f"{args[0]} vs {args[1]}: {len(diff)} differences, {len(missing)} missing; "
              f"first: {diff[:3]} {missing[:3]}")
        return 1
    print(f"{args[0]} vs {args[1]}: all {len(new)} cells equal")
    return 0


if __name__ == "__main__":
    sys.exit(main())
