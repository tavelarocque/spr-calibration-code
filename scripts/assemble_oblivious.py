"""Assemble data/oblivious_v1.csv from a `spr obl grid` CSV and its basis exports.

Usage: python3 scripts/assemble_oblivious.py grid.csv bases_dir out.csv

For every instance with eps* > 0 the closed-form step (scripts/closed_form.py) is run
on the exported basis. A closed form is recorded only if the exact checks pass: the
basis strategy is valid at the root, and the root lies between the certified rational
lower and upper bounds. eps* = 1/2 instances take the closed form 1/2 directly, since
both certified bounds equal 1/2 exactly.
"""

import csv
import math
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(__file__))
from closed_form import analyse  # noqa: E402

TIME_LIMIT = 180  # seconds per instance for the symbolic step


def _worker(path, q):
    try:
        q.put(analyse(path, verbose=False))
    except Exception as e:
        q.put({"ok": False, "reason": repr(e)[:120]})


def analyse_limited(path):
    """analyse() in a child process; gives up after TIME_LIMIT seconds."""
    import multiprocessing as mp
    q = mp.Queue()
    p = mp.Process(target=_worker, args=(path, q))
    p.start()
    p.join(TIME_LIMIT)
    if p.is_alive():
        p.terminate()
        p.join()
        return {"ok": False, "reason": f"not found within {TIME_LIMIT} s"}
    return q.get() if not q.empty() else {"ok": False, "reason": "no result"}


def main(grid, bases, out):
    if os.path.exists(out):
        print(f"refusing to overwrite {out}")
        return 1
    rows = list(csv.DictReader(open(grid)))
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["n", "s", "eps_lower", "eps_upper", "upper_certified", "eps_closed_form",
                     "closed_form_minpoly", "E", "ov_star", "ov_certified", "seconds", "note"])
        for r in rows:
            n, s = int(r["n"]), int(r["s"])
            lo, hi = Fraction(r["eps_lower"]), Fraction(r["eps_upper"])
            closed, minpoly, note = "", "", ""
            if lo == hi:
                closed = str(lo)
            elif lo > 0:
                path = os.path.join(bases, f"basis_n{n:02d}_s{s}.json")
                if os.path.exists(path):
                    res = analyse_limited(path)
                    if res.get("ok"):
                        closed, minpoly = res["closed"], res["minpoly"]
                    else:
                        note = f"closed form: {res.get('reason', 'check failed')}"
                    print(f"({n},{s}): {res.get('closed', '-')} [{res.get('minpoly', res.get('reason', ''))}] ok={res.get('ok')}", flush=True)
            E = math.log(float(lo) * n * s) / math.log(n * n * s) if lo > 0 else float("-inf")
            ov_cert = r["ov_lower"] == r["ov_upper"] and r["ov_lower"] != ""
            ov = r["ov_lower"] if ov_cert else f"[{r['ov_lower']}, {r['ov_upper']}]"
            secs = float(r["eps_seconds"]) + float(r["ov_seconds"])
            w.writerow([n, s, r["eps_lower"], r["eps_upper"], r["upper_certified"], closed, minpoly,
                        f"{E:.6f}", ov, ov_cert, f"{secs:.3f}", note])
    return 0


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.set_start_method("fork")
    sys.exit(main(*sys.argv[1:4]))
