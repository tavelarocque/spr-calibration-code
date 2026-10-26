"""Closed form of eps*(n, s) from an exported LP basis.

Usage: python3 scripts/closed_form.py basis.json
where basis.json comes from `spr obl active-set n s --out basis.json`.

Method. The export holds the LP basis whose exact basic solution gave the certified
lower bound: the basic structural columns S (reflection orbits) and the nonbasic,
tight rows R (U(h) - eps mu(h) = 0, plus the normalisation row), with |R| = |S|. This
script never forces |R| = |S| from the support: it takes the solver's final basis.

1. Solve M(eps) x = e_sum symbolically in eps (sympy), so every orbit mass is a
   rational function of eps.
2. Expand to the strategy pi(eps) over all tuples (both reflection members) and form
   every constraint as a rational function: each mass >= 0, and U(h) - eps mu(h) >= 0,
   D(h) - eps mu(h) >= 0 for every prefix of the support.
3. The basis stays feasible up to the first eps at which a constraint function
   crosses zero. That critical eps is a real root of one of their numerators; it is
   located in the certified bracket, and its minimal polynomial is reported.
4. Exact check at the root: pi(eps*) has total mass 1, all masses >= 0 and every
   constraint >= 0, with at least one equal to 0. So eps(pi(eps*)) = eps* exactly, and
   eps* is a certified lower bound in closed form. The root is also compared exactly
   with the certified rational bounds from the export.

Equality eps* = closed form is then established from below exactly and from above to
within the width of the certified rational bracket. It is not a proof of optimality.
"""

import json
import sys

import sympy as sp


def rat(text):
    if "/" in text:
        a, b = text.split("/")
        return sp.Rational(int(a), int(b))
    return sp.Rational(int(text))


def sign_at_root(num, minpoly, interval):
    """Exact sign of the polynomial num at the root of minpoly isolated in interval.

    Zero iff minpoly divides num. Otherwise the isolating interval (with rational
    endpoints) is bisected until num has no root in it; num then has constant sign on
    it, read off exactly at a rational point."""
    e = minpoly.gen
    if num.rem(minpoly).is_zero:
        return 0
    a, b = interval
    for _ in range(200):
        if num.count_roots(a, b) == 0:
            mid = (a + b) / 2
            return 1 if num.eval(mid) > 0 else -1
        mid = (a + b) / 2
        if minpoly.eval(mid) == 0:
            return 1 if num.eval(mid) > 0 else -1
        # keep the half containing the root of minpoly
        if minpoly.count_roots(a, mid) > 0:
            b = mid
        else:
            a = mid
    raise RuntimeError("could not separate roots")


def analyse(path, verbose=True):
    """Returns a dict with the closed form, its minimal polynomial and the checks."""
    from sympy.polys.matrices import DomainMatrix
    log = print if verbose else (lambda *a, **k: None)
    d = json.load(open(path))
    n, s, sum_row = d["n"], d["s"], d["sum_row"]
    lower, upper = rat(d["lower"]), rat(d["upper"])
    cols = [c for c in d["cols"] if c["basic"]]
    rows = d["nonbasic_rows"]
    k = len(cols)
    if k != len(rows):
        log(f"basis is not square ({k} columns, {len(rows)} rows); cannot proceed")
        return {"ok": False, "reason": "basis not square"}
    e = sp.Symbol("eps")
    K = sp.QQ.frac_field(e)
    pos = {r: i for i, r in enumerate(rows)}
    entries = [[K.zero] * k for _ in range(k)]
    for j, c in enumerate(cols):
        for row, ups, mem in c["entries"]:
            if row in pos:
                entries[pos[row]][j] = K.from_sympy(ups - e * mem)
        if sum_row in pos:
            entries[pos[sum_row]][j] = K.from_sympy(sp.Integer(2))
    rhs = [[K.zero] for _ in range(k)]
    rhs[pos[sum_row]][0] = K.one
    A = DomainMatrix(entries, (k, k), K)
    bvec = DomainMatrix(rhs, (k, 1), K)
    sol = A.lu_solve(bvec).to_Matrix()
    x = [sp.cancel(sol[j, 0]) for j in range(k)]

    def refl(t):
        return tuple(n - 1 - c for c in t)

    masses = {}
    for j, c in enumerate(cols):
        t = tuple(c["rep"])
        masses[t] = x[j]
        masses[refl(t)] = x[j]
    cons = [(f"mass{t}", m) for t, m in masses.items()]
    pref = {}
    for t, m in masses.items():
        for i in range(1, s):
            h, last = t[:i], t[i - 1]
            mu, u, dn = pref.get(h, (0, 0, 0))
            up = all(c > last for c in t[i:])
            down = all(c < last for c in t[i:])
            pref[h] = (mu + m, u + (m if up else 0), dn + (m if down else 0))
    for h, (mu, u, dn) in pref.items():
        cons.append((f"U{h}", sp.cancel(u - e * mu)))
        cons.append((f"D{h}", sp.cancel(dn - e * mu)))
    total = sp.cancel(sum(masses.values()))
    assert total == 1, f"masses sum to {total}"

    # Numerator and denominator polynomials of every constraint.
    polys = []
    for name, f in cons:
        num, den = sp.fraction(sp.together(f))
        polys.append((name, sp.Poly(sp.expand(num), e), sp.Poly(sp.expand(den), e)))

    lo_w, hi_w = lower - sp.Rational(1, 10**6), upper + sp.Rational(1, 10**6)
    cands = []
    for name, num, _ in polys:
        if num.is_zero or num.degree() < 1:
            continue
        for fac, _ in num.factor_list()[1]:
            for (a, b), _mult in fac.intervals(eps=sp.Rational(1, 10**12)):
                if (a <= hi_w) and (b >= lo_w):
                    cands.append(((a + b) / 2, fac, (a, b), name))
    if not cands:
        log("no constraint changes sign near the certified bracket")
        return {"ok": False, "reason": "no critical root near the bracket"}
    _, minpoly, interval, name = min(cands, key=lambda c: c[0])
    log(f"({n}, {s}): basis of {k} columns ({d.get('basis_source', '?')})")
    log(f"  critical constraint {name}; minimal polynomial {minpoly.as_expr()} = 0")
    roots = [r for r in sp.solve(minpoly.as_expr(), e) if r.is_real and interval[0] <= r <= interval[1]]
    closed = sp.simplify(roots[0]) if roots else sp.CRootOf(minpoly.as_expr(), 0)
    log(f"  eps* closed form: {closed}  ~ {sp.N(closed, 25)}")

    ok, tight = True, 0
    for name, num, den in polys:
        if num.is_zero:
            continue
        sn = sign_at_root(num, minpoly, interval)
        sd = sign_at_root(den, minpoly, interval)
        if sn * sd < 0:
            ok = False
            log(f"  constraint {name} is negative at eps*")
        if sn == 0 and not name.startswith("mass"):
            tight += 1
    in_bracket = bool(lower <= closed) and bool(closed <= upper)
    log(f"  exact check: all constraints >= 0 at eps*: {ok}; tight ratio constraints: {tight}")
    log(f"  certified bracket [{float(lower):.15f}, {float(upper):.15f}]; lower <= closed <= upper: {in_bracket}")
    log("  strategy at eps* (one member of each reflection pair, both have this mass):")
    for j, c in enumerate(cols):
        v = sp.radsimp(sp.simplify(x[j].subs(e, closed)))
        if v != 0:
            log(f"    {tuple(c['rep'])}: {v}  ~ {sp.N(v, 12)}")
    return {"ok": ok and in_bracket, "closed": str(closed), "minpoly": str(minpoly.as_expr()),
            "degree": minpoly.degree(), "in_bracket": in_bracket, "constraints_ok": ok,
            "basis_columns": k, "basis_source": d.get("basis_source", "?")}


def main(path):
    r = analyse(path)
    return 0 if r.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
