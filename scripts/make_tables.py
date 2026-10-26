"""Check derived statements against the committed data and write LaTeX tables and counts.

Usage: python3 scripts/make_tables.py OUTDIR

Reads only data/*. Every game value is an exact integer read from the CSVs and every count
is computed from them here. The script fails (non-zero exit) if the data contradict any of
the statements it checks: the binary-search bound, the two oblivious thresholds
(eps* = 0 below 2s - 1 cells, eps* = 1/2 from 2^s - 1 cells), the first-move observation,
and the staircase consistency count. It writes LaTeX tables and a file of number macros to
OUTDIR and prints the counts. Exponents (logarithms) are printed rounded; they are derived from exact values
by a closed-form expression, never fitted.
"""

import csv
import math
import os
import sys
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(name):
    path = os.path.join(ROOT, "data", name)
    with open(path) as f:
        return {(int(r["n"]), int(r["s"])): int(r["value"]) for r in csv.DictReader(f)}


def thresholds(table, n):
    """List of (v, least s <= 2n with opt(n, s) >= v) for v = 1 .. opt(n, 2n)."""
    out = []
    v = 1
    for s in range(1, 2 * n + 1):
        while table[(n, s)] >= v:
            out.append((v, s))
            v += 1
    return out


def threshold_table(table, caption, label):
    maxn = max(n for n, _ in table)
    maxv = max(table[(n, 2 * n)] for n in range(1, maxn + 1))
    cols = "r|" + "r" * maxv
    lines = [
        "\\begin{table}[ht]",
        "\\centering\\footnotesize",
        "\\setlength{\\tabcolsep}{4pt}",
        f"\\begin{{tabular}}{{{cols}}}",
        "$n$ & " + " & ".join(f"$v={v}$" for v in range(1, maxv + 1)) + " \\\\ \\hline",
    ]
    for n in range(1, maxn + 1):
        th = dict(thresholds(table, n))
        cells = [str(th[v]) if v in th else "" for v in range(1, maxv + 1)]
        lines.append(f"{n} & " + " & ".join(cells) + " \\\\")
    lines += ["\\end{tabular}", f"\\caption{{{caption}}}", f"\\label{{{label}}}", "\\end{table}", ""]
    return "\n".join(lines)


def main_grid(full, qv, maxn, maxs):
    maxdiff = max(full[(n, s)] - qv[(n, s)] for n in range(1, maxn + 1)
                  for s in range(1, min(maxs, 2 * n) + 1))
    lines = [
        "\\begin{table}[ht]",
        "\\centering\\footnotesize",
        "\\setlength{\\tabcolsep}{3.5pt}",
        "\\begin{tabular}{r|" + "r" * maxs + "}",
        "$n \\backslash s$ & " + " & ".join(str(s) for s in range(1, maxs + 1)) + " \\\\ \\hline",
    ]
    for n in range(1, maxn + 1):
        cells = []
        for s in range(1, maxs + 1):
            if s > 2 * n:
                cells.append("")
                continue
            v, q = full[(n, s)], qv[(n, s)]
            cells.append(f"\\textbf{{{v}}}" if v > q else str(v))
        lines.append(f"{n} & " + " & ".join(cells) + " \\\\")
    lines += [
        "\\end{tabular}",
        "\\caption{Certified values of $\\opt_{\\mathrm{full}}(n,s)$ for $n \\le "
        f"{maxn}$, $s \\le \\min({maxs}, 2n)$. A bold entry exceeds $\\opt_{{\\QV}}(n,s)$, "
        f"by at most {maxdiff} in this range. "
        "In this range $\\opt_{\\mathrm{forced}} = \\opt_{\\mathrm{full}}$ cell by cell. "
        "Full tables: Appendix~\\ref{app:tables}.}",
        "\\label{tab:main}",
        "\\end{table}",
        "",
    ]
    return "\n".join(lines)


def latex_closed_form(r):
    cf = r["eps_closed_form"].strip()
    mp = r["closed_form_minpoly"].strip()
    if not cf:
        return "not found" if "not found" in r["note"] else ""
    rep = {
        "1/2": "$1/2$",
        "1/3": "$1/3$",
        "3/2 - sqrt(5)/2": "$(3-\\sqrt5)/2$",
        "2 - sqrt(3)": "$2-\\sqrt3$",
        "1 - sqrt(2)/2": "$1-\\sqrt2/2$",
        "4/7 - sqrt(2)/7": "$(4-\\sqrt2)/7$",
        "5/2 - sqrt(17)/2": "$(5-\\sqrt{17})/2$",
        "5/4 - sqrt(17)/4": "$(5-\\sqrt{17})/4$",
    }
    if cf in rep:
        return rep[cf]
    deg = max(int(t.split("**")[1].split()[0]) if "**" in t else 1
              for t in mp.replace("-", "+").split("+") if "eps" in t)
    return f"root of $p_{{{r['n']},{r['s']}}}$ (degree {deg})"


def oblivious_table():
    path = os.path.join(ROOT, "data", "oblivious_v1.csv")
    rows = list(csv.DictReader(open(path)))
    lines = [
        "\\begin{table}[ht]",
        "\\centering\\footnotesize",
        "\\begin{tabular}{rrlllll}",
        "$n$ & $s$ & certified lower & certified upper & closed form & $E$ & $\\ov^*$ \\\\ \\hline",
    ]
    for r in rows:
        n, s = int(r["n"]), int(r["s"])
        lo, up = Fraction(r["eps_lower"]), Fraction(r["eps_upper"])
        if lo == 0:
            continue  # summarised in the caption
        lo_s = f"{float(lo):.10f}" if lo.denominator > 100 else f"${lo}$"
        up_s = f"{float(up):.10f}" if up.denominator > 100 else f"${up}$"
        e = r["E"]
        e_s = f"{float(e):.4f}"
        ov = r["ov_star"]
        lines.append(
            f"{n} & {s} & {lo_s} & {up_s} & {latex_closed_form(r)} & {e_s} & ${ov}$ \\\\"
        )
    zero = [(int(r["n"]), int(r["s"])) for r in rows if Fraction(r["eps_lower"]) == 0]
    zs = ", ".join(f"$({n},{s})$" for n, s in zero)
    lines += [
        "\\end{tabular}",
        "\\caption{Certified oblivious optima (Certified computation). Decimal bounds are exact "
        "rationals truncated to ten digits (lower bounds rounded down, upper bounds rounded up in the "
        "data file). $E$ is evaluated at the certified lower bound. Omitted rows, where the certified "
        f"bracket is $[0, 1.0002\\cdot 10^{{-9}}]$ and $\\eps^* = 0$ by Lemma~\\ref{{lem:positive}}: {zs}. "
        "Here $p_{10,5}(e)=2e^6-10e^5+23e^4-29e^3+22e^2-8e+1$ and $p_{12,5}(e)=3e^4-8e^3+11e^2-6e+1$; "
        "``not found'' means the symbolic step did not finish within 180~s.}",
        "\\label{tab:obl}",
        "\\end{table}",
        "",
    ]
    return "\n".join(lines), len(rows), len(zero)


def tensor_rows():
    out = []
    for t in range(3, 10):
        n = 2 ** t - 1
        al = math.log(t) / math.log(n)
        be = math.log((t + 1) / 2) / math.log(n)
        out.append((t, n, al, be, (be + 1) / (al + 2), (2 * be + 1) / (al + 2 * be + 2)))
    return out


def tensor_table():
    lines = [
        "\\begin{table}[ht]",
        "\\centering\\footnotesize",
        "\\begin{tabular}{rrcccc}",
        "$t$ & base $(2^t-1, t)$ & $\\alpha$ & $\\beta$ & $(\\beta+1)/(\\alpha+2)$ & "
        "$(2\\beta+1)/(\\alpha+2\\beta+2)$ \\\\ \\hline",
    ]
    for t, n, al, be, p1, p2 in tensor_rows():
        lines.append(f"{t} & $({n},{t})$ & {al:.5f} & {be:.5f} & {p1:.5f} & {p2:.5f} \\\\")
    lines += [
        "\\end{tabular}",
        "\\caption{Exponents obtained by tensorizing the binary-search base $\\opt_{\\QV}(2^t-1,t)=t$ "
        "(Fact~\\ref{fact:qv5}(2)). Here $\\alpha=\\log t/\\log(2^t-1)$ and "
        "$\\beta=\\log((t+1)/2)/\\log(2^t-1)$. The fifth column is the map of "
        "Fact~\\ref{fact:ddfgko-42}, the sixth that of \\cite[Theorem~4]{QV21}. "
        "Values are rounded; they are closed-form expressions in $t$, not fits.}",
        "\\label{tab:tensor}",
        "\\end{table}",
        "",
    ]
    return "\n".join(lines)


def main():
    if len(sys.argv) < 2:
        raise SystemExit("usage: python3 scripts/make_tables.py OUTDIR")
    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)
    full = load("values_full_v1.csv")
    forced = load("values_forced_v1.csv")
    qv = load("values_qv_v1.csv")

    def w(name, text):
        with open(os.path.join(outdir, name), "w") as f:
            f.write(text)

    w("main_grid.tex", main_grid(full, qv, 18, 14))
    w("thr_full.tex", threshold_table(
        full,
        "$\\opt_{\\mathrm{full}}$ in threshold form (Certified computation). Row $n$, column $v$: "
        "the least $s \\le 2n$ with $\\opt_{\\mathrm{full}}(n,s) \\ge v$; blank if $v > "
        "\\opt_{\\mathrm{full}}(n,2n)$. Since $\\opt$ is non-decreasing in $s$, this determines every "
        "value with $s \\le 2n$.",
        "tab:thr-full"))
    w("thr_forced.tex", threshold_table(
        forced,
        "$\\opt_{\\mathrm{forced}}$ in threshold form (Certified computation), as in "
        "Table~\\ref{tab:thr-full}.",
        "tab:thr-forced"))
    w("thr_qv.tex", threshold_table(
        qv,
        "$\\opt_{\\QV}$ in threshold form (Certified computation), as in Table~\\ref{tab:thr-full}.",
        "tab:thr-qv"))
    obl, n_obl, n_zero = oblivious_table()
    # Classes of the oblivious grid: eps* = 0 exactly when n < 2s - 1 (Lemma A), eps* = 1/2
    # exactly when n >= 2^s - 1 (Lemma B); the rest are certified to within a small bracket.
    # The certified data must agree with both lemmas.
    obl_zero = obl_half = 0
    nontriv_width = []
    for r in csv.DictReader(open(os.path.join(ROOT, "data", "oblivious_v1.csv"))):
        n, s_ = int(r["n"]), int(r["s"])
        lo, up = Fraction(r["eps_lower"]), Fraction(r["eps_upper"])
        if n < 2 * s_ - 1:
            if lo != 0:
                raise SystemExit(f"oblivious ({n},{s_}): positive lower bound below 2s - 1")
            obl_zero += 1
        elif n >= 2 ** s_ - 1:
            if lo != Fraction(1, 2) or up != Fraction(1, 2):
                raise SystemExit(f"oblivious ({n},{s_}): expected exactly 1/2")
            obl_half += 1
        else:
            if not (0 < lo and up < Fraction(1, 2)):
                raise SystemExit(f"oblivious ({n},{s_}): expected 0 < eps* < 1/2")
            nontriv_width.append(up - lo)
    w("oblivious.tex", obl)
    w("tensor.tex", tensor_table())

    # Counts used in the text.
    common = [k for k in full]
    reuse_adv = sum(1 for k in common if full[k] > qv[k])
    reuse_max = max(full[k] - qv[k] for k in common)
    full_eq_forced = sum(1 for k in common if full[k] == forced[k])
    forced_gt_qv = sum(1 for k in forced if forced[k] > qv[k])

    # E range check: 1 + (s - 1) eps <= ov* <= opt_QV(n, s) and eps <= 1/2, for 2 <= s <= n <= 22.
    best = (0.0, None)
    open_cells = []
    need = []
    excluded_by_lemma_a = 0
    for (n, s), c in qv.items():
        if s < 2 or s > n:
            continue
        eps_cap = min(Fraction(1, 2), Fraction(c - 1, s - 1))
        if n < 2 * s - 1:  # Lemma A: eps* = 0
            if eps_cap > 0 and math.log(eps_cap * n * s) / math.log(n * n * s) > 0.543896:
                excluded_by_lemma_a += 1
            continue
        if eps_cap <= 0:
            continue
        e = math.log(eps_cap * n * s) / math.log(n * n * s)
        if e > best[0]:
            best = (e, (n, s, c))
        if e > 0.543896:
            open_cells.append((n, s))
            need.append((n * n * s) ** 0.543896 / (n * s))
    # Binary-search bound B(n, s) = max(opt(n, s-1), 1 + opt(floor((n-1)/2), s-1)).
    ext = {}
    with open(os.path.join(ROOT, "data", "values_small_long_v1.csv")) as f:
        for r in csv.DictReader(f):
            ext[(r["variant"], int(r["n"]), int(r["s"]))] = int(r["value"])
    tables = {"full": full, "forced": forced, "qv": qv}
    for (var, n, s), v in ext.items():
        if (n, s) in tables[var] and tables[var][(n, s)] != v:
            raise SystemExit(f"values_small_long_v1 disagrees with the {var} table at {(n, s)}")

    def val(var, n, s):
        if n <= 0 or s <= 0:
            return 0
        if var == "qv":
            s = min(s, n)  # at most n rounds can be played without reuse
        t = tables[var]
        return t[(n, s)] if (n, s) in t else ext[(var, n, s)]

    def btight(var, maxn):
        # Counts only step-up cells, opt(n, s) > opt(n, s - 1): elsewhere B >= opt(n, s - 1)
        # = opt(n, s) holds trivially.
        cells, exc = 0, []
        for (n, s), v in sorted(tables[var].items()):
            if n > maxn:
                continue
            b = max(val(var, n, s - 1), 1 + val(var, (n - 1) // 2, s - 1))
            if b > v:
                raise SystemExit(f"binary-search bound violated: {var} {(n, s)}")
            if v == val(var, n, s - 1):
                if b != v:
                    raise SystemExit("plateau cell with B < opt")
                continue
            cells += 1
            if b < v:
                exc.append((n, s))
        return cells, exc

    bt = {key: btight(*key) for key in [("forced", 15), ("forced", 20), ("full", 18), ("qv", 22)]}
    fc, fe = bt[("forced", 15)]
    not_reuse = [c for c in fe if forced[c] == qv[c]]

    # Staircase recurrences f(n, s) = max(R(n, s), f(n, s - 1), 1 + f(n - c, s - t)) with the
    # binary-search value R as base: consistent (never above opt_forced) on every stored cell
    # for c = 3, t = 3, 4, 5, although Corollary cor:staircase rules them out for large n.
    def rbs(n, s):
        return 0 if n <= 0 or s <= 0 else min(s, (n + 1).bit_length() - 1)

    def stair_consistent(c, t):
        f = {}
        for n in range(0, max(k[0] for k in forced) + 1):
            for s in range(0, 2 * max(k[0] for k in forced) + 1):
                v = 0 if n <= 0 or s <= 0 else max(rbs(n, s), f[(n, s - 1)])
                if n >= c and s >= t and n > 0 and s > 0:
                    v = max(v, 1 + f[(n - c, s - t)])
                f[(n, s)] = v
        return sum(1 for k, v in forced.items() if f[k] <= v)

    stairs = [stair_consistent(3, t) for t in (3, 4, 5)]
    if any(x != len(forced) for x in stairs) or stair_consistent(2, 2) != 215:
        raise SystemExit(f"staircase consistency changed: {stairs}")

    # Refills along the principal variations (Observation obs:refills).
    ref, oth = {}, {}
    with open(os.path.join(ROOT, "data", "gap_decomposition_v1.csv")) as f:
        for r in csv.DictReader(f):
            if r["region"] == "stop" or r["j"] == "":
                continue
            left = int(r["r"]) - 1  # rounds remaining after the sign is placed
            d = ref if r["refill"] == "true" else oth
            placed, kept = d.get(left, (0, 0))
            d[left] = (placed + 1, kept + (r["sign_fate"] == "survives"))
    n_ref = sum(v[0] for v in ref.values())
    n_ref_kept = sum(v[1] for v in ref.values())
    n_moves = n_ref + sum(v[0] for v in oth.values())
    reweighted = (sum(ref[k][0] * Fraction(oth[k][1], oth[k][0]) for k in ref)
                  / sum(ref[k][0] for k in ref))
    with open(os.path.join(ROOT, "data", "gap_pv_summary_v1.csv")) as f:
        lines = list(csv.DictReader(f))
    lines_ref = sum(1 for r in lines if int(r["refills"]) > 0)
    reuse_lines = [r for r in lines if forced[(int(r["n"]), int(r["s"]))] > qv[(int(r["n"]), int(r["s"]))]]
    reuse_lines_ref = sum(1 for r in reuse_lines if int(r["refills"]) > 0)

    # First moves (Observation obs:first): every first move is worth at least opt(n, s - 1),
    # and at a jump column the optimal first moves form an interval symmetric about the centre.
    fm_pairs = fm_moves = 0
    with open(os.path.join(ROOT, "data", "first_moves_forced_v1.csv")) as f:
        for r in csv.DictReader(f):
            n, s_, v = int(r["n"]), int(r["s"]), int(r["value"])
            vals = [int(x) for x in r["pointer_values"].split()]
            opt_moves = sorted(int(x) for x in r["optimal_first_moves"].split())
            prev = forced[(n, s_ - 1)] if s_ > 1 else 0
            fm_pairs += 1
            fm_moves += len(vals)
            if len(vals) != n or min(vals) < prev or max(vals) != v:
                raise SystemExit(f"first-move observation fails at {(n, s_)}")
            if opt_moves != [j for j in range(n) if vals[j] == v]:
                raise SystemExit(f"first-moves file inconsistent at {(n, s_)}")
            if v > prev:
                lo, hi = opt_moves[0], opt_moves[-1]
                if opt_moves != list(range(lo, hi + 1)) or lo + hi != n - 1:
                    raise SystemExit(f"jump column {(n, s_)}: optimal first moves not a centred interval")

    tens = max(tensor_rows(), key=lambda r: r[4])
    q20 = qv[(20, 20)]
    beta1 = math.log((q20 + 1) / 2) / math.log(20)

    macros = {
        "NFullCells": len(full),
        "NForcedCells": len(forced),
        "NQVCells": len(qv),
        "ReuseAdv": reuse_adv,
        "ReuseMax": reuse_max,
        "FullEqForced": full_eq_forced,
        "ForcedGtQV": forced_gt_qv,
        "NOblCells": n_obl,
        "NOblZero": n_zero,
        "NOblZeroLemma": obl_zero,
        "NOblHalf": obl_half,
        "NOblNontrivial": len(nontriv_width),
        "OblMaxWidth": f"{float(max(nontriv_width)):.1e}".replace("e-0", "\\cdot 10^{-") + "}",
        "RangeE": f"{best[0]:.4f}",
        "RangeEn": best[1][0],
        "RangeEs": best[1][1],
        "RangeOpenCount": len(open_cells),
        "RangeExcludedA": excluded_by_lemma_a,
        "RangeNeedMin": f"{math.floor(min(need) * 1000) / 1000:.3f}",
        "RangeNeedMax": f"{math.ceil(max(need) * 1000) / 1000:.3f}",
        "RangeOpenList": ", ".join(f"$({a},{b})$" for a, b in sorted(open_cells)),
        "RangeOpenFirstN": min(open_cells)[0],
        "RangeOpenFirstS": min(open_cells)[1],
        "TensorBestT": tens[0],
        "TensorBestPhi": f"{tens[4]:.4f}",
        "TensorBestQV": f"{max(r[5] for r in tensor_rows()):.4f}",
        "QVTwentyTwenty": q20,
        "FMPairs": fm_pairs,
        "FMMoves": fm_moves,
        "PVLines": len(lines),
        "PVMoves": n_moves,
        "RefillCount": n_ref,
        "RefillLines": lines_ref,
        "RefillKept": n_ref_kept,
        "RefillRate": f"{n_ref_kept / n_ref:.3f}",
        "OtherRateReweighted": f"{float(reweighted):.3f}",
        "ReuseLines": len(reuse_lines),
        "ReuseLinesRefill": reuse_lines_ref,
        "StairConsistent": stairs[0],
        "BCellsFifteen": fc,
        "BTightFifteen": fc - len(fe),
        "BExcFifteen": len(fe),
        "BExcFifteenList": ", ".join(f"$({n},{s})$" for n, s in fe if (n, s) not in not_reuse),
        "BExcFifteenNoReuse": ", ".join(f"$({n},{s})$" for n, s in not_reuse),
        "BCellsForced": bt[("forced", 20)][0],
        "BTightForced": bt[("forced", 20)][0] - len(bt[("forced", 20)][1]),
        "BCellsFull": bt[("full", 18)][0],
        "BTightFull": bt[("full", 18)][0] - len(bt[("full", 18)][1]),
        "BCellsQV": bt[("qv", 22)][0],
        "BTightQV": bt[("qv", 22)][0] - len(bt[("qv", 22)][1]),
        "BExcQVList": ", ".join(f"$({n},{s})$" for n, s in bt[("qv", 22)][1]),
        "BetaAtOne": f"{beta1:.4f}",
        "CalAtOne": f"{(beta1 + 1) / 3:.4f}",
    }
    with open(os.path.join(outdir, "numbers.tex"), "w") as f:
        f.write("% Generated by scripts/make_tables.py from data/. Do not edit.\n")
        for k, v in macros.items():
            f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
    for k, v in macros.items():
        print(f"{k} = {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
