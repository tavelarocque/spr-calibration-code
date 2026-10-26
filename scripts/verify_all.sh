#!/usr/bin/env bash
# Default verification target (see REPRODUCE.md).
#
# Usage: scripts/verify_all.sh            (from the repository root or anywhere)
#
# Builds the workspace, runs the tests, re-verifies every stored certificate, and
# regenerates the small and medium tables, comparing each with the committed data. About
# six minutes and under 0.5 GB on an 8 GB laptop. Writes only to a temporary directory.
# Exits non-zero if any check fails. The long runs are listed in REPRODUCE.md.
set -u
cd "$(dirname "$0")/.."
ROOT=$(pwd)
OUT=$(mktemp -d "${TMPDIR:-/tmp}/spr-verify.XXXXXX")
SPR="$ROOT/target/release/spr"
MAXST=40000000
FAILED=()
START=$(date +%s)

step() { echo; echo "=== $1"; }
check() {  # check NAME CMD...: run CMD, record a failure under NAME
  local name=$1; shift
  if "$@"; then echo "--- ok: $name"; else echo "--- FAILED: $name"; FAILED+=("$name"); fi
}
expect() {  # expect NAME PATTERN CMD...: CMD must succeed and print a line matching PATTERN
  local name=$1 pat=$2; shift 2
  local log="$OUT/$(echo "$name" | tr ' /' '__').log"
  if "$@" > "$log" 2>&1 && grep -q -- "$pat" "$log"; then
    echo "--- ok: $name"
  else
    echo "--- FAILED: $name (log: $log)"; tail -5 "$log"; FAILED+=("$name")
  fi
}

step "toolchain"
echo "temporary output: $OUT"
check "rustc and cargo present" bash -c 'rustc --version && cargo --version'
check "python3 with sympy" python3 -c "import sys, sympy; print(sys.version.split()[0], 'sympy', sympy.__version__)"

step "build and test (cargo test --release)"
check "cargo build" cargo build --release -q
check "cargo test" cargo test --release -q

step "solver against the Python reference and the v0 tables"
expect "check-oracle" "check-oracle: PASS" "$SPR" check-oracle
expect "known values opt_QV(2^t - 1, t) = t" "t=5: opt_qv(31, 5) = 5 \[PASS\]" "$SPR" known-values --max-t 5

step "value tables: regenerate and compare with data/"
for spec in "forced 18" "full 16" "qv 19"; do
  set -- $spec
  check "tabulate $1 n <= $2" "$SPR" tabulate --variant "$1" --max-n "$2" --max-states $MAXST --out "$OUT/values_$1.csv"
  check "compare $1" python3 scripts/compare_values.py "$OUT/values_$1.csv" "data/values_$1_v1.csv"
done
for v in full forced; do
  check "tabulate $v n <= 9, s <= 40" "$SPR" tabulate --variant "$v" --max-n 9 --s-hi 40 --max-states $MAXST --out "$OUT/long_$v.csv"
  check "compare $v long" python3 scripts/compare_values.py "$OUT/long_$v.csv" data/values_small_long_v1.csv --where variant=$v
done

step "maximal erasure (Certified computation 4.6) on n <= 12; the full n <= 14 run is a long target"
expect "erasure check n <= 12" "no counterexample: 12223144 reachable" "$SPR" find-erasure-counterexample --max-n 12 --max-states $MAXST

step "gap lemma, principal variations, first moves (byte-identical outputs)"
expect "gap-check n <= 14" "no violation: 245759 boards, 4734974 checks" "$SPR" gap-check --max-n 14 --max-states $MAXST
check "gap-check tight total and slack pattern" python3 -c '
import re, sys
t, ok = 0, True
for l in open(sys.argv[1]):
    m = re.search(r"n=\s*(\d+):.*tight \(slack 0\) (\d+), max slack (\d+)", l)
    if m:
        t += int(m.group(2))
        ok &= int(m.group(3)) == (int(m.group(1)) - 1) // 2
print("tight", t, "slack pattern", ok)
sys.exit(0 if t == 2572341 and ok else 1)' "$OUT/gap-check_n_<=_14.log"
check "gap-pv n <= 16" "$SPR" gap-pv --max-n 16 --out "$OUT/gap_pv.csv" --summary "$OUT/gap_pv_sum.csv" --max-states $MAXST
check "cmp gap_decomposition" cmp "$OUT/gap_pv.csv" data/gap_decomposition_v1.csv
check "cmp gap_pv_summary" cmp "$OUT/gap_pv_sum.csv" data/gap_pv_summary_v1.csv
check "gap-pv diagonal n <= 18" "$SPR" gap-pv --diag-only --max-n 18 --out "$OUT/lg.csv" --summary "$OUT/lg_sum.csv" --max-states $MAXST
check "cmp long_game_pv" cmp "$OUT/lg.csv" data/long_game_pv_v1.csv
check "cmp long_game_pv_summary" cmp "$OUT/lg_sum.csv" data/long_game_pv_summary_v1.csv
check "gap-compare" "$SPR" gap-compare --out "$OUT/avo.csv" --max-states $MAXST
check "cmp adaptive_vs_oblivious" cmp "$OUT/avo.csv" data/adaptive_vs_oblivious_v1.csv
check "first-moves n <= 15" "$SPR" first-moves --variant forced --max-n 15 --out "$OUT/fm.csv"
check "cmp first_moves" cmp "$OUT/fm.csv" data/first_moves_forced_v1.csv
expect "QV (3,3): a corner first move is worth 1" "Labeler value of each pointer: 0:1 1:2 2:1" "$SPR" principal-variation 3 3 qv

step "fixed-policy values"
for spec in "binary-search 32 both" "binary-search 12 all" "walker-fresh 32 both" "walker-both 16 both" \
            "walker-centre 12 all" "walker-left 12 all" "walker-right 12 all"; do
  set -- $spec
  check "policy $1 forced n <= $2 ($3)" "$SPR" policy-eval --policy "$1" --max-n "$2" --s-mode "$3" --out "$OUT/pol_$1_forced_$3.csv" --max-states $MAXST
  check "compare policy $1 forced ($3)" python3 scripts/compare_values.py "$OUT/pol_$1_forced_$3.csv" data/policy_values_v1.csv --where policy=$1 --where labeler=forced
done
for spec in "binary-search 32 both" "binary-search 10 all" "walker-fresh 20 both" "walker-both 16 both" "walker-right 10 all"; do
  set -- $spec
  check "policy $1 full n <= $2 ($3)" "$SPR" policy-eval --policy "$1" --full --max-n "$2" --s-mode "$3" --out "$OUT/pol_$1_full_$3.csv" --max-states $MAXST
  check "compare policy $1 full ($3)" python3 scripts/compare_values.py "$OUT/pol_$1_full_$3.csv" data/policy_values_v1.csv --where policy=$1 --where labeler=full
done

step "oblivious strategies: exact re-verification of every stored strategy"
for f in data/strategies/eps_*.json; do
  claim=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['info']['eps_certified'])" "$f")
  expect "verify $(basename "$f")" "same eps" "$SPR" obl verify "$f" --claim "$claim" --naive
done
for dk in "3 1" "4 1" "4 2" "5 2"; do
  set -- $dk
  expect "binomial cube d=$1 k=$2" "exact eps = 1/$((2 ** $2)) " "$SPR" obl cube "$1" "$2"
done
expect "oblivious binary search s=4" "eps = 1/2" "$SPR" obl binary-search 4
expect "product (3,2) x (3,2)" "1/4" "$SPR" obl product data/strategies/eps_n03_s2.json data/strategies/eps_n03_s2.json
expect "product (3,2) x (7,3)" "1/4" "$SPR" obl product data/strategies/eps_n03_s2.json data/strategies/eps_n07_s3.json
expect "golden family and golden cube" "verify_golden: PASS" python3 scripts/verify_golden.py

step "oblivious certification: re-run on the fast instances (recorded time <= 10 s)"
expect "recertify eps* and ov*" "recertify_obl: .* PASS" python3 scripts/recertify_obl.py --max-seconds 10 --spr "$SPR"

step "closed forms (exact lower bounds inside the certified bracket)"
for ns in "5 3" "6 3" "7 4" "8 4" "9 4" "10 4" "11 4" "12 4" "13 4" "9 5"; do
  set -- $ns
  check "active-set $1 $2" "$SPR" obl active-set "$1" "$2" --out "$OUT/basis_$1_$2.json"
  expect "closed form ($1,$2)" "lower <= closed <= upper: True" python3 scripts/closed_form.py "$OUT/basis_$1_$2.json"
  if ! grep -q "all constraints >= 0 at eps\*: True" "$OUT/closed_form_($1,$2).log" 2>/dev/null; then
    echo "--- FAILED: closed form ($1,$2) constraint check"; FAILED+=("closed form ($1,$2) constraints")
  fi
done

step "statements checked against the data (counts, thresholds, first moves, staircase)"
check "make_tables" python3 scripts/make_tables.py "$OUT/tables"

echo
echo "=== verify_all finished in $(( $(date +%s) - START )) s; output in $OUT"
if [ ${#FAILED[@]} -eq 0 ]; then
  echo "verify_all: PASS"
  exit 0
else
  echo "verify_all: FAILED (${#FAILED[@]}):"
  printf '  %s\n' "${FAILED[@]}"
  exit 1
fi
