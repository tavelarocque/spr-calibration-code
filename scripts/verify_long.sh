#!/usr/bin/env bash
# Opt-in long runs (see REPRODUCE.md). One target per invocation:
#
#   scripts/verify_long.sh TARGET
#
# TARGET is one of: forced20 full18 qv22 erasure14 obl-grid closed-forms cube-compress
#
# Run one target at a time on an otherwise idle machine; the memory figures in
# REPRODUCE.md are for the solver alone. Every solver run is capped by --max-states.
set -eu
cd "$(dirname "$0")/.."
SPR="$(pwd)/target/release/spr"
OUT=$(mktemp -d "${TMPDIR:-/tmp}/spr-long.XXXXXX")
echo "output in $OUT"
cargo build --release -q

case "${1:-}" in
  forced20)
    "$SPR" tabulate --variant forced --max-n 20 --max-states 58000000 --out "$OUT/values_forced.csv"
    python3 scripts/compare_values.py "$OUT/values_forced.csv" data/values_forced_v1.csv ;;
  full18)
    "$SPR" tabulate --variant full --max-n 18 --max-states 50000000 --out "$OUT/values_full.csv"
    python3 scripts/compare_values.py "$OUT/values_full.csv" data/values_full_v1.csv ;;
  qv22)
    "$SPR" tabulate --variant qv --max-n 22 --max-states 75000000 --out "$OUT/values_qv.csv"
    python3 scripts/compare_values.py "$OUT/values_qv.csv" data/values_qv_v1.csv ;;
  erasure14)
    "$SPR" find-erasure-counterexample --max-n 14 --max-states 120000000 | tee "$OUT/erasure.log"
    grep -q "no counterexample" "$OUT/erasure.log" ;;
  obl-grid)
    python3 scripts/recertify_obl.py --max-seconds 100000 --spr "$SPR" ;;
  closed-forms)
    for ns in "10 5" "12 5" "13 5"; do
      set -- $ns
      "$SPR" obl active-set "$1" "$2" --out "$OUT/basis_$1_$2.json"
      python3 scripts/closed_form.py "$OUT/basis_$1_$2.json" | tee "$OUT/cf_$1_$2.log"
      grep -q "lower <= closed <= upper: True" "$OUT/cf_$1_$2.log"
    done ;;
  cube-compress)
    for dk in "3 2" "4 2" "4 3" "5 3"; do
      set -- $dk
      python3 scripts/cube_compress.py "$1" "$2" fair
      python3 scripts/cube_compress.py "$1" "$2" golden
    done ;;
  *)
    echo "usage: $0 {forced20|full18|qv22|erasure14|obl-grid|closed-forms|cube-compress}"; exit 2 ;;
esac
echo "verify_long $1: PASS"
