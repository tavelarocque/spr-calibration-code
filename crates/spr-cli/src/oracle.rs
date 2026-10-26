//! Validation (a): exact equality with reference/spr_ref.py on every cell.
//!
//! The oracle is imported, not run as a script: its `__main__` computes every
//! variant over the same range, and full SPR would stall at the forced ranges.
//! The oracle file itself is never modified.

use spr_core::{Forced, Full, Qv, Variant, VariantKind};
use spr_solve::{Config, Solver};
use std::path::Path;
use std::process::Command;
use std::time::Instant;

/// `(variant, n_max, s_max)` ranges required by the task for the live oracle.
const LIVE: [(VariantKind, u32, u32); 3] =
    [(VariantKind::Full, 11, 8), (VariantKind::Forced, 14, 12), (VariantKind::Qv, 11, 8)];

/// Rows `n = 1..=n_max`, each with `s = 1..=s_max`.
type Table = Vec<Vec<u32>>;

fn run_python(kind: VariantKind, n_max: u32, s_max: u32) -> Result<Table, String> {
    let (game, mode) = kind.oracle_args();
    let script = format!(
        "import sys\nsys.path.insert(0, 'reference')\nimport spr_ref\n\
         f = spr_ref.make_solver('{game}', '{mode}')\n\
         for n in range(1, {n_max} + 1):\n    print(' '.join(str(f(n, s)) for s in range(1, {s_max} + 1)), flush=True)\n"
    );
    let out = Command::new("python3").arg("-c").arg(&script).output().map_err(|e| e.to_string())?;
    if !out.status.success() {
        return Err(String::from_utf8_lossy(&out.stderr).into_owned());
    }
    let text = String::from_utf8_lossy(&out.stdout);
    let table: Table =
        text.lines().map(|l| l.split_whitespace().map(|v| v.parse().unwrap()).collect()).collect();
    if table.len() != n_max as usize || table.iter().any(|r| r.len() != s_max as usize) {
        return Err(format!("unexpected oracle output shape:\n{text}"));
    }
    Ok(table)
}

/// Parse the three tables of data/spr_values_v0.md.
fn read_v0(path: &Path) -> Result<Vec<(VariantKind, Table)>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut out: Vec<(VariantKind, Table)> = Vec::new();
    let mut current: Option<VariantKind> = None;
    for line in text.lines() {
        if let Some(h) = line.strip_prefix("## ") {
            current = if h.starts_with("opt(n, s), full SPR") {
                Some(VariantKind::Full)
            } else if h.starts_with("SPR with forced") {
                Some(VariantKind::Forced)
            } else if h.starts_with("opt_QV") {
                Some(VariantKind::Qv)
            } else {
                None
            };
            if let Some(k) = current {
                out.push((k, Vec::new()));
            }
        } else if let (Some(_), Some(rest)) = (current, line.strip_prefix("n=")) {
            let mut it = rest.split_whitespace();
            let n: usize = it.next().unwrap().parse().map_err(|_| format!("bad row {line:?}"))?;
            let row: Vec<u32> = it.map(|v| v.parse().unwrap()).collect();
            let table = &mut out.last_mut().unwrap().1;
            if table.len() + 1 != n {
                return Err(format!("rows out of order at {line:?}"));
            }
            table.push(row);
        }
    }
    if out.len() != 3 {
        return Err(format!("expected 3 tables in {}, found {}", path.display(), out.len()));
    }
    Ok(out)
}

fn compare<V: Variant>(expected: &Table, source: &str) -> bool {
    let mut solver = Solver::<u64, V>::new(Config::default());
    let mut cells = 0;
    let mut bad = Vec::new();
    for (i, row) in expected.iter().enumerate() {
        let n = i as u32 + 1;
        for (k, &e) in row.iter().enumerate() {
            let s = k as u32 + 1;
            let v = solver.opt(n, s).unwrap();
            cells += 1;
            if v != e {
                bad.push((n, s, v, e));
            }
        }
    }
    let s_max = expected.iter().map(|r| r.len()).max().unwrap_or(0);
    if bad.is_empty() {
        println!("  {:6} vs {source}: n <= {}, s <= {s_max}: all {cells} cells equal", V::KIND.name(), expected.len());
        true
    } else {
        println!("  {:6} vs {source}: {} MISMATCHES of {cells} cells", V::KIND.name(), bad.len());
        bad.sort();
        for (n, s, v, e) in bad.iter().take(10) {
            println!("    n={n} s={s}: rust {v}, expected {e}");
        }
        false
    }
}

fn compare_kind(kind: VariantKind, t: &Table, source: &str) -> bool {
    match kind {
        VariantKind::Full => compare::<Full>(t, source),
        VariantKind::Forced => compare::<Forced>(t, source),
        VariantKind::Qv => compare::<Qv>(t, source),
    }
}

pub fn check(no_python: bool) -> bool {
    let mut ok = true;
    println!("data/spr_values_v0.md:");
    match read_v0(Path::new("data/spr_values_v0.md")) {
        Ok(tables) => {
            for (kind, t) in &tables {
                ok &= compare_kind(*kind, t, "spr_values_v0.md");
            }
        }
        Err(e) => {
            println!("  could not read v0 tables: {e}");
            ok = false;
        }
    }
    if !no_python {
        println!("reference/spr_ref.py (live):");
        for (kind, n_max, s_max) in LIVE {
            let t0 = Instant::now();
            match run_python(kind, n_max, s_max) {
                Ok(t) => {
                    print!("  [python {:.1}s]", t0.elapsed().as_secs_f64());
                    ok &= compare_kind(kind, &t, "spr_ref.py");
                }
                Err(e) => {
                    println!("  {kind}: oracle failed: {e}");
                    ok = false;
                }
            }
        }
    }
    println!("{}", if ok { "check-oracle: PASS" } else { "check-oracle: FAIL" });
    ok
}
