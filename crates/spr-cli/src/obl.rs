//! `spr obl ...`: the oblivious no-reuse Pointer.

use clap::Subcommand;
use num_rational::BigRational;
use num_traits::ToPrimitive;
use serde_json::json;
use spr_oblivious::constructions::{binary_search, binomial_cube, binomial_cube_coins, product};
use spr_oblivious::exact::{verify, verify_naive};
use spr_oblivious::rational::{fmt, parse, to_f64};
use spr_oblivious::search::{EpsResult, ov_certified, solve};
use spr_oblivious::strategy::Strategy;
use spr_oblivious::tuples::Space;
use std::path::{Path, PathBuf};

#[derive(Subcommand)]
pub enum Obl {
    /// Find eps*(n, s) by float LP bisection and certify exact lower and upper bounds.
    Solve {
        n: usize,
        s: usize,
        #[arg(long, default_value_t = 1e-9)]
        tol: f64,
        /// Save the certified lower-bound strategy as JSON.
        #[arg(long)]
        save: Option<PathBuf>,
        /// Also print the support of the strategy.
        #[arg(long)]
        show: bool,
    },
    /// ov*(n, s) with exact bounds.
    Ov { n: usize, s: usize },
    /// Exactly verify a strategy file; optionally check a claimed eps.
    Verify {
        file: PathBuf,
        #[arg(long)]
        claim: Option<String>,
        /// Also run the independent naive verifier and require the same exact eps.
        #[arg(long)]
        naive: bool,
    },
    /// Build the DDFGKO binomial cube for (d, k) and verify it exactly.
    Cube {
        d: usize,
        k: usize,
        #[arg(long)]
        save: Option<PathBuf>,
    },
    /// Oblivious binary search on 2^s - 1 cells, verified exactly.
    BinarySearch {
        s: usize,
        #[arg(long)]
        save: Option<PathBuf>,
    },
    /// Product of two strategy files, verified exactly.
    Product {
        a: PathBuf,
        b: PathBuf,
        #[arg(long)]
        save: Option<PathBuf>,
    },
    /// Solve and certify eps* and ov* over a range, writing a CSV (never overwritten).
    Grid {
        /// Comma-separated `s:n_max` pairs, e.g. `2:13,3:13,4:13,5:10`; n runs from s.
        #[arg(long)]
        ranges: String,
        #[arg(long)]
        out: PathBuf,
        /// Directory for the certified strategies as JSON (optional).
        #[arg(long)]
        strategies: Option<PathBuf>,
        /// Skip ov*.
        #[arg(long)]
        no_ov: bool,
        /// Directory for basis exports (input to scripts/closed_form.py).
        #[arg(long)]
        bases: Option<PathBuf>,
    },
    /// Solve, then export the final basis system for scripts/closed_form.py.
    ActiveSet {
        n: usize,
        s: usize,
        #[arg(long)]
        out: PathBuf,
    },
}

pub fn exponent(n: usize, s: usize, eps: f64) -> f64 {
    let (n, s) = (n as f64, s as f64);
    (eps * n * s).ln() / (n * n * s).ln()
}

fn save_strategy(st: &Strategy, path: &Path, info: serde_json::Value) {
    if path.exists() {
        eprintln!("refusing to overwrite {}", path.display());
        std::process::exit(1);
    }
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).ok();
    }
    std::fs::write(path, serde_json::to_string_pretty(&st.to_json(info)).unwrap()).unwrap();
    println!("saved {}", path.display());
}

fn load(path: &Path) -> Strategy {
    let text = std::fs::read_to_string(path).expect("read strategy file");
    let v: serde_json::Value = serde_json::from_str(&text).expect("parse JSON");
    Strategy::from_json(&v).expect("strategy JSON")
}

pub fn report_eps(r: &EpsResult) {
    let sp = r.space;
    println!(
        "eps*({}, {}): LP bracket [{:.12}, {:.12}] after {} LP solves ({} cols, {} rows, {} nnz), {:.2}s",
        sp.n, sp.s, r.a, r.b, r.lp_solves, r.cols, r.rows, r.nnz, r.seconds
    );
    println!("  certified lower bound: {}  (= {:.12})", fmt(&r.lower), to_f64(&r.lower));
    println!(
        "  upper bound: {} (= {:.12}) [{}]",
        fmt(&r.upper.value),
        to_f64(&r.upper.value),
        if r.upper.certified { &r.upper.how } else { "float, uncertified" }
    );
    println!(
        "  support: {} tuples; E at lower bound = {:.6}",
        r.lower_strategy.entries.len(),
        exponent(sp.n, sp.s, to_f64(&r.lower))
    );
}

fn print_support(st: &Strategy) {
    let mut e = st.entries.clone();
    e.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    for (t, p) in e {
        println!("    {:?}  {:.10}", t, p.to_f64().unwrap_or(f64::NAN));
    }
}

pub fn run(cmd: Obl) -> bool {
    match cmd {
        Obl::Solve { n, s, tol, save, show } => {
            let r = match solve(Space::new(n, s), tol) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("error: {e}");
                    return false;
                }
            };
            report_eps(&r);
            if show {
                print_support(&r.lower_strategy);
            }
            if let Some(p) = save {
                save_strategy(
                    &r.lower_strategy,
                    &p,
                    json!({"eps_certified": fmt(&r.lower), "eps_upper": fmt(&r.upper.value),
                           "upper_certified": r.upper.certified, "source": "spr obl solve"}),
                );
            }
            true
        }
        Obl::Ov { n, s } => match ov_certified(Space::new(n, s)) {
            Ok(r) => {
                let exact = r.lower == r.upper;
                println!(
                    "ov*({n}, {s}): float {:.12}; certified [{}, {}] {} ({:.2}s)",
                    r.float,
                    fmt(&r.lower),
                    fmt(&r.upper),
                    if exact { "EXACT" } else { "interval" },
                    r.seconds
                );
                true
            }
            Err(e) => {
                eprintln!("error: {e}");
                false
            }
        },
        Obl::Verify { file, claim, naive } => {
            let st = load(&file);
            match verify(&st) {
                Ok(rep) => {
                    if naive {
                        match verify_naive(&st) {
                            Ok(e) if e == rep.eps => println!("independent naive verifier: same eps"),
                            Ok(e) => {
                                println!("REJECTED: naive verifier gives {}, main verifier {}", fmt(&e), fmt(&rep.eps));
                                return false;
                            }
                            Err(e) => {
                                println!("REJECTED by the naive verifier: {e}");
                                return false;
                            }
                        }
                    }
                    println!(
                        "{}: n={} s={} support={} eps = {} (= {:.12}), worst prefix {:?}",
                        file.display(),
                        st.n,
                        st.s,
                        st.entries.len(),
                        fmt(&rep.eps),
                        to_f64(&rep.eps),
                        rep.worst
                    );
                    if let Some(c) = claim {
                        let c: BigRational = parse(&c).expect("claim must be p/q");
                        let ok = rep.eps >= c;
                        println!("claim eps >= {}: {}", fmt(&c), if ok { "VERIFIED" } else { "REJECTED" });
                        return ok;
                    }
                    true
                }
                Err(e) => {
                    println!("REJECTED: {e}");
                    false
                }
            }
        }
        Obl::Cube { d, k, save } => {
            println!("binomial cube d={d} k={k}: {} coins", binomial_cube_coins(d, k));
            let st = match binomial_cube(d, k) {
                Ok(st) => st,
                Err(e) => {
                    eprintln!("error: {e}");
                    return false;
                }
            };
            let rep = verify(&st).expect("cube is a valid strategy");
            println!(
                "  n = {}, s = {}, support = {} tuples, exact eps = {} (expected 2^-{k}), E = {:.6}",
                st.n,
                st.s,
                st.entries.len(),
                fmt(&rep.eps),
                exponent(st.n, st.s, to_f64(&rep.eps))
            );
            if st.entries.len() <= 5000 {
                println!("  independent naive verifier: {}", fmt(&verify_naive(&st).unwrap()));
            }
            if let Some(p) = save {
                save_strategy(&st, &p, json!({"construction": format!("DDFGKO binomial cube d={d} k={k}"), "eps": fmt(&rep.eps)}));
            }
            true
        }
        Obl::BinarySearch { s, save } => {
            let st = binary_search((1 << s) - 1, s);
            let rep = verify(&st).unwrap();
            println!("oblivious binary search n={} s={s}: eps = {}", st.n, fmt(&rep.eps));
            if let Some(p) = save {
                save_strategy(&st, &p, json!({"construction": "oblivious binary search", "eps": fmt(&rep.eps)}));
            }
            true
        }
        Obl::Product { a, b, save } => {
            let (pa, pb) = (load(&a), load(&b));
            let (ea, eb) = (verify(&pa).expect("first strategy").eps, verify(&pb).expect("second strategy").eps);
            let st = product(&pa, &pb).expect("product");
            let rep = verify(&st).expect("product is a valid strategy");
            let bound = &ea * &eb;
            println!(
                "product: n={} s={} support={}; eps1 = {}, eps2 = {}, eps1*eps2 = {}; exact eps(product) = {} -> {}",
                st.n,
                st.s,
                st.entries.len(),
                fmt(&ea),
                fmt(&eb),
                fmt(&bound),
                fmt(&rep.eps),
                if rep.eps >= bound { "lemma holds" } else { "LEMMA FAILS" }
            );
            if let Some(p) = save {
                save_strategy(&st, &p, json!({"construction": "product", "eps": fmt(&rep.eps)}));
            }
            rep.eps >= bound
        }
        Obl::Grid { ranges, out, strategies, no_ov, bases } => grid(&ranges, &out, strategies.as_deref(), no_ov, bases.as_deref()),
        Obl::ActiveSet { n, s, out } => {
            let r = solve(Space::new(n, s), 1e-11).expect("solve");
            report_eps(&r);
            crate::obl_export::export(&r, &out);
            true
        }
    }
}

fn grid(ranges: &str, out: &Path, strategies: Option<&Path>, no_ov: bool, bases: Option<&Path>) -> bool {
    use std::io::Write;
    if out.exists() {
        eprintln!("refusing to overwrite {}", out.display());
        return false;
    }
    let mut f = std::fs::File::create(out).unwrap();
    writeln!(f, "n,s,lp_a,lp_b,eps_lower,eps_upper,upper_certified,E_lower,ov_lower,ov_upper,ov_exact,eps_seconds,ov_seconds,cols,nnz").unwrap();
    for part in ranges.split(',') {
        let (s, nmax) = part.split_once(':').expect("ranges are s:n_max");
        let (s, nmax): (usize, usize) = (s.trim().parse().unwrap(), nmax.trim().parse().unwrap());
        let mut prev: Option<Strategy> = None;
        for n in s.max(2)..=nmax {
            let sp = Space::new(n, s);
            let mut r = match solve(sp, 1e-10) {
                Ok(r) => r,
                Err(e) => {
                    println!("n={n} s={s}: stopped ({e})");
                    break;
                }
            };
            // Monotonicity in n: the certified (n - 1, s) strategy is a strategy on n
            // cells that never uses the last one. Re-verified exactly at n.
            if let Some(p) = prev.take() {
                let embedded = Strategy { n, s, entries: p.entries };
                if let Ok(rep) = verify(&embedded)
                    && rep.eps > r.lower
                {
                    r.lower = rep.eps;
                    r.lower_strategy = embedded;
                }
            }
            prev = Some(r.lower_strategy.clone());
            if let Some(dir) = bases {
                std::fs::create_dir_all(dir).ok();
                crate::obl_export::export(&r, &dir.join(format!("basis_n{n:02}_s{s}.json")));
            }
            let (ovl, ovu, ovs) = if no_ov {
                (String::new(), String::new(), 0.0)
            } else {
                match ov_certified(sp) {
                    Ok(o) => (fmt(&o.lower), fmt(&o.upper), o.seconds),
                    Err(e) => (format!("error: {e}"), String::new(), 0.0),
                }
            };
            let lower_f = to_f64(&r.lower);
            writeln!(
                f,
                "{n},{s},{:.15},{:.15},{},{},{},{:.6},{},{},{},{:.3},{:.3},{},{}",
                r.a,
                r.b,
                fmt(&r.lower),
                fmt(&r.upper.value),
                r.upper.certified,
                exponent(n, s, lower_f),
                ovl,
                ovu,
                !no_ov && ovl == ovu,
                r.seconds,
                ovs,
                r.cols,
                r.nnz
            )
            .unwrap();
            f.flush().unwrap();
            println!(
                "n={n:2} s={s}: eps in [{:.10}, {:.10}] ({}), E={:.4}, ov=[{ovl}, {ovu}], {:.2}s",
                lower_f,
                to_f64(&r.upper.value),
                if r.upper.certified { "certified" } else { "upper uncertified" },
                exponent(n, s, lower_f),
                r.seconds + ovs
            );
            if let Some(dir) = strategies {
                let p = dir.join(format!("eps_n{n:02}_s{s}.json"));
                if !p.exists() {
                    save_strategy(
                        &r.lower_strategy,
                        &p,
                        json!({"eps_certified": fmt(&r.lower), "eps_upper": fmt(&r.upper.value),
                               "upper_certified": r.upper.certified, "source": "spr obl grid"}),
                    );
                }
            }
        }
    }
    true
}
