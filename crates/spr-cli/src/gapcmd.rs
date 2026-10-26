//! Gap-lemma commands: `gap-check`, `gap-pv`, `gap-compare`.

use num_rational::BigRational;
use num_traits::ToPrimitive;
use spr_core::{Forced, Qv, cell_char};
use spr_solve::gap::{check_n, pv_rows};
use spr_solve::{Config, Solver};
use std::collections::HashMap;
use std::fs::File;
use std::io::Write as _;
use std::path::Path;
use std::time::Instant;

fn create(path: &Path) -> Option<File> {
    if path.exists() {
        eprintln!("refusing to overwrite {}; write a new version instead", path.display());
        return None;
    }
    Some(File::create(path).expect("create output"))
}

/// Item 1: exhaustive check of the gap lemma, n = 1 ..= max_n.
pub fn gap_check(max_n: u32, cfg: Config) -> bool {
    let mut solver = Solver::<u64, Forced>::new(cfg);
    let t0 = Instant::now();
    let (mut boards, mut checks) = (0usize, 0u64);
    for n in 1..=max_n {
        let t = Instant::now();
        match check_n(&mut solver, n) {
            None => {
                println!("# STOPPED at n={n}: not reached (memory), table limit {} states", solver.states());
                return false;
            }
            Some(Err(v)) => {
                println!("\n*** GAP LEMMA VIOLATION ***");
                println!(
                    "n = {}, r = {}, board {} : V = {}, |x| = {}, g = {}, opt_forced(g, r) = {}, so V < |x| + opt(g, r)",
                    v.n,
                    v.r,
                    v.board,
                    v.v,
                    v.board.count(),
                    v.g,
                    v.opt_g
                );
                println!("reached from 0^{} by:", v.n);
                for (b, j, sign) in &v.path {
                    println!("  {b}  point {j}  write {}", cell_char(*sign));
                }
                return false;
            }
            Some(Ok(row)) => {
                boards += row.boards;
                checks += row.checks;
                println!(
                    "n={:2}: {} reachable boards, {} (x, r) checks with |x| + r <= 2n, tight (slack 0) {}, max slack {}, {:.2}s, {} states",
                    row.n,
                    row.boards,
                    row.checks,
                    row.tight,
                    row.max_slack,
                    t.elapsed().as_secs_f64(),
                    solver.states()
                );
            }
        }
    }
    println!(
        "# no violation: {boards} boards, {checks} checks, {:.2}s, {} states",
        t0.elapsed().as_secs_f64(),
        solver.states()
    );
    true
}

/// Items 2 and 3: principal variations with the gap decomposition and refill events.
pub fn gap_pv(max_n: u32, out: &Path, summary: &Path, cfg: Config, diag_only: bool) -> bool {
    let (Some(mut f), Some(mut fs)) = (create(out), create(summary)) else { return false };
    writeln!(
        f,
        "n,s,step,r,board,V,count,g,opt_g_r,slack,region,j,refill,erased,sign,sign_fate"
    )
    .unwrap();
    writeln!(
        fs,
        "n,s,opt,moves,positive_slack_positions,moves_in_gap,moves_left,moves_right,refills,refills_surviving,refills_erased,final_count"
    )
    .unwrap();
    let mut solver = Solver::<u64, Forced>::new(cfg);
    for n in 1..=max_n {
        let t = Instant::now();
        let s_list: Vec<u32> = if diag_only { vec![n, 2 * n] } else { (1..=2 * n).collect() };
        for s in s_list {
            let Some(rows) = pv_rows(&mut solver, n, s) else {
                println!("# STOPPED at n={n} s={s}: not reached (memory)");
                return false;
            };
            let (mut pos, mut gap_m, mut left, mut right, mut refills, mut ref_surv, mut ref_er) = (0, 0, 0, 0, 0, 0, 0);
            let moves = rows.iter().filter(|r| r.j.is_some()).count();
            for r in &rows {
                if r.slack > 0 {
                    pos += 1;
                }
                match r.region {
                    "gap" => gap_m += 1,
                    "left" => left += 1,
                    "right" => right += 1,
                    _ => {}
                }
                if r.refill {
                    refills += 1;
                    if r.erased_at.is_none() { ref_surv += 1 } else { ref_er += 1 }
                }
                let fate = match (r.j, r.erased_at) {
                    (None, _) => String::new(),
                    (Some(_), None) => "survives".into(),
                    (Some(_), Some(k)) => format!("erased at step {k}"),
                };
                let erased: Vec<String> = r.erased.iter().map(|c| c.to_string()).collect();
                writeln!(
                    f,
                    "{n},{s},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                    r.step,
                    r.r,
                    r.board,
                    r.v,
                    r.count,
                    r.g,
                    r.opt_g,
                    r.slack,
                    r.region,
                    r.j.map_or(String::new(), |j| j.to_string()),
                    r.refill,
                    erased.join(" "),
                    r.sign.map_or(String::new(), |c| cell_char(c).to_string()),
                    fate
                )
                .unwrap();
            }
            let last = rows.last().unwrap();
            writeln!(
                fs,
                "{n},{s},{},{moves},{pos},{gap_m},{left},{right},{refills},{ref_surv},{ref_er},{}",
                rows[0].v,
                last.count
            )
            .unwrap();
        }
        println!("n={n:2}: all s <= {} done in {:.2}s, {} states", 2 * n, t.elapsed().as_secs_f64(), solver.states());
    }
    true
}

/// Minimal CSV reader that understands double-quoted fields.
fn read_csv(path: &str) -> Vec<HashMap<String, String>> {
    let text = std::fs::read_to_string(path).expect("read csv");
    let mut lines = text.lines();
    let split = |line: &str| -> Vec<String> {
        let mut out = Vec::new();
        let mut cur = String::new();
        let mut q = false;
        for ch in line.chars() {
            match ch {
                '"' => q = !q,
                ',' if !q => out.push(std::mem::take(&mut cur)),
                _ => cur.push(ch),
            }
        }
        out.push(cur);
        out
    };
    let header = split(lines.next().unwrap());
    lines.map(|l| header.iter().cloned().zip(split(l)).collect()).collect()
}

fn table(path: &str) -> HashMap<(u32, u32), u32> {
    read_csv(path)
        .into_iter()
        .map(|r| ((r["n"].parse().unwrap(), r["s"].parse().unwrap()), r["value"].parse().unwrap()))
        .collect()
}

/// Item 4: adaptive against oblivious, cell by cell.
pub fn gap_compare(out: &Path, cfg: Config) -> bool {
    let Some(mut f) = create(out) else { return false };
    let forced_t = table("data/values_forced_v1.csv");
    let qv_t = table("data/values_qv_v1.csv");
    let obl = read_csv("data/oblivious_v1.csv");
    let mut forced = Solver::<u64, Forced>::new(cfg);
    let mut qv = Solver::<u64, Qv>::new(cfg);
    // opt_forced: from the v1 table when stored, otherwise computed exactly.
    let opt = |m: u32, r: u32, solver: &mut Solver<u64, Forced>| -> u32 {
        if m == 0 || r == 0 {
            return 0;
        }
        match forced_t.get(&(m, r)) {
            Some(&v) => v,
            None => solver.opt(m, r).expect("table limit"),
        }
    };
    // Pure binary search through the gap lemma: R(n, s) = max(R(n, s-1), 1 + R((n-1)/2, s-1)).
    fn pure(n: u32, s: u32, memo: &mut HashMap<(u32, u32), u32>) -> u32 {
        if n == 0 || s == 0 {
            return 0;
        }
        if let Some(&v) = memo.get(&(n, s)) {
            return v;
        }
        let v = pure(n, s - 1, memo).max(1 + pure((n - 1) / 2, s - 1, memo));
        memo.insert((n, s), v);
        v
    }
    let mut memo = HashMap::new();
    writeln!(
        f,
        "n,s,opt,opt_qv,ov_star,eps_star_lower_times_s,R_binary_search,B_one_gap_step,opt_gt_ov,explanation,pv_refills"
    )
    .unwrap();
    let mut counts: HashMap<&'static str, u32> = HashMap::new();
    for r in &obl {
        let (n, s): (u32, u32) = (r["n"].parse().unwrap(), r["s"].parse().unwrap());
        let o = opt(n, s, &mut forced);
        let oq = match qv_t.get(&(n, s)) {
            Some(&v) => v,
            None => qv.opt(n, s).expect("table limit"),
        };
        let ov = spr_oblivious::rational::parse(&r["ov_star"]).expect("ov* is exact");
        let el = spr_oblivious::rational::parse(&r["eps_lower"]).unwrap();
        let es = &el * BigRational::from_integer((s as i64).into());
        let rr = pure(n, s, &mut memo);
        let b = opt(n, s - 1, &mut forced).max(1 + opt((n - 1) / 2, s - 1, &mut forced));
        let gt = BigRational::from_integer((o as i64).into()) > ov;
        let refills = spr_solve::gap::pv_rows(&mut forced, n, s)
            .map(|rows| rows.iter().filter(|x| x.refill).count())
            .expect("table limit");
        let expl: &'static str = if !gt {
            "opt <= ov*"
        } else if o == rr {
            "binary search via the gap lemma (opt = R)"
        } else if o == b {
            "one gap step then optimal play (opt = B > R)"
        } else {
            "beyond B"
        };
        *counts.entry(expl).or_default() += 1;
        writeln!(
            f,
            "{n},{s},{o},{oq},{},{},{rr},{b},{gt},{expl},{refills}",
            spr_oblivious::rational::fmt(&ov),
            es.to_f64().map_or(String::new(), |v| format!("{v:.10}"))
        )
        .unwrap();
    }
    for (k, v) in &counts {
        println!("{v:3} cells: {k}");
    }
    true
}

/// Exact fixed-policy values against a best-responding Labeler.
fn make_policy(policy: &str) -> Option<Box<dyn spr_solve::policy::Policy>> {
    use spr_solve::policy::{BinarySearch, GapRule, Walker, WalkerBoth, WalkerFresh};
    Some(match policy {
        "binary-search" => Box::new(BinarySearch),
        "walker-right" => Box::new(Walker { gap: GapRule::Right }),
        "walker-centre" => Box::new(Walker { gap: GapRule::Centre }),
        "walker-left" => Box::new(Walker { gap: GapRule::Left }),
        "walker-fresh" => Box::new(WalkerFresh),
        "walker-both" => Box::new(WalkerBoth),
        _ => return None,
    })
}

/// One game of a fixed policy against the Labeler's first minimising reply.
pub fn policy_line(policy: &str, full: bool, n: u32, s: u32) -> bool {
    use spr_core::Board;
    use spr_solve::policy::{Evaluator, worst_reply};
    let Some(pol) = make_policy(policy) else { return false };
    let mut ev = Evaluator::new(pol.as_ref(), full, 50_000_000);
    let v = ev.value_root(n, s).expect("limit");
    println!("{} ({}) n={n} s={s}: value {v}", pol.name(), if full { "full" } else { "forced" });
    let (mut x, mut r) = (Board::<u64>::empty(n), s);
    while r > 0 {
        let w = ev.value(x, r);
        let Some(j) = pol.choose(&x, r) else { break };
        if w == x.count() && ev.value(x, r) == x.count() {
            // stopping is optimal for the Pointer here
            let (_, _, _, cont) = worst_reply(&mut ev, x, r, j);
            if cont <= x.count() {
                println!("  r={r:2} {x}  |x|={} stop", x.count());
                return true;
            }
        }
        let (e, sign, y, _) = worst_reply(&mut ev, x, r, j);
        let er: Vec<u32> = spr_core::cells_of(e).collect();
        println!("  r={r:2} {x}  |x|={} W={w}  point {j}  erase {er:?} write {}", x.count(), cell_char(sign));
        x = y;
        r -= 1;
    }
    println!("  r={r:2} {x}  |x|={} end", x.count());
    true
}

pub fn policy_eval(policy: &str, full: bool, max_n: u32, s_mode: &str, out: &Path, max_states: usize) -> bool {
    use spr_solve::policy::Evaluator;
    let Some(pol) = make_policy(policy) else {
        eprintln!("unknown policy {policy}");
        return false;
    };
    let Some(mut f) = create(out) else { return false };
    writeln!(f, "policy,labeler,n,s,value,states,nodes,seconds").unwrap();
    let lab = if full { "full" } else { "forced" };
    for n in 1..=max_n {
        let s_list: Vec<u32> = match s_mode {
            "diag" => vec![n],
            "2n" => vec![2 * n],
            "both" => vec![n, 2 * n],
            _ => (1..=2 * n).collect(),
        };
        let mut line = Vec::new();
        for s in s_list {
            // A fresh evaluator per (n, s) keeps memory bounded by the largest instance.
            let mut ev = Evaluator::new(pol.as_ref(), full, max_states);
            let t = Instant::now();
            let Some(v) = ev.value_root(n, s) else {
                println!("# STOPPED at n={n} s={s}: not reached (memory), {} states", ev.states());
                return true;
            };
            writeln!(f, "{},{lab},{n},{s},{v},{},{},{:.3}", pol.name(), ev.states(), ev.nodes, t.elapsed().as_secs_f64()).unwrap();
            f.flush().unwrap();
            line.push(format!("s={s}:{v} ({} st, {:.1}s)", ev.states(), t.elapsed().as_secs_f64()));
        }
        println!("{} {lab} n={n:2}: {}", pol.name(), line.join("  "));
    }
    true
}
