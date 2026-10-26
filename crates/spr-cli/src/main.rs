//! `spr`: tabulate exact SPR values, check them against the Python oracle, search
//! for erasure-lemma counterexamples, and print principal variations.

mod gapcmd;
mod obl;
mod obl_export;
mod oracle;

use clap::{Parser, Subcommand};
use spr_core::{Board, Forced, Full, MINUS, Qv, Variant, VariantKind, Word, cell_char};
use spr_solve::erasure::ErasureSearch;
use spr_solve::strategy::{analyse, principal_variation};
use spr_solve::{Config, Solver};
use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser)]
#[command(name = "spr", about = "Exact solver for Sign Preservation with Reuse")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(clap::Args, Clone, Copy)]
struct SolverArgs {
    /// Disable the monotonicity-in-r bounds.
    #[arg(long)]
    no_monotone: bool,
    /// Abort once the table holds this many states.
    #[arg(long, default_value_t = usize::MAX)]
    max_states: usize,
    /// Allocate the table for --max-states up front.
    #[arg(long)]
    reserve: bool,
    /// QV: use plain banking only (no burned-cell deletion, no permanent-sign banking).
    #[arg(long)]
    no_qv_reduce: bool,
}

impl SolverArgs {
    fn config(self) -> Config {
        Config {
            monotone: !self.no_monotone,
            max_states: self.max_states,
            reserve: self.reserve,
            qv_reduce: !self.no_qv_reduce,
        }
    }
}

#[derive(Subcommand)]
enum Cmd {
    /// Solve opt(n, s) for n in min-n..=max-n and s in 1..=2n, with one shared table.
    Tabulate {
        #[arg(long)]
        variant: VariantKind,
        #[arg(long, default_value_t = 1)]
        min_n: u32,
        #[arg(long)]
        max_n: u32,
        /// Cap on s (default 2n).
        #[arg(long)]
        max_s: Option<u32>,
        /// Solve s = 1..=S for every n, ignoring the 2n default and --max-s (used for the
        /// out-of-table values opt(g, r) with r > 2g that the binary-search check needs).
        #[arg(long)]
        s_hi: Option<u32>,
        /// CSV output (n,s,value,states_stored,nodes,seconds). Never overwritten.
        #[arg(long)]
        out: Option<PathBuf>,
        #[command(flatten)]
        solver: SolverArgs,
    },
    /// Compare against reference/spr_ref.py (run live) and data/spr_values_v0.md.
    CheckOracle {
        /// Skip running Python; compare against the v0 tables only.
        #[arg(long)]
        no_python: bool,
    },
    /// Search the full variant's reachable states for a strictly better non-maximal erasure.
    FindErasureCounterexample {
        #[arg(long, default_value_t = 1)]
        min_n: u32,
        #[arg(long)]
        max_n: u32,
        /// Cap on s (default 2n).
        #[arg(long)]
        max_s: Option<u32>,
        #[command(flatten)]
        solver: SolverArgs,
    },
    /// All optimal first moves and one principal variation for opt(n, s).
    PrincipalVariation {
        n: u32,
        s: u32,
        variant: VariantKind,
        /// Start from this board instead of 0^n (cells 0 + - #).
        #[arg(long)]
        board: Option<String>,
    },
    /// Every optimal Pointer first move for opt(n, s), n <= max-n, s <= 2n, as CSV.
    FirstMoves {
        #[arg(long)]
        variant: VariantKind,
        #[arg(long)]
        max_n: u32,
        /// CSV output. Never overwritten.
        #[arg(long)]
        out: PathBuf,
    },
    /// A forced-variant optimal line annotated with the QV value of the same history.
    ReuseLine { n: u32, s: u32 },
    /// Exhaustive check of the gap lemma (forced variant).
    GapCheck {
        #[arg(long)]
        max_n: u32,
        #[command(flatten)]
        solver: SolverArgs,
    },
    /// Principal variations with gap decomposition and refill events.
    GapPv {
        #[arg(long)]
        max_n: u32,
        #[arg(long)]
        out: PathBuf,
        #[arg(long)]
        summary: PathBuf,
        /// Only s = n and s = 2n.
        #[arg(long)]
        diag_only: bool,
        #[command(flatten)]
        solver: SolverArgs,
    },
    /// Adaptive against oblivious, cell by cell.
    GapCompare {
        #[arg(long)]
        out: PathBuf,
        #[command(flatten)]
        solver: SolverArgs,
    },
    /// Exact value of a fixed Pointer policy against a best-responding Labeler.
    PolicyEval {
        #[arg(long)]
        policy: String,
        /// Labeler may erase any subset (lower bound on opt_full); default forced.
        #[arg(long)]
        full: bool,
        #[arg(long)]
        max_n: u32,
        /// diag (s = n), 2n, both, or all (s = 1..2n).
        #[arg(long, default_value = "both")]
        s_mode: String,
        #[arg(long)]
        out: PathBuf,
        #[arg(long)]
        max_states: usize,
    },
    /// Print one game of a fixed policy against a best-responding Labeler.
    PolicyLine {
        #[arg(long)]
        policy: String,
        #[arg(long)]
        full: bool,
        n: u32,
        s: u32,
    },
    /// The oblivious no-reuse Pointer.
    #[command(subcommand)]
    Obl(obl::Obl),
    /// Qiao-Valiant Lemma 5 check: opt_qv(2^t - 1, t) = t, and report opt_qv(2^t - 2, t).
    KnownValues {
        #[arg(long, default_value_t = 5)]
        max_t: u32,
    },
}

/// Run `$body` with `W` and `V` bound to the word type and variant chosen at run time.
macro_rules! dispatch {
    ($kind:expr, $wide:expr, $f:ident ( $($arg:expr),* )) => {
        match ($kind, $wide) {
            (VariantKind::Full, false) => $f::<u64, Full>($($arg),*),
            (VariantKind::Forced, false) => $f::<u64, Forced>($($arg),*),
            (VariantKind::Qv, false) => $f::<u64, Qv>($($arg),*),
            (VariantKind::Full, true) => $f::<u128, Full>($($arg),*),
            (VariantKind::Forced, true) => $f::<u128, Forced>($($arg),*),
            (VariantKind::Qv, true) => $f::<u128, Qv>($($arg),*),
        }
    };
}

fn main() {
    let cli = Cli::parse();
    let ok = match cli.cmd {
        Cmd::Tabulate { variant, min_n, max_n, max_s, s_hi, out, solver } => {
            dispatch!(variant, max_n > 32, tabulate(min_n, max_n, max_s, s_hi, out, solver))
        }
        Cmd::CheckOracle { no_python } => oracle::check(no_python),
        Cmd::FindErasureCounterexample { min_n, max_n, max_s, solver } => {
            if max_n > 32 {
                find_erasure::<u128>(min_n, max_n, max_s, solver)
            } else {
                find_erasure::<u64>(min_n, max_n, max_s, solver)
            }
        }
        Cmd::PrincipalVariation { n, s, variant, board } => {
            let wide = board.as_ref().map_or(n, |b| b.chars().count() as u32) > 32;
            dispatch!(variant, wide, pv(n, s, board.as_deref()))
        }
        Cmd::KnownValues { max_t } => known_values(max_t),
        Cmd::FirstMoves { variant, max_n, out } => dispatch!(variant, max_n > 32, first_moves(max_n, out)),
        Cmd::ReuseLine { n, s } => reuse_line(n, s),
        Cmd::Obl(c) => obl::run(c),
        Cmd::PolicyLine { policy, full, n, s } => gapcmd::policy_line(&policy, full, n, s),
        Cmd::PolicyEval { policy, full, max_n, s_mode, out, max_states } => {
            gapcmd::policy_eval(&policy, full, max_n, &s_mode, &out, max_states)
        }
        Cmd::GapCheck { max_n, solver } => gapcmd::gap_check(max_n, solver.config()),
        Cmd::GapPv { max_n, out, summary, diag_only, solver } => gapcmd::gap_pv(max_n, &out, &summary, solver.config(), diag_only),
        Cmd::GapCompare { out, solver } => gapcmd::gap_compare(&out, solver.config()),
    };
    if !ok {
        std::process::exit(1);
    }
}

fn tabulate<W: Word, V: Variant>(
    min_n: u32,
    max_n: u32,
    max_s: Option<u32>,
    s_hi_all: Option<u32>,
    out: Option<PathBuf>,
    args: SolverArgs,
) -> bool {
    let mut csv = match &out {
        Some(p) => {
            if p.exists() {
                eprintln!("refusing to overwrite {}; write a new version instead", p.display());
                return false;
            }
            let mut f = File::create(p).expect("create csv");
            writeln!(f, "n,s,value,states_stored,nodes,seconds").unwrap();
            Some(f)
        }
        None => None,
    };
    let mut solver = Solver::<W, V>::new(args.config());
    let t_all = Instant::now();
    println!("# {} variant, rows n, columns s = 1..2n, monotone={}", V::KIND, !args.no_monotone);
    for n in min_n..=max_n {
        let s_hi = s_hi_all.unwrap_or_else(|| max_s.map_or(2 * n, |m| m.min(2 * n)));
        let mut row = Vec::new();
        let t_row = Instant::now();
        for s in 1..=s_hi {
            let nodes0 = solver.nodes;
            let t = Instant::now();
            let Some(v) = solver.opt(n, s) else {
                println!("n={n:2}  {}", row.join(" "));
                println!(
                    "# STOPPED at n={n} s={s}: table limit of {} states reached after {:.1}s",
                    solver.states(),
                    t_all.elapsed().as_secs_f64()
                );
                return true;
            };
            let secs = t.elapsed().as_secs_f64();
            if let Some(f) = csv.as_mut() {
                writeln!(f, "{n},{s},{v},{},{},{secs:.6}", solver.states(), solver.nodes - nodes0).unwrap();
                f.flush().unwrap();
            }
            row.push(format!("{v:2}"));
        }
        println!(
            "n={n:2}  {}    [{:.2}s, {} states, {:.2} GB table]",
            row.join(" "),
            t_row.elapsed().as_secs_f64(),
            solver.states(),
            solver.table_bytes() as f64 / 1e9
        );
    }
    println!("# total {:.2}s, {} states, {} nodes", t_all.elapsed().as_secs_f64(), solver.states(), solver.nodes);
    let (by_r, by_len) = solver.histogram();
    println!("# stored states by r: {by_r:?}");
    println!("# stored states by core length: {by_len:?}");
    true
}

fn find_erasure<W: Word>(min_n: u32, max_n: u32, max_s: Option<u32>, args: SolverArgs) -> bool {
    let mut search = ErasureSearch::new(Solver::<W, Full>::new(args.config()));
    let t_all = Instant::now();
    let (mut total_pairs, mut total_checks, mut total_ptr) = (0u64, 0u64, 0u64);
    for n in min_n..=max_n {
        let s = max_s.map_or(2 * n, |m| m.min(2 * n));
        let t = Instant::now();
        let Some((st, cx)) = search.check_row(n, s) else {
            println!("# STOPPED at n={n}: solver table limit reached");
            break;
        };
        total_pairs += st.pairs;
        total_checks += st.core_checks;
        total_ptr += st.pointer_checks;
        println!(
            "n={n:2} s<={s:2}: {} reachable boards, {} reachable (x, r) pairs, {} new (core, r) checked, {} pointer comparisons, {:.2}s",
            st.boards,
            st.pairs,
            st.core_checks,
            st.pointer_checks,
            t.elapsed().as_secs_f64()
        );
        if let Some(cx) = cx {
            println!("\n*** COUNTEREXAMPLE TO THE ERASURE LEMMA ***");
            println!("n = {}, r = {}, board {} , pointer j = {}", cx.n, cx.r, cx.board, cx.j);
            println!(
                "best over all erasure sets = {} (E = {:?}, sign {}), best under maximal erasure = {}",
                cx.best_any,
                cx.best_any_option.0,
                cell_char(cx.best_any_option.1),
                cx.best_max
            );
            println!("reached from 0^{} by:", cx.n);
            for m in &cx.path {
                println!("  {}  point {}  erase {:?}  write {}", m.before, m.j, m.erased, cell_char(m.sign));
            }
            return true;
        }
    }
    println!(
        "# no counterexample: {total_pairs} reachable (x, r) pairs, {total_checks} distinct (core, r) checks, {total_ptr} pointer comparisons, {:.2}s",
        t_all.elapsed().as_secs_f64()
    );
    true
}

/// `prefix` followed by the board as `[ + 0 - ]`, with a caret under cell `mark`.
fn diagram<W: Word>(prefix: &str, x: Board<W>, mark: Option<u32>) -> String {
    let cells: String = (0..x.len).map(|i| format!(" {}", cell_char(x.cell(i)))).collect();
    let mut s = format!("{prefix}[{cells} ]");
    if let Some(j) = mark {
        s.push('\n');
        s.push_str(&" ".repeat(prefix.len() + 2 + 2 * j as usize));
        s.push_str(&format!("^ j={j}"));
    }
    s
}

fn pv<W: Word, V: Variant>(n: u32, s: u32, board: Option<&str>) -> bool {
    let x = match board {
        Some(b) => match Board::<W>::parse(b) {
            Some(x) => x,
            None => {
                eprintln!("bad board {b:?}");
                return false;
            }
        },
        None => Board::empty(n),
    };
    let mut solver = Solver::<W, V>::new(Config::default());
    let a = analyse(&mut solver, x, s).unwrap();
    println!("{} variant, board {x} (n = {}), s = {s}: value {}", V::KIND, x.len, a.value);
    if a.stop_optimal {
        println!("stopping now is optimal (|x| = {})", x.count());
    }
    let opt = a.optimal_pointers();
    println!("optimal first pointers ({} of {}): {:?}", opt.len(), a.pointer_values.len(), opt);
    let all: Vec<String> = a.pointer_values.iter().map(|(j, v)| format!("{j}:{v}")).collect();
    println!("Labeler value of each pointer: {}", all.join(" "));
    println!("\nprincipal variation:");
    for st in principal_variation(&mut solver, x, s).unwrap() {
        match st.action {
            Some((j, e, sign)) => {
                let prefix = format!("r={:2} V={} |x|={}  ", st.r, st.value, st.board.count());
                println!("{}", diagram(&prefix, st.board, Some(j)));
                let sgn = if sign == MINUS { '-' } else { '+' };
                println!("       optimal pointers {:?}; Labeler erases {e:?}, writes {sgn}", st.optimal);
            }
            None => {
                let prefix = format!("r={:2} V={} |x|={}  ", st.r, st.value, st.board.count());
                println!("{}", diagram(&prefix, st.board, None));
                println!("       Pointer stops, payoff {}", st.board.count());
            }
        }
    }
    true
}

fn known_values(max_t: u32) -> bool {
    let mut ok = true;
    for t in 1..=max_t {
        let n = (1u32 << t) - 1;
        let mut s = Solver::<u64, Qv>::new(Config::default());
        let t0 = Instant::now();
        let v = s.opt(n, t).unwrap();
        let w = s.opt(n - 1, t).unwrap();
        let pass = v == t;
        ok &= pass;
        println!(
            "t={t}: opt_qv({n}, {t}) = {v} [{}]   opt_qv({}, {t}) = {w} (expected t-1 = {}, reported only)   {:.2}s",
            if pass { "PASS" } else { "FAIL" },
            n - 1,
            t - 1,
            t0.elapsed().as_secs_f64()
        );
    }
    ok
}

fn first_moves<W: Word, V: Variant>(max_n: u32, out: PathBuf) -> bool {
    if out.exists() {
        eprintln!("refusing to overwrite {}; write a new version instead", out.display());
        return false;
    }
    let mut f = File::create(&out).expect("create csv");
    writeln!(f, "n,s,value,optimal_first_moves,pointer_values").unwrap();
    let mut solver = Solver::<W, V>::new(Config::default());
    for n in 1..=max_n {
        let t = Instant::now();
        for s in 1..=2 * n {
            let a = analyse(&mut solver, Board::empty(n), s).unwrap();
            let opt: Vec<String> = a.optimal_pointers().iter().map(|j| j.to_string()).collect();
            let pv: Vec<String> = a.pointer_values.iter().map(|(_, v)| v.to_string()).collect();
            writeln!(f, "{n},{s},{},{},{}", a.value, opt.join(" "), pv.join(" ")).unwrap();
        }
        println!("n={n:2} done in {:.2}s", t.elapsed().as_secs_f64());
    }
    true
}

fn reuse_line(n: u32, s: u32) -> bool {
    let mut forced = Solver::<u64, Forced>::new(Config::default());
    let mut qv = Solver::<u64, Qv>::new(Config::default());
    let steps = spr_solve::strategy::reuse_line(&mut forced, &mut qv, n, s).unwrap();
    println!(
        "forced line for n={n} s={s}: opt_forced = {}, opt_qv = {}",
        steps[0].v_forced,
        qv.opt(n, s).unwrap()
    );
    println!("used cells shown as u in the second row; V_qv is the QV value of the same history\n");
    for st in &steps {
        let used: String = (0..st.board.len)
            .map(|i| if st.used >> (2 * i) & 1 == 1 { " u" } else { "  " })
            .collect();
        let (fresh, reused) = st.best_fresh_and_reused();
        let fmt = |v: Option<u32>| v.map_or("-".to_string(), |v| v.to_string());
        let prefix = format!("r={:2} V={} Vqv={} |x|={}  ", st.r, st.v_forced, st.v_qv, st.board.count());
        let mark = st.action.as_ref().map(|a| a.0);
        println!("{}", diagram(&prefix, st.board, None));
        println!("{}[{used} ]", " ".repeat(prefix.len()));
        if let Some(j) = mark {
            println!("{}^ j={j}", " ".repeat(prefix.len() + 2 + 2 * j as usize));
        }
        let vals: Vec<String> = st.pointer_values.iter().map(|(j, v)| format!("{j}:{v}")).collect();
        println!(
            "       best value via a fresh cell {}, via a reused cell {}; pointer values {}",
            fmt(fresh),
            fmt(reused),
            vals.join(" ")
        );
        match &st.action {
            Some((j, e, sign)) => {
                let reuse = if st.used >> (2 * j) & 1 == 1 { "  <-- REUSED CELL" } else { "" };
                let sgn = if *sign == MINUS { '-' } else { '+' };
                println!("       Pointer {j}{reuse}; Labeler erases {e:?}, writes {sgn}\n");
            }
            None => println!("       Pointer stops, payoff {}", st.board.count()),
        }
    }
    true
}
