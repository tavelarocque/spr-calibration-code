//! Exhaustive check of the packed board against the naive Vec<u8> model on every
//! board of length n <= 9, including burned cells, for both word widths.

use spr_core::naive;
use spr_core::{Board, VariantKind, Word, cells_of};

fn all_boards(n: usize, codes: u8) -> impl Iterator<Item = Vec<u8>> {
    let total = (codes as u64).pow(n as u32);
    (0..total).map(move |mut k| {
        let mut v = vec![0u8; n];
        for c in v.iter_mut() {
            *c = (k % codes as u64) as u8;
            k /= codes as u64;
        }
        v
    })
}

fn check_board<W: Word>(x: &[u8]) {
    let b = Board::<W>::from_cells(x);
    assert_eq!(b.to_cells(), x);
    assert_eq!(b.count(), naive::count(x));
    let empt: Vec<usize> = cells_of(b.empties()).map(|j| j as usize).collect();
    assert_eq!(empt, naive::empties(x));

    // Reflection.
    assert_eq!(b.reflect().to_cells(), naive::reflect(x), "reflect {x:?}");
    let c = b.canon();
    assert!(c == b || c == b.reflect());

    // Banking: the stripped cells are exactly the maximal '+' prefix and '-' suffix.
    let (p, q) = naive::bank_lengths(x);
    let (bank, core) = b.bank();
    assert_eq!(bank as usize, p + q, "bank {x:?}");
    assert_eq!(core.to_cells(), x[p..x.len() - q].to_vec(), "core {x:?}");

    // QV reduction: same bank, same core up to reflection.
    let (qb, qcore) = naive::qv_reduce(x);
    let (pb, pcore) = b.reduce_qv();
    assert_eq!(pb, qb, "qv bank {x:?}");
    let pc = pcore.to_cells();
    assert!(pc == qcore || pc == naive::reflect(&qcore), "qv core {x:?}: {pc:?} vs {qcore:?}");
    assert_eq!(pcore, pcore.canon());

    // Move generation and successors, for every variant.
    for kind in VariantKind::ALL {
        for j in naive::empties(x) {
            let er = naive::erasable(x, j);
            let m = b.erasable(j as u32);
            let packed: Vec<usize> = cells_of(m).map(|i| i as usize).collect();
            assert_eq!(packed, er, "erasable {x:?} j={j}");
            for e in naive::erasure_options(&er, kind) {
                let mut em = W::ZERO;
                for &i in &e {
                    em = em | (W::ONE << (2 * i as u32));
                }
                for sign in [spr_core::PLUS, spr_core::MINUS] {
                    let y = naive::apply(x, j, &e, sign, kind);
                    let by = b.apply(j as u32, em, sign, kind.burn());
                    assert_eq!(by.to_cells(), y);
                }
            }
        }
    }
}

#[test]
fn packed_matches_naive_all_boards_up_to_9() {
    let mut checked = 0u64;
    for n in 0..=9 {
        for x in all_boards(n, 4) {
            check_board::<u64>(&x);
            check_board::<u128>(&x);
            checked += 1;
        }
    }
    assert_eq!(checked, (0..=9).map(|n| 4u64.pow(n)).sum::<u64>());
}

#[test]
fn full_width_boards() {
    // Edge cases at the word boundary: all '+', all '-', mixed, at 32 and 64 cells.
    for n in [31usize, 32] {
        for x in [vec![1u8; n], vec![2u8; n], (0..n).map(|i| (i % 3) as u8).collect()] {
            check_board::<u64>(&x);
            check_board::<u128>(&x);
        }
    }
    for n in [63usize, 64] {
        for x in [vec![1u8; n], vec![2u8; n], (0..n).map(|i| (i % 4) as u8).collect()] {
            let b = Board::<u128>::from_cells(&x);
            assert_eq!(b.reflect().to_cells(), naive::reflect(&x));
            let (p, q) = naive::bank_lengths(&x);
            assert_eq!(b.bank().0 as usize, p + q);
        }
    }
}
