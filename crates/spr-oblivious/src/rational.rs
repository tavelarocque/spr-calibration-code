//! Turning floats into exact rationals.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};

/// The exact value of a finite f64 as a rational (every f64 is a dyadic rational).
pub fn exact_f64(x: f64) -> BigRational {
    BigRational::from_float(x).expect("finite float")
}

/// `round(x * 2^k)` as an integer.
pub fn dyadic_numerator(x: f64, k: u32) -> BigInt {
    let scaled = exact_f64(x) * BigRational::from_integer(BigInt::one() << k);
    scaled.round().to_integer()
}

/// The best rational approximation of `x` with denominator at most `max_den`
/// (continued fractions, last convergent within the bound).
pub fn best_rational(x: f64, max_den: u64) -> BigRational {
    let target = exact_f64(x);
    let (mut h0, mut h1) = (BigInt::zero(), BigInt::one());
    let (mut k0, mut k1) = (BigInt::one(), BigInt::zero());
    let mut r = target.clone();
    let max = BigInt::from(max_den);
    for _ in 0..64 {
        let a = r.floor().to_integer();
        let h2 = &a * &h1 + &h0;
        let k2 = &a * &k1 + &k0;
        if k2 > max {
            break;
        }
        (h0, h1) = (h1, h2);
        (k0, k1) = (k1, k2);
        let frac = &r - BigRational::from_integer(a);
        if frac.is_zero() {
            break;
        }
        r = frac.recip();
    }
    if k1.is_zero() {
        return target.round();
    }
    BigRational::new(h1, k1)
}

pub fn to_f64(r: &BigRational) -> f64 {
    r.to_f64().unwrap_or(f64::NAN)
}

/// `p/q` text form.
pub fn fmt(r: &BigRational) -> String {
    if r.denom().is_one() { r.numer().to_string() } else { format!("{}/{}", r.numer(), r.denom()) }
}

pub fn parse(s: &str) -> Option<BigRational> {
    let s = s.trim();
    match s.split_once('/') {
        Some((a, b)) => {
            let (a, b): (BigInt, BigInt) = (a.trim().parse().ok()?, b.trim().parse().ok()?);
            (!b.is_zero()).then(|| BigRational::new(a, b))
        }
        None => s.parse::<BigInt>().ok().map(BigRational::from_integer),
    }
}

pub fn is_nonneg(r: &BigRational) -> bool {
    !r.is_negative()
}
