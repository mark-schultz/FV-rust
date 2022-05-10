#![allow(dead_code, unused_variables)]

mod arith;
mod encrypt;
mod utils;

// Parems once I figure them out
const SIGMA: f64 = 1000.;
// The exponent of the degree of the polynomail ring, i.e. we are working
// mod x^d+1 for d = 2^DEGREE_EXP.
const DEGREE_EXP: usize = 10; // 10;
const DEGREE: usize = 1 << DEGREE_EXP;
const DELTA_R: usize = DEGREE;
#[allow(dead_code)]
const Q_EXP: usize = 1359; // 1359;
const Q_EXP_BYTES_PLUS_1: usize = Q_EXP / 8 + 1;
// Can't construct Q in a const context, it is always Integer::One() << N though
// Same with t, although t = 2 is chosen in the paper
const T_EXP: usize = 1;
const T_EXP_BYTES_PLUS_1: usize = T_EXP / 8 + 1;

// Small RELIN terms
// const RELIN_EXP: usize = 1;
// const RELIN_TERMS: usize = Q_EXP; // floor(log_T(q))

// Large (ish) RELIN terms
const RELIN_EXP: usize = Q_EXP >> 1;
const RELIN_TERMS: usize = 2; // floor(log_T(q))

#[derive(Debug)]
pub enum Errors {
    ParseIntError,
    ParsePolyError,
    RNGError,
}
