#![feature(trace_macros)]

mod arith;
mod encrypt;
mod utils;

// Parems once I figure them out
#[allow(dead_code)]
const SIGMA: f64 = 2.;
// The exponent of the degree of the polynomail ring, i.e. we are working
// mod x^d+1 for d = 2^DEGREE_EXP.
const DEGREE_EXP: usize = 10;
const DEGREE: usize = 1 << DEGREE_EXP;
#[allow(dead_code)]
const N: usize = 1359;
// Can't construct Q in a const context, it is always Integer::One() << N though
// Same with t, although t = 2 is chosen in the paper
const T: usize = 1;

#[derive(Debug)]
pub enum Errors {
    ParseIntError,
    ParsePolyError,
    RNGError,
}
