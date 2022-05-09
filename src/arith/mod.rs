//! Cyclotomic Polynomial Arithmetic
//! In particular, arithmetic in Z[x]/(f(x)) for f(x) = x^d + 1 for d = 2^k
//!
//! Sample parameters (just before Section 7) are q = 2^n and d = 2^k for k= 10, n > 1358.

use crate::Errors;
use num_bigint::{BigInt, Sign};
use num_traits::identities::{One, Zero};
use std::{
    ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign},
    str::FromStr,
};

#[derive(Clone, Eq, PartialEq, Debug)]
struct Integer(BigInt);

impl Add<Integer> for Integer {
    type Output = Integer;
    fn add(self, rhs: Self) -> Self::Output {
        Integer(self.0 + rhs.0)
    }
}

impl AddAssign<Integer> for Integer {
    fn add_assign(&mut self, rhs: Integer) {
        self.0 += rhs.0
    }
}

impl Neg for Integer {
    type Output = Integer;
    fn neg(self) -> Self::Output {
        Integer(-self.0)
    }
}

impl Sub<Integer> for Integer {
    type Output = Integer;
    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl SubAssign<Integer> for Integer {
    fn sub_assign(&mut self, rhs: Integer) {
        self.0 -= rhs.0
    }
}

impl Mul<Integer> for Integer {
    type Output = Integer;
    fn mul(self, rhs: Self) -> Self::Output {
        Integer(self.0 * rhs.0)
    }
}

impl Zero for Integer {
    fn zero() -> Self {
        Self(BigInt::zero())
    }
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl One for Integer {
    fn one() -> Self {
        Self(BigInt::one())
    }
}

/// From a Base10 string, in Big Endian (msb first) order.
impl FromStr for Integer {
    type Err = Errors;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (sign, body) = match s.strip_prefix("-") {
            None => (Sign::Minus, s),
            Some(x) => (Sign::Plus, x),
        };
        let body: Vec<u8> = body
            .chars()
            .flat_map(|c| c.to_digit(10).ok_or(Errors::ParseIntError).map(|x| x as u8))
            .collect();
        let res = BigInt::from_radix_be(sign, &body, 10);
        res.ok_or(Errors::ParseIntError).map(|i| Integer(i))
    }
}

trait Ring: Add + AddAssign + Sub + SubAssign + Neg + Mul + Zero + One + Clone {}
impl<T> Ring for T where T: Add + AddAssign + Sub + SubAssign + Neg + Mul + Zero + One + Clone {}

/// The exponent of the degree of the polynomail ring, i.e. we are working
/// mod x^d+1 for d = 2^DEGREE_EXP.
const DEGREE_EXP: usize = 10;
const DEGREE: usize = 1 << DEGREE_EXP;

/// Index into coeffs is the degree, i.e. coeffs[0] is constant term.
#[derive(Clone, Eq, PartialEq, Debug)]
struct CycloPoly<R: Ring> {
    coeffs: Vec<R>,
}

impl<R: Ring> Zero for CycloPoly<R> {
    fn zero() -> Self::Output {
        let mut coeffs: Vec<R> = Vec::with_capacity(DEGREE);
        for _ in 0..DEGREE {
            coeffs.push(R::zero());
        }
        Self { coeffs }
    }
    fn is_zero(&self) -> bool {
        for val in self.coeffs.iter() {
            if !val.is_zero() {
                return false;
            }
        }
        true
    }
}

impl<R: Ring> One for CycloPoly<R> {
    fn one() -> Self::Output {
        let mut output = Self::zero();
        output.coeffs[0] = R::one();
        output
    }
}

impl<R: Ring> Add<CycloPoly<R>> for CycloPoly<R> {
    type Output = CycloPoly<R>;
    fn add(self, rhs: CycloPoly<R>) -> Self::Output {
        let coeffs: Vec<R> = self
            .coeffs
            .into_iter()
            .zip(rhs.coeffs.into_iter())
            .map(|(a, b)| a + b)
            .collect();
        Self { coeffs }
    }
}

impl<R: Ring + Sub<Output = R>> Sub<CycloPoly<R>> for CycloPoly<R> {
    type Output = CycloPoly<R>;
    fn sub(self, rhs: CycloPoly<R>) -> Self::Output {
        let coeffs: Vec<R> = self
            .coeffs
            .into_iter()
            .zip(rhs.coeffs.into_iter())
            .map(|(a, b)| a - b)
            .collect();
        Self { coeffs }
    }
}

/// Multiplication specialized to mod x^d+1.
/// Does (full) polynomial multiplication, then manually reduces things
impl<R: Ring> Mul<CycloPoly<R>> for CycloPoly<R> {
    type Output = CycloPoly<R>;
    fn mul(self, rhs: CycloPoly<R>) -> Self::Output {
        let mut unreduced_res = CycloPoly::<R>::zero();
        unreduced_res.coeffs.extend(CycloPoly::<R>::zero().coeffs);
        for (i, val) in unreduced_res.coeffs.iter_mut().enumerate() {
            for j in 0..=i {
                // When indexing should work
                if j < DEGREE && i < DEGREE + j {
                    *val += self.coeffs[j].clone() * rhs.coeffs[i - j].clone();
                }
            }
        }
        // reduction step
        let mut res = CycloPoly::<R>::zero();
        for (i, val) in res.coeffs.iter_mut().enumerate() {
            // Would have expected these two operations to be the other way around
            // But it matches up with my SAGE test cases.
            // Worst case this introduces a sign error (that will be consistent throughout
            // hopefully)
            *val -= unreduced_res.coeffs[i].clone();
            *val += unreduced_res.coeffs[DEGREE + i].clone();
        }
        res
    }
}

/// Input is [a, b, c, ..., n], where a,b,c all implement FromStr,
impl<R: Ring + FromStr> FromStr for CycloPoly<R> {
    type Err = Errors;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut coeffs = Vec::new();
        let stripped = s
            .strip_prefix("[")
            .ok_or(Errors::ParsePolyError)?
            .strip_suffix("]")
            .ok_or(Errors::ParsePolyError)?;
        for val in stripped.split(", ") {
            let parsed_val = R::from_str(&val).map_err(|_| Errors::ParseIntError)?;
            coeffs.push(parsed_val);
        }
        Ok(CycloPoly { coeffs })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    /// Reads test cases (generated in Sage) from a given file.
    /// Test cases are formatted as
    ///
    /// p1_coeffs
    /// p2_coeffs
    /// p3_coeffs
    /// p1'_coeffs...
    ///
    /// where p3 = p1 (some operation) p2.

    fn read_from_file(
        path: &str,
    ) -> Vec<(CycloPoly<Integer>, CycloPoly<Integer>, CycloPoly<Integer>)> {
        let contents = fs::read_to_string(path).unwrap().to_string();
        let mut results = Vec::new();
        let lines: Vec<&str> = contents.lines().collect();
        for batch in lines.as_slice().chunks(3) {
            // Why isn't chunks defined directly on `Iterator?`
            let p1 = CycloPoly::<Integer>::from_str(batch[0]).expect("p1 parses");
            let p2 = CycloPoly::<Integer>::from_str(batch[1]).expect("p2 parses");
            let p3 = CycloPoly::<Integer>::from_str(batch[2]).expect("p3 parses");

            results.push((p1, p2, p3));
        }
        results
    }

    #[test]
    fn test_addition() {
        let test_cases = read_from_file("src/arith/add_test_cases.txt");
        for (p1, p2, p3) in test_cases.into_iter() {
            assert_eq!(p1 + p2, p3);
        }
    }
    #[test]
    fn test_multiplication() {
        let test_cases = read_from_file("src/arith/mul_test_cases.txt");
        for (p1, p2, p3) in test_cases.into_iter() {
            assert_eq!(p1 * p2, p3);
        }
    }
}
