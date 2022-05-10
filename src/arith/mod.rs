//! Cyclotomic Polynomial Arithmetic
//! In particular, arithmetic in Z[x]/(f(x)) for f(x) = x^d + 1 for d = 2^k
//!
//! Sample parameters (just before Section 7) are q = 2^n and d = 2^k for k= 10, n > 1358.
#![allow(dead_code)]

use crate::{
    utils::{CDT_table, RandGenerator},
    Errors, DEGREE,
};
use num_bigint::{BigInt, Sign};
use num_traits::identities::{One, Zero};
use std::{
    ops::{
        Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Shl, ShlAssign, Shr,
        ShrAssign, Sub, SubAssign,
    },
    str::FromStr,
};

#[derive(Clone, Eq, PartialEq, Debug)]
pub(crate) struct Integer(pub(crate) BigInt);

#[macro_export]
macro_rules! impl_assign {
    ( $trait_name:ident, $fun_name:ident, $type: ty, $backing_type:ty) => {
        impl $trait_name<&$type> for $type {
            fn $fun_name(&mut self, rhs: &$type) {
                <$backing_type>::$fun_name(&mut self.0, &rhs.0);
            }
        }
    };
}
use impl_assign;

#[macro_export]
macro_rules! impl_consuming {
    ( $trait_name:ident, $assigning_fun_name:ident, $fun_name:ident, $type: ty) => {
        impl $trait_name<$type> for $type {
            type Output = $type;
            fn $fun_name(self, rhs: $type) -> Self::Output {
                let mut output = self.clone();
                output.$assigning_fun_name(&rhs);
                output
            }
        }
    };
}
use impl_consuming;

#[macro_export]
macro_rules! impl_op {
    ( $assigning_trait_name:ident, $trait_name:ident, $assigning_fun_name:ident, $fun_name:ident, $type: ty, $backing_type:ty) => {
        impl_assign!(
            $assigning_trait_name,
            $assigning_fun_name,
            $type,
            $backing_type
        );
        impl_consuming!($trait_name, $assigning_fun_name, $fun_name, $type);
    };
}
use impl_op;

impl_op!(AddAssign, Add, add_assign, add, Integer, BigInt);
impl_op!(SubAssign, Sub, sub_assign, sub, Integer, BigInt);
impl_op!(MulAssign, Mul, mul_assign, mul, Integer, BigInt);
impl_op!(DivAssign, Div, div_assign, div, Integer, BigInt);
impl_op!(RemAssign, Rem, rem_assign, rem, Integer, BigInt);

impl ShlAssign<usize> for Integer {
    fn shl_assign(&mut self, rhs: usize) {
        self.0 <<= rhs;
    }
}
impl Shl<usize> for Integer {
    type Output = Integer;
    fn shl(self, rhs: usize) -> Self::Output {
        let mut output = self.clone();
        output <<= rhs;
        output
    }
}

impl ShrAssign<usize> for Integer {
    fn shr_assign(&mut self, rhs: usize) {
        self.0 >>= rhs;
    }
}
impl Shr<usize> for Integer {
    type Output = Integer;
    fn shr(self, rhs: usize) -> Self::Output {
        let mut output = self.clone();
        output >>= rhs;
        output
    }
}

impl Integer {
    /// Centered reduction to an integer in (-q/2, q/2].
    ///
    /// Formula via case analysis of parity of q mod 2
    fn modulo_assign(&mut self, rhs: &Integer) {
        dbg!(&self);
        let mut shift = rhs.clone();
        shift.0 -= 1_i32;
        shift >>= 1;
        dbg!(&shift);
        *self += &shift;
        *self %= rhs;
        *self -= &shift;
        // hacky
        if self.0.sign() == Sign::Minus {
            *self += rhs;
            *self += &shift;
            *self %= rhs;
            *self -= &shift;
        }
    }
    fn modulo(self, rhs: Integer) -> Self {
        let mut output = self.clone();
        Integer::modulo_assign(&mut output, &rhs);
        output
    }
    fn uniform_sample<const MOD_EXP_BYTES_PLUS_ONE: usize>(
        rng: &RandGenerator,
        modulus: &Integer,
    ) -> Self {
        let mut buff = [0 as u8; MOD_EXP_BYTES_PLUS_ONE];
        rng.fill(&mut buff)
            .expect("RNG Error: Will not try to recover from");
        let mut output = Integer(BigInt::from_bytes_be(Sign::Plus, &buff));
        output.modulo_assign(modulus);
        output
    }
    /// Distance to the nearest multiple of q.
    fn modulo_norm(&self, q: &Integer) -> Integer {
        let mut res = self.clone();
        res %= &q;
        let candidate = res.clone();
        res = -res;
        res += &q;
        // Both res and candidate should be in [0, q)
        // smaller one is the size.
        if res.0 < candidate.0 {
            res
        } else {
            candidate
        }
    }
    /// Noise-tolerant Encoding/Decoding used.
    fn encode(&mut self, q: &Integer) {
        *self *= &q;
    }

    fn decode(&mut self, q: &Integer) {
        *self += &self.modulo_norm(q);
        *self /= &q;
    }
}

impl From<i32> for Integer {
    fn from(x: i32) -> Self {
        Integer(BigInt::from(x))
    }
}

impl Neg for Integer {
    type Output = Integer;
    fn neg(self) -> Self::Output {
        Integer(-self.0)
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
            None => (Sign::Plus, s),
            Some(x) => (Sign::Minus, x),
        };
        let body: Vec<u8> = body
            .chars()
            .flat_map(|c| c.to_digit(10).ok_or(Errors::ParseIntError).map(|x| x as u8))
            .collect();
        let res = BigInt::from_radix_be(sign, &body, 10);
        res.ok_or(Errors::ParseIntError).map(|i| Integer(i))
    }
}

pub(crate) trait AdditiveGroup:
    Add + Sub + Neg<Output = Self> + Zero + for<'a> AddAssign<&'a Self> + for<'a> SubAssign<&'a Self>
{
}
impl<T> AdditiveGroup for T where
    T: Add
        + Sub
        + Neg<Output = Self>
        + Zero
        + for<'a> AddAssign<&'a Self>
        + for<'a> SubAssign<&'a Self>
{
}

pub(crate) trait Ring:
    AdditiveGroup + Mul + for<'a> MulAssign<&'a Self> + One + Clone
{
}
impl<T> Ring for T where T: AdditiveGroup + Mul + for<'a> MulAssign<&'a Self> + One + Clone {}

/// Index into coeffs is the degree, i.e. coeffs[0] is constant term.
#[derive(Clone, Eq, PartialEq, Debug)]
pub(crate) struct CycloPoly<R: Ring> {
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
    fn one() -> Self {
        let mut output = Self::zero();
        output.coeffs[0] = R::one();
        output
    }
}

#[macro_export]
macro_rules! impl_assign_ring {
    ( $trait_name:ident, $fun_name:ident) => {
        impl<R: Ring> $trait_name<&CycloPoly<R>> for CycloPoly<R> {
            fn $fun_name(&mut self, rhs: &CycloPoly<R>) {
                for i in 0..DEGREE {
                    <R>::$fun_name(&mut self.coeffs[i], &rhs.coeffs[i]);
                }
            }
        }
    };
}
use impl_assign_ring;

#[macro_export]
macro_rules! impl_consuming_ring {
    ( $trait_name:ident, $assigning_fun_name:ident, $fun_name:ident) => {
        impl<R: Ring> $trait_name<CycloPoly<R>> for CycloPoly<R> {
            type Output = CycloPoly<R>;
            fn $fun_name(self, rhs: CycloPoly<R>) -> Self::Output {
                let mut output = self.clone();
                output.$assigning_fun_name(&rhs);
                output
            }
        }
    };
}
use impl_consuming_ring;

#[macro_export]
macro_rules! impl_op_ring {
    ( $assigning_trait_name:ident, $trait_name:ident, $assigning_fun_name:ident, $fun_name:ident) => {
        impl_assign_ring!($assigning_trait_name, $assigning_fun_name);
        impl_consuming_ring!($trait_name, $assigning_fun_name, $fun_name);
    };
}
use impl_op_ring;

impl_op_ring!(AddAssign, Add, add_assign, add);
impl_op_ring!(SubAssign, Sub, sub_assign, sub);

impl<R: Ring + for<'a> RemAssign<&'a R>> RemAssign<&R> for CycloPoly<R> {
    fn rem_assign(&mut self, rhs: &R) {
        for i in 0..DEGREE {
            self.coeffs[i] %= rhs;
        }
    }
}

impl<R: Ring + for<'a> DivAssign<&'a R>> DivAssign<&R> for CycloPoly<R> {
    fn div_assign(&mut self, rhs: &R) {
        for i in 0..DEGREE {
            self.coeffs[i] /= rhs;
        }
    }
}

impl<R: Ring> Neg for CycloPoly<R> {
    type Output = CycloPoly<R>;
    fn neg(self) -> Self::Output {
        let coeffs = self.coeffs.into_iter().map(|c| c.neg()).collect();
        CycloPoly { coeffs }
    }
}

impl<R: Ring> Mul<CycloPoly<R>> for CycloPoly<R> {
    type Output = CycloPoly<R>;
    fn mul(self, rhs: CycloPoly<R>) -> Self::Output {
        let mut output = self.clone();
        output *= &rhs;
        output
    }
}

/// Multiplication specialized to mod x^d+1.
/// Does (full) polynomial multiplication, then manually reduces things
impl<R: Ring> MulAssign<&CycloPoly<R>> for CycloPoly<R> {
    fn mul_assign(&mut self, rhs: &CycloPoly<R>) {
        let mut unreduced_res = CycloPoly::<R>::zero();
        unreduced_res.coeffs.extend(CycloPoly::<R>::zero().coeffs);
        for (i, val) in unreduced_res.coeffs.iter_mut().enumerate() {
            for j in 0..=i {
                // When indexing should work
                if j < DEGREE && i < DEGREE + j {
                    let mut temp = self.coeffs[j].clone();
                    temp *= &rhs.coeffs[i - j];
                    *val += &temp;
                }
            }
        }
        // reduction step
        for (i, val) in self.coeffs.iter_mut().enumerate() {
            *val = R::zero();
            *val += &unreduced_res.coeffs[i];
            *val -= &unreduced_res.coeffs[DEGREE + i];
        }
    }
}

impl MulAssign<&Integer> for CycloPoly<Integer> {
    fn mul_assign(&mut self, rhs: &Integer) {
        for i in 0..DEGREE {
            self.coeffs[i] *= rhs;
        }
    }
}

impl Mul<Integer> for CycloPoly<Integer> {
    type Output = CycloPoly<Integer>;
    fn mul(self, rhs: Integer) -> Self::Output {
        let mut output = self.clone();
        output *= &rhs;
        output
    }
}

impl CycloPoly<Integer> {
    /// Centered reduction
    pub(crate) fn modulo_assign(&mut self, rhs: &Integer) {
        for i in 0..DEGREE {
            self.coeffs[i].modulo_assign(rhs);
        }
    }
    pub(crate) fn modulo(self, rhs: Integer) -> Self {
        let mut output = self.clone();
        output.modulo_assign(&rhs);
        output
    }
    pub(crate) fn all_ones() -> Self {
        let mut output = Self::zero();
        for i in 0..DEGREE {
            output.coeffs[i] = Integer::one();
        }
        output
    }
    pub(crate) fn sample(dist: &CDT_table) -> Self {
        let mut output = Self::zero();
        for i in 0..DEGREE {
            output.coeffs[i] = dist.sample().into();
        }
        output
    }
    pub(crate) fn uniform_sample<const MOD_EXP_BYTES_PLUS_ONE: usize>(
        rng: &RandGenerator,
        modulus: &Integer,
    ) -> Self {
        let mut output = Self::zero();
        for i in 0..DEGREE {
            output.coeffs[i] = Integer::uniform_sample::<MOD_EXP_BYTES_PLUS_ONE>(rng, modulus);
        }
        output
    }
    /// L_infty form of the modulus norm
    pub(crate) fn max_norm(self, q: &Integer) -> Integer {
        let mut out = Integer::zero();
        for i in 0..DEGREE {
            let size = self.coeffs[i].modulo_norm(&q);
            if size.0 > out.0 {
                out = size;
            }
        }
        out
    }
    /// Noise-tolerant Encoding/Decoding used.
    pub(crate) fn encode(&mut self, q: &Integer) {
        for i in 0..DEGREE {
            self.coeffs[i].encode(&q);
        }
    }

    pub(crate) fn decode(&mut self, q: &Integer) {
        for i in 0..DEGREE {
            self.coeffs[i].decode(&q);
        }
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
    // non-centered reduction
    fn rem_verify(q: i32) -> bool {
        let big_q: Integer = q.into();
        let mut test: Integer = (-1).into();
        let not_expected = test.clone();
        test %= &not_expected;
        if &test == &not_expected {
            return false;
        }
        for i in 0..q {
            let expected: Integer = i.into();
            let mut test = expected.clone();
            test %= &big_q;
            if &expected != &test {
                return false;
            }
        }
        for i in q..(2 * q) {
            let expected: Integer = (i % q).into();
            let mut test = expected.clone();
            test %= &big_q;
            if &expected != &test {
                return false;
            }
        }
        true
    }
    #[test]
    fn test_rem() {
        assert!(rem_verify(11));
        assert!(rem_verify(12));
    }

    #[test]
    fn modulo_two_test() {
        let moduli: Integer = 2_i32.into();
        let reduced = Integer::zero();
        let mut test_reduction = Integer::zero();
        test_reduction.modulo_assign(&moduli);
        assert_eq!(&test_reduction, &reduced);
        let reduced = Integer::one();
        let mut test_reduction = Integer::one();
        test_reduction.modulo_assign(&moduli);
        assert_eq!(&test_reduction, &reduced);
    }
    #[test]
    fn modulo_three_test() {
        let moduli: Integer = 3_i32.into();
        for i in -1..=1 {
            let reduced: Integer = i.into();
            let mut test_reduction = reduced.clone();
            test_reduction.modulo_assign(&moduli);
            assert_eq!(&test_reduction, &reduced);
        }
        for i in 2..=4 {
            let reduced: Integer = (i - 3).into();
            let mut test_reduction: Integer = i.into();
            test_reduction.modulo_assign(&moduli);
            assert_eq!(&test_reduction, &reduced);
        }
    }

    #[test]
    #[ignore]
    fn test_addition() {
        let test_cases = read_from_file("src/arith/add_test_cases.txt");
        for (p1, p2, p3) in test_cases.into_iter() {
            assert_eq!(p1 + p2, p3);
        }
    }
    #[test]
    #[ignore]
    fn test_multiplication() {
        let test_cases = read_from_file("src/arith/mul_test_cases.txt");
        for (p1, p2, p3) in test_cases.into_iter() {
            assert_eq!(p1 * p2, p3);
        }
    }
    #[test]
    fn test_encoding_decoding() {
        let rng = RandGenerator::new();
        let mut q = Integer::one();
        q <<= 10;
        let mut m = Integer::uniform_sample::<{ 2 << 8 + 1 }>(&rng, &q);
        let ref_m = m.clone();
        m.encode(&q);
        for noise_idx in -511..512 {
            dbg!(&noise_idx);
            let noise: Integer = noise_idx.into();
            let mut noised_m = m.clone();
            noised_m += &noise;
            noised_m.decode(&q);
            assert_eq!(&noised_m, &ref_m);
        }
    }
    // Recurring bug
    // have hacky fix.
    #[test]
    fn test_mod_minus_1() {
        let mut val_to_mod = -Integer::one();
        let expected = Integer::one();
        let modulus: Integer = 2_i32.into();
        val_to_mod.modulo_assign(&modulus);
        assert_eq!(&val_to_mod, &expected);
    }
}
