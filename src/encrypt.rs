//! Encryption Scheme

use std::ops::{Add, AddAssign, Mul, Sub, SubAssign};

use crate::arith::{CycloPoly, Integer};
use crate::utils::{CDT_table, RandGenerator};
use crate::{Q_EXP, Q_EXP_BYTES_PLUS_1, RELIN_BASE, RELIN_TERMS, SIGMA, T_EXP};
use num_traits::One;

struct LPR {
    gauss_dist: CDT_table,
    sk: SecretKey,
    pk: PublicKey,
    rlk: RelinearizationKey,
    ctxt_mod: CtxtModuli,
    ptxt_mod: PtxtModuli,
}

struct SecretKey(CycloPoly<Integer>);
struct PublicKey(CycloPoly<Integer>, CycloPoly<Integer>);
struct RelinearizationKey(Vec<(CycloPoly<Integer>, CycloPoly<Integer>)>);
#[derive(Clone)]
struct CipherText(CycloPoly<Integer>, CycloPoly<Integer>);
struct PtxtModuli(Integer);
struct CtxtModuli(Integer);

impl LPR {
    fn sk_gen(dist: &CDT_table) -> SecretKey {
        // Binary secrets optimization
        SecretKey(CycloPoly::<Integer>::binary_sample(&dist.rng))
    }
    fn pk_gen(dist: &CDT_table, sk: &SecretKey, modulus: &CtxtModuli) -> PublicKey {
        let a = CycloPoly::uniform_sample::<Q_EXP_BYTES_PLUS_1>(&dist.rng, &modulus.0);
        let mut b = a.clone();
        b *= &sk.0;
        let e = CycloPoly::sample(&dist);
        b += &e;
        b = -b;
        b.modulo_assign(&modulus.0);
        PublicKey(b, a)
    }
    fn rlk_gen(dist: &CDT_table, sk: &SecretKey, modulus: &CtxtModuli) -> RelinearizationKey {
        let mut out: Vec<(CycloPoly<Integer>, CycloPoly<Integer>)> = Vec::new();
        let base: Integer = RELIN_BASE.into();
        let mut multiplier = Integer::one();
        for i in 0..=RELIN_TERMS {
            let a = CycloPoly::uniform_sample::<Q_EXP_BYTES_PLUS_1>(&dist.rng, &modulus.0);
            let mut b = a.clone();
            b *= &sk.0;
            let e = CycloPoly::sample(&dist);
            b += &e;
            b = -b;
            let mut s_squared = sk.0.clone();
            s_squared *= &sk.0;
            s_squared *= &multiplier;
            multiplier *= &base; // T^i
            b += &s_squared;
            b.modulo_assign(&modulus.0);
            out.push((b, a))
        }
        RelinearizationKey(out)
    }

    fn gen_moduli() -> (PtxtModuli, CtxtModuli) {
        let mut ptxt_mod: Integer = 1_i32.into();
        ptxt_mod <<= T_EXP;
        let mut ctxt_mod: Integer = 1_i32.into();
        ctxt_mod <<= Q_EXP;
        (PtxtModuli(ptxt_mod), CtxtModuli(ctxt_mod))
    }
    fn gen_scaling_factor() -> Integer {
        let (t, mut q) = LPR::gen_moduli();
        q.0 /= &t.0;
        q.0
    }

    fn new() -> Self {
        let rng = RandGenerator::new();
        let dist = CDT_table::new(&rng, SIGMA);
        let (t, q) = LPR::gen_moduli();
        let sk = LPR::sk_gen(&dist);
        let pk = LPR::pk_gen(&dist, &sk, &q);
        let rlk = LPR::rlk_gen(&dist, &sk, &q);
        LPR {
            gauss_dist: dist,
            sk,
            pk,
            rlk,
            ctxt_mod: q,
            ptxt_mod: t,
        }
    }
    fn enc(pk: &PublicKey, dist: &CDT_table, m: &CycloPoly<Integer>) -> CipherText {
        let mut p0 = pk.0.clone();
        let (_, q) = LPR::gen_moduli();
        let u = CycloPoly::<Integer>::sample(&dist);
        p0 *= &u;
        let e1 = CycloPoly::<Integer>::sample(&dist);
        p0 += &e1;
        let delta = LPR::gen_scaling_factor();
        let mut m = m.clone();
        m.encode(&delta);
        p0 += &m;
        p0.modulo_assign(&q.0);

        let mut p1 = pk.1.clone();
        p1 *= &u;
        let e2 = CycloPoly::<Integer>::sample(&dist);
        p1 += &e2;
        p1.modulo_assign(&q.0);

        CipherText(p0, p1)
    }
    fn dec(&self, ctxt: &CipherText) -> CycloPoly<Integer> {
        let t = &self.ptxt_mod;
        let q = &self.ctxt_mod;
        let mut c = ctxt.1.clone();
        c *= &self.sk.0;
        c += &ctxt.0;
        c.modulo_assign(&q.0);
        let delta = LPR::gen_scaling_factor();
        c.decode(&delta);
        c.modulo_assign(&t.0);
        c
    }
}

// Can't implement `zero` as we can't implement `is_zero`.

impl AddAssign<&CipherText> for CipherText {
    fn add_assign(&mut self, rhs: &CipherText) {
        self.0 += &rhs.0;
        self.1 += &rhs.1;
    }
}

impl Add<CipherText> for CipherText {
    type Output = CipherText;
    fn add(self, rhs: Self) -> Self::Output {
        let mut out = self.clone();
        out += &rhs;
        out
    }
}

impl SubAssign<&CipherText> for CipherText {
    fn sub_assign(&mut self, rhs: &CipherText) {
        self.0 -= &rhs.0;
        self.1 -= &rhs.1;
    }
}

impl Sub<CipherText> for CipherText {
    type Output = CipherText;
    fn sub(self, rhs: Self) -> Self::Output {
        let mut out = self.clone();
        out -= &rhs;
        out
    }
}

struct UnlinearizedCiphertext(CycloPoly<Integer>, CycloPoly<Integer>, CycloPoly<Integer>);

impl Mul<CipherText> for CipherText {
    type Output = UnlinearizedCiphertext;
    fn mul(self, rhs: CipherText) -> Self::Output {
        let mut c0 = self.0.clone();
        c0 *= &rhs.0;

        let mut c1 = self.0.clone();
        c1 *= &rhs.1;
        let mut tmp = self.1.clone();
        tmp *= &rhs.0;
        c1 += &tmp;

        let mut c2 = self.1.clone();
        c2 *= &rhs.1.clone();

        let (t, q) = LPR::gen_moduli();
        c0 *= &t.0;
        c1 *= &t.0;
        c2 *= &t.0;

        c0 /= &q.0;
        c1 /= &q.0;
        c2 /= &q.0;

        UnlinearizedCiphertext(c0, c1, c2)
    }
}

impl UnlinearizedCiphertext {
    fn linearize(self, rlk: &RelinearizationKey) -> CipherText {
        let (_, q) = LPR::gen_moduli();
        let mut c0_prime = self.0.clone();
        let mut c1_prime = self.1.clone();
        let c2_decomp = self.2.clone().base_decompose();
        for i in 0..=RELIN_TERMS {
            let mut tmp = rlk.0[i].0.clone();
            tmp *= &c2_decomp[i];
            c0_prime += &tmp;
        }
        c0_prime.modulo_assign(&q.0);
        for i in 0..=RELIN_TERMS {
            let mut tmp = rlk.0[i].1.clone();
            tmp *= &c2_decomp[i];
            c1_prime += &tmp;
        }
        c1_prime.modulo_assign(&q.0);
        CipherText(c0_prime, c1_prime)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{DELTA_R, T_EXP_BYTES_PLUS_1};
    use num_traits::Zero;

    #[test]
    fn test_modulo_relation() {
        let (t, q) = LPR::gen_moduli();
        let mut delta = LPR::gen_scaling_factor();
        let mut expected_delta = q.0.clone();
        expected_delta >>= 1;
        assert_eq!(&expected_delta, &delta);
        let mut r_t = q.0.clone();
        r_t %= &t.0;
        delta *= &t.0;
        r_t += &delta;
        assert_eq!(&q.0, &r_t);
    }
    #[test]
    fn test_pub_key_relation() {
        let keys = LPR::new();
        let mut res: CycloPoly<Integer> = keys.sk.0.clone();
        res *= &keys.pk.1;
        res += &keys.pk.0;
        res.modulo_assign(&keys.ctxt_mod.0);
        // Should be B-bounded now.
        let norm = res.max_norm(&keys.ctxt_mod.0);
        assert!(norm.0 < keys.gauss_dist.bound.into());
    }

    #[test]
    fn test_noise_size() {
        let keys = LPR::new();
        let mut rand_message = CycloPoly::<Integer>::uniform_sample::<T_EXP_BYTES_PLUS_1>(
            &keys.gauss_dist.rng,
            &keys.ptxt_mod.0,
        );
        rand_message.modulo_assign(&keys.ptxt_mod.0);
        let ctxt = LPR::enc(&keys.pk, &keys.gauss_dist, &rand_message);
        let (mut c0, mut c1) = (ctxt.0, ctxt.1);
        c1 *= &keys.sk.0;
        c0 += &c1;
        c0.modulo_assign(&keys.ctxt_mod.0);
        let mut delta = (&keys.ctxt_mod.0).clone();
        delta /= &keys.ptxt_mod.0;
        rand_message *= &delta;
        c0 -= &rand_message;
        let norm = c0.max_norm(&keys.ctxt_mod.0);
        let dist_bound: usize = (&keys.gauss_dist.bound)
            .abs()
            .try_into()
            .expect("dist_bound is on the order of 10sigma << 2^32");
        let bound = 2 * DELTA_R * dist_bound * dist_bound + dist_bound;
        assert!(norm.0 < bound.into());
    }
    #[test]
    fn test_enc_zero() {
        let keys = LPR::new();
        let m = CycloPoly::<Integer>::zero();
        let ctxt = LPR::enc(&keys.pk, &keys.gauss_dist, &m);
        let dec_m = keys.dec(&ctxt);
        assert_eq!(&dec_m, &m);
    }

    #[test]
    fn test_enc() {
        let keys = LPR::new();
        let mut rand_message = CycloPoly::<Integer>::uniform_sample::<T_EXP_BYTES_PLUS_1>(
            &keys.gauss_dist.rng,
            &keys.ptxt_mod.0,
        );
        rand_message.modulo_assign(&keys.ptxt_mod.0);
        let ctxt = LPR::enc(&keys.pk, &keys.gauss_dist, &rand_message);
        let dec_m = keys.dec(&ctxt);
        assert_eq!(dec_m, rand_message);
    }
    #[test]
    fn test_add() {
        let keys = LPR::new();
        let m1 = CycloPoly::<Integer>::uniform_sample::<T_EXP_BYTES_PLUS_1>(
            &keys.gauss_dist.rng,
            &keys.ptxt_mod.0,
        );
        let m2 = CycloPoly::<Integer>::uniform_sample::<T_EXP_BYTES_PLUS_1>(
            &keys.gauss_dist.rng,
            &keys.ptxt_mod.0,
        );
        let mut m3 = m1.clone();
        m3 += &m2;
        let c1 = LPR::enc(&keys.pk, &keys.gauss_dist, &m1);
        let c2 = LPR::enc(&keys.pk, &keys.gauss_dist, &m2);
        let mut c3 = c2;
        c3 += &c1;
        let dec_m3 = keys.dec(&c3);

        // Hom.ops only valid mod t.
        m3.modulo_assign(&keys.ptxt_mod.0);
        assert_eq!(&dec_m3, &m3);
    }
    #[test]
    fn test_mul() {
        let keys = LPR::new();
        let m1 = CycloPoly::<Integer>::uniform_sample::<T_EXP_BYTES_PLUS_1>(
            &keys.gauss_dist.rng,
            &keys.ptxt_mod.0,
        );
        let m2 = CycloPoly::<Integer>::uniform_sample::<T_EXP_BYTES_PLUS_1>(
            &keys.gauss_dist.rng,
            &keys.ptxt_mod.0,
        );
        let mut m3 = m1.clone();
        m3 *= &m2;
        let c1 = LPR::enc(&keys.pk, &keys.gauss_dist, &m1);
        let c2 = LPR::enc(&keys.pk, &keys.gauss_dist, &m2);
        let c3 = c1 * c2;
        let c3 = c3.linearize(&keys.rlk);
        let dec_m3 = keys.dec(&c3);

        // Hom.ops only valid mod t.
        m3.modulo_assign(&keys.ptxt_mod.0);
        assert_eq!(&dec_m3, &m3);
    }
}
