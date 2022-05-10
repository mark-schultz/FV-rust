//! Encryption Scheme

use std::ops::AddAssign;

use crate::arith::{CycloPoly, Integer};
use crate::utils::{CDT_table, RandGenerator};
use crate::{DELTA_R, Q_EXP, Q_EXP_BYTES_PLUS_1, SIGMA, T_EXP, T_EXP_BYTES_PLUS_1};
use num_traits::identities::{Zero};

struct LPR {
    gauss_dist: CDT_table,
    sk: SecretKey,
    pk: PublicKey,
    ctxt_mod: CtxtModuli,
    ptxt_mod: PtxtModuli,
}

struct SecretKey(CycloPoly<Integer>);
struct PublicKey(CycloPoly<Integer>, CycloPoly<Integer>);
struct CipherText(CycloPoly<Integer>, CycloPoly<Integer>);
struct PtxtModuli(Integer);
struct CtxtModuli(Integer);

impl AddAssign<&CipherText> for CipherText {
    fn add_assign(&mut self, rhs: &CipherText) {
        self.0 += &rhs.0;
        self.1 += &rhs.1;
    }
}

impl LPR {
    fn sk_gen(dist: &CDT_table) -> SecretKey {
        SecretKey(CycloPoly::<Integer>::sample(dist))
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
        LPR {
            gauss_dist: dist,
            sk,
            pk,
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
        dbg!(&p0);

        let mut p1 = pk.1.clone();
        p1 *= &u;
        let e2 = CycloPoly::<Integer>::sample(&dist);
        p1 += &e2;
        p1.modulo_assign(&q.0);
        dbg!(&p1);

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

#[cfg(test)]
mod tests {
    use super::*;

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
        res.modulo_assign(dbg!(&keys.ctxt_mod.0));
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
}
