//! Encryption Scheme

use crate::arith::{CycloPoly, Integer};
use crate::utils::{CDT_table, RandGenerator};
use crate::{N, SIGMA, T};

struct LPR {
    gauss_dist: CDT_table,
    sk: SecretKey,
    pk: PublicKey,
    ctxt_mod: Integer,
    ptxt_mod: Integer,
}

struct SecretKey(CycloPoly<Integer>);
struct PublicKey(CycloPoly<Integer>, CycloPoly<Integer>);
struct CipherText(CycloPoly<Integer>, CycloPoly<Integer>);

impl LPR {
    fn sk_gen(dist: &CDT_table) -> SecretKey {
        SecretKey(CycloPoly::<Integer>::sample(dist))
    }
    fn pk_gen(dist: &CDT_table, sk: &SecretKey, modulus: &Integer) -> PublicKey {
        let mut a = CycloPoly::uniform_sample(&dist.rng, modulus);
        let e = CycloPoly::sample(dist);
        a *= &sk.0;
        a += &e;
        a = -a;
        let mut b = a.clone();
        b.modulo_assign(modulus);
        PublicKey(b, a)
    }
    fn gen_moduli() -> (Integer, Integer) {
        let mut ptxt_mod: Integer = 1_i32.into();
        ptxt_mod <<= T;
        let mut ctxt_mod: Integer = 1_i32.into();
        ctxt_mod <<= N;
        (ptxt_mod, ctxt_mod)
    }
    fn new() -> Self {
        let rng = RandGenerator::new();
        let dist = CDT_table::new(&rng, SIGMA);
        let sk = LPR::sk_gen(&dist);
        let (ptxt_mod, ctxt_mod) = LPR::gen_moduli();
        let pk = LPR::pk_gen(&dist, &sk, &ctxt_mod);
        LPR {
            gauss_dist: dist,
            sk,
            pk,
            ctxt_mod,
            ptxt_mod,
        }
    }
    fn enc(pk: &PublicKey, dist: &CDT_table, m: &CycloPoly<Integer>) -> CipherText {
        let mut p0 = pk.0.clone();
        let mut p1 = pk.1.clone();
        let (t, q) = LPR::gen_moduli();
        let mut delta = q.clone();
        delta /= &t;
        let u = CycloPoly::<Integer>::sample(&dist);
        p0 *= &u;
        p1 *= &u;
        let e1 = CycloPoly::<Integer>::sample(&dist);
        let e2 = CycloPoly::<Integer>::sample(&dist);
        p1 += &e2;
        p0 += &e1;
        let mut m = m.clone();
        m *= &delta;
        p0 += &m;
        p0.modulo_assign(&q);
        p1.modulo_assign(&q);
        CipherText(p0, p1)
    }
    fn dec(&self, ctxt: &CipherText) -> CycloPoly<Integer> {
        let (mut c0, mut c1) = (ctxt.0.clone(), ctxt.1.clone());
        c1 *= &self.sk.0;
        c0 += &c1;
        let q = &self.ctxt_mod;
        let t = &self.ptxt_mod;
        c0.modulo_assign(q);
        c0 *= t;
        // Simulating centered rounding
        // Normally, add 1/2
        // Here, add (q/2)/q.
        let mut half_q: Integer = 1_i32.into();
        half_q <<= N - 1;
        let mut to_add = CycloPoly::<Integer>::all_ones();
        to_add *= &half_q;
        c0 += &to_add;
        c0 /= q;
        c0.modulo_assign(t);
        c0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_modulo_relation() {
        let (ptxt, ctxt) = LPR::gen_moduli();
        let mut delta = ctxt.clone();
        delta /= &ptxt;
        let mut expected_delta: Integer = 1_i32.into();
        expected_delta <<= N - 1;
        assert_eq!(&expected_delta, &delta);
        let mut r_t = ctxt.clone();
        r_t %= &ptxt;
        delta *= &ptxt;
        r_t += &delta;
        assert_eq!(&ctxt, &r_t);
    }

    #[test]
    fn test_enc() {
        let keys = LPR::new();
        let mut rand_message =
            CycloPoly::<Integer>::uniform_sample(&keys.gauss_dist.rng, &keys.ptxt_mod);
        rand_message.modulo_assign(&keys.ptxt_mod);
        let ctxt = LPR::enc(&keys.pk, &keys.gauss_dist, &rand_message);
        let dec_m = keys.dec(&ctxt);
        assert_eq!(dec_m, rand_message);
    }
}
