//! Utility Functions
#![allow(dead_code)]

use crate::Errors;
use ring::rand::{SecureRandom, SystemRandom};
use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
};

#[derive(Clone, Debug)]
pub(crate) struct RandGenerator {
    pub(crate) rng: Arc<Mutex<dyn SecureRandom>>,
}

impl RandGenerator {
    pub(crate) fn fill(&self, dest: &mut [u8]) -> Result<(), Errors> {
        self.rng
            .lock()
            .map_err(|_| Errors::RNGError)?
            .fill(dest)
            .map_err(|_| Errors::RNGError)
    }
    pub(crate) fn new() -> Self {
        let rng = Arc::new(Mutex::new(SystemRandom::new()));
        RandGenerator { rng }
    }
}

/// Sampling from a Discrete Gaussian using CDT sampling/"inversion sampling".
/// Many limitations given the time constraint, most obvious is the (perhaps too low) precision
/// the probabilities are computed to.
/// It is plausibly fine, see https://ieeexplore.ieee.org/document/8295226 , but I have not
/// verified.

#[allow(non_camel_case_types)]
#[derive(Debug)]
pub(crate) struct CDT_table {
    pub(crate) rng: RandGenerator,
    probs: HashMap<i32, f64>,
    pub(crate) bound: i32,
}

impl CDT_table {
    pub(crate) fn new(rng: &RandGenerator, sigma: f64) -> Self {
        // Tail cut
        let bound = (10. * sigma.ceil()) as i32;
        let mut probs: HashMap<i32, f64> = HashMap::new();
        for idx in -bound..=bound {
            let exponent = -std::f64::consts::PI * ((idx as f64 / sigma).powi(2));
            probs.insert(idx, f64::exp(exponent));
        }
        let total_prob_mass: f64 = probs.values().sum();
        let probs = probs
            .iter_mut()
            .map(|(k, v)| (*k, *v / total_prob_mass))
            .collect();
        let rng = rng.clone();
        CDT_table { rng, probs, bound }
    }
    pub(crate) fn sample(&self) -> i32 {
        // Typical way of generating a uniformly random (in [0,1]) double.
        // 53 as f64's have 53 bit exponents
        let mut buff = [0 as u8; 64 / 8];
        self.rng
            .fill(&mut buff)
            .expect("RNG Error: Will not try to recover from");
        let sampled = u64::from_be_bytes(buff) % (1 << 53);
        let sampled: f64 = sampled as f64 / ((1_u64 << 53) as f64);
        let mut cumulative = 0.;
        for i in -self.bound..self.bound {
            let upper_bound = cumulative + self.probs.get(&i).unwrap();
            if cumulative <= sampled && sampled < upper_bound {
                return i;
            } else {
                cumulative = upper_bound;
            }
        }
        self.bound
    }
}

#[cfg(test)]

mod tests {

    // Testing could be greatly improved using the GLITCH test suite
    // https://eprint.iacr.org/2017/438
    // Additionally, many other sampling techniques are better for large sigma, for example
    // convolution sampling.
    use super::*;
    #[test]
    fn test_table_gen() {
        let rng = RandGenerator::new();
        let table = CDT_table::new(&rng, 10.);
        let total_prob_mass: f64 = table.probs.values().sum();
        assert!((1. - total_prob_mass).abs() < 0.00001);
    }

    #[test]
    fn test_sampling() {
        const N: u32 = 1000;
        let rng = RandGenerator::new();
        let table = CDT_table::new(&rng, 10.);
        let mut samples = Vec::new();
        for _ in 0..N {
            samples.push(table.sample());
        }
        // Just checking if Pr[-1] ~ Pr[1], and the mean is zero.
        // Very basic checks
        let mut neg_counts = 0;
        let mut pos_counts = 0;
        let mut total = 0;
        for sample in samples.iter() {
            total += sample;
            if *sample < 0 {
                neg_counts += 1;
            } else if *sample > 0 {
                pos_counts += 1;
            }
        }
        // Test really pulled out of a hat
        let diff: i32 = pos_counts - neg_counts;
        assert!(diff.abs() < 100);
        let mean = (total as f64) / (N as f64);
        assert!(mean.abs() < 1.);
    }
}
