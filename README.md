# FV-rust
A Rust implementation of the [Fan-Vercauteren (FV) FHE Scheme](https://eprint.iacr.org/2012/144.pdf).
The goal of this repository is to get *as much* of FV implemented as possible,
given a fairly severe timing constraint of

* ~16 hours total available, where
* a few of the hours will be used to prepare a presentation.

This is to say that the focus will be on the speed of implementing things.
Broadly, I will proceed in the implementation in the following order

1. Cyclotomic Polynomial Arithmetic
2. Encryption/Decryption
3. Homomorphic Addition
4. Homomorphic Multiplication
5. Bootstrapping?
