mod arith;
mod utils;

#[derive(Debug)]
pub enum Errors {
    ParseIntError,
    ParsePolyError,
    RNGError,
}
