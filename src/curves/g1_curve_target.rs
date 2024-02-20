use plonky2::{field::extension::Extendable, hash::hash_types::RichField};
use subtle::Choice;

use crate::fields::fq_target::FqTarget;

#[derive(Clone, Debug)]
pub struct G1AffineTarget<F: RichField + Extendable<D>, const D: usize> {
    pub x: FqTarget<F, D>,
    pub y: FqTarget<F, D>,
    pub infinity: Choice,
}
