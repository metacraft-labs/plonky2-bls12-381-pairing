use plonky2::{field::extension::Extendable, hash::hash_types::RichField};
use subtle::Choice;

use crate::fields::fq2_target::Fq2Target;

#[derive(Clone, Debug)]
pub struct G2ProjectiveTarget<F: RichField + Extendable<D>, const D: usize> {
    pub x: Fq2Target<F, D>,
    pub y: Fq2Target<F, D>,
    pub z: Fq2Target<F, D>,
}

#[derive(Clone, Debug)]
pub struct G2AffineTarget<F: RichField + Extendable<D>, const D: usize> {
    pub x: Fq2Target<F, D>,
    pub y: Fq2Target<F, D>,
    pub infinity: Choice, // not sure if it's going to work in circuit
}
