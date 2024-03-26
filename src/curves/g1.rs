use plonky2::{field::extension::Extendable, hash::hash_types::RichField};

use crate::fields::fq_target::FqTarget;

pub struct G1AffineTarget<F: RichField + Extendable<D>, const D: usize> {
    pub x: FqTarget<F, D>,
    pub y: FqTarget<F, D>,
    pub infinity: bool,
}

impl<F: RichField + Extendable<D>, const D: usize> G1AffineTarget<F, D> {
    pub fn xy(&self) -> Option<(&self::FqTarget<F, D>, &self::FqTarget<F, D>)> {
        (!self.infinity).then(|| (&self.x, &self.y))
    }
}

pub struct G1Prepared<F: RichField + Extendable<D>, const D: usize>(pub G1AffineTarget<F, D>);
