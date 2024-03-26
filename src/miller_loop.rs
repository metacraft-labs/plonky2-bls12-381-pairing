use plonky2::{
    field::extension::Extendable, hash::hash_types::RichField,
    plonk::circuit_builder::CircuitBuilder,
};

use crate::{
    curves::{g1::G1AffineTarget, g2::EllCoeffTarget},
    fields::fq12_target::Fq12Target,
};

fn ell<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    f: &mut Fq12Target<F, D>,
    g2_coeffs: EllCoeffTarget<F, D>,
    p: G1AffineTarget<F, D>,
) {
    let c0 = g2_coeffs.0;
    let mut c1 = g2_coeffs.1;
    let mut c2 = g2_coeffs.2;
    let (px, py) = p.xy().unwrap();

    let c2 = c2.mul_assign_by_fp(builder, py.clone());
    let c1 = c1.mul_assign_by_fp(builder, px.clone());
    f.mul_by_014(builder, &c0, &c1, &c2);
}
