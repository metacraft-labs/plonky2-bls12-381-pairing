use ark_bls12_381::Fq;
use num::BigInt;
use num::{One, Zero};
use plonky2::{
    field::extension::Extendable,
    hash::hash_types::RichField,
    iop::target::BoolTarget,
    plonk::{circuit_builder::CircuitBuilder, circuit_data::CircuitConfig},
};

use crate::fields::bls12_381base::Bls12_381Base;
use crate::fields::fq_target::FqTarget;

#[derive(Clone, Debug)]
pub struct G1AffineTarget<F: RichField + Extendable<D>, const D: usize> {
    pub x: FqTarget<F, D>,
    pub y: FqTarget<F, D>,
    infinity: BoolTarget,
}

impl<F: RichField + Extendable<D>, const D: usize> G1AffineTarget<F, D> {
    pub fn is_identity(&self) -> BoolTarget {
        self.infinity
    }

    /// Returns the identity of the group: the point at infinity.
    pub fn identity() -> Self {
        let config = CircuitConfig::pairing_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        Self {
            x: FqTarget::constant(&mut builder, Fq::zero()),
            y: FqTarget::constant(&mut builder, Fq::one()),
            infinity: BoolTarget::new_unsafe(builder.one()),
        }
    }

    pub fn generator() -> Self {
        let config = CircuitConfig::pairing_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        Self {
            x: FqTarget::fp_constant(
                &mut builder,
                Bls12_381Base([
                    0x5cb3_8790_fd53_0c16,
                    0x7817_fc67_9976_fff5,
                    0x154f_95c7_143b_a1c1,
                    0xf0ae_6acd_f3d0_e747,
                    0xedce_6ecc_21db_f440,
                    0x1201_7741_9e0b_fb75,
                ]),
            ),
            y: FqTarget::fp_constant(
                &mut builder,
                Bls12_381Base([
                    0xbaac_93d5_0ce7_2271,
                    0x8c22_631a_7918_fd8e,
                    0xdd59_5f13_5707_25ce,
                    0x51ac_5829_5040_5194,
                    0x0e1c_8c3f_ad00_59c0,
                    0x0bbc_3efc_5008_a26a,
                ]),
            ),
            infinity: BoolTarget::new_unsafe(builder.zero()),
        }
    }

    pub fn is_point_equal_to(&self, builder: &mut CircuitBuilder<F, D>, rhs: &Self) -> BoolTarget {
        // The only cases in which two points are equal are
        // 1. infinity is set on both
        // 2. infinity is not set on both, and their coordinates are equal

        let inf_set_on_both = builder.and(self.infinity, rhs.infinity);
        let inf_not_set_on_self = builder.not(self.infinity);
        let inf_not_set_on_rhs = builder.not(rhs.infinity);
        let inf_not_set_on_both = builder.and(inf_not_set_on_self, inf_not_set_on_rhs);
        let inf_not_set_on_both = builder.not(inf_not_set_on_both);

        let x_eq_x = self.x.is_equal(builder, &rhs.x);
        let y_eq_y = self.y.is_equal(builder, &rhs.y);

        let x_y_are_eq = builder.and(x_eq_x, y_eq_y);
        let second_pred = builder.and(x_y_are_eq, inf_not_set_on_both);

        builder.or(inf_set_on_both, second_pred)
    }
}

#[cfg(test)]
mod tests {
    use super::G1AffineTarget;
    use plonky2::{
        field::{goldilocks_field::GoldilocksField, types::Field},
        iop::witness::{PartialWitness, WitnessWrite},
        plonk::{
            circuit_builder::CircuitBuilder, circuit_data::CircuitConfig,
            config::PoseidonGoldilocksConfig,
        },
    };
    type F = GoldilocksField;
    type C = PoseidonGoldilocksConfig;
    const D: usize = 2;

    #[test]
    fn test_g1affine_point_equality() {
        let a: G1AffineTarget<F, D> = G1AffineTarget::generator();
        let b: G1AffineTarget<F, D> = G1AffineTarget::identity();

        let config = CircuitConfig::pairing_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let a_point_eq_a_point = a.is_point_equal_to(&mut builder, &a);

        let mut pw = PartialWitness::new();
        pw.set_target(a_point_eq_a_point.target, F::ONE);
        let data = builder.build::<C>();
        dbg!(data.common.degree_bits());
        let _proof = data.prove(pw);
    }
}
