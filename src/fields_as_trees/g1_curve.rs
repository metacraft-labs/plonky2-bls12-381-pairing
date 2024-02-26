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

    pub fn is_point_equal_to(&self, builder: &mut CircuitBuilder<F, D>, rhs: &Self) -> bool {
        // The only cases in which two points are equal are
        // 1. infinity is set on both
        // 2. infinity is not set on both, and their coordinates are equal

        let one = builder.constant(F::from_canonical_u8(1));
        let is_self_inf_true = self.infinity.target.eq(&one);
        let is_rhs_inf_true = rhs.infinity.target.eq(&one);

        let self_x_eq_rhs_x = (self.x.is_equal(builder, &rhs.x)).target.eq(&one);
        let self_y_eq_rhs_y = (self.y.is_equal(builder, &rhs.y)).target.eq(&one);

        let onee = FqTarget::constant(builder, Fq::one());
        let oneee = FqTarget::constant(builder, Fq::zero());

        let test_purposes = onee.is_equal(builder, &oneee).target.eq(&one);

        println!("test_purposes: {:?}", test_purposes);
        println!("is_self_inf_true: {:?}", !is_self_inf_true);
        println!("is_rhs_inf_true: {:?}", !is_rhs_inf_true);
        println!("self_x_eq_rhs_x: {:?}", self_x_eq_rhs_x);
        println!("self_y_eq_rhs_y: {:?}", self_y_eq_rhs_y);
        (is_self_inf_true & is_rhs_inf_true) | ((!is_self_inf_true) & (!is_rhs_inf_true) & self_x_eq_rhs_x & self_y_eq_rhs_y)
    }
}

#[cfg(test)]
mod tests {
    use plonky2::{
        field::goldilocks_field::GoldilocksField,
        iop::witness::PartialWitness,
        plonk::{
            circuit_builder::CircuitBuilder, circuit_data::CircuitConfig,
            config::PoseidonGoldilocksConfig,
        },
    };
    use super::G1AffineTarget;
    type F = GoldilocksField;
    type C = PoseidonGoldilocksConfig;
    const D: usize = 2;


    #[test]
    fn test_affine_point_equality() {
        let a: G1AffineTarget<F, D> = G1AffineTarget::generator();
        let b: G1AffineTarget<F, D> = G1AffineTarget::identity();

        let config = CircuitConfig::pairing_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        assert!(a.is_point_equal_to(&mut builder, &a));

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        dbg!(data.common.degree_bits());
        let _proof = data.prove(pw);
    }
}
