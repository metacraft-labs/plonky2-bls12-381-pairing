use ark_bls12_381::{Fq, G1Affine};
use ark_ec::AffineRepr;
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
    pub fn identity(builder: &mut CircuitBuilder<F, D>) -> Self {
        Self {
            x: FqTarget::constant(builder, Fq::zero()),
            y: FqTarget::constant(builder, Fq::one()),
            infinity: builder._true(),
        }
    }

    pub fn experimental_generator(builder: &mut CircuitBuilder<F, D>) -> Self {
        Self {
            x: FqTarget::constant(builder, *G1Affine::generator().x().unwrap()),
            y: FqTarget::constant(builder, *G1Affine::generator().y().unwrap()),
            infinity: builder._false(),
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

        let x_eq_x = self.x.is_equal(builder, &rhs.x);
        let y_eq_y = self.y.is_equal(builder, &rhs.y);

        let x_y_are_eq = builder.and(x_eq_x, y_eq_y);
        let second_pred = builder.and(x_y_are_eq, inf_not_set_on_both);

        builder.or(inf_set_on_both, second_pred)
    }

    pub fn conditional_select(
        builder: &mut CircuitBuilder<F, D>,
        a: Self,
        b: Self,
        flag: BoolTarget,
    ) -> Self {
        Self {
            x: FqTarget::select(builder, &a.x, &b.x, &flag),
            y: FqTarget::select(builder, &a.y, &b.y, &flag),
            infinity: builder.or(a.infinity, b.infinity),
        }
    }

    pub fn neg(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
        let one = FqTarget::constant(builder, Fq::one());
        let y_neg = self.y.neg(builder);
        Self {
            x: self.x.clone(),
            y: FqTarget::select(builder, &y_neg, &one, &self.infinity),
            infinity: self.infinity,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::fields::fq_target::FqTarget;

    use super::G1AffineTarget;
    use ark_bls12_381::{Fq, G1Affine};
    use ark_ec::AffineRepr;
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
        let config = CircuitConfig::pairing_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let g1_affine_generator = G1AffineTarget::experimental_generator(&mut builder);

        let const_one = FqTarget::constant(&mut builder, Fq::from(1));
        let g1_affine_x = g1_affine_generator.x;
        let g1_affine_x_1 = g1_affine_x.add(&mut builder, &const_one);
        let g1_affine_y = g1_affine_generator.y;
        let g1_affine_y_1 = g1_affine_y.add(&mut builder, &const_one);

        let a: G1AffineTarget<F, D> = G1AffineTarget {
            x: g1_affine_x.clone(),
            y: g1_affine_y.clone(),
            infinity: builder._false(),
        };

        let b: G1AffineTarget<F, D> = G1AffineTarget {
            x: g1_affine_x,
            y: g1_affine_y.clone(),
            infinity: builder._false(),
        };
        let case_1 = a.is_point_equal_to(&mut builder, &b);

        let b: G1AffineTarget<F, D> = G1AffineTarget {
            x: g1_affine_x_1.clone(),
            y: g1_affine_y_1.clone(),
            infinity: builder._false(),
        };
        let case_2 = a.is_point_equal_to(&mut builder, &b);

        let b: G1AffineTarget<F, D> = G1AffineTarget {
            x: g1_affine_x_1.clone(),
            y: g1_affine_y_1,
            infinity: builder._true(),
        };
        let case_3 = a.is_point_equal_to(&mut builder, &b);

        let a: G1AffineTarget<F, D> = G1AffineTarget {
            x: g1_affine_x_1,
            y: g1_affine_y,
            infinity: builder._true(),
        };
        let case_4 = a.is_point_equal_to(&mut builder, &b);

        let mut pw = PartialWitness::new();
        pw.set_target(case_1.target, F::ONE);
        pw.set_target(case_2.target, F::ZERO);
        pw.set_target(case_3.target, F::ZERO);
        pw.set_target(case_4.target, F::ONE);
        let data = builder.build::<C>();
        dbg!(data.common.degree_bits());
        let _proof = data.prove(pw);
    }
}
