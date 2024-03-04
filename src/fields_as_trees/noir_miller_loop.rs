use ark_bls12_381::Fq;
use plonky2::{
    field::extension::Extendable, hash::hash_types::RichField,
    plonk::circuit_builder::CircuitBuilder,
};
use ark_std::One;

use crate::{fields::{bls12_381base::Bls12_381Base, fq_target::FqTarget}, utils::constants::NAF_DIGIT};

use super::{fq12_target_tree::Fq12Target, fq2_target_tree::Fq2Target, fq6_target_tree::Fq6Target};

#[derive(Clone, Debug)]
struct LineEval<F: RichField + Extendable<D>, const D: usize> {
    pub l: Fq12Target<F, D>,
    pub T: G2PointTarget<F, D>,
}

#[derive(Clone, Debug)]
pub struct G2PointTarget<F: RichField + Extendable<D>, const D: usize> {
    pub x: Fq2Target<F, D>,
    pub y: Fq2Target<F, D>,
    pub z: Fq2Target<F, D>,
}

impl<F: RichField + Extendable<D>, const D: usize> G2PointTarget<F, D> {
    pub fn negate(self, builder: &mut CircuitBuilder<F, D>) -> Self {
        Self {
            x: self.x,
            y: self.y.neg(builder),
            z: self.z,
        }
    }

    pub fn from_affine(
        builder: &mut CircuitBuilder<F, D>,
        x: Fq2Target<F, D>,
        y: Fq2Target<F, D>,
    ) -> Self {
        Self {
            x,
            y,
            z: Fq2Target::one(builder),
        }
    }

    pub fn generator(builder: &mut CircuitBuilder<F, D>) -> Self {
        let x = Fq2Target {
            c0: FqTarget::fp_constant(
                builder,
                Bls12_381Base([
                    0xb8bd_21c1_c856_80d4,
                    0xefbb_05a8_2603_ac0b,
                    0x77d1_e37a_640b_51b4,
                    0x023b_40fa_d47a_e4c6,
                    0x5110_c52d_2705_0826,
                    0x910a_8ff0_b2a2_4a02
                ]),
            ),
            c1: FqTarget::fp_constant(
                builder,
                Bls12_381Base([
                    0x7e2b_045d_057d_ace5,
                    0x575d_9413_12f1_4c33,
                    0x4950_7fdc_bb61_dab5,
                    0x1ab6_2099_d0d0_6b59,
                    0x654f_2788_a0d3_ac7d,
                    0x609f_7152_602b_e013
                ]),
            ),
        };

        let y = Fq2Target {
            c0: FqTarget::fp_constant(
                builder,
                Bls12_381Base([
                    0x0128_b808_8654_93e1,
                    0x89a2_ac3b_ccc9_3a92,
                    0x2cd1_6051_699a_426d,
                    0xa7d3_bd8c_aa9b_fdad,
                    0x1a35_2eda_c6cd_c98c,
                    0x116e_7d72_27d5_e50c,
                ]),
            ),
            c1: FqTarget::fp_constant(
                builder,
                Bls12_381Base([
                    0xbe79_5ff0_5f07_a9aa,
                    0xa11d_ec5c_270d_373f,
                    0xab99_2e57_ab92_7426,
                    0xaf63_a785_7e28_3ecb,
                    0x998b_c22b_b0d2_ac32,
                    0xcc34_a72e_a0c4_0606
                ]),
            ),
        };

        G2PointTarget::from_affine(
            builder,
            x,
            y,
        )
    }    
}

#[derive(Clone, Debug)]
pub struct G1PointTarget<F: RichField + Extendable<D>, const D: usize> {
    pub x: FqTarget<F, D>,
    pub y: FqTarget<F, D>,
    pub z: FqTarget<F, D>,
}

impl<F: RichField + Extendable<D>, const D: usize> G1PointTarget<F, D> {
    pub fn from_affine(
        builder: &mut CircuitBuilder<F, D>,
        x: FqTarget<F, D>,
        y: FqTarget<F, D>,
    ) -> Self {
        let one = FqTarget::constant(builder, Fq::one());
        Self {
            x,
            y,
            z: one,
        }
    }

    pub fn generator(builder: &mut CircuitBuilder<F, D>) -> Self {
        let x = FqTarget::fp_constant(
            builder,
            Bls12_381Base([
                0xbbc6_22db_0af0_3afb,
                0xef1a_7af9_3fe8_556c,
                0x58ac_1b17_3f3a_4ea1,
                0x05b9_7497_4f8c_68c3,
                0x0fac_a94f_8c63_9526,
                0x94d7_9731_a7d3_f117
            ]),
        );

        let y = FqTarget::fp_constant(
            builder,
            Bls12_381Base([
                0xe1e7_c546_2923_aa0c,
                0xe48a_88a2_44c7_3cd0,
                0xedb3_042c_cb18_db00,
                0xf60a_d0d5_95e0_f5fc,
                0xe48a_1d74_ed30_9ea0,
                0xf1a0_aae3_81f4_b308
            ]),
        );

        G1PointTarget::from_affine(
            builder,
            x,
            y,
        )
    }
}

// Algorithm 26 from: https://eprint.iacr.org/2010/354.pdf
fn line_function_double_point<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    q: G2PointTarget<F, D>,
    p: G1PointTarget<F, D>,
) -> LineEval<F, D> {
    let tmp0 = q.x.square(builder);
    let tmp1 = q.y.square(builder);
    let tmp2 = tmp1.square(builder);
    let x = tmp1.add(builder, &q.x);
    let x = x.square(builder);
    let x = x.sub(builder, &tmp0);
    let tmp3 = x.sub(builder, &tmp2);
    let tmp3 = tmp3.double(builder);
    let tmp0_doubled = tmp0.double(builder);
    let tmp4 = tmp0.add(builder, &tmp0_doubled);
    let tmp6 = q.x.add(builder, &tmp4);
    let tmp5 = tmp4.square(builder);
    let tmp3_doubled = tmp3.double(builder);
    let X_T = tmp5.sub(builder, &tmp3_doubled);
    let Z_T = q.y.add(builder, &q.z);
    let Z_T = Z_T.square(builder);
    let Z_T = Z_T.sub(builder, &tmp1);
    let q_z_squared = q.z.square(builder);
    let Z_T = Z_T.sub(builder, &q_z_squared);
    let Y_T = tmp3.sub(builder, &X_T);
    let Y_T = Y_T.mul(builder, &tmp4);
    let tmp2_8 = tmp2.double(builder);
    let tmp2_8 = tmp2_8.double(builder);
    let tmp2_8 = tmp2_8.double(builder);
    let Y_T = Y_T.sub(builder, &tmp2_8);
    let q_z_doubled = q.z.double(builder);
    let tmp3 = tmp4.mul(builder, &q_z_doubled);
    let tmp3 = tmp3.double(builder);
    let tmp3 = tmp3.neg(builder);
    let tmp3 = tmp3.mul_by_b0(builder, p.x);
    let tmp6 = tmp6.square(builder);
    let tmp6 = tmp6.sub(builder, &tmp0);
    let tmp6 = tmp6.sub(builder, &tmp5);
    let tmp1_4 = tmp1.double(builder);
    let tmp1_4 = tmp1_4.double(builder);
    let tmp6 = tmp6.sub(builder, &tmp1_4);
    let q_z_squared = q.z.square(builder);
    let tmp0 = Z_T.mul(builder, &q_z_squared);
    let tmp0 = tmp0.double(builder);
    let tmp0 = tmp0.mul_by_b0(builder, p.y);
    let a0 = Fq6Target {
        c0: tmp0,
        c1: Fq2Target::zero(builder),
        c2: Fq2Target::zero(builder),
    };

    let a1 = Fq6Target {
        c0: tmp3,
        c1: tmp6,
        c2: Fq2Target::zero(builder),
    };

    let l = Fq12Target { c0: a0, c1: a1 };

    let T = G2PointTarget {
        x: X_T,
        y: Y_T,
        z: Z_T,
    };

    LineEval { l, T }
}

// Algorithm 27 from: https://eprint.iacr.org/2010/354.pdf
fn line_function_add_point<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    q: G2PointTarget<F, D>,
    r: G2PointTarget<F, D>,
    p: G1PointTarget<F, D>,
) -> LineEval<F, D> {
    let r_z_square = r.z.square(builder);
    let t0 = q.x.mul(builder, &r_z_square);
    let t1 = q.y.add(builder, &r.z);
    let t1 = t1.square(builder);
    let q_y_square = q.y.square(builder);
    let t1 = t1.sub(builder, &q_y_square);
    let r_z_square = r.z.square(builder);
    let t1 = t1.sub(builder, &r_z_square);
    let t1 = t1.mul(builder, &r_z_square);
    let t2 = t0.sub(builder, &r.x);
    let t3 = t2.square(builder);
    let t4 = t3.double(builder);
    let t4 = t4.double(builder);
    let t5 = t4.mul(builder, &t2);
    let r_y_double = r.y.double(builder);
    let t6 = t1.sub(builder, &r_y_double);
    let t9 = t6.mul(builder, &q.x);
    let t7 = r.x.mul(builder, &t4);
    let X_T = t6.square(builder);
    let X_T = X_T.sub(builder, &t5);
    let t_7_double = t7.double(builder);
    let X_T = X_T.sub(builder, &t_7_double);
    let Z_T = r.z.add(builder, &t2);
    let Z_T = Z_T.square(builder);
    let r_z_square = r.z.square(builder);
    let Z_T = Z_T.sub(builder, &r_z_square);
    let Z_T = Z_T.sub(builder, &t3);
    let t10 = q.y.add(builder, &Z_T);
    let t8 = t7.sub(builder, &X_T);
    let t8 = t8.mul(builder, &t6);
    let t0 = r.y.mul(builder, &t5);
    let t0 = t0.double(builder);
    let Y_T = t8.sub(builder, &t0);
    let t10 = t10.square(builder);
    let q_y_square = q.y.square(builder);
    let t10 = t10.sub(builder, &q_y_square);
    let z_t_square = Z_T.square(builder);
    let t10 = t10.sub(builder, &z_t_square);
    let t9 = t9.double(builder);
    let t9 = t9.sub(builder, &t10);
    let t10 = Z_T.mul_by_b0(builder, p.y);
    let t10 = t10.double(builder);
    let t6 = t6.neg(builder);
    let t1 = t6.mul_by_b0(builder, p.x);
    let t1 = t1.double(builder);

    let l0 = Fq6Target {
        c0: t10,
        c1: Fq2Target::zero(builder),
        c2: Fq2Target::zero(builder),
    };

    let l1 = Fq6Target {
        c0: t1,
        c1: t9,
        c2: Fq2Target::zero(builder),
    };

    let l = Fq12Target { c0: l0, c1: l1 };

    let T = G2PointTarget {
        x: X_T,
        y: Y_T,
        z: Z_T,
    };

    LineEval { l, T }
}

pub fn noir_miller_loop_implementation<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    Q: G2PointTarget<F, D>,
    P: G1PointTarget<F, D>,
) -> Fq12Target<F, D> {
    let mut T = Q.clone();
    let mut f = Fq12Target::one(builder);

    for i in 0..12 {
        let line_evaluated = line_function_double_point(builder, T.clone(), P.clone());
        f = f.square(builder);
        f = f.mul(builder, &line_evaluated.l);
        T = line_evaluated.T;
        if NAF_DIGIT[i] == 2 {
            let Q_NEG = Q.clone().negate(builder);
            let line_evaluated = line_function_add_point(builder, T.clone(), Q_NEG, P.clone());
            f = f.mul(builder, &line_evaluated.l);
            T = line_evaluated.T;
        } else if NAF_DIGIT[i] == 1 {
            let line_evaluated = line_function_add_point(builder, T.clone(), Q.clone(), P.clone());
            f = f.mul(builder, &line_evaluated.l);
            T = line_evaluated.T;
        }
    }

    // Q1 <- pi_p(Q)
    let q1x = Q.x.conjugate(builder);
    let q1y = Q.y.conjugate(builder);
    let q1x = q1x.mul_by_non_residue_1_power_2(builder);
    let q1y = q1y.mul_by_non_residue_1_power_3(builder);
    let Q1 = G2PointTarget::from_affine(builder, q1x, q1y);

    // Q2 <- pi_p_square(Q);
    let q2x = Q.x.mul_by_non_residue_2_power_2(builder);
    let q2y = Q.y.mul_by_non_residue_2_power_3(builder);
    let q2y = q2y.neg(builder);
    // Q2 is negated above to use directly in the line_function
    let Q2 = G2PointTarget::from_affine(builder, q2x, q2y);

    // Line eval with Q1
    let line_evaluated = line_function_add_point(builder, T.clone(), Q1, P.clone());
    f = f.mul(builder, &line_evaluated.l);
    T = line_evaluated.T;

    // Line eval with Q2
    let line_evaluated = line_function_add_point(builder, T, Q2, P);
    let f = f.mul(builder, &line_evaluated.l);
    f
}

mod tests {
    use ark_bls12_381::{G1Affine, G2Affine};
    use ark_ec::AffineRepr;
    use ark_ff::UniformRand;
    use plonky2::{field::{goldilocks_field::GoldilocksField, types::Field}, iop::witness::{PartialWitness, WitnessWrite}, plonk::{circuit_builder::CircuitBuilder, circuit_data::CircuitConfig, config::PoseidonGoldilocksConfig}};

    use crate::{fields::fq_target::FqTarget, fields_as_trees::{fq2_target_tree::Fq2Target, noir_miller_loop::{noir_miller_loop_implementation, G1PointTarget, G2PointTarget}}};

    type F = GoldilocksField;
    type C = PoseidonGoldilocksConfig;
    const D: usize = 2;

    #[test]
    pub fn noir_miller_loop() {
        let config = CircuitConfig::pairing_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);

        println!("====================================================================================");
        let rng = &mut rand::thread_rng();
        let p0 = G1Affine::rand(rng);
        let p1 = G1Affine::rand(rng);
        let q0 = G2Affine::rand(rng);
        let q1 = G2Affine::rand(rng);

        let zero = FqTarget::zero(&mut builder);
        let zero_fq2 = Fq2Target::zero(&mut builder);
        let p0_x = FqTarget::constant(&mut builder, *p0.x().unwrap());
        let p0_y = FqTarget::constant(&mut builder, *p0.y().unwrap());
        let p0: G1PointTarget<F, D> = G1PointTarget {
            x: p0_x,
            y: p0_y,
            z: zero.clone(),
        };

        let q0_x = Fq2Target {
            c0: FqTarget::constant(&mut builder, q0.x().unwrap().c0),
            c1: FqTarget::constant(&mut builder, q0.x().unwrap().c1),
        };
        let q0_y = Fq2Target {
            c0: FqTarget::constant(&mut builder, q0.y().unwrap().c0),
            c1: FqTarget::constant(&mut builder, q0.y().unwrap().c1),
        };

        let q0: G2PointTarget<F, D> = G2PointTarget {
            x: q0_x,
            y: q0_y,
            z: zero_fq2.clone()
        };

        let p1_x = FqTarget::constant(&mut builder, *p1.x().unwrap());
        let p1_y = FqTarget::constant(&mut builder, *p1.y().unwrap());

        let p1: G1PointTarget<F, D> = G1PointTarget {
            x: p1_x,
            y: p1_y,
            z: zero,
        };

        let q1_x = Fq2Target {
            c0: FqTarget::constant(&mut builder, q1.x().unwrap().c0),
            c1: FqTarget::constant(&mut builder, q1.x().unwrap().c1),
        };
        let q1_y = Fq2Target {
            c0: FqTarget::constant(&mut builder, q1.y().unwrap().c0),
            c1: FqTarget::constant(&mut builder, q1.y().unwrap().c1),
        };

        let q1: G2PointTarget<F, D> = G2PointTarget {
            x: q1_x,
            y: q1_y,
            z: zero_fq2
        };

        let q0 = G2PointTarget::generator(&mut builder);
        let p0 = G1PointTarget::generator(&mut builder);

        let x = noir_miller_loop_implementation(&mut builder, q0, p0);
        println!("x is: {:?}", x);

        let mut pw = PartialWitness::new();
        let data = builder.build::<C>();
        dbg!(data.common.degree_bits());
        let _proof = data.prove(pw);
    }
}