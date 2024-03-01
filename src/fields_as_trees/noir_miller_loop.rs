use plonky2::{
    field::extension::Extendable, hash::hash_types::RichField,
    plonk::circuit_builder::CircuitBuilder,
};

use crate::{fields::fq_target::FqTarget, utils::constants::NAF_DIGIT};

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
}

#[derive(Clone, Debug)]
pub struct G1PointTarget<F: RichField + Extendable<D>, const D: usize> {
    pub x: FqTarget<F, D>,
    pub y: FqTarget<F, D>,
    pub z: FqTarget<F, D>,
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
    let mut T = Q;
    let mut f = Fq12Target::one(builder);

    for i in 0..12 {
        let line_evaluated = line_function_double_point(builder, T.clone(), P.clone());
        f = f.square(builder);
        f = f.mul(builder, &line_evaluated.l);
        T = line_evaluated.T;
        if NAF_DIGIT[i] == 2 {
            let Q_NEG = Q.negate(builder);
            let line_evaluated = line_function_add_point(builder, T, Q_NEG, P);
            f = f.mul(builder, &line_evaluated.l);
            T = line_evaluated.T;
        } else if NAF_DIGIT[i] == 1 {
            let line_evaluated = line_function_add_point(builder, T, Q, P);
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
    let line_evaluated = line_function_add_point(builder, T, Q1, P);
    f = f.mul(builder, &line_evaluated.l);
    T = line_evaluated.T;

    // Line eval with Q2
    let line_evaluated = line_function_add_point(builder, T, Q2, P);
    let f = f.mul(builder, &line_evaluated.l);
    f
}
