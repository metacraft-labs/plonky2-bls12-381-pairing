use ark_bls12_381::{Fq, Fq2};
use ark_ec::short_weierstrass::SWCurveConfig;
use ark_ff::BitIteratorBE;
use num::One;
use plonky2::{
    field::extension::Extendable, hash::hash_types::RichField,
    plonk::circuit_builder::CircuitBuilder,
};

use crate::{
    fields::{fq2_target::Fq2Target, fq_target::FqTarget},
    utils::constants::BLS_X,
};

#[derive(Clone, Debug)]
pub struct G2AffineTarget<F: RichField + Extendable<D>, const D: usize> {
    pub x: Fq2Target<F, D>,
    pub y: Fq2Target<F, D>,
    pub infinity: bool,
}

impl<F: RichField + Extendable<D>, const D: usize> G2AffineTarget<F, D> {
    fn xy(&self) -> Option<(&self::Fq2Target<F, D>, &self::Fq2Target<F, D>)> {
        (!self.infinity).then(|| (&self.x, &self.y))
    }
}

#[derive(Clone, Debug)]
pub struct G2ProjectiveTarget<F: RichField + Extendable<D>, const D: usize> {
    pub x: Fq2Target<F, D>,
    pub y: Fq2Target<F, D>,
    pub z: Fq2Target<F, D>,
}

impl<F: RichField + Extendable<D>, const D: usize> G2ProjectiveTarget<F, D> {
    fn double_in_place(
        &mut self,
        builder: &mut CircuitBuilder<F, D>,
        two_inv: &FqTarget<F, D>,
    ) -> EllCoeffTarget<F, D> {
        let a = self.x.mul(builder, &self.y);
        let a = a.mul_assign_by_fp(builder, two_inv.clone());
        let b = self.y.simple_square(builder);
        let c = self.z.simple_square(builder);
        let coeff_b = ark_bls12_381::g2::Config::COEFF_B;
        let coeff_b = Fq2Target::constant(builder, coeff_b);
        let c_double = c.double(builder);
        let c_double = c_double.add(builder, &c);
        let e = coeff_b.mul(builder, &c_double);
        let e_double = e.double(builder);
        let f = e_double.add(builder, &e);
        let g = b.add(builder, &f);
        let g = g.mul_assign_by_fp(builder, two_inv.clone());
        let y_add_z = self.y.add(builder, &self.z);
        let y_z_squared = y_add_z.simple_square(builder);
        let b_add_c = b.add(builder, &c);
        let h = y_z_squared.sub(builder, &b_add_c);
        let i = e.sub(builder, &b);
        let j = self.x.simple_square(builder);
        let e_square = e.simple_square(builder);

        let b_sub_f = b.sub(builder, &f);

        self.x = a.mul(builder, &b_sub_f);
        let g_square = g.simple_square(builder);
        let e_sq_double = e_square.double(builder);
        let e_sq_double = e_sq_double.add(builder, &e_square);
        self.y = g_square.sub(builder, &e_sq_double);
        self.z = b.mul(builder, &h);

        let j_double = j.double(builder);
        let j_triple = j_double.add(builder, &j);

        (i, j_triple, h.neg(builder))
    }

    fn add_in_place(
        &mut self,
        builder: &mut CircuitBuilder<F, D>,
        q: &G2AffineTarget<F, D>,
    ) -> EllCoeffTarget<F, D> {
        let (qx, qy) = q.xy().unwrap();
        let qy_z = qy.mul(builder, &self.z);
        let theta = self.y.sub(builder, &qy_z);
        let qx_z = qx.mul(builder, &self.z);
        let lambda = self.x.sub(builder, &qx_z);
        let c = theta.simple_square(builder);
        let d = lambda.simple_square(builder);
        let e = lambda.mul(builder, &d);
        let f = self.z.mul(builder, &c);
        let g = self.x.mul(builder, &d);
        let g_double = g.double(builder);
        let e_add_f = e.add(builder, &f);
        let h = e_add_f.add(builder, &g_double);
        self.x = lambda.mul(builder, &h);
        let g_sub_h = g.sub(builder, &h);
        let e_mul_y = e.mul(builder, &self.y);
        let theta_mul_g_sub_h = theta.mul(builder, &g_sub_h);
        self.y = theta_mul_g_sub_h.sub(builder, &e_mul_y);
        self.z = self.z.mul(builder, &e);
        let lambda_qy = lambda.mul(builder, &qy);
        let theta_qx = theta.mul(builder, &qx);
        let j = theta_qx.sub(builder, &lambda_qy);

        (j, theta.neg(builder), lambda)
    }
}

#[derive(Clone, Debug)]
pub struct G2PreparedTarget<F: RichField + Extendable<D>, const D: usize> {
    /// Stores the coefficients of the line evaluations as calculated in
    /// <https://eprint.iacr.org/2013/722.pdf>
    pub ell_coeffs: Vec<EllCoeffTarget<F, D>>,
    pub infinity: bool,
}

pub(crate) type EllCoeffTarget<F, const D: usize> =
    (Fq2Target<F, D>, Fq2Target<F, D>, Fq2Target<F, D>);

impl<F: RichField + Extendable<D>, const D: usize> G2PreparedTarget<F, D> {
    fn is_zero(&self) -> bool {
        self.infinity
    }

    fn from(builder: &mut CircuitBuilder<F, D>, q: G2AffineTarget<F, D>) -> Self {
        let one = FqTarget::constant(builder, Fq::one()); // Fq::two
        let two = one.add(builder, &one);
        let two_inv = two.inv(builder);
        let zero = G2PreparedTarget {
            ell_coeffs: vec![],
            infinity: true,
        };

        q.xy().map_or(zero, |(q_x, q_y)| {
            let mut ell_coeffs = vec![];
            let mut r = G2ProjectiveTarget {
                x: q_x.clone(),
                y: q_y.clone(),
                z: Fq2Target::constant(builder, Fq2::one()),
            };

            for i in BitIteratorBE::new([BLS_X]).skip(1) {
                ell_coeffs.push(r.double_in_place(builder, &two_inv));

                if i {
                    ell_coeffs.push(r.add_in_place(builder, &q));
                }
            }

            Self {
                ell_coeffs,
                infinity: false,
            }
        })
    }
}
