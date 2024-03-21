use plonky2::{
    field::extension::Extendable,
    hash::hash_types::RichField,
    iop::target::BoolTarget,
    plonk::{circuit_builder::CircuitBuilder, circuit_data::CircuitConfig},
};

use crate::utils::constants::{BLS_X, BLS_X_IS_NEGATIVE};

use super::{
    fq12_target_tree::Fq12Target,
    fq2_target_tree::Fq2Target,
    fq6_target_tree::Fq6Target,
    g1_curve::G1AffineTarget,
    g2_curve::{G2AffineTarget, G2ProjectiveTarget},
};

#[derive(Clone, Debug)]
pub struct Gt<F: RichField + Extendable<D>, const D: usize>(pub Fq12Target<F, D>);

impl<F: RichField + Extendable<D>, const D: usize> Gt<F, D> {
    // comparison
    pub fn compare(&self, builder: &mut CircuitBuilder<F, D>, rhs: &Self) {
        Fq12Target::connect(builder, &self.0, &rhs.0)
    }
}

#[derive(Clone, Debug)]
pub struct MillerLoopResult<F: RichField + Extendable<D>, const D: usize>(pub Fq12Target<F, D>);

impl<F: RichField + Extendable<D>, const D: usize> MillerLoopResult<F, D> {
    pub fn default(builder: &mut CircuitBuilder<F, D>) -> Self {
        MillerLoopResult(Fq12Target::one(builder))
    }

    pub fn fq4_square(
        builder: &mut CircuitBuilder<F, D>,
        a: Fq2Target<F, D>,
        b: Fq2Target<F, D>,
    ) -> (Fq2Target<F, D>, Fq2Target<F, D>) {
        let t0 = a.square(builder);
        let t1 = b.square(builder);
        let mut t2 = t1.mul_by_nonresidue(builder);
        let c0 = t2.add(builder, &t0);
        t2 = a.add(builder, &b);
        t2 = t2.square(builder);
        t2 = t2.sub(builder, &t0);
        let c1 = t2.sub(builder, &t1);

        (c0, c1)
    }

    pub fn cyclotomic_square(
        builder: &mut CircuitBuilder<F, D>,
        f: Fq12Target<F, D>,
    ) -> Fq12Target<F, D> {
        let mut z0 = f.c0.c0;
        let mut z4 = f.c0.c1;
        let mut z3 = f.c0.c2;
        let mut z2 = f.c1.c0;
        let mut z1 = f.c1.c1;
        let mut z5 = f.c1.c2;

        let (t0, t1) = Self::fq4_square(builder, z0.clone(), z1.clone());

        // For A
        z0 = t0.sub(builder, &z0);
        z0 = z0.add(builder, &z0);
        z0 = z0.add(builder, &t0);

        z1 = t1.add(builder, &z1);
        z1 = z1.add(builder, &z1);
        z1 = z1.add(builder, &t1);

        let (mut t0, t1) = Self::fq4_square(builder, z2.clone(), z3.clone());
        let (t2, t3) = Self::fq4_square(builder, z4.clone(), z5.clone());

        // For C
        z4 = t0.sub(builder, &z4);
        z4 = z4.add(builder, &z4);
        z4 = z4.add(builder, &t0);

        z5 = t1.add(builder, &z5);
        z5 = z5.add(builder, &z5);
        z5 = z5.add(builder, &t1);

        // For B

        t0 = t3.mul_by_nonresidue(builder);
        z2 = t0.add(builder, &z2);
        z2 = z2.add(builder, &z2);
        z2 = z2.add(builder, &t0);

        z3 = t2.sub(builder, &z3);
        z3 = z3.add(builder, &z3);
        z3 = z3.add(builder, &t2);

        Fq12Target {
            c0: Fq6Target {
                c0: z0,
                c1: z4,
                c2: z3,
            },

            c1: Fq6Target {
                c0: z2,
                c1: z1,
                c2: z5,
            },
        }
    }

    pub fn cycolotomic_exp(
        builder: &mut CircuitBuilder<F, D>,
        f: Fq12Target<F, D>,
    ) -> Fq12Target<F, D> {
        let x = BLS_X;
        let mut tmp = Fq12Target::one(builder);
        let mut found_one = false;
        for i in (0..64).rev().map(|b| ((x >> b) & 1) == 1) {
            if found_one {
                tmp = Self::cyclotomic_square(builder, tmp)
            } else {
                found_one = i;
            }

            if i {
                tmp.mul(builder, &f);
            }
        }

        tmp.conjugate(builder)
    }

    pub fn f_conversion(
        builder: &mut CircuitBuilder<F, D>,
        mut t0: Fq12Target<F, D>,
        mut t1: Fq12Target<F, D>,
    ) -> Fq12Target<F, D> {
        let mut t2 = t0.mul(builder, &t1);
        t1 = t2.clone();
        t2 = t2.frobenius_map(builder).frobenius_map(builder);
        t2 = t2.mul(builder, &t1);
        t1 = Self::cyclotomic_square(builder, t2.clone()).conjugate(builder);
        let mut t3 = Self::cycolotomic_exp(builder, t2.clone());
        let mut t4 = Self::cyclotomic_square(builder, t3.clone());
        let mut t5 = t1.mul(builder, &t3);
        t1 = Self::cycolotomic_exp(builder, t5.clone());
        t0 = Self::cycolotomic_exp(builder, t1.clone());
        let mut t6 = Self::cycolotomic_exp(builder, t0.clone());
        t6 = t6.mul(builder, &t4);
        t4 = Self::cycolotomic_exp(builder, t6.clone());
        t5 = t5.conjugate(builder);
        let xx = t5.mul(builder, &t2);
        t4 = t4.mul(builder, &xx);
        t5 = t2.conjugate(builder);
        t1 = t1.mul(builder, &t2);
        t1 = t1
            .frobenius_map(builder)
            .frobenius_map(builder)
            .frobenius_map(builder);
        t6 = t6.mul(builder, &t5);
        t6 = t6.frobenius_map(builder);
        t3 = t3.mul(builder, &t0);
        t3 = t3.frobenius_map(builder).frobenius_map(builder);
        t3 = t3.mul(builder, &t1);
        t3 = t3.mul(builder, &t6);
        let f = t3.mul(builder, &t4);

        f
    }

    pub fn final_exponentiation(&self, builder: &mut CircuitBuilder<F, D>) -> Gt<F, D> {
        let f = &self.0;
        let t0 = f
            .frobenius_map(builder)
            .frobenius_map(builder)
            .frobenius_map(builder)
            .frobenius_map(builder)
            .frobenius_map(builder)
            .frobenius_map(builder);
        let f_inverted = f.inv(builder);

        Gt(Self::f_conversion(builder, t0, f_inverted))
    }
}

#[derive(Clone, Debug)]
pub struct G2PreparedTarget<F: RichField + Extendable<D>, const D: usize> {
    pub infinity: BoolTarget,
    pub coeffs: Vec<(Fq2Target<F, D>, Fq2Target<F, D>, Fq2Target<F, D>)>,
}

impl<F: RichField + Extendable<D>, const D: usize> From<G2AffineTarget<F, D>>
    for G2PreparedTarget<F, D>
{
    fn from(q: G2AffineTarget<F, D>) -> Self {
        struct Adder<A: RichField + Extendable<B>, const B: usize> {
            cur: G2ProjectiveTarget<A, B>,
            base: G2AffineTarget<A, B>,
            coeffs: Vec<(Fq2Target<A, B>, Fq2Target<A, B>, Fq2Target<A, B>)>,
        }

        impl<A: RichField + Extendable<B>, const B: usize> MillerLoopDriver for Adder<A, B> {
            type Output = ();

            fn point_doubling_and_line_evaluation(&mut self, _: Self::Output) -> Self::Output {
                let config = CircuitConfig::pairing_config();
                let mut builder = CircuitBuilder::<A, B>::new(config);
                let coeffs = point_doubling_and_line_evaluation(&mut builder, &mut self.cur);
                self.coeffs.push(coeffs);
            }
            fn point_addition_and_line_evaluation(&mut self, _: Self::Output) -> Self::Output {
                let config = CircuitConfig::pairing_config();
                let mut builder = CircuitBuilder::<A, B>::new(config);
                let coeffs =
                    point_addition_and_line_evaluation(&mut builder, &mut self.cur, &mut self.base);
                self.coeffs.push(coeffs);
            }
            fn square_output(_: Self::Output) -> Self::Output {}
            fn conjugate(_: Self::Output) -> Self::Output {}
            fn one() -> Self::Output {}
        }

        let config = CircuitConfig::pairing_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);

        let is_identity = q.is_identity();
        let g2_generator = G2AffineTarget::experimental_generator(&mut builder);
        let q = G2AffineTarget::conditional_select(&mut builder, q, g2_generator, is_identity);
        let mut adder = Adder {
            cur: G2ProjectiveTarget::from(q.clone()),
            base: q,
            coeffs: Vec::with_capacity(68), //WHY?
        };

        miller_loop::<F, D, _>(&mut adder);

        assert_eq!(adder.coeffs.len(), 68);

        G2PreparedTarget {
            infinity: is_identity,
            coeffs: adder.coeffs,
        }
    }
}

trait MillerLoopDriver {
    type Output;

    fn point_doubling_and_line_evaluation(&mut self, f: Self::Output) -> Self::Output;
    fn point_addition_and_line_evaluation(&mut self, f: Self::Output) -> Self::Output;
    fn square_output(f: Self::Output) -> Self::Output;
    fn conjugate(f: Self::Output) -> Self::Output;
    fn one() -> Self::Output;
}

pub fn multi_miller_loop<F: RichField + Extendable<D>, const D: usize>(
    terms: &[(&G1AffineTarget<F, D>, &G2PreparedTarget<F, D>)],
) -> MillerLoopResult<F, D> {
    struct Adder<'a, 'b, 'c, A: RichField + Extendable<B>, const B: usize> {
        terms: &'c [(&'a G1AffineTarget<A, B>, &'b G2PreparedTarget<A, { B }>)],
        index: usize,
    }

    impl<'a, 'b, 'c, A: RichField + Extendable<B>, const B: usize> MillerLoopDriver
        for Adder<'a, 'b, 'c, A, B>
    {
        type Output = Fq12Target<A, B>;

        fn point_doubling_and_line_evaluation(&mut self, mut f: Self::Output) -> Self::Output {
            let config = CircuitConfig::pairing_config();
            let mut builder = CircuitBuilder::<A, B>::new(config);
            let index = self.index;
            for term in self.terms {
                let either_identity = builder.or(term.0.is_identity(), term.1.infinity);

                let new_f = ell(&mut builder, f.clone(), &term.1.coeffs[index], term.0);
                // conditional_select ->
                // f = Fq12Target::select(&mut builder, &new_f, &f, &either_identity);
                f = Fq12Target::conditional_select(&mut builder, &new_f, &f, either_identity);
            }
            self.index += 1;

            f
        }

        fn point_addition_and_line_evaluation(&mut self, mut f: Self::Output) -> Self::Output {
            let config = CircuitConfig::pairing_config();
            let mut builder = CircuitBuilder::<A, B>::new(config);
            let index = self.index;
            for term in self.terms {
                let either_identity = builder.or(term.0.is_identity(), term.1.infinity);

                let new_f = ell(&mut builder, f.clone(), &term.1.coeffs[index], term.0);
                // conditional_select ->
                // f = Fq12Target::select(&mut builder, &new_f, &f, &either_identity);
                f = Fq12Target::conditional_select(&mut builder, &new_f, &f, either_identity);
            }
            self.index += 1;

            f
        }

        fn square_output(f: Self::Output) -> Self::Output {
            let config = CircuitConfig::pairing_config();
            let mut builder = CircuitBuilder::<A, B>::new(config);
            f.square(&mut builder)
        }

        fn conjugate(f: Self::Output) -> Self::Output {
            let config = CircuitConfig::pairing_config();
            let mut builder = CircuitBuilder::<A, B>::new(config);
            f.conjugate(&mut builder)
        }

        fn one() -> Self::Output {
            let config = CircuitConfig::pairing_config();
            let mut builder = CircuitBuilder::<A, B>::new(config);
            Fq12Target::one(&mut builder)
        }
    }

    let mut adder = Adder { terms, index: 0 };

    let tmp = miller_loop::<F, D, _>(&mut adder);

    MillerLoopResult(tmp)
}

fn miller_loop<F: RichField + Extendable<D>, const D: usize, M: MillerLoopDriver>(
    driver: &mut M,
) -> M::Output {
    let mut f = M::one();
    let mut found_one = false;

    for i in (0..64).rev().map(|b| (((BLS_X >> 1) >> b) & 1) == 1) {
        if !found_one {
            found_one = i;
            continue;
        }

        f = driver.point_doubling_and_line_evaluation(f);

        if i {
            f = driver.point_addition_and_line_evaluation(f);
        }

        f = M::square_output(f);
    }

    f = driver.point_doubling_and_line_evaluation(f);

    if BLS_X_IS_NEGATIVE {
        f = M::conjugate(f);
    }

    f
}

// VERIFIED I TIME WHILE FOCUSED
// Adaptation of Algorithm 26, https://eprint.iacr.org/2010/354.pdf
fn point_doubling_and_line_evaluation<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    r: &mut G2ProjectiveTarget<F, D>,
) -> (Fq2Target<F, D>, Fq2Target<F, D>, Fq2Target<F, D>) {
    let tmp0 = r.x.square(builder);
    let tmp1 = r.y.square(builder);
    let tmp2 = tmp1.square(builder);
    let tmp3 = tmp1.add(builder, &r.x);
    let tmp3 = tmp3.square(builder);
    let tmp3 = tmp3.sub(builder, &tmp0);
    let tmp3 = tmp3.sub(builder, &tmp2);
    let tmp3 = tmp3.add(builder, &tmp3);
    let tmp0_double = tmp0.add(builder, &tmp0);
    let tmp4 = tmp0.add(builder, &tmp0_double);
    let tmp6 = r.x.add(builder, &tmp4);
    let tmp5 = tmp4.square(builder);
    let zsquared = r.z.square(builder);
    r.x = tmp5.sub(builder, &tmp3);
    r.x = r.x.sub(builder, &tmp3);
    let rz_ry = r.z.add(builder, &r.y);
    r.z = rz_ry.square(builder);
    r.z = r.z.sub(builder, &tmp1);
    r.z = r.z.sub(builder, &zsquared);
    let tmp3_min_rx = tmp3.sub(builder, &r.x);
    r.y = tmp3_min_rx.mul(builder, &tmp4);
    let tmp2 = tmp2.add(builder, &tmp2);
    let tmp2 = tmp2.add(builder, &tmp2);
    let tmp2 = tmp2.add(builder, &tmp2);
    r.y = r.y.sub(builder, &tmp2);
    let tmp3 = tmp4.mul(builder, &zsquared);
    let tmp3 = tmp3.add(builder, &tmp3);
    let tmp3 = tmp3.conjugate(builder);
    let tmp6 = tmp6.square(builder);
    let tmp6 = tmp6.sub(builder, &tmp0);
    let tmp6 = tmp6.sub(builder, &tmp5);
    let tmp1 = tmp1.add(builder, &tmp1);
    let tmp1 = tmp1.add(builder, &tmp1);
    let tmp6 = tmp6.sub(builder, &tmp1);
    let tmp0 = r.z.mul(builder, &zsquared);
    let tmp0 = tmp0.add(builder, &tmp0);

    (tmp0, tmp3, tmp6)
}

// VERIFIED I TIME WHILE FOCUSED
// Adaptation of Algorithm 27, https://eprint.iacr.org/2010/354.pdf
fn point_addition_and_line_evaluation<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    r: &mut G2ProjectiveTarget<F, D>,
    q: &G2AffineTarget<F, D>,
) -> (Fq2Target<F, D>, Fq2Target<F, D>, Fq2Target<F, D>) {
    let zsquared = r.z.square(builder);
    let ysquared = q.y.square(builder);
    let t0 = zsquared.mul(builder, &q.x);
    let qy_add_rz = q.y.add(builder, &r.z);
    let qy_rz_squared = qy_add_rz.square(builder);
    let qy_rz_squared = qy_rz_squared.sub(builder, &ysquared);
    let qy_rz_squared = qy_rz_squared.sub(builder, &zsquared);
    let t1 = qy_rz_squared.mul(builder, &zsquared);
    let t2 = t0.sub(builder, &r.x);
    let t3 = t2.square(builder);
    let t4 = t3.add(builder, &t3);
    let t4 = t4.add(builder, &t4);
    let t5 = t4.mul(builder, &t2);
    let t6 = t1.sub(builder, &r.y);
    let t6 = t6.sub(builder, &r.y);
    let t9 = t6.mul(builder, &q.x);
    let t7 = t4.mul(builder, &r.x);
    r.x = t6.square(builder);
    r.x = r.x.sub(builder, &t5);
    r.x = r.x.sub(builder, &t7);
    r.x = r.x.sub(builder, &t7);
    let rz_add_t2 = r.z.add(builder, &t2);
    r.z = rz_add_t2.square(builder);
    r.z = r.z.sub(builder, &zsquared);
    r.z = r.z.sub(builder, &t3);
    let t10 = q.y.add(builder, &r.z);
    let t7_sub_rx = t7.sub(builder, &r.x);
    let t8 = t7_sub_rx.mul(builder, &t6);
    let t0 = r.y.mul(builder, &t5);
    let t0 = t0.add(builder, &t0);
    r.y = t8.sub(builder, &t0);
    let t10 = t10.square(builder);
    let t10 = t10.sub(builder, &ysquared);
    let ztsquared = r.z.square(builder);
    let t10 = t10.sub(builder, &ztsquared);
    let t9 = t9.add(builder, &t9);
    let t9 = t9.sub(builder, &t10);
    let t10 = r.z.add(builder, &r.z);
    let t6 = t6.conjugate(builder);
    let t1 = t6.add(builder, &t6);

    (t10, t1, t9)
}

// VERIFIED I TIME
fn ell<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    f: Fq12Target<F, D>,
    coeffs: &(Fq2Target<F, D>, Fq2Target<F, D>, Fq2Target<F, D>),
    p: &G1AffineTarget<F, D>,
) -> Fq12Target<F, D> {
    let c0 = &coeffs.0;
    let c1 = &coeffs.1;

    c0.c0.mul(builder, &p.y);
    c0.c1.mul(builder, &p.y);

    c1.c0.mul(builder, &p.x);
    c1.c1.mul(builder, &p.x);

    f.mul_by_014(builder, &coeffs.2, &c1, &c0)
}

pub fn pairing<F: RichField + Extendable<D>, const D: usize>(
    p: &G1AffineTarget<F, D>,
    q: &G2AffineTarget<F, D>,
) -> Gt<F, D> {
    struct Adder<A: RichField + Extendable<B>, const B: usize> {
        cur: G2ProjectiveTarget<A, B>,
        base: G2AffineTarget<A, B>,
        p: G1AffineTarget<A, B>,
    }

    impl<A: RichField + Extendable<B>, const B: usize> MillerLoopDriver for Adder<A, B> {
        type Output = Fq12Target<A, B>;

        fn point_doubling_and_line_evaluation(&mut self, f: Self::Output) -> Self::Output {
            let config = CircuitConfig::pairing_config();
            let mut builder = CircuitBuilder::<A, B>::new(config);
            let coeffs = point_doubling_and_line_evaluation(&mut builder, &mut self.cur);
            ell(&mut builder, f, &coeffs, &self.p)
        }

        fn point_addition_and_line_evaluation(&mut self, f: Self::Output) -> Self::Output {
            let config = CircuitConfig::pairing_config();
            let mut builder = CircuitBuilder::<A, B>::new(config);
            let coeffs =
                point_addition_and_line_evaluation(&mut builder, &mut self.cur, &self.base);
            ell(&mut builder, f, &coeffs, &self.p)
        }

        fn square_output(f: Self::Output) -> Self::Output {
            let config = CircuitConfig::pairing_config();
            let mut builder = CircuitBuilder::<A, B>::new(config);
            f.square(&mut builder)
        }

        fn conjugate(f: Self::Output) -> Self::Output {
            let config = CircuitConfig::pairing_config();
            let mut builder = CircuitBuilder::<A, B>::new(config);
            f.conjugate(&mut builder)
        }

        fn one() -> Self::Output {
            let config = CircuitConfig::pairing_config();
            let mut builder = CircuitBuilder::<A, B>::new(config);
            Fq12Target::one(&mut builder)
        }
    }

    let config = CircuitConfig::pairing_config();
    let mut builder = CircuitBuilder::<F, D>::new(config);
    let either_identity = builder.or(p.is_identity(), q.infinity);

    let g1_generator = G1AffineTarget::experimental_generator(&mut builder);
    let p =
        G1AffineTarget::conditional_select(&mut builder, p.clone(), g1_generator, either_identity);

    let g2_generator = G2AffineTarget::experimental_generator(&mut builder);
    let q =
        G2AffineTarget::conditional_select(&mut builder, q.clone(), g2_generator, either_identity);

    let mut adder = Adder {
        cur: G2ProjectiveTarget::from(q.clone()),
        base: q,
        p,
    };

    let tmp = miller_loop::<F, D, _>(&mut adder);

    let one = Fq12Target::one(&mut builder);
    let tmp = MillerLoopResult(Fq12Target::conditional_select(
        &mut builder,
        &tmp,
        &one,
        either_identity,
    ));
    tmp.final_exponentiation(&mut builder)
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Fq, Fq12, Fr, G1Affine, G2Affine};
    use ark_ec::AffineRepr;
    use ark_ff::{Field, UniformRand};
    use num::BigUint;
    use plonky2::{
        field::{goldilocks_field::GoldilocksField, types::Field as OtherField},
        iop::witness::{PartialWitness, WitnessWrite},
        plonk::{
            circuit_builder::CircuitBuilder, circuit_data::CircuitConfig,
            config::PoseidonGoldilocksConfig,
        },
    };

    use crate::{
        fields::fq_target::FqTarget,
        fields_as_trees::{
            fq12_target_tree::Fq12Target,
            fq2_target_tree::Fq2Target,
            fq6_target_tree::Fq6Target,
            g1_curve::G1AffineTarget,
            g2_curve::G2AffineTarget,
            miller_loop::{multi_miller_loop, MillerLoopResult},
        },
    };

    type F = GoldilocksField;
    type C = PoseidonGoldilocksConfig;
    const D: usize = 2;

    #[test]
    fn test_gt_generator() {
        let config = CircuitConfig::pairing_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let rng = &mut rand::thread_rng();
        let p0 = G1Affine::rand(rng);
        // println!("G1Affine: {}");
        println!(
            "===================================================================================="
        );
        let p1 = G1Affine::rand(rng);
        let q0 = G2Affine::rand(rng);
        let q1 = G2Affine::rand(rng);

        let q0_x = Fq2Target {
            c0: FqTarget::constant(&mut builder, q0.x().unwrap().c0),
            c1: FqTarget::constant(&mut builder, q0.x().unwrap().c1),
        };
        let q0_y = Fq2Target {
            c0: FqTarget::constant(&mut builder, q0.y().unwrap().c0),
            c1: FqTarget::constant(&mut builder, q0.y().unwrap().c1),
        };
        let q0: G2AffineTarget<F, D> = G2AffineTarget {
            x: q0_x.clone(),
            y: q0_y.clone(),
            infinity: builder._false(),
        };

        let p0_x = FqTarget::constant(&mut builder, *p0.x().unwrap());
        let p0_y = FqTarget::constant(&mut builder, *p0.y().unwrap());
        let p0: G1AffineTarget<F, D> = G1AffineTarget {
            x: p0_x.clone(),
            y: p0_y.clone(),
            infinity: builder._false(),
        };

        let p1_x = FqTarget::constant(&mut builder, *p1.x().unwrap());
        let p1_y = FqTarget::constant(&mut builder, *p1.y().unwrap());
        let p1: G1AffineTarget<F, D> = G1AffineTarget {
            x: p1_x.clone(),
            y: p1_y.clone(),
            infinity: builder._false(),
        };

        let q1_x = Fq2Target {
            c0: FqTarget::constant(&mut builder, q1.x().unwrap().c0),
            c1: FqTarget::constant(&mut builder, q1.x().unwrap().c1),
        };
        let q1_y = Fq2Target {
            c0: FqTarget::constant(&mut builder, q1.y().unwrap().c0),
            c1: FqTarget::constant(&mut builder, q1.y().unwrap().c1),
        };

        let q1: G2AffineTarget<F, D> = G2AffineTarget {
            x: q1_x.clone(),
            y: q1_y.clone(),
            infinity: builder._false(),
        };

        let r0 = multi_miller_loop(&[(&p0.clone(), &q0.clone().into())]);
        let r1 = multi_miller_loop(&[(&p1.clone(), &q1.clone().into())]);

        println!("r0: {:?}", &r0);
        println!(
            "===================================================================================="
        );
        // let r_expected = r0.0.mul(&mut builder, &r1.0);
        // let r = multi_miller_loop(&[(
        //     &P0, &Q0.into()), (&P1, &Q1.into()
        // )]).0;

        println!("r1: {:?}", &r1);

        let one = Fq12Target::one(&mut builder);
        let _one = Fq12Target::one(&mut builder);
        let x = one.is_equal(&mut builder, &_one);
        // Fq12Target::connect(&mut builder, &_one, &one);

        let mut pw = PartialWitness::new();
        pw.set_target(x.target, F::ONE);
        let data = builder.build::<C>();
        dbg!(data.common.degree_bits());
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_miller_loop() {
        let config = CircuitConfig::pairing_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);

        let result_mml = multi_miller_loop(&[(
            &G1AffineTarget::<F, D>::generator(),
            &G2AffineTarget::<F, D>::identity(&mut builder).into(),
        )])
        .0;

        let one = Fq12Target::one(&mut builder);

        Fq12Target::connect(&mut builder, &result_mml, &one);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        dbg!(data.common.degree_bits());
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_final_exponentiation() {
        let config = CircuitConfig::pairing_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let rng = &mut rand::thread_rng();
        let x = Fq12::rand(rng);
        let x_target = Fq12Target {
            c0: Fq6Target {
                c0: Fq2Target {
                    c0: FqTarget::constant(&mut builder, x.c0.c0.c0),
                    c1: FqTarget::constant(&mut builder, x.c0.c0.c1),
                },
                c1: Fq2Target {
                    c0: FqTarget::constant(&mut builder, x.c0.c1.c0),
                    c1: FqTarget::constant(&mut builder, x.c0.c1.c1),
                },
                c2: Fq2Target {
                    c0: FqTarget::constant(&mut builder, x.c0.c2.c0),
                    c1: FqTarget::constant(&mut builder, x.c0.c2.c1),
                },
            },
            c1: Fq6Target {
                c0: Fq2Target {
                    c0: FqTarget::constant(&mut builder, x.c1.c0.c0),
                    c1: FqTarget::constant(&mut builder, x.c1.c0.c1),
                },
                c1: Fq2Target {
                    c0: FqTarget::constant(&mut builder, x.c1.c1.c0),
                    c1: FqTarget::constant(&mut builder, x.c1.c1.c1),
                },
                c2: Fq2Target {
                    c0: FqTarget::constant(&mut builder, x.c1.c2.c0),
                    c1: FqTarget::constant(&mut builder, x.c1.c2.c1),
                },
            },
        };
        let ops = MillerLoopResult(x_target);
        let final_x = MillerLoopResult::final_exponentiation(&ops, &mut builder);

        use ark_ff::PrimeField;
        let p: BigUint = Fq::MODULUS.into();
        let r: BigUint = Fr::MODULUS.into();
        let exp = (p.pow(12) - 1u32) / r;
        let final_x2 = x.pow(&exp.to_u64_digits());
        let final_x2 = Fq12Target {
            c0: Fq6Target {
                c0: Fq2Target {
                    c0: FqTarget::constant(&mut builder, final_x2.c0.c0.c0),
                    c1: FqTarget::constant(&mut builder, final_x2.c0.c0.c1),
                },
                c1: Fq2Target {
                    c0: FqTarget::constant(&mut builder, final_x2.c0.c1.c0),
                    c1: FqTarget::constant(&mut builder, final_x2.c0.c1.c1),
                },
                c2: Fq2Target {
                    c0: FqTarget::constant(&mut builder, final_x2.c0.c2.c0),
                    c1: FqTarget::constant(&mut builder, final_x2.c0.c2.c1),
                },
            },
            c1: Fq6Target {
                c0: Fq2Target {
                    c0: FqTarget::constant(&mut builder, final_x2.c1.c0.c0),
                    c1: FqTarget::constant(&mut builder, final_x2.c1.c0.c1),
                },
                c1: Fq2Target {
                    c0: FqTarget::constant(&mut builder, final_x2.c1.c1.c0),
                    c1: FqTarget::constant(&mut builder, final_x2.c1.c1.c1),
                },
                c2: Fq2Target {
                    c0: FqTarget::constant(&mut builder, final_x2.c1.c2.c0),
                    c1: FqTarget::constant(&mut builder, final_x2.c1.c2.c1),
                },
            },
        };

        Fq12Target::connect(&mut builder, &final_x.0, &final_x2);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        dbg!(data.common.degree_bits());
        let _proof = data.prove(pw);
    }
}
