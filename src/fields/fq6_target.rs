use itertools::Itertools;
use plonky2::{
    field::extension::Extendable, hash::hash_types::RichField, iop::target::BoolTarget,
    plonk::circuit_builder::CircuitBuilder,
};
use plonky2_ecdsa::gadgets::nonnative::CircuitBuilderNonNative;

use super::{fq2_target::Fq2Target, fq_target::FqTarget};

#[derive(Debug, Clone)]
pub struct Fq6Target<F: RichField + Extendable<D>, const D: usize> {
    pub coeffs: [FqTarget<F, D>; 6],
}

impl<F: RichField + Extendable<D>, const D: usize> Fq6Target<F, D> {
    pub fn empty(builder: &mut CircuitBuilder<F, D>) -> Self {
        let coeffs = [(); 6]
            .iter()
            .map(|_| FqTarget::empty(builder))
            .collect_vec()
            .try_into()
            .unwrap();
        Fq6Target { coeffs }
    }

    pub fn new(coeffs: Vec<FqTarget<F, D>>) -> Self {
        Fq6Target {
            coeffs: coeffs.try_into().unwrap(),
        }
    }

    pub fn connect(builder: &mut CircuitBuilder<F, D>, lhs: &Self, rhs: &Self) {
        for i in 0..6 {
            builder.connect_nonnative(&lhs.coeffs[i].target, &rhs.coeffs[i].target);
        }
    }

    pub fn select(
        builder: &mut CircuitBuilder<F, D>,
        a: &Self,
        b: &Self,
        flag: &BoolTarget,
    ) -> Self {
        let selected = a
            .coeffs
            .iter()
            .zip(b.coeffs.iter())
            .map(|(a, b)| FqTarget::select(builder, a, b, flag))
            .collect_vec();

        Self {
            coeffs: selected.try_into().unwrap(),
        }
    }

    // pub fn constant(builder: &mut CircuitBuilder<F, D>, c: Fq6) -> Self {
    //     let coeffs = [c.c0, c.c1]
    //         .iter()
    //         .map(|x| FqTarget::constant(builder, x.clone()))
    //         .collect_vec()
    //         .try_into()
    //         .unwrap();
    //     Self { coeffs }
    // }

    // pub fn constant(builder: &mut CircuitBuilder<F, D>, c: Fq6) -> Self {
    //     let c: MyFq12 = c.into();
    //     let coeffs = c
    //         .coeffs
    //         .iter()
    //         .map(|x| FqTarget::constant(builder, x.clone()))
    //         .collect_vec()
    //         .try_into()
    //         .unwrap();
    //     Self { coeffs }
    // }

    pub fn add(&self, builder: &mut CircuitBuilder<F, D>, rhs: &Self) -> Self {
        let coeffs = self
            .coeffs
            .iter()
            .enumerate()
            .map(|(i, x)| x.add(builder, &rhs.coeffs[i]))
            .collect_vec()
            .try_into()
            .unwrap();
        Fq6Target { coeffs }
    }

    pub fn neg(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
        let coeffs = self
            .coeffs
            .iter()
            .map(|x| x.neg(builder))
            .collect_vec()
            .try_into()
            .unwrap();
        Fq6Target { coeffs }
    }

    pub fn sub(&self, builder: &mut CircuitBuilder<F, D>, rhs: &Self) -> Self {
        let coeffs = self
            .coeffs
            .iter()
            .enumerate()
            .map(|(i, x)| x.sub(builder, &rhs.coeffs[i]))
            .collect_vec()
            .try_into()
            .unwrap();
        Fq6Target { coeffs }
    }

    // Using https://eprint.iacr.org/2022/367.pdf#section.4 - should be reviewed again!!!
    pub fn mul(&self, builder: &mut CircuitBuilder<F, D>, rhs: &Self) -> Self {
        let a_c0 = Fq2Target::new(self.coeffs[..2].to_vec());
        let a_c1 = Fq2Target::new(self.coeffs[2..4].to_vec());
        let a_c2 = Fq2Target::new(self.coeffs[4..6].to_vec());

        let b_c0 = Fq2Target::new(rhs.coeffs[..2].to_vec());
        let b_c1 = Fq2Target::new(rhs.coeffs[2..4].to_vec());
        let b_c2 = Fq2Target::new(rhs.coeffs[4..6].to_vec());

        let b10_p_b11 = b_c1.coeffs[0].add(builder, &b_c1.coeffs[1]);
        let b10_m_b11 = b_c1.coeffs[0].sub(builder, &b_c1.coeffs[1]);
        let b20_p_b21 = b_c2.coeffs[0].add(builder, &b_c2.coeffs[1]);
        let b20_m_b21 = b_c2.coeffs[0].sub(builder, &b_c2.coeffs[1]);

        self.clone()
    }

    // pub fn mul(&self, builder: &mut CircuitBuilder<F, D>, rhs: &Self) -> Self {
    //     let a0 = self.coeffs[0].clone();
    //     let a1 = self.coeffs[1].clone();

    //     let b0 = rhs.coeffs[0].clone();
    //     let b1 = rhs.coeffs[1].clone();

    //     let a0_b0 = a0.mul(builder, &b0);
    //     let a1_b1 = a1.mul(builder, &b1);

    //     let c0 = a0_b0.sub(builder, &a1_b1);

    //     let a0_b1 = a0.mul(builder, &b1);
    //     let a1_b0 = a1.mul(builder, &b0);

    //     let c1 = a0_b1.add(builder, &a1_b0);

    //     Fq2Target { coeffs: [c0, c1] }

    //     // Fp2 {
    //     //     c0: Fp::sum_of_products([self.c0, -self.c1], [rhs.c0, rhs.c1]),
    //     //     c1: Fp::sum_of_products([self.c0, self.c1], [rhs.c1, rhs.c0]),
    //     // }
    // }

    pub fn mul_by_01(
        self,
        builder: &mut CircuitBuilder<F, D>,
        c0: &Fq2Target<F, D>,
        c1: &Fq2Target<F, D>,
    ) -> Self {
        let fq6_c0 = Fq2Target::new(self.coeffs[..2].to_vec());
        let fq6_c1 = Fq2Target::new(self.coeffs[2..4].to_vec());
        let fq6_c2 = Fq2Target::new(self.coeffs[4..6].to_vec());

        let a_a = fq6_c0.mul(builder, c0);
        let b_b = fq6_c1.mul(builder, c1);

        let t1 = fq6_c2.mul(builder, c1);
        let t1 = t1.mul_by_nonresidue(builder);
        let t1 = t1.add(builder, &a_a);

        let c0_add_c1 = c0.add(builder, c1);
        let fq6_c0_add_fq6_c1 = fq6_c0.add(builder, &fq6_c1);
        let t2 = c0_add_c1.mul(builder, &fq6_c0_add_fq6_c1);
        let t2 = t2.sub(builder, &a_a);
        let t2 = t2.sub(builder, &b_b);

        let t3 = fq6_c2.mul(builder, c0);
        let t3 = t3.add(builder, &b_b);

        Self::new(vec![
            t1.coeffs[0].clone(),
            t1.coeffs[1].clone(),
            t2.coeffs[0].clone(),
            t2.coeffs[1].clone(),
            t3.coeffs[0].clone(),
            t3.coeffs[1].clone(),
        ])
    }

    pub fn mul_by_1(self, builder: &mut CircuitBuilder<F, D>, c1: &Fq2Target<F, D>) -> Self {
        let fq6_c0 = Fq2Target::new(self.coeffs[..2].to_vec());
        let fq6_c1 = Fq2Target::new(self.coeffs[2..4].to_vec());
        let fq6_c2 = Fq2Target::new(self.coeffs[4..6].to_vec());

        let c0 = fq6_c2.mul(builder, c1);
        let c0 = c0.mul_by_nonresidue(builder);
        let c1 = fq6_c0.mul(builder, c1);
        let c2 = fq6_c1.mul(builder, &c1);

        Self::new(vec![
            c0.coeffs[0].clone(),
            c0.coeffs[1].clone(),
            c1.coeffs[0].clone(),
            c1.coeffs[1].clone(),
            c2.coeffs[0].clone(),
            c2.coeffs[1].clone(),
        ])
    }

    /// Multiply by quadratic nonresidue v.
    pub fn mul_by_nonresidue(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
        // Given a + bv + cv^2, this produces
        //     av + bv^2 + cv^3
        // but because v^3 = u + 1, we have
        //     c(u + 1) + av + v^2

        let fq6_c0 = Fq2Target::new(self.coeffs[..2].to_vec());
        let fq6_c1 = Fq2Target::new(self.coeffs[2..4].to_vec());
        let fq6_c2 = Fq2Target::new(self.coeffs[4..6].to_vec());

        let c0 = fq6_c2.mul_by_nonresidue(builder);
        let c1 = fq6_c0;
        let c2 = fq6_c1;

        Self::new(vec![
            c0.coeffs[0].clone(),
            c0.coeffs[1].clone(),
            c1.coeffs[0].clone(),
            c1.coeffs[1].clone(),
            c2.coeffs[0].clone(),
            c2.coeffs[1].clone(),
        ])
    }

    pub fn conditional_mul(
        &self,
        builder: &mut CircuitBuilder<F, D>,
        x: &Self,
        flag: &BoolTarget,
    ) -> Self {
        let muled = self.mul(builder, x);
        Self::select(builder, &muled, &self, flag)
    }

    // pub fn div(&self, builder: &mut CircuitBuilder<F, D>, other: &Self) -> Self {
    //     let inv = other.inv(builder);
    //     self.mul(builder, &inv)
    // }

    // pub fn inv(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
    //     let inv = Self::empty(builder);
    //     builder.add_simple_generator(Fq12InverseGenerator::<F, D> {
    //         x: self.clone(),
    //         inv: inv.clone(),
    //     });
    //     let one = Self::constant(builder, Fq12::ONE);
    //     let x_mul_inv = self.mul(builder, &inv);
    //     Self::connect(builder, &x_mul_inv, &one);
    //     inv
    // }

    // pub fn confugate(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
    //     let mut coeffs = self.coeffs.clone();
    //     coeffs[1] = coeffs[1].neg(builder);
    //     coeffs[3] = coeffs[3].neg(builder);
    //     coeffs[5] = coeffs[5].neg(builder);
    //     coeffs[7] = coeffs[7].neg(builder);
    //     coeffs[9] = coeffs[9].neg(builder);
    //     coeffs[11] = coeffs[11].neg(builder);
    //     Self { coeffs }
    // }
}
