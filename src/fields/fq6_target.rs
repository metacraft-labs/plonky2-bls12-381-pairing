use ark_bls12_381::{Fq, Fq6};
use ark_ff::Field;
use itertools::Itertools;
use num_bigint::BigUint;
use plonky2::{
    field::extension::Extendable,
    hash::hash_types::RichField,
    iop::{
        generator::{GeneratedValues, SimpleGenerator},
        target::{BoolTarget, Target},
        witness::{PartitionWitness, WitnessWrite},
    },
    plonk::circuit_builder::CircuitBuilder,
    util::serialization::{Buffer, IoError},
};
use plonky2_ecdsa::gadgets::{
    biguint::{GeneratedValuesBigUint, WitnessBigUint},
    nonnative::CircuitBuilderNonNative,
};

use super::{fq2_target::Fq2Target, fq_target::FqTarget, helpers::{mul_by_01, mul_by_1}};

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

    pub fn mul(&self, builder: &mut CircuitBuilder<F, D>, rhs: &Self) -> Self {
        let a = self;
        let b = rhs;
        let mut a0b0_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(11);
        let mut a0b1_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(11);
        let mut a1b0_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(11);
        let mut a1b1_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(11);
        for i in 0..6 {
            for j in 0..6 {
                let coeff00 = a.coeffs[i].mul(builder, &b.coeffs[j]);
                let coeff01 = a.coeffs[i].mul(builder, &b.coeffs[j + 6]);
                let coeff10 = a.coeffs[i + 6].mul(builder, &b.coeffs[j]);
                let coeff11 = a.coeffs[i + 6].mul(builder, &b.coeffs[j + 6]);
                if i + j < a0b0_coeffs.len() {
                    a0b0_coeffs[i + j] = a0b0_coeffs[i + j].add(builder, &coeff00);
                    a0b1_coeffs[i + j] = a0b1_coeffs[i + j].add(builder, &coeff01);
                    a1b0_coeffs[i + j] = a1b0_coeffs[i + j].add(builder, &coeff10);
                    a1b1_coeffs[i + j] = a1b1_coeffs[i + j].add(builder, &coeff11);
                } else {
                    a0b0_coeffs.push(coeff00);
                    a0b1_coeffs.push(coeff01);
                    a1b0_coeffs.push(coeff10);
                    a1b1_coeffs.push(coeff11);
                }
            }
        }

        let mut a0b0_minus_a1b1: Vec<FqTarget<F, D>> = Vec::with_capacity(11);
        let mut a0b1_plus_a1b0: Vec<FqTarget<F, D>> = Vec::with_capacity(11);
        for i in 0..11 {
            let a0b0_minus_a1b1_entry = a0b0_coeffs[i].sub(builder, &a1b1_coeffs[i]);
            let a0b1_plus_a1b0_entry = a0b1_coeffs[i].add(builder, &a1b0_coeffs[i]);
            a0b0_minus_a1b1.push(a0b0_minus_a1b1_entry);
            a0b1_plus_a1b0.push(a0b1_plus_a1b0_entry);
        }

        let const_one = FqTarget::constant(builder, Fq::from(1));
        let mut out_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(12);
        for i in 0..6 {
            if i < 5 {
                // let coeff: Fq = a0b0_minus_a1b1[i] + Fq::from(1) * a0b0_minus_a1b1[i + 6]
                //     - a0b1_plus_a1b0[i + 6];
                let term0 = a0b0_minus_a1b1[i].clone();
                let term1 = a0b0_minus_a1b1[i + 6].mul(builder, &const_one);
                let term2 = a0b1_plus_a1b0[i + 6].neg(builder);
                let term0_plus_term1 = term0.add(builder, &term1);
                let coeff = term0_plus_term1.add(builder, &term2);
                out_coeffs.push(coeff);
            } else {
                out_coeffs.push(a0b0_minus_a1b1[i].clone());
            }
        }
        for i in 0..6 {
            if i < 5 {
                // let coeff: Fq = a0b1_plus_a1b0[i]
                //     + a0b0_minus_a1b1[i + 6]
                //     + Fq::from(1) * a0b1_plus_a1b0[i + 6];
                let term0 = a0b1_plus_a1b0[i].clone();
                let term1 = a0b0_minus_a1b1[i + 6].clone();
                let term2 = a0b1_plus_a1b0[i + 6].mul(builder, &const_one);
                let term0_plus_term1 = term0.add(builder, &term1);
                let coeff = term0_plus_term1.add(builder, &term2);
                out_coeffs.push(coeff);
            } else {
                out_coeffs.push(a0b1_plus_a1b0[i].clone());
            }
        }
        Self {
            coeffs: out_coeffs.try_into().unwrap(),
        }
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

    pub fn conditional_mul(
        &self,
        builder: &mut CircuitBuilder<F, D>,
        x: &Self,
        flag: &BoolTarget,
    ) -> Self {
        let muled = self.mul(builder, x);
        Self::select(builder, &muled, &self, flag)
    }
}