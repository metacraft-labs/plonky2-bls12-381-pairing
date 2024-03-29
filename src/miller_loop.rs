use ark_bls12_381::Fq12;
use ark_ff::vec::IntoIter;
use ark_ff::BitIteratorBE;
use ark_std::cfg_chunks_mut;
use num::One;
use plonky2::{
    field::extension::Extendable, hash::hash_types::RichField,
    plonk::circuit_builder::CircuitBuilder,
};

use crate::{
    curves::{
        g1::{G1AffineTarget, G1PreparedTarget},
        g2::{EllCoeffTarget, G2PreparedTarget},
    },
    fields::{fq12_target::Fq12Target, fq2_target::Fq2Target},
    utils::constants::{BLS_X, BLS_X_IS_NEGATIVE},
};

pub fn multi_miller_loop<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    a: impl IntoIterator<Item = impl Into<G1PreparedTarget<F, D>>>,
    b: impl IntoIterator<Item = impl Into<G2PreparedTarget<F, D>>>,
) -> Fq12Target<F, D> {
    use itertools::Itertools;

    println!("[1]");

    let mut pairs = a
        .into_iter()
        .zip_eq(b)
        .filter_map(|(p, q)| {
            let (p, q) = (p.into(), q.into());
            match !p.0.is_zero() && !q.is_zero() {
                true => Some((p, q.ell_coeffs.into_iter())),
                false => None,
            }
        })
        .collect::<Vec<_>>();

    println!("[2]");
    let mut pairs_f_storage: Vec<Fq12Target<F, D>> = Vec::new();

    for pairs in cfg_chunks_mut!(pairs, 4) {
        let mut f_to = Fq12Target::constant(builder, Fq12::one());
        for i in BitIteratorBE::without_leading_zeros([BLS_X]).skip(1) {
            f_to = f_to.mul(builder, &f_to); // square in place
            for (p, coeffs) in pairs.iter_mut() {
                f_to = ell_target(builder, &f_to, coeffs.next().unwrap(), p.0.clone());
            }
            if i {
                for (p, coeffs) in pairs.iter_mut() {
                    f_to = ell_target(builder, &f_to, coeffs.next().unwrap(), p.0.clone());
                }
            }
        }
        pairs_f_storage.push(f_to)
    }

    let f = Fq12Target::multiply_elements(builder, pairs_f_storage.into_iter()).unwrap();

    // let f: Fq12Target<F, D> = cfg_chunks_mut!(pairs, 4)
    //     .map(|pairs| {
    //         let mut f = Fq12Target::constant(builder, Fq12::one());
    //         // for i in BitIteratorBE::without_leading_zeros([BLS_X]).skip(1) {
    //         //     f = f.mul(builder, &f); // square in place
    //         //     for (p, coeffs) in pairs.iter_mut() {
    //         //         f = ell_target(builder, &f, coeffs.next().unwrap(), p.0.clone());
    //         //     }
    //         //     if i {
    //         //         for (p, coeffs) in pairs.iter_mut() {
    //         //             f = ell_target(builder, &f, coeffs.next().unwrap(), p.0.clone());
    //         //         }
    //         //     }
    //         // }
    //         f
    //     })
    //     .product();

    println!("[3]");

    // let mut f = Fq12Target::constant(builder, Fq12::one());
    if BLS_X_IS_NEGATIVE {
        f.conjugate(builder);
    }

    println!("[4]");

    f
}

fn ell_target<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    f: &Fq12Target<F, D>,
    g2_coeffs: EllCoeffTarget<F, D>,
    p: G1AffineTarget<F, D>,
) -> Fq12Target<F, D> {
    let c0 = g2_coeffs.0;
    let c1 = g2_coeffs.1;
    let c2 = g2_coeffs.2;
    let (px, py) = p.xy().unwrap();

    let c2 = c2.mul_assign_by_fp(builder, py.clone());
    let c1 = c1.mul_assign_by_fp(builder, px.clone());
    let f = f.mul_by_014(builder, &c0, &c1, &c2);

    f
}

pub fn getter_of_prepared_pairs_target<F: RichField + Extendable<D>, const D: usize>(
    a: impl IntoIterator<Item = impl Into<G1PreparedTarget<F, D>>>,
    b: impl IntoIterator<Item = impl Into<G2PreparedTarget<F, D>>>,
) -> Vec<(
    G1PreparedTarget<F, D>,
    IntoIter<(Fq2Target<F, D>, Fq2Target<F, D>, Fq2Target<F, D>)>,
)> {
    use itertools::Itertools;

    let pairs = a
        .into_iter()
        .zip_eq(b)
        .filter_map(|(p, q)| {
            let (p, q) = (p.into(), q.into());
            match !p.0.is_zero() && !q.is_zero() {
                true => Some((p, q.ell_coeffs.into_iter())),
                false => None,
            }
        })
        .collect::<Vec<_>>();

    pairs
}

pub fn connect_pairs<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    lhs: Vec<(
        G1PreparedTarget<F, D>,
        IntoIter<(Fq2Target<F, D>, Fq2Target<F, D>, Fq2Target<F, D>)>,
    )>,
    rhs: Vec<(
        G1PreparedTarget<F, D>,
        IntoIter<(Fq2Target<F, D>, Fq2Target<F, D>, Fq2Target<F, D>)>,
    )>,
) {
    G1AffineTarget::connect(builder, &lhs[0].0 .0, &rhs[0].0 .0);
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Fq12, Fq2, G1Affine, G2Affine};
    use ark_ec::pairing::Pairing;
    use ark_ff::UniformRand;
    use plonky2::{
        field::goldilocks_field::GoldilocksField,
        iop::witness::PartialWitness,
        plonk::{
            circuit_builder::CircuitBuilder, circuit_data::CircuitConfig,
            config::PoseidonGoldilocksConfig,
        },
    };

    use crate::{
        curves::{
            g1::{G1AffineTarget, G1PreparedTarget},
            g2::{G2AffineTarget, G2PreparedTarget},
        },
        fields::{fq12_target::Fq12Target, fq2_target::Fq2Target},
        miller_loop::multi_miller_loop,
        native::miller_loop::{ell, getter_of_prepared_pairs, G1Prepared, G2Prepared},
    };

    use super::{ell_target, getter_of_prepared_pairs_target};

    type F = GoldilocksField;
    type C = PoseidonGoldilocksConfig;
    const D: usize = 2;

    #[test]
    fn test_miller_loop_circuit() {
        let config = CircuitConfig::pairing_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let rng = &mut rand::thread_rng();
        let p = G1Affine::rand(rng);
        let q = G2Affine::rand(rng);
        let r_expected = ark_bls12_381::Bls12_381::miller_loop(p, q).0;

        let p_prepared_t = [G1PreparedTarget(G1AffineTarget::constant(&mut builder, p))];
        let q_t = G2AffineTarget::constant(&mut builder, q);
        let q_prepared_t = [G2PreparedTarget::from(&mut builder, q_t)];

        let r_t = multi_miller_loop(&mut builder, p_prepared_t, q_prepared_t);

        let r_expected_t = Fq12Target::constant(&mut builder, r_expected);

        println!("[5]");

        Fq12Target::connect(&mut builder, &r_t, &r_expected_t);

        let pw = PartialWitness::<F>::new();
        let data = builder.build::<C>();
        dbg!(data.common.degree_bits());
        println!("[6]");
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_ell_target() {
        let config = CircuitConfig::pairing_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let rng = &mut rand::thread_rng();
        let p = G1Affine::rand(rng);
        let g2_coeff_c0: Fq2 = Fq2::rand(rng);
        let g2_coeff_c1: Fq2 = Fq2::rand(rng);
        let g2_coeff_c2: Fq2 = Fq2::rand(rng);
        let f: Fq12 = Fq12::rand(rng);

        let g2_coeff_c0_t = Fq2Target::constant(&mut builder, g2_coeff_c0);
        let g2_coeff_c1_t = Fq2Target::constant(&mut builder, g2_coeff_c1);
        let g2_coeff_c2_t = Fq2Target::constant(&mut builder, g2_coeff_c2);
        let p_t = G1AffineTarget::constant(&mut builder, p);
        let f_t = Fq12Target::constant(&mut builder, f);

        let mut f = f;
        ell(&mut f, (g2_coeff_c0, g2_coeff_c1, g2_coeff_c2), p);
        let f_t = ell_target(
            &mut builder,
            &f_t,
            (g2_coeff_c0_t, g2_coeff_c1_t, g2_coeff_c2_t),
            p_t,
        );

        let f_expected = Fq12Target::constant(&mut builder, f);
        Fq12Target::connect(&mut builder, &f_t, &f_expected);

        let pw = PartialWitness::<F>::new();
        let data = builder.build::<C>();
        dbg!(data.common.degree_bits());
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_getter_of_prepared_pairs() {
        let config = CircuitConfig::pairing_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let rng = &mut rand::thread_rng();
        let p = G1Affine::rand(rng);
        let q = G2Affine::rand(rng);
        let native_pairs = getter_of_prepared_pairs([G1Prepared(p)], [G2Prepared::from(q)]);
        let native_pair_g1 =
            G1PreparedTarget(G1AffineTarget::constant(&mut builder, native_pairs[0].0 .0));
        let native_pair_g2_coeffs = native_pairs[0].1.clone().next().unwrap();

        let p_prepared_t = [G1PreparedTarget(G1AffineTarget::constant(&mut builder, p))];
        let q_t = G2AffineTarget::constant(&mut builder, q);
        let q_prepared_t = [G2PreparedTarget::from(&mut builder, q_t)];

        let target_pairs = getter_of_prepared_pairs_target(p_prepared_t, q_prepared_t);

        // connect_pairs(&mut builder, &native_pairs, &target_pairs);

        let pw = PartialWitness::<F>::new();
        let data = builder.build::<C>();
        dbg!(data.common.degree_bits());
        let _proof = data.prove(pw);
    }
}
