use ark_bls12_381::Fq12;

use super::{
    final_exponentiation::final_exponentiation,
    miller_loop::{multi_miller_loop_native, G1Prepared, G2Prepared},
};

pub fn native_pairing(
    a: impl IntoIterator<Item = impl Into<G1Prepared>>,
    b: impl IntoIterator<Item = impl Into<G2Prepared>>,
) -> Fq12 {
    final_exponentiation(multi_miller_loop_native(a, b).into())
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{G1Affine, G2Affine};
    use ark_ec::pairing::Pairing;
    use ark_ff::UniformRand;

    use crate::native::{
        miller_loop::{G1Prepared, G2Prepared},
        pairing::native_pairing,
    };

    #[test]
    fn test_native_pairing() {
        let rng = &mut rand::thread_rng();
        let p0 = G1Affine::rand(rng);
        let q0 = G2Affine::rand(rng);
        let ark_pairing = ark_bls12_381::Bls12_381::pairing(p0, q0).0;
        let native_pairing = native_pairing([G1Prepared(p0)], [G2Prepared::from(q0)]);
        assert_eq!(ark_pairing, native_pairing);
    }
}
