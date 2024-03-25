use crate::utils::{constants::BLS_X, helpers::MyFq12};
use ark_bls12_381::Fq12;
use ark_ff::Field;

use super::fin_exp_helpers::{conjugate_fp12, experimental_pow, frobenius_map_native};

// out = in^{ (q^6 - 1)*(q^2 + 1) }
pub fn easy_part_native<'v>(a: MyFq12) -> MyFq12 {
    let f1 = conjugate_fp12(a);
    let f2 = {
        let f1_fp12: Fq12 = f1.into();
        let a_fp12: Fq12 = a.into();
        let divided = f1_fp12 / a_fp12;
        divided.into()
    };
    let f3 = frobenius_map_native(f2, 2);
    let f = f3 * f2;
    f
}

pub fn hard_part_native(r: Fq12) -> Fq12 {
    let mut y0 = r.square();
    let mut y1: Fq12 = experimental_pow(r.into(), vec![BLS_X]).into();
    y1 = *y1.conjugate_in_place();
    let mut y2 = r;
    y2 = *y2.conjugate_in_place();

    y1 = y1 * y2;
    y2 = experimental_pow(y1.into(), vec![BLS_X].into()).into();
    y2 = *y2.conjugate_in_place();
    y1 = *y1.conjugate_in_place();
    y1 = y1 * y2;
    y2 = experimental_pow(y1.into(), vec![BLS_X].into()).into();
    y2 = *y2.conjugate_in_place();
    y1 = frobenius_map_native(y1.into(), 1).into();
    y1 = y1 * y2;
    let r = r * y0;
    y0 = experimental_pow(y1.into(), vec![BLS_X].into()).into();
    y0 = *y0.conjugate_in_place();
    y2 = experimental_pow(y0.into(), vec![BLS_X].into()).into();
    y2 = *y2.conjugate_in_place();
    y0 = y1;
    y0 = frobenius_map_native(y0.into(), 2).into();
    y1 = *y1.conjugate_in_place();
    y1 *= y2;
    y1 *= y0;
    let r = r * y1;
    r
}

pub fn final_exponentiation(a: MyFq12) -> Fq12 {
    let f0 = easy_part_native(a);
    let f = hard_part_native(f0.into());
    f
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fq12;
    use ark_ec::pairing::{MillerLoopOutput, Pairing};
    use ark_ff::UniformRand;

    use crate::native::final_exponentiation::final_exponentiation;

    #[test]
    fn test_native_final_exponentiation() {
        let rng = &mut rand::thread_rng();
        let x = Fq12::rand(rng);
        let ark_final_exp_result =
            ark_bls12_381::Bls12_381::final_exponentiation(MillerLoopOutput(x))
                .unwrap()
                .0;
        let final_exp_native_result = final_exponentiation(x.into());

        assert_eq!(ark_final_exp_result, final_exp_native_result);
    }
}
