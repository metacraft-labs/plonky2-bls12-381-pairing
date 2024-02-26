use plonky2::{
    field::extension::Extendable, hash::hash_types::RichField,
    plonk::circuit_builder::CircuitBuilder,
};
use plonky2_ecdsa::gadgets::nonnative::CircuitBuilderNonNative;

use super::{fq2_target_tree::Fq2Target, fq6_target_tree::Fq6Target};

#[derive(Clone, Debug)]
pub struct Fq12Target<F: RichField + Extendable<D>, const D: usize> {
    pub c0: Fq6Target<F, D>,
    pub c1: Fq6Target<F, D>,
}

impl<F: RichField + Extendable<D>, const D: usize> Fq12Target<F, D> {
    pub fn zero(builder: &mut CircuitBuilder<F, D>) -> Self {
        Self {
            c0: Fq6Target::zero(builder),
            c1: Fq6Target::zero(builder),
        }
    }

    pub fn one(builder: &mut CircuitBuilder<F, D>) -> Self {
        Self {
            c0: Fq6Target::zero(builder),
            c1: Fq6Target::one(builder),
        }
    }

    pub fn add(&self, builder: &mut CircuitBuilder<F, D>, rhs: Self) -> Self {
        Self {
            c0: self.c0.add(builder, rhs.c0),
            c1: self.c1.add(builder, rhs.c1),
        }
    }

    pub fn sub(&self, builder: &mut CircuitBuilder<F, D>, rhs: Self) -> Self {
        Self {
            c0: self.c0.sub(builder, rhs.c0),
            c1: self.c1.sub(builder, rhs.c1),
        }
    }

    pub fn neg(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
        Self {
            c0: self.c0.neg(builder),
            c1: self.c1.neg(builder),
        }
    }

    pub fn mul(&self, builder: &mut CircuitBuilder<F, D>, rhs: &Self) -> Self {
        let aa = self.c0.mul(builder, &rhs.c0);
        let bb = self.c1.mul(builder, &rhs.c1);
        let o = rhs.c0.add(builder, rhs.c1.clone());
        let c1 = self.c1.add(builder, self.c0.clone());
        let c1 = c1.mul(builder, &o);
        let c1 = c1.sub(builder, aa.clone());
        let c1 = c1.sub(builder, bb.clone());
        let c0 = bb.mul_by_nonresidue(builder);
        let c0 = c0.add(builder, aa);
        Self { c0, c1 }
    }

    pub fn square(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
        let ab = self.c0.mul(builder, &self.c1);
        let c0c1 = self.c0.add(builder, self.c1.clone());
        let c0 = self.c1.mul_by_nonresidue(builder);
        let c0 = c0.add(builder, self.c0.clone());
        let c0 = c0.mul(builder, &c0c1);
        let c0 = c0.sub(builder, ab.clone());

        let c1 = ab.add(builder, ab.clone());
        let tmp = ab.mul_by_nonresidue(builder);
        let c0 = c0.sub(builder, tmp);
        Self { c0, c1 }
    }

    pub fn mul_by_014(
        &self,
        builder: &mut CircuitBuilder<F, D>,
        c0: &Fq2Target<F, D>,
        c1: &Fq2Target<F, D>,
        c4: &Fq2Target<F, D>,
    ) -> Self {
        let aa = self.c0.mul_by_01(builder, c0, c1);
        let bb = self.c1.mul_by_1(builder, c4);
        let o = c1.add(builder, c4.clone());
        let c1 = self.c1.add(builder, self.c0.clone());
        let c1 = c1.mul_by_01(builder, c0, &o);
        let c1 = c1.sub(builder, aa.clone());
        let c1 = c1.sub(builder, bb.clone());
        let c0 = bb;
        let c0 = c0.mul_by_nonresidue(builder);
        let c0 = c0.add(builder, aa);

        Self { c0, c1 }
    }

    pub fn connect(builder: &mut CircuitBuilder<F, D>, lhs: &Self, rhs: &Self) {
        builder.connect_nonnative(&lhs.c0.c0.c0.target, &rhs.c0.c0.c0.target);
        builder.connect_nonnative(&lhs.c0.c0.c1.target, &rhs.c0.c0.c1.target);
        builder.connect_nonnative(&lhs.c0.c1.c0.target, &rhs.c0.c1.c0.target);
        builder.connect_nonnative(&lhs.c0.c1.c1.target, &rhs.c0.c1.c1.target);
        builder.connect_nonnative(&lhs.c0.c2.c0.target, &rhs.c0.c2.c0.target);
        builder.connect_nonnative(&lhs.c0.c2.c1.target, &rhs.c0.c2.c1.target);
        builder.connect_nonnative(&lhs.c1.c0.c0.target, &rhs.c1.c0.c0.target);
        builder.connect_nonnative(&lhs.c1.c0.c1.target, &rhs.c1.c0.c1.target);
        builder.connect_nonnative(&lhs.c1.c1.c0.target, &rhs.c1.c1.c0.target);
        builder.connect_nonnative(&lhs.c1.c1.c1.target, &rhs.c1.c1.c1.target);
        builder.connect_nonnative(&lhs.c1.c2.c0.target, &rhs.c1.c2.c0.target);
        builder.connect_nonnative(&lhs.c1.c2.c1.target, &rhs.c1.c2.c1.target);
    }
}

#[cfg(test)]
mod tests {
    use plonky2::{
        field::goldilocks_field::GoldilocksField,
        iop::witness::PartialWitness,
        plonk::{
            circuit_builder::CircuitBuilder, circuit_data::CircuitConfig,
            config::PoseidonGoldilocksConfig,
        },
    };

    use crate::{
        fields::{bls12_381base::Bls12_381Base, fq_target::FqTarget},
        fields_as_trees::{fq2_target_tree::Fq2Target, fq6_target_tree::Fq6Target},
    };

    use super::Fq12Target;

    type F = GoldilocksField;
    type C = PoseidonGoldilocksConfig;
    const D: usize = 2;

    #[test]
    fn test_fq12_arithmetic() {
        let config = CircuitConfig::pairing_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let a = Fq12Target {
            c0: Fq6Target {
                c0: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x47f9_cb98_b1b8_2d58,
                            0x5fe9_11eb_a3aa_1d9d,
                            0x96bf_1b5f_4dd8_1db3,
                            0x8100_d27c_c925_9f5b,
                            0xafa2_0b96_7464_0eab,
                            0x09bb_cea7_d8d9_497d,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x0303_cb98_b166_2daa,
                            0xd931_10aa_0a62_1d5a,
                            0xbfa9_820c_5be4_a468,
                            0x0ba3_643e_cb05_a348,
                            0xdc35_34bb_1f1c_25a6,
                            0x06c3_05bb_19c0_e1c1,
                        ]),
                    ),
                },
                c1: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x46f9_cb98_b162_d858,
                            0x0be9_109c_f7aa_1d57,
                            0xc791_bc55_fece_41d2,
                            0xf84c_5770_4e38_5ec2,
                            0xcb49_c1d9_c010_e60f,
                            0x0acd_b8e1_58bf_e3c8,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x8aef_cb98_b15f_8306,
                            0x3ea1_108f_e4f2_1d54,
                            0xcf79_f69f_a1b7_df3b,
                            0xe4f5_4aa1_d16b_1a3c,
                            0xba5e_4ef8_6105_a679,
                            0x0ed8_6c07_97be_e5cf,
                        ]),
                    ),
                },
                c2: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0xcee5_cb98_b15c_2db4,
                            0x7159_1082_d23a_1d51,
                            0xd762_30e9_44a1_7ca4,
                            0xd19e_3dd3_549d_d5b6,
                            0xa972_dc17_01fa_66e3,
                            0x12e3_1f2d_d6bd_e7d6,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0xad2a_cb98_b173_2d9d,
                            0x2cfd_10dd_0696_1d64,
                            0x0739_6b86_c6ef_24e8,
                            0xbd76_e2fd_b1bf_c820,
                            0x6afe_a7f6_de94_d0d5,
                            0x1099_4b0c_5744_c040,
                        ]),
                    ),
                },
            },
            c1: Fq6Target {
                c0: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x47f9_cb98_b1b8_2d58,
                            0x5fe9_11eb_a3aa_1d9d,
                            0x96bf_1b5f_4dd8_1db3,
                            0x8100_d27c_c925_9f5b,
                            0xafa2_0b96_7464_0eab,
                            0x09bb_cea7_d8d9_497d,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x0303_cb98_b166_2daa,
                            0xd931_10aa_0a62_1d5a,
                            0xbfa9_820c_5be4_a468,
                            0x0ba3_643e_cb05_a348,
                            0xdc35_34bb_1f1c_25a6,
                            0x06c3_05bb_19c0_e1c1,
                        ]),
                    ),
                },
                c1: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x46f9_cb98_b162_d858,
                            0x0be9_109c_f7aa_1d57,
                            0xc791_bc55_fece_41d2,
                            0xf84c_5770_4e38_5ec2,
                            0xcb49_c1d9_c010_e60f,
                            0x0acd_b8e1_58bf_e3c8,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x8aef_cb98_b15f_8306,
                            0x3ea1_108f_e4f2_1d54,
                            0xcf79_f69f_a1b7_df3b,
                            0xe4f5_4aa1_d16b_1a3c,
                            0xba5e_4ef8_6105_a679,
                            0x0ed8_6c07_97be_e5cf,
                        ]),
                    ),
                },
                c2: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0xcee5_cb98_b15c_2db4,
                            0x7159_1082_d23a_1d51,
                            0xd762_30e9_44a1_7ca4,
                            0xd19e_3dd3_549d_d5b6,
                            0xa972_dc17_01fa_66e3,
                            0x12e3_1f2d_d6bd_e7d6,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0xad2a_cb98_b173_2d9d,
                            0x2cfd_10dd_0696_1d64,
                            0x0739_6b86_c6ef_24e8,
                            0xbd76_e2fd_b1bf_c820,
                            0x6afe_a7f6_de94_d0d5,
                            0x1099_4b0c_5744_c040,
                        ]),
                    ),
                },
            },
        };

        let b = Fq12Target {
            c0: Fq6Target {
                c0: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x47f9_cb98_b1b8_2d58,
                            0x5fe9_11eb_a3aa_1d9d,
                            0x96bf_1b5f_4dd8_1db3,
                            0x8100_d272_c925_9f5b,
                            0xafa2_0b96_7464_0eab,
                            0x09bb_cea7_d8d9_497d,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x0303_cb98_b166_2daa,
                            0xd931_10aa_0a62_1d5a,
                            0xbfa9_820c_5be4_a468,
                            0x0ba3_643e_cb05_a348,
                            0xdc35_34bb_1f1c_25a6,
                            0x06c3_05bb_19c0_e1c1,
                        ]),
                    ),
                },
                c1: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x46f9_cb98_b162_d858,
                            0x0be9_109c_f7aa_1d57,
                            0xc791_bc55_fece_41d2,
                            0xf84c_5770_4e38_5ec2,
                            0xcb49_c1d9_c010_e60f,
                            0x0acd_b8e1_58bf_e348,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x8aef_cb98_b15f_8306,
                            0x3ea1_108f_e4f2_1d54,
                            0xcf79_f69f_a1b7_df3b,
                            0xe4f5_4aa1_d16b_1a3c,
                            0xba5e_4ef8_6105_a679,
                            0x0ed8_6c07_97be_e5cf,
                        ]),
                    ),
                },
                c2: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0xcee5_cb98_b15c_2db4,
                            0x7159_1082_d23a_1d51,
                            0xd762_30e9_44a1_7ca4,
                            0xd19e_3dd3_549d_d5b6,
                            0xa972_dc17_01fa_66e3,
                            0x12e3_1f2d_d6bd_e7d6,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0xad2a_cb98_b173_2d9d,
                            0x2cfd_10dd_0696_1d64,
                            0x0739_6b86_c6ef_24e8,
                            0xbd76_e2fd_b1bf_c820,
                            0x6afe_a7f6_de94_d0d5,
                            0x1099_4b0c_5744_c040,
                        ]),
                    ),
                },
            },
            c1: Fq6Target {
                c0: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x47f9_cb98_b1b8_2d58,
                            0x5fe9_11eb_a3aa_1d9d,
                            0x96bf_1b5f_4dd2_1db3,
                            0x8100_d27c_c925_9f5b,
                            0xafa2_0b96_7464_0eab,
                            0x09bb_cea7_d8d9_497d,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x0303_cb98_b166_2daa,
                            0xd931_10aa_0a62_1d5a,
                            0xbfa9_820c_5be4_a468,
                            0x0ba3_643e_cb05_a348,
                            0xdc35_34bb_1f1c_25a6,
                            0x06c3_05bb_19c0_e1c1,
                        ]),
                    ),
                },
                c1: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x46f9_cb98_b162_d858,
                            0x0be9_109c_f7aa_1d57,
                            0xc791_bc55_fece_41d2,
                            0xf84c_5770_4e38_5ec2,
                            0xcb49_c1d9_c010_e60f,
                            0x0acd_b8e1_58bf_e3c8,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x8aef_cb98_b15f_8306,
                            0x3ea1_108f_e4f2_1d54,
                            0xcf79_f69f_a117_df3b,
                            0xe4f5_4aa1_d16b_1a3c,
                            0xba5e_4ef8_6105_a679,
                            0x0ed8_6c07_97be_e5cf,
                        ]),
                    ),
                },
                c2: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0xcee5_cb98_b15c_2db4,
                            0x7159_1082_d23a_1d51,
                            0xd762_30e9_44a1_7ca4,
                            0xd19e_3dd3_549d_d5b6,
                            0xa972_dc17_01fa_66e3,
                            0x12e3_1f2d_d6bd_e7d6,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0xad2a_cb98_b173_2d9d,
                            0x2cfd_10dd_0696_1d64,
                            0x0739_6b86_c6ef_24e8,
                            0xbd76_e2fd_b1bf_c820,
                            0x6afe_a7f6_de94_d0d5,
                            0x1099_4b0c_5744_c040,
                        ]),
                    ),
                },
            },
        };

        let c = Fq12Target {
            c0: Fq6Target {
                c0: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x47f9_cb98_71b8_2d58,
                            0x5fe9_11eb_a3aa_1d9d,
                            0x96bf_1b5f_4dd8_1db3,
                            0x8100_d27c_c925_9f5b,
                            0xafa2_0b96_7464_0eab,
                            0x09bb_cea7_d8d9_497d,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x0303_cb98_b166_2daa,
                            0xd931_10aa_0a62_1d5a,
                            0xbfa9_820c_5be4_a468,
                            0x0ba3_643e_cb05_a348,
                            0xdc35_34bb_1f1c_25a6,
                            0x06c3_05bb_19c0_e1c1,
                        ]),
                    ),
                },
                c1: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x46f9_cb98_b162_d858,
                            0x0be9_109c_f7aa_1d57,
                            0x7791_bc55_fece_41d2,
                            0xf84c_5770_4e38_5ec2,
                            0xcb49_c1d9_c010_e60f,
                            0x0acd_b8e1_58bf_e3c8,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x8aef_cb98_b15f_8306,
                            0x3ea1_108f_e4f2_1d54,
                            0xcf79_f69f_a1b7_df3b,
                            0xe4f5_4aa1_d16b_133c,
                            0xba5e_4ef8_6105_a679,
                            0x0ed8_6c07_97be_e5cf,
                        ]),
                    ),
                },
                c2: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0xcee5_cb98_b15c_2db4,
                            0x7159_1082_d23a_1d51,
                            0xd762_40e9_44a1_7ca4,
                            0xd19e_3dd3_549d_d5b6,
                            0xa972_dc17_01fa_66e3,
                            0x12e3_1f2d_d6bd_e7d6,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0xad2a_cb98_b173_2d9d,
                            0x2cfd_10dd_0696_1d64,
                            0x0739_6b86_c6ef_24e8,
                            0xbd76_e2fd_b1bf_c820,
                            0x6afe_a7f6_de94_d0d5,
                            0x1099_4b0c_1744_c040,
                        ]),
                    ),
                },
            },
            c1: Fq6Target {
                c0: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x47f9_cb98_b1b8_2d58,
                            0x5fe9_11eb_a3aa_1d9d,
                            0x96bf_1b5f_4dd8_1db3,
                            0x8100_d27c_c925_9f5b,
                            0xafa2_0b96_7464_0eab,
                            0x09bb_cea7_d8d9_497d,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x0303_cb98_b166_2daa,
                            0xd931_10aa_0a62_1d5a,
                            0xbfa9_820c_5be4_a468,
                            0x0ba3_643e_cb05_a348,
                            0xdc35_34bb_1f1c_25a6,
                            0x06c3_05bb_19c0_e1c1,
                        ]),
                    ),
                },
                c1: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x46f9_cb98_b162_d858,
                            0x0be9_109c_f7aa_1d57,
                            0xc791_bc55_fece_41d2,
                            0xf84c_5770_4e38_5ec2,
                            0xcb49_c1d3_c010_e60f,
                            0x0acd_b8e1_58bf_e3c8,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0x8aef_cb98_b15f_8306,
                            0x3ea1_108f_e4f2_1d54,
                            0xcf79_f69f_a1b7_df3b,
                            0xe4f5_4aa1_d16b_1a3c,
                            0xba5e_4ef8_6105_a679,
                            0x0ed8_6c07_97be_e5cf,
                        ]),
                    ),
                },
                c2: Fq2Target {
                    c0: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0xcee5_cb98_b15c_2db4,
                            0x7159_1082_d23a_1d51,
                            0xd762_30e9_44a1_7ca4,
                            0xd19e_3dd3_549d_d5b6,
                            0xa972_dc17_01fa_66e3,
                            0x12e3_1f2d_d6bd_e7d6,
                        ]),
                    ),
                    c1: FqTarget::fp_constant(
                        &mut builder,
                        Bls12_381Base([
                            0xad2a_cb98_b173_2d9d,
                            0x2cfd_10dd_0696_1d64,
                            0x0739_6b86_c6ef_24e8,
                            0xbd76_e2fd_b1bf_c820,
                            0x6afe_a7f6_de94_d0d5,
                            0x1099_4b0c_5744_1040,
                        ]),
                    ),
                },
            },
        };

        let a = a.square(&mut builder);
        let a = a.square(&mut builder);
        let a = a.add(&mut builder, c.clone());

        let b = b.square(&mut builder);
        let b = b.square(&mut builder);
        let b = b.add(&mut builder, a.clone());

        let c = c.square(&mut builder);
        let c = c.square(&mut builder);
        let c = c.add(&mut builder, b.clone());

        let a_squared = a.square(&mut builder);
        let a_mul_a = a.mul(&mut builder, &a);

        let b_squared = b.square(&mut builder);
        let b_mul_b = b.mul(&mut builder, &b);

        let c_squared = c.square(&mut builder);
        let c_mul_c = c.mul(&mut builder, &c);

        Fq12Target::connect(&mut builder, &a_squared, &a_mul_a);
        Fq12Target::connect(&mut builder, &b_squared, &b_mul_b);
        Fq12Target::connect(&mut builder, &c_squared, &c_mul_c);

        let a_plus_b = a.add(&mut builder, b.clone());
        let a_plus_b_mul_c_squared = a_plus_b.mul(&mut builder, &c_squared);

        let c_mul_c_mul_a = c_mul_c.mul(&mut builder, &a);
        let c_mul_c_mul_b = c_mul_c.mul(&mut builder, &b);
        let c_c_a_plus_c_c_b = c_mul_c_mul_a.add(&mut builder, c_mul_c_mul_b);

        Fq12Target::connect(&mut builder, &c_c_a_plus_c_c_b, &a_plus_b_mul_c_squared);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        dbg!(data.common.degree_bits());
        let _proof = data.prove(pw);
    }
}
