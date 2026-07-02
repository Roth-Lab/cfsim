import pytest

import pandas as pd 

import argparse

from cfsim.simulate import simulate

from cfsim.tests.build_data import main

from dataclasses import dataclass, asdict


@dataclass
class Args:
    hapclone_data_file: str
    hapclone_results_file: str
    snp_file: str
    clone_prevalence_file: str
    out_file: str
    clone_prevalence_prior: float
    tumour_content: float
    coverage: float
    read_length: int
    seed: int


def test_implementation():

    # RUN WITH NEW CODE 

    args = Args(
        hapclone_data_file="/home/matteo/projects/cfdna/data/hapclone/TFRI004/hapclone_refit/data.h5",
        hapclone_results_file="/home/matteo/projects/cfdna/data/hapclone/TFRI004/hapclone_refit/merged_results.tsv.gz",
        snp_file="/home/matteo/projects/cfdna/data/hapclone/TFRI004/hapclone_refit/rephased_snps.bcf",
        clone_prevalence_file="/home/matteo/projects/cfdna/wfs/configs/hpc-configs/local/clone-prevs/TFRI004/clone_prevs_clone_00.tsv",
        out_file="/home/matteo/projects/lrn/tutorials/pixi-play-ground/cfmonorepo/packages/cfsim/tests/test_simulate/results/test.tsv",
        clone_prevalence_prior=0.1,
        tumour_content=0.1,
        coverage=1.,
        read_length=150,
        seed=0
    )

    cfdna = simulate(**asdict(args))

    # RUN WITH OLD CODE 

    args.out_file = "/home/matteo/projects/lrn/tutorials/pixi-play-ground/cfmonorepo/packages/cfsim/tests/test_simulate/results/test_old.tsv"

    parsed_args = argparse.Namespace(**asdict(args))

    cfdna_ols = main(parsed_args)

    # CHECK RESULTS ARE EQUAL 

    df_new = pd.read_csv("/home/matteo/projects/lrn/tutorials/pixi-play-ground/cfmonorepo/packages/cfsim/tests/test_simulate/results/test.tsv", sep='\t')

    df_old = pd.read_csv("/home/matteo/projects/lrn/tutorials/pixi-play-ground/cfmonorepo/packages/cfsim/tests/test_simulate/results/test_old.tsv", sep='\t')

    is_equal = df_new.equals(df_old)

    assert is_equal