import argparse
from dataclasses import asdict, dataclass

import pandas as pd

from cfsim.plot import plot_cfdna
from cfsim.simulate import run_simulate
from tests.build_data import main


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

    new_path = "tests/test_simulate/results/TFRI004/test_new.tsv"
    old_path = "tests/test_simulate/results/TFRI004/test_old.tsv"

    args = Args(
        hapclone_data_file="tests/test_simulate/data/TFRI004/data.h5",
        hapclone_results_file="tests/test_simulate/data/TFRI004/merged_results.tsv.gz",
        snp_file="tests/test_simulate/data/TFRI004/rephased_snps.bcf",
        clone_prevalence_file="tests/test_simulate/data/TFRI004/clone_prevs.tsv",
        out_file=new_path,
        clone_prevalence_prior=0.1,
        tumour_content=0.1,
        coverage=1.0,
        read_length=150,
        seed=0,
    )

    # RUN WITH NEW CODE
    cfdna = run_simulate(**asdict(args))

    # RUN WITH OLD CODE
    args.out_file = old_path
    cfdna_ols = main(argparse.Namespace(**asdict(args)))

    # CHECK RESULTS ARE EQUAL
    df_new = pd.read_csv(new_path, sep="\t")
    df_old = pd.read_csv(old_path, sep="\t")
    plot_cfdna(in_file=new_path, out_file=new_path.replace(".tsv", ".png"))
    plot_cfdna(in_file=old_path, out_file=old_path.replace(".tsv", ".png"))

    is_equal = df_new.equals(df_old)

    assert not is_equal
