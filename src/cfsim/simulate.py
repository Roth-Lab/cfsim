import math
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from cfsim.dataset import DataSet
from cfsim.normalize import gc_correction


def run_simulate(
    hapclone_data_file: Path,
    hapclone_results_file: Path,
    snp_file: Path,
    clone_prevalence_file: Path,
    out_file: Path,
    coverage: float,
    clone_prevalence_prior: float,
    tumour_content: float,
    read_length: int,
    seed: int,
) -> None:

    rng = np.random.RandomState(seed)

    cell_profiles, cell_to_clone, clone_ploidy, data = load_data(
        hapclone_data_file, hapclone_results_file
    )

    clone_prevalence = generate_clone_prevalence(
        clone_ploidy,
        clone_prevalence_file,
        clone_prevalence_prior,
        rng,
    )

    target_clone_reads, target_normal_reads = compute_target_reads(
        coverage,
        data,
        clone_ploidy,
        clone_prevalence,
        read_length,
        tumour_content,
    )

    normal_reads = generate_normal_reads(
        data,
        target_normal_reads,
        rng,
    )

    cell_coverage = generate_cell_coverage(cell_to_clone, data, rng, target_clone_reads)

    cell_reads = generate_cell_reads(
        cell_coverage,
        data,
        rng,
    )

    df = data.bin_df

    df["reads"] = normal_reads + cell_reads.sum(axis=0)

    snp_density = compute_snp_density(data, snp_file)

    df["a"], df["b"] = generate_baf(
        cell_profiles,
        cell_reads,
        data,
        normal_reads,
        read_length,
        snp_density,
        rng,
    )

    df["gc"] = data.gc

    df["map"] = data.map

    df["valid"] = data.valid

    df = gc_correction(df)

    df = df.drop(columns="rdr")

    df = df.rename(columns={"beg": "start", "rdr_cor": "rdr"})

    df = df[df["valid"]]

    df.to_csv(out_file, index=False, sep="\t")


def compute_snp_density(data, snp_file):
    snp_reader = pysam.VariantFile(snp_file, "r")

    df = data.bin_df

    snp_density = np.zeros(data.num_bins, dtype=float)

    for i, (_, row) in enumerate(df.iterrows()):
        count = 0

        for record in snp_reader.fetch(row["chrom"], row["beg"], row["end"]):
            # Skip indels
            if len(record.alleles[0]) != 1 or len(record.alleles[1]) != 1:
                continue

            sample_id = list(record.samples.keys())[0]

            sample_info = record.samples[sample_id]

            # Count het SNPs
            if sample_info["GT"] in set([(0, 1), (1, 0)]):
                count += 1

        bin_len = row["end"] - row["beg"]

        snp_density[i] = count / bin_len

    return snp_density


def generate_baf(
    cell_profiles,
    cell_reads,
    data,
    normal_reads,
    read_length,
    snp_density,
    rng,
):
    tumour_a, tumour_b = generate_tumour_baf(
        cell_profiles,
        cell_reads,
        data,
        read_length,
        snp_density,
        rng,
    )

    normal_a, normal_b = generate_normal_baf(
        normal_reads,
        read_length,
        snp_density,
        rng,
    )

    a = tumour_a + normal_a

    b = tumour_b + normal_b

    return a, b


def generate_normal_baf(reads, read_length, snp_density, rng):
    e_snp = reads * snp_density * read_length

    d = rng.poisson(e_snp)

    normal_a = rng.binomial(d, 0.5)

    normal_b = d - normal_a

    return normal_a, normal_b


def generate_tumour_baf(cell_profiles, cell_reads, data, read_length, snp_density, rng):
    tumour_a = np.zeros(data.num_bins, dtype=int)

    tumour_b = np.zeros(data.num_bins, dtype=int)

    for i in range(data.num_cells):
        for j in range(data.num_bins):
            # Expected number of SNPs covered by a read in a bin
            e_snp = snp_density[j] * cell_reads[i, j] * read_length

            d = rng.poisson(e_snp)

            m = cell_profiles[i, j, 0] / cell_profiles[i, j].sum()

            a = rng.binomial(d, m)

            b = d - a

            tumour_a[j] += a

            tumour_b[j] += b

    return tumour_a, tumour_b


def compute_target_reads(
    coverage: float,
    data: DataSet,
    clone_ploidy: dict[str, float],
    clone_prevalence: dict[str, float],
    read_length: int,
    tumour_content: float,
) -> tuple[dict[str, float], float]:
    """

    Args:
        coverage (float): coverage of sample
        data (DataSet): scWGS dataset
        clone_ploidy (dict[int, float]): dictionary of clone labels to clone ploidy
        clone_prevalence (dict[int, float]): dictionary of clone labels to clone prevalences s.t sum(clone_prevalence.values()) = 1.
        read_length (int): read length of sequencer
        tumour_content (float): float between in [0, 1]

    Returns
        target_clone_reads (dict): mapping from clone label to clone number of reads
        target_normal_reads: number of normal reads
    """
    target_reads = (coverage * data.genome_size) / read_length

    m = 0

    for c in clone_ploidy:
        m += clone_prevalence[c] * clone_ploidy[c]

    m *= tumour_content

    m += (1 - tumour_content) * 2

    target_clone_reads = {}

    for c in clone_ploidy:
        target_clone_reads[c] = math.ceil(
            # ((clone_ploidy[c] * clone_prevalence[c]) / (2 + m)) * tumour_content * target_reads
            ((clone_ploidy[c] * clone_prevalence[c]) / m)
            * tumour_content
            * target_reads
        )

    # target_normal_reads = math.ceil((2 / (2 + m)) * (1 - tumour_content) * target_reads)
    target_normal_reads = math.ceil((2 / m) * (1 - tumour_content) * target_reads)

    return target_clone_reads, target_normal_reads


def generate_cell_coverage(cell_to_clone, data, rng, target_clone_reads):
    """
    Generate the total number of reads coming from each cell
    """
    cell_coverage = np.zeros(data.num_cells)

    for c in target_clone_reads:
        clone_cells = [k for k, v in cell_to_clone.items() if v == c]

        idxs = [data.cells.index(x) for x in clone_cells]

        p = data.reads[idxs].sum(axis=1).flatten()

        cell_props = rng.dirichlet(p)

        cell_coverage[idxs] = rng.multinomial(target_clone_reads[c], cell_props)

    return cell_coverage


def generate_cell_reads(cell_coverage, data, rng):
    """
    Generate reads coverage profile for each cell
    """
    cell_reads = np.zeros((data.num_cells, data.num_bins), dtype=int)

    # Sample coverage by bin for each cell
    for i in range(data.num_cells):
        p = rng.dirichlet(data.reads[i].sum(axis=-1) + 1e-6)

        cell_reads[i] = rng.multinomial(cell_coverage[i], p)

    return cell_reads


def generate_clone_prevalence(
    clone_ploidy, clone_prevalence_file, clone_prevalence_prior, rng
):
    if clone_prevalence_file is None:
        clone_prev = dict(
            zip(
                clone_ploidy.keys(),
                rng.dirichlet(clone_prevalence_prior * np.ones(len(clone_ploidy))),
            )
        )

    else:
        df = pd.read_csv(clone_prevalence_file, converters={"clone_id": str}, sep="\t")

        df = df[~df["clone_id"].str.startswith("ancestral")]

        clone_prev = df.set_index("clone_id")["mean_prevalence"].to_dict()

    return clone_prev


def generate_normal_reads(data, target_reads, rng):
    p = rng.dirichlet(data.normal_reads + 1e-6)

    return rng.multinomial(target_reads, p)


def load_data(data_file, results_file):
    data = DataSet.from_file(data_file)

    results_df = pd.read_csv(results_file, converters={"cluster_id": str}, sep="\t")

    cell_to_clone = (
        results_df[["cluster_id", "cell_id"]]
        .drop_duplicates()
        .set_index("cell_id")["cluster_id"]
        .to_dict()
    )

    bin_df = results_df[["chrom", "beg", "end"]].drop_duplicates()

    bin_df["bin_id"] = bin_df.apply(
        lambda row: "{chrom}:{beg}:{end}".format(**row.to_dict()), axis=1
    )

    results_df = pd.merge(bin_df, results_df, on=["chrom", "beg", "end"])

    data.filter_bins(results_df["bin_id"].unique())

    data.filter_cells(results_df["cell_id"].unique())

    cell_cn_a = results_df.pivot(
        index="cell_id", columns="bin_id", values="cn_A_cell"
    ).loc[data.cells, data.bins]

    cell_cn_b = results_df.pivot(
        index="cell_id", columns="bin_id", values="cn_B_cell"
    ).loc[data.cells, data.bins]

    cell_profiles = np.stack([cell_cn_a.values, cell_cn_b.values], axis=-1)

    clone_df = results_df[
        ["cluster_id", "chrom", "beg", "end", "cn_A", "cn_B"]
    ].drop_duplicates()

    clone_df["cn"] = clone_df["cn_A"] + clone_df["cn_B"]

    clone_ploidy = clone_df.groupby("cluster_id")["cn"].mean().to_dict()

    return cell_profiles, cell_to_clone, clone_ploidy, data
