from pathlib import Path

def simulate(
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
    print(hapclone_data_file)
    print("Simulated cfDNA")