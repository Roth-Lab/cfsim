import click

from cfsim.plot import plot_cfdna
from cfsim.simulate import run_simulate
from cfsim.wig import run_generate_wig


@click.command(name="simulate")
@click.option(
    "-d",
    "--hapclone-data-file",
    required=True,
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-r",
    "--hapclone-results-file",
    required=True,
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-s",
    "--snp-file",
    required=True,
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-o",
    "--out-file",
    required=True,
    type=click.Path(exists=False, resolve_path=True),
)
@click.option(
    "--clone-prevalence-file",
    required=True,
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-c",
    "--coverage",
    required=True,
    type=float,
)
@click.option(
    "-p",
    "--clone-prevalence-prior",
    required=True,
    type=float,
)
@click.option(
    "-t",
    "--tumour-content",
    required=True,
    type=float,
)
@click.option(
    "--read-length",
    required=True,
    type=click.IntRange(0),
    help="",
)
@click.option(
    "--seed",
    required=True,
    type=click.IntRange(0),
    help="",
)
def simulate(**kwargs):
    run_simulate(**kwargs)


@click.command(name="plot-cfdna")
@click.option(
    "-i",
    "--in-file",
    required=True,
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-o",
    "--out-file",
    required=True,
    type=click.Path(exists=False, resolve_path=True),
)
def plot(**kwargs):
    plot_cfdna(**kwargs)


@click.command(name="generate-wig")
@click.option(
    "-i",
    "--in-file",
    required=True,
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-o",
    "--out-file",
    required=True,
    type=click.Path(exists=False, resolve_path=True),
)
def generate_wig(**kwargs):
    run_generate_wig(**kwargs)


@click.group(name="cfsim")
@click.version_option()
def main():
    pass


main.add_command(simulate)
main.add_command(plot)
main.add_command(generate_wig)
