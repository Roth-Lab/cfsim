import click

from cfsim.simulate import run_simulate

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

@click.group(name="cfsim")
@click.version_option()
def main():
    pass


main.add_command(simulate)
