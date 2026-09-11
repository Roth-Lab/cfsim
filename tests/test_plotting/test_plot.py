from pathlib import Path

from cfsim.plot import plot_cfdna


def test_plot_cfdna():
    in_file = "tests/test_plotting/data/ctdna.tsv.gz"
    out_file = Path("tests/test_plotting/results/test.png")
    out_file.parent.mkdir(exist_ok=True, parents=True)
    plot_cfdna(in_file=in_file, out_file=out_file)
