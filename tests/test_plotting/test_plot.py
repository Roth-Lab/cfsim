import pytest

import pandas as pd 

from pathlib import Path

from cfsim.plot import plot_cfdna

def test_plot_cfdna():
    in_file = "/home/matteo/projects/lrn/tutorials/pixi-play-ground/cfmonorepo/packages/cfsim/tests/test_plotting/data/ctdna.tsv.gz"
    out_file = Path("/home/matteo/projects/lrn/tutorials/pixi-play-ground/cfmonorepo/packages/cfsim/tests/test_plotting/results/test.png")
    out_file.parent.mkdir(exist_ok=True, parents=True)
    plot_cfdna(in_file=in_file, out_file=out_file)
