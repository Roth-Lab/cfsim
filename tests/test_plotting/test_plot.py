import pytest

import pandas as pd 

from pathlib import Path

from cfsim.plot import plot_cfdna

def test_plot_cfdna():
    in_file = "/home/matteo/projects/lrn/tutorials/pixi-play-ground/cfmonorepo/packages/cfsim/tests/test_plotting/data/test.tsv"
    out_file = Path("/home/matteo/projects/lrn/tutorials/pixi-play-ground/cfmonorepo/packages/cfsim/tests/test_plotting/results/test.png")
    out_file.parent.mkdir(exist_ok=True, parents=True)
    df = pd.read_csv(in_file, sep='\t')
    plot_cfdna(cfdna=df, out_file=out_file)
