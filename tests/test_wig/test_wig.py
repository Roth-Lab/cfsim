import pytest

import subprocess

from cfsim.wig import run_generate_wig

def test_wig():
    tsv_file = "tests/test_wig/data/test.tsv.gz"
    out_file = "tests/test_wig/results/test.wig"
    test_file = "tests/test_wig/data/test.wig"
    run_generate_wig(tsv_file, out_file)
    result = subprocess.run(['diff', test_file, out_file], capture_output=True, text=True)
    assert result.returncode == 0


if __name__ == "__main__":
    test_wig()


