import pytest

from cfsim.cli import plot

from tests.helpers import load_config, check_cli_from_config_dict


def test_cli_simulate():

    config_file = "tests/test_cli/test_cli_plot.yaml"

    config = load_config(config_file)

    check_cli_from_config_dict(command=plot, config_dict=config)

