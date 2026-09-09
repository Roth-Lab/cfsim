import pytest

from cfsim.cli import generate_wig

from tests.helpers import load_config, check_cli_from_config_dict


def test_cli_simulate():

    config_file = "tests/test_cli/test_cli_generate_wig.yaml"

    config = load_config(config_file)

    check_cli_from_config_dict(command=generate_wig, config_dict=config)