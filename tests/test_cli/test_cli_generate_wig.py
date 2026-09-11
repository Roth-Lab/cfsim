from cfsim.cli import generate_wig
from tests.helpers import check_cli_from_config_dict, load_config


def test_cli_simulate():

    config_file = "tests/test_cli/test_cli_generate_wig.yaml"

    config = load_config(config_file)

    check_cli_from_config_dict(command=generate_wig, config_dict=config)
