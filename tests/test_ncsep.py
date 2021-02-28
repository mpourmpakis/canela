import click
from click.testing import CliRunner
from canela.bin.ncsep import ncsep
import canela.data
import importlib.resources


def test_ncsep__successful_exit_code():
    # use this to test if click script worked
    with importlib.resources.path(canela.data, 'au25_pet18_-1.xyz') as p:
        runner = CliRunner()
        result = runner.invoke(ncsep, [str(p)])
        assert result.exit_code == 0
        # make sure output is correct
        for test_str in ['N-core: 13', 'N-shellint: 12', 'dimer: 6']:
            assert test_str in result.output

