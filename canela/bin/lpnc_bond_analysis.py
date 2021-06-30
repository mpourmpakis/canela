import canela
import click
from canela import __version__
from canela.lpnc.plotting import visualize_bond_distance_analysis


@click.command(name='lpnc_bond_analysis',
               context_settings={'help_option_names': ['-h', '--help'],
                                 'show_default': True})
@click.version_option(__version__)
@click.argument('nc_path', type=str)
@click.option('-s', '--save-path', type=str, metavar='<s>',
              help='path to save bond analysis plot')
@click.option('--title', type=str, default=None, metavar='<s>',
              help='add title to plot')
@click.option('--scale', type=float, default=1.1, metavar='<f>',
              help='bond distance = covalent_radii x scale')
def main(nc_path, save_path, scale, title):
    """Creates bond analysis plot to compare bond lengths based on atom ids

    nc_path: path to nc geometry file (.xyz, .pdb, etc.)
    """
    # use function in canela.lpnc.plotting to create plot
    visualize_bond_distance_analysis(nc_path,
                                     save_path=save_path,
                                     title=title,
                                     scale=scale)
