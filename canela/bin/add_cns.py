import click
import ase.io
from canela import __version__
from canela.lpnc import Bonds


@click.command(name='add_cns',
               context_settings={'help_option_names': ['-h', '--help'],
                                 'show_default': True})
@click.version_option(__version__)
@click.argument('geom_path', type=str, metavar='geom_path')
@click.option('--scale', type=float, default=1.1, metavar='<f>',
              help='bond distance = covalent_radii x scale')
@click.option('--save_path', type=str, default=None, metavar='<s>',
              help='if given, used <save_path>.xyz as file name')
def main(geom_path: str, scale: float, save_path: str) -> None:
    """Adds CNs to XYZ geometry file for visualization in ASE GUI.

    \b
    - Sets CNs as "inital_charges" to leverage ASE's vis tools
    - Saves xyz as "<geom file name>-WITH-CNs.xyz by default"
    - Else uses <save_path>.xyz as file name
    - NOTE: Currently uses covalent_radii * scale for bonding.
            Future versions should allow for custom radii.

    \b
    To vis CNs:
    1. Open new xyz with ase (ase gui <path 2 xyz>)
    2. Select 'View' -> 'Show Labels' -> 'Initial Charges'

    geom_path: path to geometry file (.xyz, .pdb, etc.)
    """
    # create ase.Atoms object
    atoms = ase.io.read(geom_path)

    # get bond details
    bonds = Bonds(atoms, scale)

    # set CNs as initial_charges of atoms object
    atoms.set_initial_charges(bonds.cns)

    # if no save_path given, create file name from geom_path
    if save_path is None:
        # get filename by removing the extension
        name = geom_path.rsplit('.', 1)[0]

        # add the standard extension
        save_path = name + '_WITH-CNs.xyz'
    else:
        # ensure save_path has xyz extension
        if not save_path.endswith('.xyz'):
            save_path += '.xyz'

    # save atoms object with CNs
    atoms.write(save_path)
    print(f'Saved {save_path}.')
