import ase
import ase.io


class LPNCException(Exception):
    """LPNC Exceptions"""


def get_nc_label(n, m, ligand='SCH_3'):
    return f'$\\rm M_{{{n}}}({ligand})_{{{m}}}$'


def make_atoms_obj(atoms):
    """checks to see if path can be read in as an ase.Atoms object
    - returns path back if it is already an ase.Atoms object

    Args:
    - path (str): path to geometry file

    Returns:
    - ase.Atoms: atoms object
    """
    if isinstance(atoms, ase.Atoms):
        return atoms
    try:
        atoms = ase.io.read(atoms)
        return atoms
    except FileNotFoundError:
        raise LPNCException("path does not lead to a "
                            "supported LPNC geometry file")


def make_parity(ax, add_parity_line=True):
    """Modifies axes to make a parity plot"""
    lims = list(ax.get_xlim()) + list(ax.get_ylim())
    a, b = min(lims), max(lims)
    ax.set_xlim(a, b)
    ax.set_ylim(a, b)
    if add_parity_line:
        ax.plot([a, b], [a, b], color='k', zorder=-1000, label='Parity')
    ax.set_aspect(1)


def flatten_ls(val, _tot=[]):
    """flattens an arbitrary list of values (ints, floats, str, etc.) and lists
    """
    if not isinstance(val, list):
        return _tot + [val]
    else:
        for i in val:
            _tot = flatten_ls(i, _tot)
        return _tot
