import ase
import ase.io


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
    except:
        raise ValueError("path does not lead to a "
                         "supported geometry file")


def flatten_ls(val, _tot=[]):
    """flattens an arbitrary list of values (ints, floats, str, etc.) and lists
    """
    if not isinstance(val, list):
        return _tot + [val]
    else:
        for i in val:
            _tot = flatten_ls(i, _tot)
        return _tot

