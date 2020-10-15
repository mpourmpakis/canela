import ase.io.cube


def read_cube(path):
    """
    Read .cube format structure and voxel data via ASE
    """
    with open(path, 'r') as handle:
        cube = ase.io.cube.read_cube(
          handle, read_data=True, program=None, verbose=False)
    return cube
