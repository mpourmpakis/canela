import numpy as np


def calc_grid(cube):
    """
    Takes ASE cube object with data as input.
    Calculates spatial coordinates of voxels
    """
    grid_dims = cube['data'].shape
    grid_dims = np.array(grid_dims)

    cell_vecs = cube['atoms'].cell.tolist()  # one list for each lattice vector
    cell_vecs = np.array(cell_vecs)

    # Spacing vectors of grid. Grid points on cell boundaries. Hence, "-1" term
    d_a = cell_vecs[0] / (grid_dims[0] - 1)
    d_b = cell_vecs[1] / (grid_dims[1] - 1)
    d_c = cell_vecs[2] / (grid_dims[2] - 1)

    grid_coords_flat = []
    for ix in range(grid_dims[0]):
        for iy in range(grid_dims[1]):
            for iz in range(grid_dims[2]):
                coord = ix*d_a + iy*d_b + iz*d_c
                grid_coords_flat.append(coord)
    return np.array(grid_coords_flat)  # convert to np array
