import pytest
import importlib.resources

from ase import Atoms

from canela.io import cube, pickle
import canela.data


@pytest.fixture(scope='module')
def au13rho():
    # get path to Au13 cube file with rho data on 2x2x2 grid
    with importlib.resources.path(canela.data, 'au13rho.cube') as f:
        return f

    
def test_read_cube__au13_imports_as_dict(au13rho):
    # ASE should import .cube files as dict 
    assert isinstance(cube.read_cube(au13rho), dict)

def test_read_cube__au13_dict_has_atoms_object(au13rho):
    # basic check of ASE expected import behavior
    cube_dict = cube.read_cube(au13rho)
    atoms = cube_dict['atoms']
    assert isinstance(atoms, Atoms)

def test_read_cube__au13_atoms_formula_is_Au13(au13rho):
    # lame check of an Atoms method to see if ASE is broken
    # Should refine based on needs of dscribe.soap
    atoms = cube.read_cube(au13rho)['atoms']
    assert atoms.get_chemical_formula() == 'Au13'

def test_read_cube__au13_data_shape_is_2x2x2(au13rho):
    # make sure rho data is numpy array interpretable as voxels
    cube_dict = cube.read_cube(au13rho)
    data = cube_dict['data']
    assert data.shape == (2, 2, 2)
