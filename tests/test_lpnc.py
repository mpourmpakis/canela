import pytest
import canela.lpnc as lpnc
import canela.data
import os
import importlib.resources
import ase
import ase.io

@pytest.fixture(scope='module')
def au25():
    # get Au25(PET)18 nanocluster ase.Atoms object
    with importlib.resources.path(canela.data, 'au25_pet18_-1.xyz') as p:
        return ase.io.read(p)

@pytest.fixture(scope='module')
def au25cs_dets(au25):
    # get au25 core-shell details
    return lpnc.get_core_shell(au25)


def test_coreshell__au25_core_is_13_au_atoms(au25, au25cs_dets):
    assert au25[au25cs_dets['core']].get_chemical_formula() == 'Au13'

def test_coreshell__au25_shell_is_all_atoms_not_in_core(au25, au25cs_dets):
    assert set(au25cs_dets['shell']) == set(range(len(au25))) - set(au25cs_dets['core'])


def test_count_motifs__au25_shell_only_has_6_dimers(au25):
    # get motifs dict
    motifs = lpnc.count_motifs(au25, full_cluster=True)

    # au25 should have 6 dimers
    assert 2 in motifs
    assert motifs[2].shape == (6, 5)
    assert au25[motifs[2].flatten()].get_chemical_formula() == 'Au12S18'

    # there should only be the 6 dimers in motifs dict
    assert set(motifs) == {2}

