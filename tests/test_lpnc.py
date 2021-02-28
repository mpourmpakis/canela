import pytest
import canela.lpnc as lpnc
import canela.data
import os
import importlib.resources
import ase
import ase.io
from ase.build import molecule

@pytest.fixture(scope='module')
def au25():
    # get Au25(PET)18 nanocluster ase.Atoms object
    with importlib.resources.path(canela.data, 'au25_pet18_-1.xyz') as p:
        return ase.io.read(p)

@pytest.fixture(scope='module')
def au25cs_dets(au25):
    # get au25 core-shell details
    return lpnc.get_core_shell(au25)


@pytest.fixture(scope='module')
def ncs_dict():
    all_ncs = {'au18_schex14': {
                    'n_core': 8,
                    'n_shellint': 10,
                    'motifs': {1: 2, 4: 2}},
               'au20_tbbt16': {
                    'n_core': 7,
                    # actual n_shellint = 10, but algorithm just sums
                    # Au's in the shell
                    'n_shellint': 12,
                    'motifs': {-8: 1, 1: 2, 3: 1}},
               'au25_pet18_-1': {
                   'n_core': 13,
                   'n_shellint': 12,
                   'motifs': {2: 6}},
               'au60s6_pmt36': {
                   'n_core': 20,
                   'n_shellint': 32,
                   'motifs': {-7: 1, -6: 1, -3: 3,
                              -2: 1, -1: 6, 8: 2}}}

    for nc in all_ncs:
        with importlib.resources.path(canela.data, nc + '.xyz') as p:
            all_ncs[nc]['atoms'] = ase.io.read(p)
    return all_ncs


def test_get_core_shell__ncs_match_correct_core_shell_info(ncs_dict):
    for nc in ncs_dict:
        # get core shell info
        cs_info = lpnc.get_core_shell(ncs_dict[nc]['atoms'])
      
        assert [nc, len(cs_info['core'])] == [nc, ncs_dict[nc]['n_core']]
        assert [nc, cs_info['nshellint']] == [nc, ncs_dict[nc]['n_shellint']]


def test_count_motifs__ncs_match_correct_motif_counts(ncs_dict):
    for nc in ncs_dict:
        # calculate motif info
        motifs = lpnc.count_motifs(ncs_dict[nc]['atoms'], full_cluster=True)

        # the same motif types are found (same dict keys)
        assert set(motifs) == set(ncs_dict[nc]['motifs'])

        # the amount of each motif type matches
        assert {k: len(v) for k, v
                in motifs.items()} == ncs_dict[nc]['motifs']


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


def test_flatten_ls__flattens_iterable_of_vals_and_iterables():
    assert lpnc.flatten_ls([1, 2, 'a', 4.4]) == [1, 2, 'a', 4.4]
    assert lpnc.flatten_ls([[1], [2], [3], [4]]) == [1, 2, 3, 4]
    assert lpnc.flatten_ls([0, [1, 2], [3, [4, 5, [6, 7]], 8], [9, 10, 11], 12]) == list(range(13))


def test_flatten_ls__turns_nonlist_into_list():
    assert lpnc.flatten_ls((1, 2, 3)) == [(1, 2, 3)]
    assert lpnc.flatten_ls(10) == [10]


def test_get_bonds__find_all_methane_bonds():
    methane = molecule('CH4')
    bonds = lpnc.get_bonds(methane)
    assert bonds.tolist() == [[0, 1], [0, 2], [0, 3], [0, 4]]


def test_get_bonds__0scale_return_empty_list():
    methane = molecule('CH4')
    assert lpnc.get_bonds(methane, scale=0) == []

