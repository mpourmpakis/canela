import pytest
import canela.lpnc as lpnc
import canela.lpnc.utils as utils
import canela.data
import os
import importlib.resources
import ase
import ase.io
from ase.build import molecule


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
                   'motifs': {-1: 6, 2: 4, 3: 8, 4: 2}}
               }

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
        motifs = lpnc.count_motifs(ncs_dict[nc]['atoms'])

        # the same motif types are found (same dict keys)
        assert set(motifs) == set(ncs_dict[nc]['motifs'])

        # the amount of each motif type matches
        assert {k: len(v) for k, v
                in motifs.items()} == ncs_dict[nc]['motifs']


def test_flatten_ls__flattens_iterable_of_vals_and_iterables():
    assert utils.flatten_ls([1, 2, 'a', 4.4]) == [1, 2, 'a', 4.4]
    assert utils.flatten_ls([[1], [2], [3], [4]]) == [1, 2, 3, 4]
    assert utils.flatten_ls([0, [1, 2], [3, [4, 5, [6, 7]], 8], [9, 10, 11], 12]) == list(range(13))


def test_flatten_ls__turns_nonlist_into_list():
    assert utils.flatten_ls((1, 2, 3)) == [(1, 2, 3)]
    assert utils.flatten_ls(10) == [10]


def test_Bonds__find_all_methane_bonds():
    methane = molecule('CH4')
    bonds = lpnc.Bonds(methane)
    assert bonds.bond_arr.tolist() == [[0, 1], [0, 2], [0, 3], [0, 4]]


def test_get_bonds__0scale_return_empty_list():
    methane = molecule('CH4')
    bonds = lpnc.Bonds(methane, scale=0)
    assert bonds.bond_arr.tolist() == []

