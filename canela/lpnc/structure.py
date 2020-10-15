from __future__ import division
import sys
import os
import numpy as np
import ase
import ase.io
import ase.neighborlist
from ase.data import covalent_radii
from ase.data import chemical_symbols
import ase.visualize


# default bonding scale (SCALE * covalent_radii)
SCALE = 1.0

# center output
CEN = 40

# center for value outputs
VCEN = int(CEN // 1.8)

# set of transition metals
METALS = set()
[METALS.add(sym) for r in [21, 39, 71, 103] for sym in chemical_symbols[r:r+10]]

# set of metals and S
MS = METALS.copy()
MS.add('S')

# names of motifs
MOTIFNAMES = ['bridging S', 'monomer', 'dimer', 'trimer', 'tetramer',
              'pentamer', 'hexamer', 'heptamer', 'octamer'] + \
              ['%imer' % i for i in range(9, 2000)]


def build_neighborlist(atom, scale=SCALE):
    """creates NeighborList object based on atom

    Arguments:
        atom {ase.Atoms} -- atoms object

    Keyword Arguments:
        scale {float} -- scales covalent radii of each atom
                         (default: {1.0})

    Returns:
        {NeighborList}: neighborlist object that can calculate all neighbors
                        of a given atom
    """
    radii = covalent_radii[atom.numbers] * scale

    n = ase.neighborlist.NeighborList(
            cutoffs=radii,
            self_interaction=False)
    n.update(atom)
    return n


def flatten_ls(val, _tot=[]):
    """
    Flattens an arbitrary list of values (ints, floats, str, etc.) and lists
    """
    if not isinstance(val, list):
        return _tot + [val]
    else:
        for i in val:
            _tot = flatten_ls(i, _tot)
        return _tot


def get_bonds(atom, scale=SCALE, neighbor_list=None,
              return_neighbor_list=False):
    """finds bonds between atoms based on bonding radii

    Arguments:
        atom {ase.Atoms} -- atoms object

    Keyword Arguments:
        scale {float} -- scales covalent bonding radii of atoms
                         (default: {1.0})
        neighbor_list {NeighborList} -- a neighborlist that was already
                                            built for atoms object
                                            (default: {None})
        return_neighbor_list {bool} -- if True, neighbor_list is returned with
                                       bonds list
                                       (default: {False})

    Returns:
        bonds list {list of lists} -- [[bond-1-atom-1, bond-1-atom-2],
                                       [bond-2-atom-1, bond-2-atom-2],
                                       [bond-3-atom-1, bond-3-atom-2],
                                       ...etc]
        neighbor_list -- if return_neighbor_list is True
    """
    if neighbor_list is None:
        n = build_neighborlist(atom, scale=scale)
    else:
        n = neighbor_list
        n.update(atom)

    if not n.nneighbors:
        return []

    bonds = np.zeros((n.nneighbors, 2), int)
    spot1 = 0
    for atomi in range(len(atom)):
        # get neighbors of atomi
        neighs = n.get_neighbors(atomi)[0]

        # find second cutoff in matrix
        spot2 = spot1 + len(neighs)

        # add bonds to matrix
        bonds[spot1:spot2, 0] = atomi
        bonds[spot1:spot2, 1] = neighs

        # shift down matrix
        spot1 = spot2

        # once all bonds have been found break loop
        if spot1 == n.nneighbors:
            break

    if return_neighbor_list:
        return bonds, n
    else:
        return bonds


def get_core_shell(atom, neighbor_list=None, scale=SCALE,
                   sulfido_in_core=False, show=False):
    """
    Separates metal NC into core and shell based on "divide and protect" theory

    Arguments:
        atom {ase.Atoms} -- metal NC atoms object

    Keyword Arguments:
        neighbor_list {NeighborList} -- ase NeighborList object for metal NC
                                        (default: {None})
        scale {float} -- scale covalent radii of each atom - for determining
                         neighbors when no neighborlist is passed in
                         (default: {1.0})
        sulfido_in_core {bool} -- True: sulfido atoms were included in core
                                  False: sulfido atoms were included in shell
                                  (default: {False})
        show {bool} -- prints details of core and shell if True
                       (defauls: {False})

    Returns:
        core shell info {dict} -- {core: core atom indices,
                                   shell: shell atom indices,
                                   sulfido: sulfido atom indices,
                                   bridge: bridging S indices,
                                   nshellint: num shell interactions,
                                   corecnavg: avg CN of core atoms
                                              (includes bonds to shell),
                                   justcorecnavg: avg CN of core excluding
                                                  bonds to shell
                                   }
    """
    # if atom isn't ase object, assume it is a path to NC structure file
    if not isinstance(atom, ase.Atoms):
        atom = ase.io.read(atom)

    # determine if NC has R group (check for C's and H's)
    hasr = (atom.numbers == 6).any() or (atom.numbers == 1).any()

    # calculate bonds list (and neighborlist if necessary)
    if neighbor_list is None:
        bonds, neighbor_list = get_bonds(atom, scale=scale,
                                         return_neighbor_list=True)
    else:
        bonds = get_bonds(atom, neighbor_list=neighbor_list)

    # first, find sulfido atoms (if any)
    sulfido = []
    # can only find sulfido's if R group is present
    # with no R groups, cannot differentiate b/n sulfido and bridge
    if hasr:
        sulfurs = np.array([i.index for i in atom if i.symbol == 'S'])
        for s in sulfurs:
            bonded_to = atom[[i for i in np.unique(bonds[
                                                   np.where(bonds == s)[0]
                                                   ]) if i != s]]
            mets = sum([1 for i in bonded_to if i.symbol in METALS])
            if mets == len(bonded_to) and mets > 2:
                sulfido.append(s)

    # core atom indices (include sulfido atoms if <sulfido_in_core>)
    core = sulfido[:] if sulfido_in_core else []
    # calc avg CN of core atoms
    corecnavg = 0
    for a in atom:
        if a.symbol in METALS:
            # find S neighbors that aren't already in core
            # (ignores sulfido atoms if <sulfido_in_core>)
            neighs = [i for i
                      in np.unique(bonds[np.where(bonds == a.index)[0]])
                      if i != a.index]
            s_neighs = sum([1 for s in neighs
                            if atom[s].symbol == 'S' and s not in core])

            # less than two S neighbors = core atom
            if s_neighs < 2:
                # add CN
                cn = len(neighs)
                corecnavg += cn

                # add index to list of core atoms
                core.append(a.index)

    # get avg core CN
    corecnavg /= len(core)

    # number of atoms in core
    num_core = len(core)

    # get core CN avg excluding core-shell bonds
    b = get_bonds(atom[core])
    justcorecnavg = np.unique(b, return_counts=True)[1].mean()

    # get shell atoms
    shell = [k.index for k in atom if k.index not in core]

    if num_core:
        # calculate min shell-to-core distance for Au shell atoms
        dists = [min([atom.get_distance(sh, c) for c in core])
                 for sh in shell if atom[sh].symbol in METALS]

        # add Au atoms to nshellint if their distance is < 5 to core
        nshellint = sum([1 for m in dists if m < 5.0])
    else:
        nshellint = 0

    # find bridging motifs
    # if no R group, bridging S will have no bonds
    # else they'll have 1 bond
    match = 1 if hasr else 0

    # create matrix only containing shell bonds
    shell_bonds = np.array([b for b in bonds if np.isin(b, shell).all()])

    # find bridging S's by matching <match> CN
    bridges = [s for s in shell
               if (shell_bonds == s).sum() == match and
               atom[s].symbol == 'S']

    # add in bridging S's to nshellint
    nshellint += len(bridges)

    # create info dict
    info = {'core': core,
            'shell': shell,
            'sulfido': sulfido,
            'bridge': bridges,
            'nshellint': nshellint,
            'corecnavg': corecnavg,
            'justcorecnavg': justcorecnavg}

    # print summary of info (if show)
    if show:
        print('')
        print(atom.get_chemical_formula('metal').center(CEN, '-'))

        if not hasr:
            print('Unable to find sulfidos (no Rs in NC)'.rjust(CEN))

        print('----- Sep. Info -----'.center(CEN))
        print('N-core'.rjust(VCEN) + ': %i' % (len(info['core'])))
        print('N-shellint'.rjust(VCEN) + ': %i' % (info['nshellint']))
        print('Core CN Avg'.rjust(VCEN) + ': %.3f' % (info['corecnavg']))
        jccn = 'Just Core CN Avg'.rjust(VCEN)
        print(jccn + ': %.3f\n' % (info['justcorecnavg']))

    return info


def count_motifs(atom, full_cluster=False, scale=SCALE,
                 show=False, sulfido=[], sulfido_in_core=False):
    """
    Algorithmically determine motif types and counts of metal NC

    Arguments:
        atom {ase.Atoms} -- metal NC atoms object

    Keyword Arguments:
        full_cluster {bool} -- if False, atoms object only contains shell
                               (default: {False})
        scale {float} -- scales covalent radii when calculating neighborlist
                         (default: {1.0})
        show {bool} -- if True, motif info is printed
                       (default: {False})
        sulfido {list} -- list of sulfido atom indices found from
                          get_core_shell function
                          (default: {[]})
        sulfido_in_core {bool} -- True: sulfido atoms were included in core
                                  False: sulfido atoms were included in shell
                                  (default: {False})

    Returns:
        motif info {dict}: {-1: [sulfido indices],
                             0: [bridging S indices],
                             1: [monomer (S-M-S) indices],
                             2: [dimer (S-M-S-M-S) indices],
                             ...etc}
    """

    if isinstance(atom, str):
        atom = ase.io.read(atom)

    fc_atom = atom.copy()

    # dictionary of motifs
    # key: motif type (1 - monomer, 3 - trimer, etc.)
    # negative values: -1 = sulfido, -n = "n-meric ring"
    # values: lists of Au and S indices for motif
    all_motifs = {}

    # separate into shell if full cluster
    if full_cluster:
        cs_res = get_core_shell(atom, scale=scale, show=False)
        atom = ase.Atoms(atom[cs_res['shell']])
        shell_i = cs_res['shell']
        sulfido = cs_res['sulfido']
    else:
        shell_i = np.arange(len(atom))

    # create list to map aus_indices back to orig atom_indices
    # finds metal and S atoms that are in shell
    mapping_i = np.array([i.index
                          for i in fc_atom
                          if i.symbol in MS and
                          i.index in shell_i])

    # make atoms obj of just S and metals (no R)
    aus = ase.Atoms(fc_atom[mapping_i])

    # add sulfido atoms to motifs (if any)
    if len(sulfido):
        all_motifs[-1] = sulfido

    # get bonds list
    bonds = get_bonds(aus, scale=scale)

    # S-M-S-M-...-S motif building"""
    aus_i = set(range(len(aus)))
    motif = []
    used = set()
    ends_found = [0, 0]
    while aus_i:
        if not motif:
            i = aus_i.pop()
            motif = [i]
            used.add(i)
            # print(aus[i].symbol, end='\r')

        # find atoms bonded to i
        for i, last in zip([motif[-1], motif[0]], [True, False]):
            if ends_found[last]:
                continue

            bonded2 = np.unique(bonds[np.where(bonds == i)[0]])
            for b in bonded2:
                # find next link in motif
                # - cannot be same atom type
                # - cannot be already included in motif
                if (b != i and (aus[b].symbol != aus[i].symbol)
                   and b not in used):
                    motif.insert(len(motif) * last, b)
                    # print('-'.join([aus[z].symbol for z in motif]), end='\r')

                    # remove b from aus_i
                    aus_i.remove(b)

                    # add b to used
                    used.add(b)
                    done = True
                    break
            else:
                ends_found[last] = 1

        # once both motif ends found, add it to all_motifs
        if sum(ends_found) == 2 or not aus_i:
            # else all of motif has been found
            if len(motif) == 1:
                mtype = 0
            elif len(motif) % 2:
                mtype = int((len(motif) - 1) / 2)
            else:
                mtype = -int(len(motif) / 2)
            atom_indices = mapping_i[motif].tolist()
            if mtype not in all_motifs:
                all_motifs[mtype] = [atom_indices]
            else:
                all_motifs[mtype].append(atom_indices)

            # reset motif list
            motif = []
            ends_found = [0, 0]
            # print('')

    for m in all_motifs:
        all_motifs[m] = np.array(all_motifs[m])

    if show:
        print_motifs(all_motifs, sulfido_in_core=sulfido_in_core)

    return all_motifs


def print_motifs(motifs_dict, sulfido_in_core=False):
    """prints motif types and counts of dict from count_motifs function

    Arguments:
        motifs_dict {dict} -- motifs dictionary returned from count_motifs
                              function

    Keyword Arguments:
        sulfido_in_core {bool} -- True: sulfido atoms were included in core
                                  False: sulfido atoms were included in shell
                                  (default: {False})
    """
    print('---- Motifs Info ----'.center(CEN))
    if -1 in motifs_dict:
        if sulfido_in_core:
            print('(sulfidos in core)'.center(CEN))
        else:
            print('(sulfidos in shell)'.center(CEN))
    for m in sorted(motifs_dict):
        if m == -1:
            name = 'sulfido'
        elif m < -1:
            name = MOTIFNAMES[-m] + 'ic-ring'
        else:
            name = MOTIFNAMES[m]
        print(name.rjust(VCEN) + ': %i' % (len(motifs_dict[m])))


def save_view_atom(baseatom, options, args, action='save', ne_core=False):
    """creates atom object of args passed in and saves or visualizes it

    Arguments:
        baseatom {ase.Atoms} -- full atoms object to take
        options {dict} -- sections of atoms object that can be pieced together
                          to make temp atom
        args {list} -- list of args that should match options keys to make
                       temp atom

    Keyword Arguments:
        action {str} -- either save or vis temp atom
                        (default: {'save'})
        ne_core {bool} -- convert core metal atoms to Ne
                          (default: {True})
    """
    # build atoms object based on args
    showme = []
    for arg in args:
        arg = arg.lower()
        if arg in options:
            add = options[arg]
            if arg == 'core' and (action == 'vis' or ne_core):
                for a in add:
                    if a.symbol in METALS:
                        a.symbol = 'Ne'
            showme += list(add)
        else:
            # quit if incorrect option given
            print('ERROR'.center(CEN))
            cant = 'cannot %s "%s"' % (action, arg)
            print(cant.center(CEN))
            title = '%s OPTIONS:' % action.upper()
            print(title.center(CEN))
            for o in options:
                print(o.center(CEN))
            return

    # remove duplicate atoms
    compare = []
    final = []
    for s in showme:
        # remove duplicate atoms (based on symbol in position)
        a_id = [s.symbol, s.x, s.y, s.z]
        if a_id not in compare:
            final.append(s)
            compare.append(a_id)

    final = ase.Atoms(final)
    name = '-'.join(map(str.lower, args))

    # save/view atoms obj based on action
    if action == 'save':
        final.write(name + '.xyz')
    elif action == 'vis':
        ase.visualize.view(final)

    # successful action is printed to output
    outp = '%s: %s' % (action, name.replace('-', ', '))
    print(outp.center(CEN))


def summ_nc_dir(dirpath, scale=SCALE, sulfido_in_core=False):
    """
    Calculates core shell info and motifs of all XYZ files in a given directory

    Arguments:
        dirpath {str} -- path to a directory containing NC .xyz files

    Keywork Arguments:
        scale {float} -- scale covalent radii of each atom - for determining
                         neighbors when no neighborlist is passed in
                         (default: {1.0})
        sulfido_in_core {bool} -- True: sulfido atoms were included in core
                                  False: sulfido atoms were included in shell
                                  (default: {False})
    """
    for f in os.listdir(dirpath):
        # read in xyz path
        atom = ase.io.read(os.path.join(dirpath, f))

        # split atoms into core and shell
        info = get_core_shell(atom,
                              scale=scale)

        shell = atom[info['shell']].copy()
        all_motifs = count_motifs(shell,
                                  scale=scale,
                                  show=True,
                                  sulfido=info['sulfido'])

        if info['sulfido']:
            info_sic = get_core_shell(atom, scale=scale, sulfido_in_core=True,
                                      show=False)
            shell_sic = atom[info_sic['shell']].copy()
            all_motifs_sic = count_motifs(shell_sic,
                                          scale=scale,
                                          show=True,
                                          sulfido=info_sic['sulfido'],
                                          sulfido_in_core=True)
        print('-' * CEN)

