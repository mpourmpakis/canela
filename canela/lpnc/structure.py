from __future__ import division
import canela.lpnc.utils as utils
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
MOTIFNAMES = {0: 'bridging S', 1: 'monomer', 2: 'dimer', 3: 'trimer',
              4: 'tetramer', 5: 'pentamer', 6: 'hexamer', 7: 'heptamer',
              8: 'octamer'}


def build_neighborlist(atom, scale=SCALE):
    """creates NeighborList object based on atom

    Arguments:
        atom (ase.Atoms): atoms object

    Keyword Arguments:
        scale (float): scales covalent radii of each atom
                       (default: 1.0)

    Returns:
        (NeighborList): ase neighborlist object that can calculate all
                        neighbors of a given atom
    """
    radii = covalent_radii[atom.numbers] * scale

    n = ase.neighborlist.NeighborList(
            cutoffs=radii,
            self_interaction=False)
    n.update(atom)
    return n


def build_coord_dict(bonds):
    """creates a dict of atoms' 1st neighbors using a bond list

    Arguments:
        bonds (np.ndarray | list): matrix of atom indices involved in bonds

    Returns:
        (dict): {atom_index: [1st neighbor indices]}
    """
    coords = {k: set() for k in np.unique(bonds)}
    for i, j in bonds:
        coords[i].add(j)
        coords[j].add(i)
    # convert value sets to sorted lists
    return {k: sorted(v) for k, v in coords.items()}


def flatten_ls(val, _tot=[]):
    """flattens an arbitrary list of values (ints, floats, str, etc.) and lists
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
        atom (ase.Atoms): atoms object

    Keyword Arguments:
        scale (float): scales covalent bonding radii of atoms
                         (default: 1.0)
        neighbor_list (NeighborList): a neighborlist that was already
                                            built for atoms object
                                            (default: {None})
        return_neighbor_list (bool): if True, neighbor_list is returned with
                                       bonds list
                                       (default: False)

    Returns:
        bonds (list of lists): [[bond-1-atom-1, bond-1-atom-2],
                                [bond-2-atom-1, bond-2-atom-2],
                                [bond-3-atom-1, bond-3-atom-2],
                                     ...etc]
        neighbor_list: if return_neighbor_list is True
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


def get_core_shell(atom, neighbor_list=None, scale=SCALE, show=False):
    """separates LPNC into core and shell based on "divide and protect" theory

    Arguments:
        atom (ase.Atoms): metal NC atoms object

    Keyword Arguments:
        neighbor_list (NeighborList): ase NeighborList object for metal NC
                                        (default: None)
        scale (float): scale covalent radii of each atom - for determining
                         neighbors when no neighborlist is passed in
                         (default: 1.0)
        show (bool): prints details of core and shell if True
                       (defauls: False)

    Returns:
        core shell info (dict): {core: core atom indices,
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
    # create atoms object from path or do nothing if already atoms object
    atom = utils.make_atoms_obj(atom)

    # determine if NC has R group (check for C's and H's)
    hasr = (atom.numbers == 6).any() or (atom.numbers == 1).any()

    # calculate bonds list (and neighborlist if necessary)
    if neighbor_list is None:
        bonds, neighbor_list = get_bonds(atom, scale=scale,
                                         return_neighbor_list=True)
    else:
        bonds = get_bonds(atom, neighbor_list=neighbor_list)

    # create coord dict
    coords = build_coord_dict(bonds)

    # first, find sulfido atoms (if any)
    sulfido = []
    # can only find sulfido's if R group is present
    # with no R groups, cannot differentiate b/n sulfido and bridge
    if hasr:
        # iterate over sulfur atoms
        for s in np.where(atom.symbols == 'S')[0]:
            # get indices of neighbors
            bonded_to = coords[s]

            # count number of metal neighbors
            mets = sum(a.symbol in METALS for a in atom[bonded_to])

            # S is a sulfido if all neighbors are metals and > 2 neighbors
            if mets == len(bonded_to) > 2:
                sulfido.append(s)

    # initialize list of core atom indices
    core = []

    # calc avg CN of core atoms
    corecnavg = 0
    for a in atom:
        if a.symbol in METALS:
            # find S neighbors that aren't already in core
            neighs = coords[a.index]
            s_neighs = sum(atom[neighs].symbols == 'S')

            # less than two S neighbors = core atom
            if s_neighs < 2:
                # add CN
                cn = len(neighs)
                corecnavg += cn

                # add index to list of core atoms
                core.append(a.index)

    # get shell atoms
    shell = np.array(list(set(range(len(atom))) - set(core)))

    # get core CN avg excluding core-shell bonds
    justcorecnavg = 0

    # only make these calcs if there is a core
    if len(core):
        # get core CN avg excluding core-shell bonds
        b = get_bonds(atom[core])
        justcorecnavg = np.bincount(b.flatten()).mean()

        # get avg core CN
        corecnavg /= len(core)

        # calculate min shell-to-core distance for Au shell atoms
        metal_shell = [sh for sh in shell if atom[sh].symbol in METALS]

        all_dists = atom.get_all_distances()
        dists = all_dists[np.vstack(metal_shell), core].min(1)

        # add Au atoms to nshellint if their distance is < 5 to core
        nshellint = sum(dists < 5)
    else:
        nshellint = 0

    # find bridging motifs
    # if no R group, bridging S will have no bonds
    # else they'll have 1 bond
    match = int(hasr)

    # create matrix only containing shell bonds
    shell_bonds = bonds[np.isin(bonds, shell).all(1)]

    # find bridging S's by matching <match> CN
    s_shell = shell[atom[shell].numbers == 16]

    bridges = s_shell[np.bincount(shell_bonds.flatten())[s_shell] == match]

    # add in bridging S's to nshellint
    nshellint += len(bridges)

    # create info dict
    info = {'core': core,
            'shell': shell.tolist(),
            'sulfido': sulfido,
            'bridge': bridges.tolist(),
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
        # create list of sep details
        sep_dets = [
            'N-core'.rjust(VCEN) + f': {len(core)}',
            'N-shellint'.rjust(VCEN) + f': {nshellint}',
            'Core CN Avg'.rjust(VCEN) + f': {corecnavg:.3f}',
            'Just Core CN Avg'.rjust(VCEN) + f': {justcorecnavg:.3f}\n']

        # if no core, only print shell interactions
        if not len(core):
            print('(NO CORE FOUND)'.center(CEN))
            print(sep_dets[1])

        # else print all core details
        else:
            print('\n'.join(sep_dets))

    return info


def count_motifs(atom, scale=SCALE, show=False, sulfido=[]):
    """algorithmically determine motif types and counts of metal NC

    Arguments:
        atom (ase.Atoms): metal NC atoms object

    Keyword Arguments:
        full_cluster (bool): if False, atoms object only contains shell
                               (default: False)
        scale (float): scales covalent radii when calculating neighborlist
                         (default: 1.0)
        show (bool): if True, motif info is printed
                       (default: False)
        sulfido (list): list of sulfido atom indices found from
                          get_core_shell function
                          (default: [])

    Returns:
        motif info (dict): {-1: [sulfido indices],
                             0: [bridging S indices],
                             1: [monomer (S-M-S) indices],
                             2: [dimer (S-M-S-M-S) indices],
                             ...etc}
    """
    # create atoms object from path or do nothing if already atoms object
    atom = utils.make_atoms_obj(atom)

    fc_atom = atom.copy()

    # dictionary of motifs
    # key: motif type (1 - monomer, 3 - trimer, etc.)
    # negative values: -1 = sulfido, -n = "n-meric ring"
    # values: lists of Au and S indices for motif
    all_motifs = {}


    # determine if atoms object is full NC or just shell
    # if # metal >= # S, it must be the full NC
    ns = nm = 0
    for a in atom:
        if a.symbol == 'S':
            ns += 1
        elif a.symbol in METALS:
            nm += 1
    full_cluster = nm >= ns

    # separate into shell if full cluster
    if full_cluster:
        cs_res = get_core_shell(atom, scale=scale, show=False)
        atom = ase.Atoms(atom[cs_res['shell']])
        shell_i = cs_res['shell']
        sulfido = cs_res['sulfido']
    else:
        shell_i = np.arange(len(atom))

    # create list to map ms_indices back to orig atom_indices
    # finds metal and S atoms that are in shell
    mapping_i = np.array([i.index
                          for i in fc_atom
                          if i.symbol in MS and
                          i.index in shell_i])

    # make atoms obj of just S and metals (no R)
    ms = ase.Atoms(fc_atom[mapping_i])

    # get mapped sulfido atoms (if none, set to empty list)
    ms_sulfido = []
    if len(sulfido):
        ms_sulfido = np.where(np.vstack(sulfido) == mapping_i)[1]

    # add sulfido atoms to motifs (if any)
    if len(sulfido):
        all_motifs[-1] = sulfido

    # get bonds list
    bonds = get_bonds(ms, scale=scale)

    # create coordination dict
    coords = build_coord_dict(bonds)

    # S-M-S-M-...-S motif building"""
    ms_i = set(range(len(ms))) - set(ms_sulfido)
    motif = []
    used = set()
    ends_found = [0, 0]

    # set max iterations to avoid endless while loop
    max_iter = 1000
    for _ in range(max_iter):
        if not motif:
            # if no M S atoms left, terminate loop (all motifs found)
            if not len(ms_i):
                break
            i = ms_i.pop()
            motif = [i]
            used.add(i)

        for i, last in zip([motif[-1], motif[0]], [1, 0]):
            if ends_found[last]:
                continue

            bonded_to = coords.get(i, [])
            for b in bonded_to:
                # only look at new atoms
                if b in used:
                    continue

                # find next link in motif
                # - motif must have a sulfur everyother atom!
                if sum(ms[[b, i]].symbols == 'S') == 1:
                    motif.insert(len(motif) * last, b)

                    # add b to used and remove from ms_i
                    # iff it is not a sulfido atom
                    if b not in ms_sulfido:
                        used.add(b)
                        ms_i.remove(b)
                    else:
                        ends_found[last] = 1
    
                    break
            else:
                ends_found[last] = 1

        # once both motif ends found, add it to all_motifs
        if sum(ends_found) == 2:
            # use number of atoms in motif to determine integer name (mtype)
            # S-M-S-M-S: 5 atoms // 2 = 2: dimer
            mtype = len(motif) // 2

            # if len(motif) is even, negate mtype to indicate ring
            # -S-M-S-M-S-M-S-M-: (8 atom ring // 2)* - 1 = -4: tetrameric ring
            if len(motif) % 2 == 0:
                # ring motif should have bound ends
                # assert sorted([motif[0], motif[-1]]) in bonds.tolist()
                mtype *= -1

            # get the correct indices that map back to atoms obj passed in
            atom_indices = mapping_i[motif].tolist()

            if mtype not in all_motifs:
                all_motifs[mtype] = [atom_indices]
            else:
                all_motifs[mtype].append(atom_indices)

            # reset motif list
            motif = []
            ends_found = [0, 0]

    # raise ValueError if unable to classify all M S atoms into motifs
    # within <max_iter>
    else:
        raise ValueError(f"Motif algorithm exceeded {max_iter:,} iterations.")

    for m in all_motifs:
        all_motifs[m] = np.array(all_motifs[m])

    # if show, print motif types and counts of dict in easy-to-read format
    if show:
        print('---- Motifs Info ----'.center(CEN))
        for m in sorted(all_motifs):
            if m == -1:
                name = 'sulfido'
            else:
                name = MOTIFNAMES.get(abs(m), '%imer' % abs(m))
                if m < -1:
                    name += 'ic-ring'
            print(name.rjust(VCEN) + ': %i' % (len(all_motifs[m])))

    return all_motifs


def save_view_atom(baseatom, options, args, action='save', neon_core=False):
    """creates atom object of args passed in and saves or visualizes it

    Arguments:
        baseatom (ase.Atoms): full atoms object to take
        options (dict): sections of atoms object that can be pieced together
                          to make temp atom
        args (list): list of args that should match options keys to make
                       temp atom

    Keyword Arguments:
        action (str): either save or vis temp atom
                        (default: 'save')
        neon_core (bool): convert core metal atoms to Ne
                          (default: True)
    """
    # if visualizing sulfidos, convert them to P
    # build atoms object based on args
    showme = []
    for arg in args:
        arg = arg.lower()
        if arg in options:
            add = options[arg]
            if arg == 'core' and (action == 'vis' or neon_core):
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
    compare = set()
    final = []
    for s in showme:
        # remove duplicate atoms (based on symbol in position)
        a_id = (s.symbol, s.x, s.y, s.z)
        if a_id not in compare:
            final.append(s)
            compare.add(a_id)

    final = ase.Atoms(final)
    name = '-'.join(map(str.lower, args))

    # save/view atoms obj based on action
    if action == 'save':
        final.write(name + '.xyz')
    elif action == 'vis':
        # if sulfidos are present, convert them to P before visualizing
        if 'sulfido' in options:
            sulf_tags = options['sulfido'].get_tags()
            s0_i = np.where(np.vstack(sulf_tags) == final.get_tags())[1]
            final.symbols[s0_i] = 'P'

        ase.visualize.view(final)

    # successful action is printed to output
    outp = '%s: %s' % (action, name.replace('-', ', '))
    print(outp.center(CEN))


def summ_nc_dir(dirpath, scale=SCALE):
    """calculates core shell info and motifs of all XYZ files in a given directory

    Arguments:
        dirpath (str): path to a directory containing NC .xyz files

    Keywork Arguments:
        scale (float): scale covalent radii of each atom - for determining
                         neighbors when no neighborlist is passed in
                         (default: 1.0)
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
        print('-' * CEN)

