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


class Bonds(object):
    def __init__(self, atoms, scale=1.0):
        # create atoms object from path or do nothing if already atoms object
        self.atoms = utils.make_atoms_obj(atoms)

        self.scale = scale

        # radii array
        self.radii = covalent_radii[self.atoms.numbers] * self.scale

        # create neighborlist
        self.neighborlist = ase.neighborlist.NeighborList(
                                    cutoffs=self.radii,
                                    self_interaction=False)
        self.neighborlist.update(self.atoms)

        # array of i-j atoms involved in a bond
        # shape = (n_bonds, 2)
        self.bond_arr = None

        # dict of atoms' 1st neighbors using a bond list
        self.coord_dict = {a: set() for a in range(len(self.atoms))}

        # array of coordination numbers (CNs)
        self.cns = None

        # dynamically build half bond array if getter is called
        self._halfbond_arr = None

        # calculate bonding, coordination dict, and cns
        self._get_bonds()

    @property
    def halfbond_arr(self):
        """array of all halfbonds
        - includes i-j and j-i in array
        - shape = (2 * n_bonds, 2)
        """
        if self._halfbond_arr is None:
            self._halfbond_arr = np.concatenate(
                                    (self.bond_arr, self.bond_arr[:, ::-1]))
        return self._halfbond_arr

    def _get_bonds(self):
        """finds bonds between atoms based on bonding radii
    
        Sets:
        self.bond_arr (2d array of ints): [[bond-1-atom-1, bond-1-atom-2],
                                          [bond-2-atom-1, bond-2-atom-2],
                                          [bond-3-atom-1, bond-3-atom-2],
                                           ...etc]
        self.coord_dict (dict: list): {0: [2, 3, 4]}

        self.cns (1d array of ints): [12, 12, 6, 6, etc]
        """
        # shorten to "n" for simplicity
        n = self.neighborlist
        bonds = np.zeros((n.nneighbors, 2), int)
        spot1 = 0
        for atomi in range(len(self.atoms)):
            # get neighbors of atomi
            neighs = n.get_neighbors(atomi)[0]

            # find second cutoff in matrix
            spot2 = spot1 + len(neighs)

            # add bonds to matrix
            bonds[spot1:spot2, 0] = atomi
            bonds[spot1:spot2, 1] = neighs

            # add neighbors and atomi details to coord dict
            self.coord_dict[atomi] |= set(neighs)
            [self.coord_dict[n].add(atomi) for n in neighs]

            # shift down matrix
            spot1 = spot2

            # once all bonds have been found break loop
            if spot1 == n.nneighbors:
                break

        # set bond list attribute
        self.bond_arr = bonds

        # convert coord_dict values from sets to sorted lists
        self.coord_dict = {k: sorted(v) for k, v in self.coord_dict.items()}

        # compute coordination numbers (CNs)
        self.cns = np.bincount(bonds.flatten())


def get_core_shell(atom, bonds=None, scale=SCALE, show=False):
    """separates LPNC into core and shell based on "divide and protect" theory

    Arguments:
        atom (ase.Atoms): metal NC atoms object

    Keyword Arguments:
        bonds (Bonds): bond object containing all bond details
                       (default: None: Bonds obj will be built)
        scale (float): scale covalent radii of each atom - for determining
                       bonds when no Bonds object is passed in
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
    # create Bonds object if not given
    if bonds is None:
        bonds = Bonds(atom, scale=scale)

    # first, find sulfido atoms (if any)
    sulfido = []
    # can only find sulfido's if R group is present
    # with no R groups, cannot differentiate b/n sulfido and bridge
    if hasr:
        # iterate over sulfur atoms
        for s in np.where(atom.symbols == 'S')[0]:
            # get indices of neighbors
            bonded_to = bonds.coord_dict[s]

            # count number of metal neighbors
            mets = sum(a.symbol in METALS for a in atom[bonded_to])

            # S is a sulfido if all neighbors are metals and > 2 neighbors
            if mets == len(bonded_to) > 2:
                sulfido.append(s)

    # initialize list of core atom indices
    core = []

    for a in atom:
        if a.symbol in METALS:
            # find S neighbors that aren't already in core
            neighs = bonds.coord_dict[a.index]
            s_neighs = sum(atom[neighs].symbols == 'S')

            # less than two S neighbors = core atom
            if s_neighs < 2:
                # add index to list of core atoms
                core.append(a.index)

    # get shell atoms
    shell = np.array(list(set(range(len(atom))) - set(core)))

    # get core CN avg
    corecnavg = 0

    # get core CN avg excluding core-shell bonds
    justcorecnavg = 0

    # get number of shell-core interactions
    nshellint = 0

    # only make these calcs if there is a core
    if len(core):
        # calc avg CN of core atoms
        corecnavg = bonds.cns[core].mean()

        # calc core CN avg excluding core-shell bonds
        core_set = set(core)
        justcore_cns = map(len, (set(bonds.coord_dict[c]) & core_set
                                 for c in core))
        justcorecnavg = np.mean(list(justcore_cns))

        # calculate min shell-to-core distance for M shell atoms
        metal_shell = [sh for sh in shell if atom[sh].symbol in METALS]

        all_dists = atom.get_all_distances()
        dists = all_dists[np.vstack(metal_shell), core].min(1)

        # add M atoms to nshellint if their distance is < 5 to core
        nshellint = sum(dists < 5)

    # find bridging motifs
    # if no R group, bridging S will have no bonds
    # else they'll have 1 bond
    match = int(hasr)

    # create matrix only containing shell bonds
    shell_bonds = bonds.bond_arr[np.isin(bonds.bond_arr, shell).all(1)]

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
    temp_ms_sulfido = None
    sulfido_counts = {}
    if len(sulfido):
        ms_sulfido = np.where(np.vstack(sulfido) == mapping_i)[1]

        # track use of sulfidos in motifs (3 uses per sulfido atom)
        sulfido_counts = {s: 0 for s in ms_sulfido}

        # add sulfido atoms to motifs (if any)
        all_motifs[-1] = sulfido
    
    # create Bonds object
    bonds = Bonds(ms, scale=scale)

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

            # use sulfidos first (each sulfido should be involved in
            # three different motifs)
            for i in sulfido_counts:
                if sulfido_counts[i] < 3:
                    # if starting with a sulfido, we've already
                    # found the starting end
                    ends_found[0] = 1
                    break

            # if no sulfidos or no sulfido uses left, start from a random atom
            else:
                i = ms_i.pop()
                used.add(i)

            # initialize new motif with atom i
            motif = [i]

        for i, last in zip([motif[-1], motif[0]], [1, 0]):
            if ends_found[last]:
                continue

            bonded_to = bonds.coord_dict.get(i, [])
            for b in bonded_to:
                # only look at new atoms
                if b in used or b in motif:
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
        if all(ends_found):
            assert len(set(motif)) == len(motif)

            # increment sulfido usage
            if len(ms_sulfido):
                for m in motif:
                    if m in sulfido_counts:
                        sulfido_counts[m] += 1

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

    # convert motifs to arrays
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

