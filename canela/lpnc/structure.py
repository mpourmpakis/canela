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
SCALE = 1.1

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
    """get all bond information about molecular system

    Attributes:

    Arguments:
    atoms (ase.Atoms | str): atoms object or path to molecule

    Keyword Arguments:
    scale (float): scales the bonding criteria
                   - bonding is determined by:
                   (covalent_radii * scale)
                   (DEFAULT: 1.1)

    """
    def __init__(self, atoms, scale=SCALE):
        # create atoms object from path or do nothing if already atoms object
        self.atoms = utils.make_atoms_obj(atoms)

        self.scale = scale

        # radii array
        self.radii = covalent_radii[self.atoms.numbers] * self.scale

        # create neighborlist
        self.neighborlist = ase.neighborlist.NeighborList(
                                    cutoffs=self.radii,
                                    skin=0,
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


class LPNC(object):
    """dissect structure of ligand-protected metal nanocluster (LPNC)
    - divide-and-protect theory to distinguish core and shell atoms
    - get all shell motif information, including motif types and counts
    - get atom ids, which contains fine details on its structural properties
        - core atom layer (surface, etc.)
        - shell atom's motif type and position (end S vs. middle S in motif)

    Arguments:
        atoms (ase.Atoms | str): LPNC atoms object or path to geometry file

    Keyward Arguments:
        scale (float): scale covalent radii of each atom - for determining
                       bonds when no Bonds object is passed in

    Attributes:
    atoms (ase.Atoms): LPNC atoms object
    info (dict): {'core': core atom indices, 'shell': shell atom indices, ...}
    core (ase.Atoms): atoms object containing LPNC core atoms
    shell (ase.Atoms): atoms object containing LPNC shell atoms
    ligands (list[ase.Atoms]): list of each ligand as an ase.Atoms object
    motifs (dict): {motif type: np.ndarray(motif indices)}
    ids (np.ndarray[str]): atom ids

    """
    def __init__(self, atoms, scale=SCALE):
        # create atoms object from path or do nothing if already atoms object
        self.atoms = utils.make_atoms_obj(atoms)
        self.atoms = self.atoms.copy()

        # set tag indices to atoms object
        # that way, any new atoms objects will have a mapping index
        # back to original atoms object
        self.atoms.set_tags(range(len(self.atoms)))

        self.scale = scale
        self.bonds = Bonds(self.atoms, self.scale)

        self.info = get_core_shell(self.atoms, bonds=self.bonds)

        # create core and shell atoms objects
        self.core = self.atoms[self.info['core']]
        self.shell = self.atoms[self.info['shell']]

        # get number of metal atoms, number of sulfur atoms
        self.n_m = np.isin(self.atoms.symbols, list(METALS)).sum()
        self.n_s = sum(self.atoms.symbols == 'S')


        # number of core atoms
        self.n_core = len(self.core)

        # get average core CN
        self.core_cn_avg = self.info['corecnavg']
        self.just_core_cn_avg = self.info['justcorecnavg']

        # create ase.Atoms objects of each ligand
        self.ligands = None
        self._get_ligands()

        # get number of sulfido atoms
        self.n_sulfido = len(self.info['sulfido'])

        # calculate number of ligands
        self.n_ligand = self.n_s - self.n_sulfido

        # get number of C and H per ligand
        self.n_c_per_lig = sum(self.atoms.symbols == 'C') / self.n_ligand
        self.n_h_per_lig = sum(self.atoms.symbols == 'H') / self.n_ligand

        # get motif details
        self.motifs = None
        self._get_motifs()

        # get atom structural ids
        self.ids = id_atoms(self.atoms, self.info, self.motifs, self.scale)

        # feature vector (fingerprint, fp)
        # n metals, n sulfurs, n core atoms, average CN of core atoms
        self.fp = [self.n_m, self.n_s, self.n_core, self.core_cn_avg]

        # add number of C and H per ligand to fingerprint
        self.fp += [self.n_c_per_lig, self.n_h_per_lig]

        # add motif counts to fingerprint
        # n sulfido, n bridge, n monomer, n dimer, n trimer,
        self.fp += [len(list(self.motifs.get(i, []))) for i in range(-2, 4)]

        # convert fingerprint to array
        self.fp = np.array(self.fp)

    def __len__(self):
        return len(self.atoms)

    def __repr__(self):
        return self.atoms.get_chemical_formula('metal')

    def _get_ligands(self):
        self.ligands = []
        # find all ligand atoms (atoms that are NOT metals)
        all_ligands = self.atoms[np.isin(self.atoms.symbols,
                                         list(METALS), invert=True)]
        all_inds = set(range(len(all_ligands)))
        lig_bonds = Bonds(all_ligands, self.scale)
        for s in np.where(all_ligands.symbols == 'S')[0]:
            lig = {s}

            for _ in range(1000):
                last = lig.copy()
                for i in last:
                    lig |= set(lig_bonds.coord_dict[i])
                if len(lig) == len(last):
                    self.ligands.append(all_ligands[list(lig)])
                    all_inds -= lig
                    break
            else:
                raise ValueError("Unable to correctly make ligands")

        if all_inds:
            raise ValueError("Leftover atoms when trying to make ligands")

            self.ligands.append(all_ligands[list(lig)])

    def _get_motifs(self):
        shelli = np.array(self.info['shell'])

        map_s0 = []
        if self.info['sulfido']:
            map_s0 = np.where(np.vstack(self.info['sulfido']) == shelli)[1]

        # only pass in shell atoms to avoid running get_core_shell again
        shell_motifs = count_motifs(self.shell, scale=self.scale,
                                    sulfido=map_s0)

        # map shell_motif indices back to original atoms object indices
        self.motifs = {k: shelli[v] for k, v in shell_motifs.items()}


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
        all_motifs[-1] = np.vstack(sulfido)

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
        print_motifs(all_motifs)

    return all_motifs


def get_motif_name(mot_id):
    """convert motif id number to motif name

    Arguments:
    mot_id (int): id of specific motif types
                  1: monomer, 2: dimer, 3...
                  -8: octomeric ring, -6 hexameric ring...
                  -1: sulfido atom, 0: bridging thiolate
    Returns
    (str): motif name
    """
    if mot_id == -1:
        name = 'sulfido'
    else:
        name = MOTIFNAMES.get(abs(mot_id), '%imer' % abs(mot_id))
    if mot_id < -1:
        name += 'ic-ring'

    return name


def id_atoms(atoms, cs_details=None, motifs=None, scale=SCALE):
    """encode structural type of each atom in LPNC
    - (C[core] | S[shell], ...
    - For core:
        - (C, S[surface of core] | B[bulk], int[core layer])
        - core layer = 0[center], 1[layer 1], ..., n[surface]
    - For shell:
        - (S, R[R group of ligand])
        - (S, M[metal] | S[sulfur] , int[motif type], E[end] | M[middle])
        - end = sulfur atom that terminates a motif (at least one bond to core)

    Returns:
    (np.ndarray): 1D array of atom id's ordered by atoms object indices
    """
    # initialize id dict
    ids = {}

    # set tags of atoms to match indices of original atoms object
    # NOTE: tags remain with each atom even when creating new atoms objects
    atoms.set_tags(range(len(atoms)))

    # get Bonds object
    bonds = Bonds(atoms, scale=scale)

    # get core shell details dict
    if cs_details is None:
        cs_details = get_core_shell(atoms, bonds=bonds)

    # DEFINE CORE IDs
    # use iterative apprach to find each core layer
    core_atoms = set(cs_details['core'])
    toremove = set()
    layers = {}
    layer = 0
    for _ in range(1000):
        if not core_atoms:
            break
        for c in core_atoms:
            # if c is not bonded to all other core atoms,
            # it is in the current layer
            if not all(i in core_atoms for i in bonds.coord_dict[c]):
                layers[layer] = layers.get(layer, []) + [c]
                toremove.add(c)
        # remove all atoms found in current layer and shift one layer down
        layer -= 1
        core_atoms -= toremove
        toremove = set()
    else:
        raise ValueError("Unable to identify core atom layers")

    # add back to layer count such that
    # 0: central core atom(s)
    # 1: layer 1
    # <toadd>: surface of core
    # final id: C_(B|S)_<layer number>_x
    toadd = abs(min(layers))
    for key in layers:
        bc = 'B' if key != 0 else 'S'
        for c in layers[key]:
            ids[c] = f'C_{bc}_{key+toadd:02d}_{atoms[c].symbol:x>2}'

    # create array of shell atom indices
    shell = np.array(cs_details['shell'])

    # DEFINE LIGAND IDs (R groups)
    r_group = shell[~np.isin(atoms[shell].symbols, list(METALS) + ['S'])]
    for r in r_group:
        ids[r] = 'S_R_xx_xx'

    # get motifs dict
    if motifs is None:
        map_s0 = []
        if cs_details['sulfido']:
            map_s0 = np.where(np.vstack(cs_details['sulfido']) == shell)[1]

        # only pass in shell atoms to avoid running get_core_shell again
        shell_motifs = count_motifs(atoms[shell], scale=scale, sulfido=map_s0)

        # map shell_motif indices back to original atoms object indices
        motifs = {k: shell[v] for k, v in shell_motifs.items()}

    # DEFINE SHELL MOTIF IDs
    for mtype in motifs:
        mtype_str = f'{mtype:02d}'

        for mot in motifs[mtype]:
            names = np.array(['x_x_xx_xx'] * len(mot))
            # handle rings
            if mtype < -1:
                inds = np.arange(len(names))
                if atoms[mot[0]].symbol in METALS:
                    mets = inds[::2]
                    sulf = inds[1::2]
                else:
                    mets = inds[1::2]
                    sulf = inds[::2]
                # define sulfur ids
                names[sulf] = f'S_S_{mtype_str}_xM'

                for m in mets:
                    sym = atoms[mot[m]].symbol
                    names[m] = f'S_M_{mtype_str}_{sym:x>2}'

            # handle bridge and sulfido S's
            elif mtype < 1:
                names[:] = f'S_S_{mtype_str}_xE'

            # handle all other "typical" motifs
            else:
                # define end and middle S's
                names[0] = f'S_S_{mtype_str}_xE'
                names[2:-1:2] = f'S_S_{mtype_str}_xM'
                names[-1] = f'S_S_{mtype_str}_xE'

                # define M's
                for i in range(1, len(names), 2):
                    sym = atoms[mot[i]].symbol
                    names[i] = f'S_M_{mtype_str}_{sym:x>2}'

            for i in range(len(names)):
                # make sure it isn't already filled in
                # this will avoid overwriting sulfido ids
                if mot[i] not in ids:
                    ids[mot[i]] = names[i]

    id_arr = np.array([ids[i] for i in sorted(ids)])

    if len(id_arr) != len(atoms):
        raise ValueError("unable to map to all atoms in LPNC")

    return id_arr


def make_nc_fingerprint(atoms):
    """
    [n_m, n_s, n_core, core_cn]
    """
    lpnc = LPNC(atoms)


def print_motifs(motifs):
    print('---- Motifs Info ----'.center(CEN))
    for m in sorted(motifs):
        name = get_motif_name(m)
        print(name.rjust(VCEN) + ': %i' % (len(motifs[m])))


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
