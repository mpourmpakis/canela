import canela.lpnc as lpnc
from canela import __version__
import argparse
import sys
import ase.io
from ase.data import chemical_symbols


# default bonding scale (SCALE * covalent_radii)
SCALE = 1.0

# center output
CEN = 40

# center for value outputs
VCEN = int(CEN // 1.8)

# set of transition metals
METALS = set()
[METALS.add(sym) for at_num in [21, 39, 71, 103]
 for sym in chemical_symbols[at_num:at_num+10]]

# set of metals and S
MS = METALS.copy()
MS.add('S')

# names of motifs
MOTIFNAMES = ['bridging S', 'monomer', 'dimer', 'trimer', 'tetramer',
              'pentamer', 'hexamer', 'heptamer', 'octamer'] + \
              ['%imer' % i for i in range(9, 2000)]


def ncsep_script():
    """ncsep script"""
    parser = argparse.ArgumentParser(
                description='Separate NC into core and shell - '
                            'can provide motif info and visualize '
                            'core and shell')

    # REQUIRED
    parser.add_argument(
        'nc_path',
        type=str,
        help='path to NC geometry file')

    # OPTIONAL
    parser.add_argument(
        '--save',
        nargs='*',
        help='saves .xyz files based on arg passed in - no argument:'
             ' core and shell files')

    parser.add_argument(
        '--scale',
        type=float,
        default=SCALE,
        help='scales covalent_radii when calculating bonds (Default: %s)'
             % (str(SCALE)))

    parser.add_argument(
        '--nomotifs',
        action='store_true',
        help='do not print motif info on NC')

    parser.add_argument(
        '--sulfidocore',
        action='store_true',
        help='sulfidos are added to core and not '
             'considered when determining motifs (Default included in shell)')

    parser.add_argument(
        '--vis',
        nargs='+',
        default=None,
        help='can visualize core, shell, ligands, nc, or motifs')

    parser.add_argument(
        '--ne_core',
        action='store_true',
        help='saves core atoms as Ne atoms')

    parser.add_argument(
        '--version',
        action='store_true',
        help='print current version of script')

    # PARSE ARGUMENTS
    args = parser.parse_args()

    # print current version
    if args.version:
        print('ncsep, ' + __version__)

    # read in NC structure
    atom = ase.io.read(args.nc_path)

    # split atoms into core and shell and print core shell info
    info = lpnc.get_core_shell(atom, scale=args.scale,
                               sulfido_in_core=args.sulfidocore,
                               show=True)

    # create core atoms obj
    core = atom[info['core']]
    core.info['cnavg'] = info['corecnavg']
    core.info['cnavg_justcore'] = info['justcorecnavg']

    if info['sulfido'] and args.sulfidocore:
        core.info['nsulfido'] = len(info['sulfido'])

    # create shell atoms obj
    shell = atom[info['shell']]
    shell.info['nshellint'] = info['nshellint']

    # MOTIF INFO
    if not args.nomotifs:
        all_motifs = lpnc.count_motifs(shell, scale=args.scale, show=True,
                                       sulfido=info['sulfido'],
                                       sulfido_in_core=args.sulfidocore)

    # create options dict for saving and visualizing
    options = {'nc': atom,
               'core': core,
               'shell': shell,
               # S and Metals forming motifs
               'motifs': ase.Atoms([s for s in shell
                                    if s.symbol in MS]),
               # just ligands
               'ligands': ase.Atoms([a for a in atom
                                     if a.symbol not in METALS])}

    # add specific motifs found to visualization options
    if not args.nomotifs:
        # atoms object that motif indices are mapped to
        for mot in all_motifs:
            useatom = shell
            if mot == -1:
                useatom = atom
                name = 'sulfido'
            elif mot == 0:
                name = 'bridge'
            elif mot < -1:
                name = MOTIFNAMES[-mot] + 'ic-ring'
            else:
                name = MOTIFNAMES[mot]

            # create atoms object for each motif type
            # (flatten lists of lists of any size)
            f = lpnc.flatten_ls(list(all_motifs[mot].flatten()))
            options[name] = useatom[f]

    # SAVE INFO
    if args.save is not None:
        print('')
        print('---- Saving XYZs ----'.center(CEN))

        # save core.xyz and shell.xyz
        if args.save == []:
            core.write('core.xyz')
            shell.write('shell.xyz')
            print('save: core, shell'.center(CEN))
        else:
            lpnc.save_view_atom(atom, options, args.save, 'save', args.ne_core)

    # VIS INFO
    if args.vis is not None:
        print('')
        print('----- Vis. Info -----'.center(CEN))
        lpnc.save_view_atom(atom, options, args.vis, 'vis', args.ne_core)

    print('-' * CEN)


def main():
    # print version if that is only arg
    if len(sys.argv) > 1 and sys.argv[1] == '--version':
        print('ncsep, ' + __version__)
    else:
        ncsep_script()

