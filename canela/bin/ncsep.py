import canela.lpnc as lpnc
import canela.lpnc.utils as utils
from canela import __version__
import click
import sys
import ase.io
from ase.data import chemical_symbols
import numpy as np


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

@click.command(name='ncsep',
               context_settings={'help_option_names':['-h', '--help'],
                                 'show_default': True})
@click.version_option(__version__)
@click.argument('nc_path', type=str)
@click.option('-s', '--save', type=str, metavar='<s>', multiple=True,
              help='saves .xyz files based on arg passed in - "cs":'
              ' saves core and shell files')
@click.option('--scale', type=float, default=1.1, metavar='<f>',
              help='bond distance = covalent_radii x scale')
@click.option('-v', '--vis', type=str, metavar='<section>', multiple=True,
              help='can visualize core, shell, ligands, motifs, and/or nc' +
                   '\n visualizes core as Ne and sulfidos as P')
@click.option('--save-neon-core', is_flag=True,
              help='saves core atoms as Ne atoms')
def ncsep(nc_path, save, scale, vis, save_neon_core):
    """Dissects LPNC structure to determine: core, shell, ligands, and motifs

    nc_path: path to nc geometry file (.xyz, .pdb, etc.)
    """
    # read in NC structure
    atom = ase.io.read(nc_path)
    atom.set_tags(range(len(atom)))

    # split atoms into core and shell and print core shell info
    info = lpnc.get_core_shell(atom, scale=scale, show=True)

    # create core atoms obj
    core = []
    if info['core']:
        core = atom[info['core']]
        core.info['cnavg'] = info['corecnavg']
        core.info['cnavg_justcore'] = info['justcorecnavg']

    # create shell atoms obj
    shell = atom[info['shell']]
    shell.info['nshellint'] = info['nshellint']

    if info['sulfido']:
        shell.info['nsulfido'] = len(info['sulfido'])

    # MOTIF INFO
    # need to map sulfidos to shell indices
    map_s0 = info['sulfido']
    if info['sulfido']:
        map_s0 = np.where(np.vstack(info['sulfido']) == info['shell'])[1]

    # get motif info dict
    all_motifs = lpnc.count_motifs(shell, scale=scale, show=True,
                                   sulfido=map_s0)

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

    if not len(options['core']):
        options.pop('core')

    # add specific motifs found to visualization options
    # atoms object that motif indices are mapped to
    for mot in all_motifs:
        if mot == -1:
            name = 'sulfido'
        elif mot == 0:
            name = 'bridge'
        elif mot < -1:
            name = MOTIFNAMES[-mot] + 'ic-ring'
        else:
            name = MOTIFNAMES[mot]

        # create atoms object for each motif type
        # (flatten lists of lists of any size)
        f = utils.flatten_ls(list(all_motifs[mot].flatten()))
        options[name] = shell[f]

    # SAVE INFO
    if save:
        click.echo('')
        click.echo('---- Saving XYZs ----'.center(CEN))
        save_view_atom(atom, options, save, 'save', save_neon_core)

    # VIS INFO
    if vis:
        click.echo('')
        click.echo('----- Vis. Info -----'.center(CEN))
        save_view_atom(atom, options, vis, 'vis', save_neon_core)

    click.echo('-' * CEN)


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
