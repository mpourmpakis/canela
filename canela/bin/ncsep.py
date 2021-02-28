import canela.lpnc as lpnc
from canela import __version__
import click
import sys
import ase.io
from ase.data import chemical_symbols


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
@click.option('--scale', type=float, default=1.0, metavar='<f>',
              help='bond distance = covalent_radii x scale')
@click.option('--no-motifs', is_flag=True, help="don't print motif info on nc")
@click.option('--sulfido-core', is_flag=True,
              help='sulfidos are added to core and not '
              'considered when determining motifs (default included in shell)')
@click.option('-v', '--vis', type=str, metavar='<section>', multiple=True,
              help='can visualize core, shell, ligands, motifs, and/or nc')
@click.option('--save-core-as-ne', is_flag=True,
              help='saves core atoms as ne atoms')
def ncsep(nc_path, save, scale, no_motifs, sulfido_core, vis, save_core_as_ne):
    """Dissects LPNC structure to determine: core, shell, ligands, and motifs

    nc_path: path to nc geometry file (.xyz, .pdb, etc.)
    """
    # read in NC structure
    atom = ase.io.read(nc_path)

    # split atoms into core and shell and print core shell info
    info = lpnc.get_core_shell(atom, scale=scale,
                               sulfido_in_core=sulfido_core,
                               show=True)

    # create core atoms obj
    core = atom[info['core']]
    core.info['cnavg'] = info['corecnavg']
    core.info['cnavg_justcore'] = info['justcorecnavg']

    if info['sulfido'] and sulfido_core:
        core.info['nsulfido'] = len(info['sulfido'])

    # create shell atoms obj
    shell = atom[info['shell']]
    shell.info['nshellint'] = info['nshellint']

    # MOTIF INFO
    if not no_motifs:
        all_motifs = lpnc.count_motifs(shell, scale=scale, show=True,
                                       sulfido=info['sulfido'],
                                       sulfido_in_core=sulfido_core)

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
    if not no_motifs:
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
    if save:
        click.echo('')
        click.echo('---- Saving XYZs ----'.center(CEN))
        lpnc.save_view_atom(atom, options, save, 'save', save_core_as_ne)

    # VIS INFO
    if vis:
        click.echo('')
        click.echo('----- Vis. Info -----'.center(CEN))
        lpnc.save_view_atom(atom, options, vis, 'vis', save_core_as_ne)

    click.echo('-' * CEN)

