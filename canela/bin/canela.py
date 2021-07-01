import click
import yaml
import importlib.resources
from canela import __version__
import canela.bin

@click.command(name='canela',
               context_settings={'help_option_names': ['-h', '--help'],
                                 'show_default': True})
@click.version_option(__version__)
@click.option('--scale', type=float, default=1.1, metavar='<f>',
              help='bond distance = covalent_radii x scale')
def main(scale):
    # separation line for cleaner format
    line = '    ' + '-' * 43
    print(line)

    header = f"""
    The Canela package contains various tools
    for comp chem analysis and visualization.

    Current Scripts:"""
    with importlib.resources.path(canela.bin, 'script_info.yml') as path:
        with open(path, 'r') as fid:
            script_info = yaml.safe_load(fid)


    print(header)
    print(line)
    for script in script_info:
        print('   ', script['name'])
        txt = script['info'].split()
        lengths = [len(t) for t in txt]
        last = 0

        # width limit for description lines
        lim = 40

        for i in range(1, len(txt)):
            if sum(lengths[last:i+1]) + (i - last) > lim:
                print('      ', ' '.join(txt[last:i]))
                last = i
        print('      ', ' '.join(txt[last:]))
    print(line)

    print('    To learn more about a script, type:')
    print('       <script_name> -h')
    print(line)
