import sys
from warnings import warn

import click
from canela import __version__

import pandas as pd
import numpy as np

"""
This script takes reaction energy information as an input.
It will return a TOF or multiple TOF(s) according to the Kozuch energy span model[1].
You can give one or multiple mechanisms simultaneously.

You may run this script in a terminal as follows:

    e_span.py path/to/reaction/file.csv temperature

Both .csv and .xlsx files are OK.
Temperature must be specified in K.

[1] Kozuch, S., & Shaik, S. (2011). Accounts of chemical research, 44(2), 101-110.
"""

# setting up constants
k_b = 8.617333262145E-05 # Boltzmann constant in eV/K
h = 4.135667696E-15 # Planck constant in eV*Hz-1

@click.command(name='e_span',
               context_settings={'help_option_names': ['-h', '--help'],
                                 'show_default': True})
@click.version_option(__version__)
@click.argument('filepath', type=str)
@click.option('-T', '--temperature', type=float, default=298.15, metavar='<f>',
              help='temperature of reaction')
def e_span(filepath,temperature):
    if filepath[-3:] == 'csv':
        data = pd.read_csv(filepath)
    elif filepath[-3:] == 'lsx':
        data = pd.read_excel(filepath)
    else:
        warn('Nothing was done. Input file must be .csv or .xlsx.')

    working = data.to_numpy()

    e_span = np.zeros(working.shape[1]-1)
    rdi_loc = np.zeros(working.shape[1]-1)
    rdts_loc = np.zeros(working.shape[1]-1)
    memory = None

    for i in range(working.shape[0]):
        post_diffs = working[i:,1:]-working[i,1:]
        pre_diffs = working[:i,1:]-working[i,1:]+working[-1,1:]-working[0,1:]
        all_diffs = np.vstack([pre_diffs,post_diffs])
        largest = np.max(all_diffs,axis=0)
        e_span = np.max([e_span,largest],axis=0)
        checks = largest >= e_span    
        
        for pos, rxn in enumerate(checks):
            if rxn == True:
                if memory == None:
                    rdi_loc[pos] = i
                    memory = i
                elif memory != None:
                    rdi_loc[pos] = memory
                rdts_loc[pos] = np.argmax(all_diffs[:,pos])

    TOF = (k_b*temperature/h)*np.exp(-e_span.astype(float)/(k_b*temperature))
    center = 54
    line = '-' * center
    print('')
    print(line)
    print('RDI/RDTS identities and TOF values for {} pathway(s)'.format(int(working.shape[1]-1)))
    print(line,'\n')

    for mech in range(int(working.shape[1]-1)):
        space = center - 15
        mech_label = data.keys()[int(mech+1)]
        rdi_label = 'RDI'
        rdi = data['state'][rdi_loc[mech]]
        rdts_label = 'RDTS'
        rdts = data['state'][rdts_loc[mech]]
        tof_label = 'TOF'
        tof = TOF[mech]
        print(f'{mech_label}{rdi_label:>{space-len(mech_label)-2}}{rdi:>{space-22}}')
        print(f'{rdts_label:>{space-2}}{rdts:>{space-22}}')
        print(f'{tof_label:>{space-2}}{tof:>{space-25}.5e} Hz\n')
    print(line,'\n')
