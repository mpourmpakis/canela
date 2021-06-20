#! /usr/bin/env python

from canela import __version__
import click
import numpy as np
from ase.thermochemistry import HarmonicThermo
import ase.units

"""
Thermochemistry calculations using the Harmonic Oscillator approximation
- Calculates ZPE, Cv, U, S, and T*S at a given temperature
- Also prints G_correction term in eV

Equations taken from:
(1) Thermochemistry-ASE documentation-Campos-wiki pages:
    https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html
(2) C.J. Cramer. Essentials of Computational Chemistry, 2nd Ed. Wiley, 2004.
(3) Thermochemistry in Gaussian - ELTE:
    http://organ.chem.elte.hu/farkas/teach/thermo.pdf
"""


@click.command(name='thermo',
               context_settings={'help_option_names': ['-h', '--help'],
                                 'show_default': True})
@click.version_option(__version__)
@click.argument('vib_path', type=str)
@click.option('-t', '--temp', type=float, metavar='<T>', default=298.15,
              help='temperature [K] used for calculations')
def thermo(vib_path, temp):
    """Calculate vibrational contributions to U, ZPE, Cv, and S

    vib_path (str): path to vibrational frequency data file (cm^-1 units)
    """
    # Read in vibrational frequencies (assumes cm^-1 units)
    vib = np.loadtxt(vib_path)

    # convert from cm-1 to eV
    vib_energies = vib * ase.units.invcm

    # initialize Harmonic-oscillator-approximated thermochemistry object
    thermo = HarmonicThermo(vib_energies)

    # Calculate U, ZPE, Cv, S, and ST.
    # Note that when E_pot = 0, U = ZPE + Cv
    u = thermo.get_internal_energy(temp, verbose=False)
    zpe = thermo.get_ZPE_correction()
    cv = u - zpe
    s = thermo.get_entropy(temp, verbose=False)
    ts = s * temp

    # calculate the G correction term
    g_correction = zpe + cv - ts

    # convert the results into a dictionary
    results = {'E_ZPE': zpe,
               'Cv_harm (0->T)': cv,
               'U': u,
               'S_harm': s,
               'T*S': ts,
               f'G_correction @ {temp} K': g_correction}

    # print the results in a nice format
    center = 50
    line = '=' * center
    print(line)
    print(f'Harmonic Thermochem Results at T = {temp} K'.center(center))
    print(line)

    # iteratively print each result in table format
    for name in results:
        lbl = 'eV/K' if name == 'S_harm' else 'eV'
        space = center - len(name) - 5
        val = results[name]
        print(f'{name}{val:>{space}.9f} {lbl}')
    print(line)
