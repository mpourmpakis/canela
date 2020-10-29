import time
import sys

import click
from dscribe.descriptors import SOAP

from canela import __version__
from canela.io.cube import read_cube
from canela.io.pickle import write_pickle
from canela.fingerprints.grid import calc_grid


@click.command()
@click.option('-filepath', default=None, help='Input cube file with data')
@click.option('-nmax', default=2, help='Number of radial basis functions')
@click.option('-lmax', default=2, help='Max degree of spherical harmonics')
@click.option('-rcut', default=5.0, help='Cutoff radius in Angstrom')
@click.option('-procs', default=1, help='Processes for parallel execution')
def create_soap(filepath, nmax, lmax, rcut, procs):
    """
    Script to calculate SOAP fingerprints on spatial grid points of
    data from input .cube file.

    Assumes Au is only species and system is nonperiodic.
    """
    species = ['Au']
    soap = SOAP(
        species=species,
        periodic=False,
        rcut=rcut,
        nmax=nmax,
        lmax=lmax,
    )
    print('Reading cube...')
    start = time.time()
    cube = read_cube(filepath)
    elapsed = time.time() - start
    print('DONE in ', elapsed, '\n')

    print('Calculating grid...')
    start = time.time()
    grid = calc_grid(cube)
    elapsed = time.time() - start
    print('DONE in ', elapsed, '\n')
    print('Shape of grid:', grid.shape)

    struc = cube['atoms']
    rho = cube['data']

    outpath = filepath.strip('cube')
    
    print('Writing grid to disk...')
    start = time.time()
    write_pickle(grid, outpath + 'grid')
    elapsed = time.time() - start
    print('DONE in ', elapsed, '\n')

    print('Flattening rho...')
    start = time.time()
    rho = rho.flatten(order='C')  # flatten in row-major order to match grid 
    rho = rho.reshape(-1, 1)
    elapsed = time.time() - start
    print('DONE in ', elapsed, '\n')
    print('Shape of rho:', rho.shape)

    print('Writing rho to disk...')
    start = time.time()
    write_pickle(rho, outpath + 'rho')
    elapsed = time.time() - start
    print('DONE in ', elapsed, '\n')
    
    print('Calculating SOAPs on', procs, 'procs')
    print('nmax, lmax, rcut =', nmax, lmax, rcut)
    start = time.time()
    soap_struc = soap.create(struc, grid, n_jobs=procs)
    elapsed = time.time() - start
    print('DONE in ', elapsed, '\n')
    print('Shape of SOAPs:', soap_struc.shape)

    print('Writing SOAPs to disk...')
    start = time.time()
    write_pickle(soap_struc, outpath
                 + 'soap_n_' + str(nmax)
                 + '_l_' + str(lmax)
                 + '_r_' + str(rcut)
                 + '_p_' + str(procs))
    elapsed = time.time() - start
    print('DONE in ', elapsed, '\n')
    print('ALL DONE')


def main():
    # print version if that is only arg
    if len(sys.argv) > 1 and sys.argv[1] == '--version':
        print('soap, ' + __version__)
    else:
        create_soap()
