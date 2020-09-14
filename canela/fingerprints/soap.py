import time

import click
from dscribe.descriptors import SOAP

from canela.io.cube import read_cube
from canela.io.pickle import write_pickle
from canela.fingerprints.grid import calc_grid


path = './'
file_cube = 'aiida-rho-ELECTRON_DENSITY-1_0.cube'

print('Reading cube...')
start = time.time()
cube = read_cube(path + file_cube)
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

print('Writing grid to disk...')
start = time.time()
write_pickle(grid, path + 'grid')
elapsed = time.time() - start
print('DONE in ', elapsed, '\n')

print('Flattening rho...')
start = time.time()
rho = rho.flatten(order='C')  # flatten in row-major order to match grid order
rho = rho.reshape(-1, 1)
elapsed = time.time() - start
print('DONE in ', elapsed, '\n')
print('Shape of rho:', rho.shape)

print('Writing rho to disk...')
start = time.time()
write_pickle(rho, path + 'rho')
elapsed = time.time() - start
print('DONE in ', elapsed, '\n')


@click.command()
@click.option('-nmax', default=2, help='Number of radial basis functions')
@click.option('-lmax', default=2, help='Max degree of spherical harmonics')
@click.option('-rcut', default=5.0, help='Cutoff radius in Angstrom')
@click.option('-procs', default=1, help='Processes for parallel execution')
def create_soap(nmax, lmax, rcut, procs):
    species = ['Au']
    soap = SOAP(
        species=species,
        periodic=False,
        rcut=rcut,
        nmax=nmax,
        lmax=lmax,
    )
    print('Calculating SOAPs on', procs, 'procs')
    print('nmax, lmax, rcut =', nmax, lmax, rcut)
    start = time.time()
    soap_struc = soap.create(struc, grid, n_jobs=procs)
    elapsed = time.time() - start
    print('DONE in ', elapsed, '\n')
    print('Shape of SOAPs:', soap_struc.shape)

    print('Writing SOAPs to disk...')
    start = time.time()
    write_pickle(soap_struc, path
                 + 'soap_n_' + str(nmax)
                 + '_l_' + str(lmax)
                 + '_r_' + str(rcut)
                 + '_p_' + str(procs))
    elapsed = time.time() - start
    print('DONE in ', elapsed, '\n')
    print('ALL DONE')


if __name__ == '__main__':
    create_soap()
