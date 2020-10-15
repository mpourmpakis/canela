from canela import __version__
import sys
import os
import shutil
from subprocess import PIPE, Popen
import argparse
import ase.io
import ase.units
from ase.utils import basestring
from ase.calculators.turbomole import Turbomole


def execute(args, input_str=None, error_test=True,
            stdout_tofile=True):
    """executes a turbomole executable and process the outputs"""
    if isinstance(args, basestring):
        args = args.split()

    if stdout_tofile:
        stdout_file = 'ASE.TM.' + args[0] + '.out'
        stdout = open(stdout_file, 'w')
    else:
        stdout = PIPE

    if input_str:
        stdin = input_str.encode()
    else:
        stdin = None

    message = 'TM command "' + args[0] + '" execution failed'
    try:
        proc = Popen(args, stdin=PIPE, stderr=PIPE, stdout=stdout)
        res = proc.communicate(input=stdin)
        if error_test:
            error = res[1].decode()
            if 'abnormally' in error or 'ended normally' not in error:
                message += ' with error:\n' + error
                message += '\nSee file ' + stdout_file + ' for details.\n'
                raise RuntimeError(message)
    except RuntimeError as err:
        raise err
    except OSError as err:
        raise OSError(err.args[1] + '\n' + message)
    else:
        print('TM command: "' + args[0] + '" successfully executed')

    if not stdout_tofile:
        return res[0].decode()


class CANELaTurbomole(Turbomole):
    def __init__(self, *args, **kwargs):
        super(CANELaTurbomole, self).__init__(*args, **kwargs)

    def initialize(self, d3=False):
        """
        extended method to make slight changes to control file
        """
        """prepare turbomole control file by running module 'define'"""
        if self.initialized:
            return
        self.verify_parameters()
        if not self.atoms:
            raise RuntimeError('atoms missing during initialization')
        if not os.path.isfile('coord'):
            raise IOError('file coord not found')

        if self.define_str is not None:
            define_str = self.define_str
        else:
            define_str = self.get_define_str()


        # missing a 'y' in define string!
        define_str = define_str.replace('eht\ny\n', 'eht\ny\ny\n')

        # run define
        execute('define', input_str=define_str)

        # process non-default initial guess
        iguess = self.parameters['initial guess']
        if isinstance(iguess, dict) and 'use' in iguess.keys():
            # "use" initial guess
            if self.parameters['multiplicity'] != 1 or self.parameters['uhf']:
                define_str = '\n\n\ny\nuse ' + iguess['use'] + '\nn\nn\nq\n'
            else:
                define_str = '\n\n\ny\nuse ' + iguess['use'] + '\nn\nq\n'
            execute('define', input_str=define_str)
        elif self.parameters['initial guess'] == 'hcore':
            # "hcore" initial guess
            if self.parameters['multiplicity'] != 1 or self.parameters['uhf']:
                delete_data_group('uhfmo_alpha')
                delete_data_group('uhfmo_beta')
                add_data_group('uhfmo_alpha', 'none file=alpha')
                add_data_group('uhfmo_beta', 'none file=beta')
            else:
                delete_data_group('scfmo')
                add_data_group('scfmo', 'none file=mos')

        self._set_post_define()

        self.initialized = True
        self.converged = False

        # APPENDED SECTIONS - CUSTOM CANELa THINGS
        # always use automatic for orbital shift
        subprocess.call(
                "sed -i 's/closedshell=.05/automatic=0.1/' control",
                shell=True)

        # read in control file
        with open('control', 'r') as fidr:
            lines = fidr.readlines()

        # turn on $marij if ri is used
        if self.parameters['use resolution of identity']:
            print('Adding $marij...')

            # add $marij if not in control file
            if '$marij\n' not in lines:
                key = 'file=auxbasis'
                subprocess.call(
                    "sed -i 's/{0}/{0}\\n$marij/' control".format(key),
                    shell=True)

        # turn on dispersion if d3
        if d3:
            print('Adding D3 correction...')
            subprocess.call("sed -i 's/{0}/{0}\\n$disp3/' control".format('$rij'),
                            shell=True)



def make_tmole_slurm(runtime=48, title=None,
                     mem=64, send_email=True,
                     graceful_end=False, fname='turbomole.sl',
                     overwrite=True):
    """Creates turbomole slurm file

    KArgs:
    - runtime (int): max wall time in hours
                     (DEFAULT: 48)

    - title (str): title of run
                   (DEFAULT: 'tmole run')

    - mem (int): max memory allocated for job (in GB)
                 (DEFAULT: 64 GB)

    - send_email (bool): send an email when job completes or fails
                         (DEFAULT: True)

    - graceful_end (bool): if True, a 'STOP' file is added to directory
                           3 hours before max wall time
                           (DEFAULT: False)

    - fname (str): name of slurm file
                   (DEFAULT: turbomole.sl)

    - overwrite (bool): if True, function will overwrite <fname> if it exists
                        (DEFAULT: True)
    """
    # fname must have .sl extension
    if not fname.endswith('.sl'):
        fname += '.sl'

    # get base file name
    basename = fname[:-3]

    # set default title if None is given
    if title is None:
        title = 'tmole run'

    # start building raw text of slurm file
    text = '#!/bin/bash\n'

    # job name
    text += '#SBATCH -J ' + '"' + title + '"\n'

    # always 1 node and 8 core/node for turbomole calcs
    text += '#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=8\n'

    # job memory (in MB)
    text += '#SBATCH --mem=%i000\n' % mem

    # run time
    text += '#SBATCH -t %02i:00:00\n' % runtime

    # always SMP for turbomole calcs
    text += '#SBATCH --cluster=smp\n'

    # add email alerts
    if send_email:
        email = os.environ['USER'] + '@pitt.edu'
        text += '#SBATCH --mail-type=END,FAIL\n'
        text += '#SBATCH --mail-user=%s\n' % email

    # spacer
    text += '\n'

    # add ulimit lines
    text += 'ulimit -s unlimited\n'
    text += 'ulimit -l unlimited\n'

    # spacer
    text += '\n'

    # only load turbomole module
    text += 'module purge\nmodule load turbomole/7.02-smp\n'

    # spacer
    text += '\n'

    # needed to run turbomole on SMP
    text += '# needed to run turbomole on SMP\n'
    text += 'export TM_PAR_FORK=2\n'

    # spacer
    text += '\n'

    # if a "graceful" end is desired
    # (turbomole is more likely to stop in between steps)
    # minimum of 6 hours
    if graceful_end:
        text += '# process should end gracefully\n'
        text += 'sleep %ih && touch stop &\n' % max(runtime - 3, 6)
        # spacer
        text += '\n'

    # add turbomole submission lines
    text += '# TURBOMOLE SUBMISSION LINES (based on runtype)\n'

    # geo opt (500 steps) is uncommented by default
    text += '# geometry optimization\n'
    text += 'jobex -c 500 -ri >& jobex.out\n'

    text += '# single point energy calculation\n'
    text += '#ridft >& ridft.out\n'

    text += '# vibrational calculation\n'
    text += '#aoforce >& force.out\n'

    text += '# excited state calculation (TDDFT)\n'
    text += '#escf >& escf.out\n'

    # spacer
    text += '\n'

    # print jobstats (custom CRC script)
    text += 'crc-job-stats.py > CRC-JOB-STATS.txt\n'

    # write the slurm file in current working directory
    # ensure no overwriting will be done (if specified)
    if not overwrite:
        i = 1
        fname = '%s_%02i.sl' % (basename, i)
        while os.path.isfile(fname):
            i += 1
        
    with open(fname, 'w') as fid:
        fid.write(text)
 

def sj_tmole(args):
    # create atoms object
    atoms = ase.io.read(args.geom_path)

    # is it a restart file?
    # TODO: IMPLEMENT RESTART KArg
    restart = False

    # calc default multiplicity
    n_electrons = atoms.get_atomic_numbers().sum() - args.charge
    default_mult = n_electrons % 2 + 1

    # ensure specified multiplicity (if given) is allowed
    if args.multiplicity:
        if args.multiplicity % 2 != default_mult % 2:
            print('Invalid multiplicity specified. Using mult=%i.'
                  % default_mult)
            args.multiplicity = default_mult
    else:
        args.multiplicity = default_mult

    # turn on uhf if not singlet
    if bool((default_mult - 1) % 2) and not args.uhf:
        print('Open shell system. Using UHF.')
        args.uhf = True

    # calc SCF convergence value
    scf = float('1E-%i' % args.scfconv) * ase.units.Ha

    # create parameters dictionary
    params = {
        'task': 'energy',
        'total charge': args.charge,
        'restart': restart,
        'multiplicity': args.multiplicity,
        'use dft': True,
        'density functional': args.xcfunc,
        'basis set name': args.basis,
        'scf iterations': args.scfiterlimit,
        'scf energy convergence': scf,
        'use resolution of identity': True,
        'ri memory': 20000,
        'force convergence': 0.02,
        'initial damping': args.idamp,
        'initial guess': 'eht',
        'uhf': args.uhf,
        'use fermi smearing': args.fermi
    }


    def print_dict(d, just=35):
        for key in sorted(d):
            print(str(key).rjust(just) + ': ' + str(d[key]))



    # initialize calculator and attach it to atoms object
    calc = CANELaTurbomole(**params)
    atoms.set_calculator(calc)

    # print final params
    if args.verbose:
        print_dict(calc.parameters)

    # create control file
    calc.initialize(d3=args.d3)

    # create slurm file
    make_tmole_slurm(runtime=args.runtime, title=args.title,
                     send_email=not args.noemail, mem=args.mem)

    # convert slurm to run ridft if jobtype is energy
    if args.jobtype.lower() == 'energy':
        subprocess.call(['sed', '-i', 's/jobex/#jobex/', 'turbomole.sl'])
        subprocess.call(['sed', '-i', 's/#ridft/ridft/', 'turbomole.sl'])
        print('turning on ridft and turning off jobex in slurm')


def parse_arguments():
    parser = argparse.ArgumentParser(prog='sj_tmole',
                                     description='Setup Turbomole jobs')

    # REQUIRED
    parser.add_argument(
        'geom_path',
        type=str,
        help='path to molecular geometry file')

    # OPTIONAL
    parser.add_argument(
        '-q', '--charge',
        type=int,
        default=0,
        metavar='I',
        help='charge of the system')

    parser.add_argument(
        '-m', '--multiplicity',
        type=int,
        metavar='I',
        help='multiplicity of system '
        '(DEFAULT: 1 for closed shell, 2 for open shell systems)')

    parser.add_argument(
        '--scfconv',
        type=int,
        default=6,
        metavar='I',
        help='SCF convergence: 1 x 10^(-I) Ha, allowed values: [4-9]')

    parser.add_argument(
        '--idamp',
        type=float,
        default=0.7,
        metavar='F',
        help='initial damping (increase if having convergence issues)')

    parser.add_argument(
        '-t', '--runtime',
        type=int,
        default=48,
        metavar='HRS',
        help='runtime of turbomole calculation (in hours) (DEFAULT: 48)')

    parser.add_argument(
        '-f', '--xcfunc',
        type=str,
        default='pbe',
        metavar='XC',
        help='specify xc-functional to use (DEFAULT: PBE)')

    parser.add_argument(
        '-b', '--basis',
        type=str,
        default='def2-SV(P)',
        metavar='BS',
        help='specify basis set to use (DEFAULT: def2-SV(P))')

    parser.add_argument(
        '-j', '--jobtype',
        type=str,
        default='GEO_OPT',
        metavar='JT',
        help='specify job type (ENERGY or the default: GEO_OPT)')

    parser.add_argument(
        '--scfiterlimit',
        type=int,
        default=300,
        metavar='I',
        help='max iterations allowed in scf calculation (DEFAULT: 300)')

    parser.add_argument(
        '--mem',
        type=int,
        default=64,
        help='specified memory (in GB) for job - default = 64GB')

    parser.add_argument(
        '--title',
        type=str,
        default=None,
        help='title of calculation')

    # FLAGS (no argument required)
    parser.add_argument(
        '--d3',
        action='store_true',
        help='turns on D3 correction')

    parser.add_argument(
        '--uhf',
        action='store_true',
        help='run spin polarized calculation (UKS/UHF)')

    parser.add_argument(
        '--fermi',
        action='store_true',
        help='use fermi smearing (at 300K)')

    parser.add_argument(
        '--noemail',
        action='store_true',
        help='slurm will not email you when job completes or fails')

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='print out more explicit details while setting up job')

    parser.add_argument(
        '--start',
        action='store_true',
        help='automatically submits the job to run')

    parser.add_argument(
        '--version',
        action='store_true',
        help='prints current version of script')

    # parse and return arguments
    args = parser.parse_args()

    # print version if necessary
    if args.version:
        print('sj_tmole, ' + __version__)

    return args


def main():
    # only print version if no geom path is given
    if len(sys.argv[1]) > 1 and sys.argv[1] == '--version':
        print('sj_tmole, ' + __version__)
    else:
        args = parse_arguments()
        sj_tmole(args)

        # submit the job
        if args.start:
            out = subprocess.check_output(['sbatch', 'turbomole.sl']).decode()
            print(out)

