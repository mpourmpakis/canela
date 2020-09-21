import setuptools

# read in version of canela package
with open('canela/_version.py', 'r') as fid:
    # this will create the __version__ variable to be used below
    exec(fid.read())

# use README file to create long description of package
# ignore images (lines that start with '![')
with open('README.md', 'r') as readme:
    description = ''.join([i for i in readme.readlines()
                           if not i.startswith('![')])

# get required packages from requirements.txt
with open('requirements.txt', 'r') as fid_reqs:
    # ignores commented lines, blank lines, and pytest
    reqs = [r.strip() for r in fid_reqs.read().splitlines()
            if r and not r.strip().startswith('#')
            and 'pytest' not in r]

setuptools.setup(
    name='canela',
    version=__version__,
    author='CANELa (Mpourmpakis Lab)',
    url='https://www.github.com/mpourmpakis/canela',
    description="CANELa tools to setup, track, and analyze comp chem calcs",
    long_description=description,
    long_description_content_type='text/markdown',
    packages=['canela'],
    entry_points={'console_scripts': ['ncsep=canela.bin.ncsep:main',
                                      'sj_tmole=canela.bin.sj_tmole:main']},
    include_package_data=True,  # include data files
    exclude_package_data={'': ['README.md']},
    python_requires='>=3.5',
    install_requires=reqs)

