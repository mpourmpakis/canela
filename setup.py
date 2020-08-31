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
    # ignores commented lines and blank lines
    reqs = [r.strip() for r in fid_reqs.read().splitlines()
            if r and not r.strip().startswith('#')]

setuptools.setup(name='canela',
                 version=__version__,
                 author='CANELa (Mpourmpakis Lab)',
                 url='https://www.github.com/mpourmpakis/canela',
                 description="CANELa tools to setup, track, and analyze comp chem calculations",
                 long_description=description,
                 long_description_content_type='text/markdown',
                 packages=setuptools.find_packages(),
                 python_requires='>=3.5',
                 install_requires=reqs)

