from setuptools import setup, find_packages


name = 'fdest'

setup(
    name=name,
    packages=find_packages(exclude=('tests')),
    description='FGCM DES Transmission curves',
    author='Eli Rykoff',
    author_email='erykoff@stanford.edu',
    url='https://github.com/erykoff/fdest',
    install_requires=['numpy', 'fitsio', 'scipy'],
    use_scm_version=True,
    setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
)
