from setuptools import setup, find_packages

exec(open('supreme/_version.py').read())

name = 'supreme'

scripts = ['scripts/supreme_map_patch.py',
           'scripts/supreme_map_tract.py',
           'scripts/supreme_mapper.py']

setup(
    name=name,
    packages=find_packages(exclude=('tests')),
    version=__version__, # noqa
    description='Survey property maps for DESC using healsparse',
    author='Eli Rykoff and others',
    author_email='erykoff@stanford.edu',
    url='https://github.com/lsstdesc/supreme',
    scripts=scripts,
)
