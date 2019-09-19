from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name = 'metax',
    version = '1.0.0',
    description='Processing reports and outputs of metagenomics taxonomic classifiers',
    long_description=readme,
    author='Simon Ye',
    author_email='mail@yesimon.com',
    url='https://github.com/yesimon/metax',
    license=license,
    package_data={
        'metax.db': ['*.sh']
    },
    packages=find_packages(exclude=('tests', 'docs')),
    install_requires=[
        'cigar',
        'ncbitax',
   ],
    dependency_links=[
        'http://github.com/yesimon/ncbitax/tarball/master#egg=ncbitax-1.0'
    ],
    entry_points = {
        'console_scripts': [
            'metax = metax.__main__:main'
        ]
    })
