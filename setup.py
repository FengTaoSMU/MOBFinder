from setuptools import setup, find_packages

setup(
    name='mobfinder',
    version='1.0.0',
    author="FengTao",
    author_email="fengtaosmu@foxmail.com",
    description="MOB typing for plasmid metagenomic fragments",
    packages=find_packages(),
    package_data={'MOBFinder': ['MOBFinder.R', 'MOBFinder_vector.w2v', 'testdata/*']},
    url='https://github.com/FengTaoSMU/MOBFinder',
    entry_points={
        'console_scripts': [
            'mobfinder = MOBFinder.MOBFinder:main',
            'mobfinder_bin = MOBFinder.MOBFinder_bin:main'
        ]
    },
    install_requires=[
        'biopython>=1.78',
        'ete3>=3.1.2',
        'numpy>=1.24.3',
        'pycurl>=7.45.1'
        'scipy>=1.10.1',
        'six>=1.16.0',
    ],
)

