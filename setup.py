from setuptools import setup, find_packages

setup(
    name='snapper-av',
    version='0.1.4',
    packages=find_packages(),
    install_requires=[
        'biopython',
        'intervaltree',
        'numpy',
        'pyfaidx',
        'pysam',
        'setuptools'
    ],
    llong_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'snapper=snapper.core:main',
        ],
    },
)