from setuptools import setup, find_packages

setup(
    name='snapper',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'pysam',
        'biopython',
        'pyfaidx',
        'intervaltree',
        'numpy',
        'subprocess32',
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