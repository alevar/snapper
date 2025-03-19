from setuptools import setup, find_packages

setup(
    name='snapper-av',
    version='1.0.2',
    author='Ales Varabyou',
    author_email='ales.varabyou@jhu.edu',
    url='https://github.com/alevar/snapper',
    description='SNAPPER: Correct Intron Shifts in Alignments via Reference Annotation',
    license='GPLv3',
    packages=find_packages(),
    install_requires=[
        'biopython',
        'intervaltree',
        'numpy',
        'pyfaidx',
        'pysam',
        'setuptools'
    ],
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
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