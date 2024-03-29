from setuptools import setup

import piquant


setup(
    name='piquant',
    version=piquant.__version__,
    url='https://github.com/lweasel/piquant',
    license='MIT License',
    author='Owen Dando',
    author_email='owen.dando@ed.ac.uk',
    packages=['piquant'],
    install_requires=[
        'Cython==0.20.1',
        'Jinja2>=2.10.1',
        'MarkupSafe==0.23',
        'Pygments==2.7.4',
        'Sphinx==1.2.2',
        'argparse==1.2.1',
        'docopt==0.6.1',
        'docutils==0.11',
        'matplotlib==1.4.1',
        'numpy==1.21.0',
        'pandas==0.13.0',
        'py==1.10.0',
        'pyparsing==2.0.2',
        'pytest==2.5.2',
        'python-dateutil==2.2',
        'pytz==2014.2',
        'schema==0.3.1',
        'scipy==0.13.3',
        'seaborn==0.3.1',
        'six==1.6.1',
    ],
    scripts=[
        'bin/analyse_quantification_run',
        'bin/assemble_quantification_data',
        'bin/calculate_reads_for_depth',
        'bin/calculate_unique_transcript_sequence',
        'bin/count_transcripts_for_genes',
        'bin/fix_antisense_reads',
        'bin/piquant',
        'bin/randomise_read_strands',
        'bin/simulate_read_bias'
    ],
    package_data={
        'piquant': ['bias_motif.pwm'],
    },
)
