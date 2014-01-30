from setuptools import setup

setup(
    name='assess_isoform_quantification',
    version=0.1,
    author="Owen Dando",
    install_requires=['docopt>=0.6.1',
                      'ez-setup>=0.9',
                      'gtf-to-genes>=1.09',
                      'py>=1.4.20',
                      'pytest>=2.5.2',
                      'schema>=0.2.0'],
    author_email='owen.dando@ed.ac.uk',
    description='Assessing performance of RNA-seq " +\
        "isoform quantification tools.',
    packages=['assess_isoform_quantification'],
    include_package_data=True,
    platforms='any'
)
