[![Build Status](https://travis-ci.org/lweasel/piquant.svg?branch=develop)](https://travis-ci.org/lweasel/piquant)

[![Coverage Status](https://coveralls.io/repos/lweasel/piquant/badge.png?branch=develop)](https://coveralls.io/r/lweasel/piquant?branch=develop)

piquant
=======

*piquant* is a pipeline of Python scripts to help assess the accuracy of quantification of [transcripts](http://en.wikipedia.org/wiki/Transcriptome) from [RNA-sequencing](http://en.wikipedia.org/wiki/RNA-Seq) data. It has three stages:

* Simulation: RNA-seq reads are simulated from a starting set of transcripts with known abundances, under specified combinations of RNA-sequencing parameters.
* Quantification: A number of transcriptome quantification tools estimate transcript abundances for each set of simulated reads.
* Analysis: The comparative accuracy of expression estimates calculated by each tool can be assessed via a range of automatically generated statistics and graphs.

The latest *piquant* documentation can be found [here](http://piquant.readthedocs.org/en/latest/).

Installation
============

Note: as *piquant* has a number of dependencies on other Python packages, it is **strongly** recommended to install in an isolated environment using the [virtualenv](http://virtualenv.readthedocs.org/en/latest/index.html>) tool. The [virtualenvwrapper](http://virtualenvwrapper.readthedocs.org/en/latest/install.html>) tool makes managing multiple virtual environments easier.

After setting up ``virtualenv`` and ``virtualenvwrapper``, create and work in a virtual environment for *piquant* using the ``virtualenvwrapper`` tool:

```
mkproject piquant
```

Clone the *piquant* GitHub repository into this environment:

```
git clone https://github.com/lweasel/piquant.git .
```

Install the Python packages required by *piquant* into the virtual environment by running:

```
pip install -r requirements.txt
```

in the tool's top-level directory. Note that it may take some time to install and build these dependencies.

Testing
=======

Run unit tests for *piquant* using the command:

```
PYTHONPATH=. py.test test
```

in the tool's top-level directory.

Changelog
=========

* 1.0: First full release (28/10/14)
* 0.1: Initial beta development.
