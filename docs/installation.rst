Installation
============

.. attention:: As *piquant* has a number of dependencies on other Python packages, it is **strongly** recommended to install in an isolated environment using the `virtualenv <http://virtualenv.readthedocs.org/en/latest/index.html>`_ tool. The `virtualenvwrapper <http://virtualenvwrapper.readthedocs.org/en/latest/install.html>`_ tool makes managing multiple virtual environments easier.

Create and work in a virtual environment for *piquant* using the ``virtualenvwrapper`` tool::

    mkproject piquant

Clone the *piquant* GitHub repository into this environment::

    git clone https://github.com/lweasel/piquant.git .

Install the *piquant* package and scripts, and their Python package dependencies, into the virtual environment by running (note that it may take some time to install and build the dependencies)::

    pip install -r requirements.txt

Run unit tests for *piquant* using the command::

    py.test test
