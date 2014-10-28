Installation
============

.. attention:: As *piquant* has a number of dependencies on other Python packages, it is **strongly** recommended to install in an isolated environment using the `virtualenv <http://virtualenv.readthedocs.org/en/latest/index.html>`_ tool. The `virtualenvwrapper <http://virtualenvwrapper.readthedocs.org/en/latest/install.html>`_ tool makes managing multiple virtual environments easier.

Create and work in a virtual environment for *piquant* using the ``virtualenvwrapper`` tool::

    mkproject piquant

Clone the *piquant* GitHub repository into this environment::

    git clone https://github.com/lweasel/piquant.git .

Install the Python packages required by *piquant* into the virtual environment (note that it may take some time to install and build these dependencies)::

    pip install -r requirements.txt

Run unit tests for *piquant* using the command::

    PYTHONPATH=. py.test test
