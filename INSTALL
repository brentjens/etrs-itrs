Installation
====================


The ``etrsitrs`` package uses the standard Python distutils for
building and installation. The package is written in Python
version 2. Use::

    user@host: ~/etrs-itrs/$ python2 setup.py install 

for a default installation. Although on most systems as of 2012,
python 2 is the default Python interpreter, use the ``python2``
interpreter explicitly to ensure you don't accidentally use python
version 3 on systems that have a dual Python installation. Use
``--prefix`` to install at an alternative location, for example::

    user@host: ~/etrs-itrs/$ python2 setup.py install --prefix=/opt/python


setup.py help
-------------

.. code::

   Common commands: (see '--help-commands' for more)

     setup.py build      will build the package underneath 'build/'
     setup.py install    will install the package

   Global options:
     --verbose (-v)  run verbosely (default)
     --quiet (-q)    run quietly (turns verbosity off)
     --dry-run (-n)  don't actually do anything
     --help (-h)     show detailed help message
     --no-user-cfg   ignore pydistutils.cfg in your home directory

   Options for 'install' command:
     --prefix            installation prefix
     --exec-prefix       (Unix only) prefix for platform-specific files
     --home              (Unix only) home directory to install under
     --user              install in user site-package
                         '/home/brentjens/.local/lib/python2.7/site-packages'
     --install-base      base installation directory (instead of --prefix or --
                         home)
     --install-platbase  base installation directory for platform-specific files
                         (instead of --exec-prefix or --home)
     --root              install everything relative to this alternate root
                         directory
     --install-purelib   installation directory for pure Python module
                         distributions
     --install-platlib   installation directory for non-pure module distributions
     --install-lib       installation directory for all module distributions
                         (overrides --install-purelib and --install-platlib)
     --install-headers   installation directory for C/C++ headers
     --install-scripts   installation directory for Python scripts
     --install-data      installation directory for data files
     --compile (-c)      compile .py to .pyc [default]
     --no-compile        don't compile .py files
     --optimize (-O)     also compile with optimization: -O1 for "python -O", -O2
                         for "python -OO", and -O0 to disable [default: -O0]
     --force (-f)        force installation (overwrite any existing files)
     --skip-build        skip rebuilding everything (for testing/debugging)
     --record            filename in which to record list of installed files

   usage: setup.py [global_opts] cmd1 [cmd1_opts] [cmd2 [cmd2_opts] ...]
      or: setup.py --help [cmd1 cmd2 ...]
      or: setup.py --help-commands
      or: setup.py cmd --help
