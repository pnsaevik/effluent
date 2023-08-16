Contributing
============

There are many ways of contributing to the package, all of which are greatly
appreciated. See below for examples of how you can make a contribution.



Bug reports
-----------

Report bugs at https://github.com/pnsaevik/effluent/issues.

If you are reporting a bug, please include

* Version number of Effluent
* Steps to reproduce bug
* Expected behaviour
* Actual behaviour
* Any details about your local setup that might be helpful in troubleshooting,
  such as operating system, python version, conda/pip environment etc.


Propose features or improvements
--------------------------------

Feature proposals can be submitted at https://github.com/pnsaevik/effluent/issues.
But please beware that Effluent is a small-scope project with a limited
budget. Only features which require little effort can be expected to be
implemented by the maintainers.


Contribute code
---------------

If you have the capacity to fix bugs yourself, or implement new features, this
is of course very welcome. In this case, the preferred approach is as follows:

1.  Open an issue on https://github.com/pnsaevik/effluent/issues, where you
    describe the bug/feature and your proposed implementation idea in broad
    strokes.

2.  Get feedback from the development team and agree on the way forward.

3.  Fork the Effluent repo on GitHub and clone your fork locally::

     $ cd name_of_dev_dir
     $ git clone git@github.com:your_name_here/effluent.git .

4.  Use ``conda`` or ``virtualenv`` to create an isolated development
    environment. For instance, assuming that ``conda`` is installed, use::

     $ conda create -n effluent python xarray netcdf4 holoviews scipy dask
     $ conda activate effluent
     $ cd name_of_dev_dir
     $ pip install pytest
     $ pip install -e .

5.  Use ``pytest`` to create one or more (failing) tests that demonstrate the
    bug/feature. Make local changes to the code until the test passes. Running
    tests is simple::

     $ cd name_of_dev_dir
     $ pytest

6.  Bump the minor version number (in ``__init__.py``) and make a descriptive
    entry in ``CHANGELOG.md``. Make changes to the documentation if necessary.
    New features should be documented by creating a new example, or editing
    an old one.

7.  Commit your changes and push your branch to GitHub::

     $ git add .
     $ git commit -m "Your detailed description of your changes."
     $ git push origin name-of-your-bugfix-or-feature

8.  Submit a pull request through the GitHub website. The developers may
    suggest changes to the code before the request is ultimately accepted.
