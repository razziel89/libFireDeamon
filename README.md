# GENERAL

This is a set of C++-based Python extensions for the ManipulateAggregates
software package (see <https://github.com/razziel89/ManipulateAggregates>).

Please see the documentation for all detailed descriptions of libFireDeamon on
<https://razziel89.github.io/libFireDeamon/> (provided via GitHub pages).
Since we lack a good tool to document both Python and C++ code, both are
documented separately.
Doxygen is used for C++ whereas Sphinx is used for Python.
The Python docs are available at:
<https://razziel89.github.io/libFireDeamon/>.
They link to the C++ documentation.

# DOWNLOAD

You can get libFireDeamon easily if you have `git` installed by executing

> git clone git://github.com/razziel89/libfiredeamon.git

In the future, the package will be made available via `pip`.
Please check the more extensive online documentation for more details.

# Contributing

Contributions are very welcome!
Please simply open a pull request.
If you would like to make large-ish contributions, it might be prudent to first
contact the maintainer to better co-ordinate those efforts.

This project uses the following auto-formatters:
* Python code: black (the uncompromising Python code formatter)
  <https://github.com/psf/black>
* C++ code: clang-format <https://clang.llvm.org/docs/ClangFormat.html> with the
  options `--style="{BasedOnStyle: llvm, ColumnLimit: 100}"`

Please make sure to auto-format your pull request with those options.
Furthermore, please document any code you add.
Happy contributing!
