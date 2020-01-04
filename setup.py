import os
import itertools
from subprocess import Popen, PIPE
from setuptools import setup, Extension, find_packages
from distutils.command.build import build as build_orig


# Package data
NAME = "FireDeamon"
VERSION = 0.1
DESCRIPTION = "Misc routines for the ManipulateAggregates package"
AUTHOR = "Torsten Sachse"
EMAIL = "torsten.sachse@gmail.com"
LONG_DESCRIPTION = DESCRIPTION
URL = "https://github.com/razziel89/libFireDeamon"

# Executables
UNAME = "uname"


class SetupError(RuntimeError):
    """Base error class for this setup file"""

    pass


class CommandNotFoundError(SetupError):
    """Raised when running a command cannot be found"""

    pass


class CannotRunCommandError(SetupError):
    """Raised when running a command returns with an error"""

    pass


def call_prog(prog, args):
    try:
        process = Popen([prog] + args, stdout=PIPE)
    except FileNotFoundError:
        raise CommandNotFoundError("Command '{}' not found".format(prog))
    (output, err) = process.communicate()
    exit_code = process.wait()
    if exit_code != 0:
        raise CannotRunCommandError(
            "Error running '{} {}', error: {}".format(prog, " ".join(args), err)
        )
    else:
        return output


def get_sys_lib_dirs():
    dirs = ["/usr/lib/"]
    try:
        arch = call_prog(UNAME, ["-m"])
    except CommandNotFoundError:
        arch = None
    if arch is not None:
        dirs.append("/usr/lib/{}-linux-gnu/".format(arch))
    return dirs


def get_sys_include_dirs():
    dirs = ["/usr/include/"]
    try:
        arch = call_prog(UNAME, ["-m"])
    except CommandNotFoundError:
        arch = None
    if arch is not None:
        dirs.append("/usr/include/{}-linux-gnu/".format(arch))
    return dirs


# Re-order the commands build_py and build_ext so that build_ext is run first. That way,
# swig-generated Python files will be available for the build_py step.
class build(build_orig):
    def finalize_options(self):
        super().finalize_options()
        condition = lambda l: l[0] == "build_ext"
        t1, t2 = itertools.tee(self.sub_commands)
        rest, sub_build_ext = (
            itertools.filterfalse(condition, t1),
            filter(condition, t2),
        )
        self.sub_commands[:] = list(sub_build_ext) + list(rest)


# Minimum functionality extension
fd_min_ext = Extension(
    NAME + "._cpp_min",
    sources=[
        "FireDeamon/cpp_min.i",
        "src/core/parallel_generic.cpp",
        "src/core/irregular_grid_interpolation.cpp",
        "src/core/arbitrary_grid_local_minima.cpp",
        "src/core/set_procname.cpp",
    ],
    include_dirs=get_sys_include_dirs() + ["include"],
    library_dirs=get_sys_lib_dirs(),
    libraries=["pthread"],
    language="c++",
    swig_opts=["-c++", "-Iinclude"],
    extra_compile_args=[
        "-std=c++0x",
        "-v",
        "-pedantic",
        "-Wall",
        "-Wextra",
        "-fPIC",
        "-O3",
    ],
    # extra_link_args=["-Wl,-z,defs"],
)

# Maximum functionality extension
fd_max_ext = Extension(
    NAME + "._cpp_max",
    sources=[
        "FireDeamon/cpp_max.i",
        "src/core/skin_surface_deamon.cpp",
        "src/core/electrostatic_potential_charges.cpp",
        "src/core/electron_density.cpp",
        "src/core/isosurface.cpp",
        "src/core/halfnum/angular_integral.cpp",
        "src/core/halfnum/radial_integral.cpp",
        "src/core/electrostatic_potential_orbitals.cpp",
        "src/core/constants.cpp",
        "src/core/orbital_overlap.cpp",
    ],
    include_dirs=get_sys_include_dirs() + ["include"],
    library_dirs=get_sys_lib_dirs(),
    libraries=["pthread", "mpfr", "gmp", "boost_thread", "CGAL"],
    language="c++",
    swig_opts=["-c++", "-Iinclude"],
    extra_compile_args=[
        "-std=c++0x",
        "-v",
        "-pedantic",
        "-Wall",
        "-Wextra",
        "-fPIC",
        "-O3",
        "-frounding-math",
    ],
    # extra_link_args=["-Wl,-z,defs"],
)

# Visualization extension
fd_vis_ext = Extension(
    NAME + ".visualize._cpp",
    sources=["FireDeamon/visualize/cpp.i", "src/visualize/interface.cpp"],
    include_dirs=get_sys_include_dirs() + ["include"],
    library_dirs=get_sys_lib_dirs(),
    libraries=["GL"],
    language="c++",
    swig_opts=["-c++", "-Iinclude"],
    extra_compile_args=[
        "-std=c++0x",
        "-v",
        "-pedantic",
        "-Wall",
        "-Wextra",
        "-fPIC",
        "-O3",
    ],
    # extra_link_args=["-Wl,-z,defs"],
)

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    author=AUTHOR,
    author_email=EMAIL,
    maintainer=AUTHOR,
    maintainer_email=EMAIL,
    long_description=LONG_DESCRIPTION,
    ext_modules=[fd_min_ext, fd_max_ext, fd_vis_ext],
    packages=find_packages(),
    cmdclass={"build": build},
    url=URL,
)
