import os
import itertools
from subprocess import Popen, PIPE
from setuptools import setup, Extension, find_packages
from distutils.command.build import build as build_orig

# Misc executables
UNAME = "uname"


def get_long_description():
    """Load the contents of README.md to use as a long description
    
    Returns:
        The content of README.md as a string.
    """
    this_directory = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
        long_description = f.read()
    return long_description


# Error classes
class SetupError(RuntimeError):
    """Base error class for this setup file"""

    pass


class CommandNotFoundError(SetupError):
    """Raised when running a command cannot be found"""

    pass


class CannotRunCommandError(SetupError):
    """Raised when running a command returns with an error"""

    pass


# Helper functions
def call_prog(prog, args):
    """Call an executable and return its standard output
    Args:
        prog: (str) - the executable's name
        args: (list of str) - arguments for the executable

    Returns:
        The standard output of the executable called with the given arguments.

    Raises:
        CommandNotFoundError if the executable cannot be found.
        CannotRunCommandError if the executable did not yield an exit code of zero.
    """
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
    """Obtain system library dirs depending on the architecture.

    This may not be necessary depending on the compiler used. If the program "uname" is
    not installed, the system's architecture cannot be determined and a standard set of
    directories is returned.

    Returns:
        A list of strings containing the paths of the system libraries.
    """
    dirs = ["/usr/lib/"]
    try:
        arch = call_prog(UNAME, ["-m"])
    except CommandNotFoundError:
        arch = None
    if arch is not None:
        dirs.append("/usr/lib/{}-linux-gnu/".format(arch))
    return dirs


def get_sys_include_dirs():
    """Obtain system include dirs depending on the architecture.

    This may not be necessary depending on the compiler used. If the program "uname" is
    not installed, the system's architecture cannot be determined and a standard set of
    directories is returned.

    Returns:
        A list of strings containing the system include paths.
    """
    dirs = ["/usr/include/"]
    try:
        arch = call_prog(UNAME, ["-m"])
    except CommandNotFoundError:
        arch = None
    if arch is not None:
        dirs.append("/usr/include/{}-linux-gnu/".format(arch))
    return dirs


class build(build_orig):
    """A class that allows re-ordering the build process.

    Re-order the commands build_py and build_ext so that build_ext is run first. That
    way, swig-generated Python files will be available for the build_py step.

    Taken from https://stackoverflow.com/questions/50239473/building-a-module-with-setuptools-and-swig
    """

    def finalize_options(self):
        super().finalize_options()
        condition = lambda l: l[0] == "build_ext"
        t1, t2 = itertools.tee(self.sub_commands)
        rest, sub_build_ext = (
            itertools.filterfalse(condition, t1),
            filter(condition, t2),
        )
        self.sub_commands[:] = list(sub_build_ext) + list(rest)


def get_ext_modules(basename):
    """Get a list of external modules depending on the desired level of functionality.

    To determine the desired level of functionality, the two environment variables
    FDINST_FULL_SURFACE_SUPPORT and FDINST_FULL_VIS_SUPPORT are evaluated. By default,
    everything is installed.
    
    If FDINST_FULL_VIS_SUPPORT is not 1, visualization will use the pyopengl package to
    use OpenGL. Effectively, this removes libGL.so as a hard dependency. 

    If FDINST_FULL_SURFACE_SUPPORT is not 1, the bare minimum required to run the
    energyscan tool will be installed. Parts of manipagg might also works but there are
    no guarantees.

    Args:
        basename: (str) - the name of the main module, used to ensure that any
            extensions are sub-modules

    Returns:
        A list of extensions to be passed to the ext_modules argument of the setup
        function.
    """

    # Minimum functionality values
    sources = [
        "FireDeamon/cpp.i",
        "src/core/parallel_generic.cpp",
        "src/core/irregular_grid_interpolation.cpp",
        "src/core/arbitrary_grid_local_minima.cpp",
        "src/core/set_procname.cpp",
    ]
    libraries = ["pthread"]
    extra_compile_args = [
        "-std=c++14",
        "-v",
        "-pedantic",
        "-Wall",
        "-Wextra",
        "-fPIC",
        "-O3",
    ]
    swig_opts = ["-c++", "-Iinclude"]

    # Maximum functionality values if so desired
    if os.environ.get("FDINST_FULL_SURFACE_SUPPORT", "1") == "1":
        sources += [
            "src/core/skin_surface_deamon.cpp",
            "src/core/electrostatic_potential_charges.cpp",
            "src/core/electron_density.cpp",
            "src/core/isosurface.cpp",
            "src/core/halfnum/angular_integral.cpp",
            "src/core/halfnum/radial_integral.cpp",
            "src/core/electrostatic_potential_orbitals.cpp",
            "src/core/constants.cpp",
            "src/core/orbital_overlap.cpp",
        ]
        libraries += ["pthread", "mpfr", "gmp", "CGAL"]
        extra_compile_args += ["-frounding-math"]
        # Define the symbol FD_FULL_SUPPORT to have the swig pre-processor include
        # everything.
        swig_opts += ["-DFD_FULL_SUPPORT"]

    # Create the main C++ extension with either min or max funtionality
    fd_cpp_ext = Extension(
        basename + "._cpp",
        sources=sources,
        include_dirs=get_sys_include_dirs() + ["include"],
        library_dirs=get_sys_lib_dirs(),
        libraries=libraries,
        language="c++",
        swig_opts=swig_opts,
        extra_compile_args=extra_compile_args,
    )
    ext_modules = [fd_cpp_ext]

    # Visualization extension if so desired
    if os.environ.get("FDINST_FULL_VIS_SUPPORT", "1") == "1":
        fd_vis_ext = Extension(
            basename + ".visualize._cpp",
            sources=["FireDeamon/visualize/cpp.i", "src/visualize/interface.cpp"],
            include_dirs=get_sys_include_dirs() + ["include"],
            library_dirs=get_sys_lib_dirs(),
            libraries=["GL"],
            language="c++",
            swig_opts=["-c++", "-Iinclude"],
            extra_compile_args=[
                "-std=c++14",
                "-v",
                "-pedantic",
                "-Wall",
                "-Wextra",
                "-fPIC",
                "-O3",
            ],
        )
        ext_modules.append(fd_vis_ext)

    return ext_modules


def main():
    # Package data
    PACKAGE_DATA = {
        "name": "FireDeamon",
        "version": 0.1,
        "description": "Misc routines for the ManipulateAggregates package",
        "author": "Torsten Sachse",
        "mail": "torsten.sachse@gmail.com",
        "url": "https://github.com/razziel89/libFireDeamon",
        "long_description": get_long_description(),
        "long_description_content_type": "text/markdown",
    }

    setup(
        ext_modules=get_ext_modules(PACKAGE_DATA["name"]),
        packages=find_packages(),
        cmdclass={"build": build},
        **PACKAGE_DATA,
    )


if __name__ == "__main__":
    main()
