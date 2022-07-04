"""
    Setup file for pybiotk.
    Use setup.cfg to configure your project.

    This file was generated with PyScaffold 4.2.3.
    PyScaffold helps you to put up the scaffold of your new Python project.
    Learn more under: https://pyscaffold.org/
"""
from setuptools import setup, Extension
from Cython.Build import cythonize
from glob import glob


ext_modules = [
    Extension("pybiotk.bx.bitset", ["src/pybiotk/bx/bitset.pyx", "src/pybiotk/bx/binBits.c", "src/pybiotk/bx/bits.c", "src/pybiotk/bx/common.c"]),
    Extension("pybiotk.bx.cluster", ["src/pybiotk/bx/cluster.c", "src/pybiotk/bx/cluster.pyx"]),
    Extension("pybiotk.bx.intersection", ["src/pybiotk/bx/intersection.pyx"])
]

if __name__ == "__main__":
    try:
        setup(name="pybiotk",
              version="1.0.1",
              use_scm_version={"version_scheme": "no-guess-dev"},
              ext_modules=cythonize(ext_modules, build_dir="build"),
              scripts=glob("rscripts/*.R") + glob("scripts/*sh"),)
    except:  # noqa
        print(
            "\n\nAn error occurred while building the project, "
            "please ensure you have the most updated version of setuptools, "
            "setuptools_scm and wheel with:\n"
            "   pip install -U setuptools setuptools_scm wheel\n\n"
        )
        raise
