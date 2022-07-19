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
    Extension("pybiotk.bx.bitset", ["src/bx/bitset.pyx", "src/bx/binBits.c", "src/bx/bits.c", "src/bx/common.c"]),
    Extension("pybiotk.bx.cluster", ["src/bx/cluster.c", "src/bx/cluster.pyx"]),
    Extension("pybiotk.bx.intersection", ["src/bx/intersection.pyx"])
]

if __name__ == "__main__":
    try:
        setup(use_scm_version={"version_scheme": "no-guess-dev", "local_scheme": "no-local-version"},
              ext_modules=cythonize(ext_modules, build_dir="build"),
              scripts=glob("rscripts/*.R") + glob("scripts/*sh"),)
    except:
        print(
            "\n\nAn error occurred while building the project, "
            "please ensure you have the most updated version of setuptools, "
            "setuptools_scm and wheel with:\n"
            "   pip install -U setuptools setuptools_scm wheel\n\n"
        )
        raise
