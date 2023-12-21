import setuptools
from setuptools import find_packages
from numpy.distutils.core import setup, Extension


fgeqdsk = Extension(
    name="fgeqdsk",
    sources=["fgeqdsk.f90"],
    f2py_options=["--quiet"],
)


setup(
    name="pgeqdsk",
    install_requires=[
        "numpy",
        "scipy",
    ],
    py_modules=["pgeqdsk"],
    ext_modules=[fgeqdsk],
)
