from setuptools import setup
from Cython.Build import cythonize

setup(
    name="cftpy",
    ext_modules=cythonize("cftpy/cythonft.pyx"),
    packages=["cftpy"],
)