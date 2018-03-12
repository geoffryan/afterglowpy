from setuptools import setup, find_packages, Extension
import numpy as np
import imp

version = imp.load_source('grbpy.version', 'version.py')

inc = [np.get_include()]
libs = []
libdirs = []

sources = ["grbpy/jetmodule.c", "grbpy/offaxis_struct_funcs.c",
                "grbpy/integrate.c"]

depends = ["grbpy/offaxis_struct_funcs.h"]

module = Extension('grbpy.jet', sources=sources, include_dirs=inc)
                        #libraries=libs, library_dirs=libdirs, depends=depends)

setup(
    name='grbpy',
    version=version.version,
    description='GRBAfterglowModels',
    packages=['grbpy'],
    ext_modules = [module],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: C",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy"],
    install_requires=['numpy>=1.10'],
    extras_require={
        'docs': ['numpydoc']
        }
    )
