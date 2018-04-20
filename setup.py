from setuptools import setup, find_packages, Extension
import numpy as np
import imp

version = imp.load_source('grbpy.version', 'version.py')

inc = [np.get_include()]
libs = []
libdirs = []

jetsources = ["grbpy/jetmodule.c", "grbpy/offaxis_struct_funcs.c",
                "grbpy/integrate.c", "grbpy/shockEvolution.c"]
jetdepends = ["grbpy/offaxis_struct_funcs.h", "grbpy/shockEvolution.h"]

shocksources = ["grbpy/shockmodule.c", "grbpy/shockEvolution.c"]
shockdepends = ["grbpy/shockEvolution.h", "grbpy/offaxis_struct_funcs.h"]

jetmodule = Extension('grbpy.jet', sources=jetsources, include_dirs=inc,
                                    depends=jetdepends)
shockmodule = Extension('grbpy.shock', sources=shocksources, include_dirs=inc,
                                    depends=shockdepends)

setup(
    name='grbpy',
    version=version.version,
    description='GRBAfterglowModels',
    packages=['grbpy'],
    ext_modules = [jetmodule, shockmodule],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: C",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy"],
    install_requires=['numpy>=1.10', 'scipy>=0.14'],
    extras_require={
        'docs': ['numpydoc']
        }
    )
