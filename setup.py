from setuptools import setup, Extension
import numpy as np
import imp

version = imp.load_source('afterglowpy.version', 'version.py')

inc = [np.get_include()]
libs = []
libdirs = []

jetsources = ["afterglowpy/jetmodule.c", "afterglowpy/offaxis_struct_funcs.c",
              "afterglowpy/integrate.c", "afterglowpy/shockEvolution.c"]
jetdepends = ["afterglowpy/offaxis_struct_funcs.h",
              "afterglowpy/shockEvolution.h"]

shocksources = ["afterglowpy/shockmodule.c", "afterglowpy/shockEvolution.c"]
shockdepends = ["afterglowpy/shockEvolution.h",
                "afterglowpy/offaxis_struct_funcs.h"]

jetmodule = Extension('afterglowpy.jet', sources=jetsources, include_dirs=inc,
                      depends=jetdepends)
shockmodule = Extension('afterglowpy.shock', sources=shocksources,
                        include_dirs=inc, depends=shockdepends)

setup(
    name='afterglowpy',
    version=version.version,
    description='GRBAfterglowModels',
    packages=['afterglowpy'],
    ext_modules=[jetmodule, shockmodule],
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
