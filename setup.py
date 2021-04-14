#!/usr/bin/env python3
import site
from setuptools import setup, Extension
import numpy as np
# import imp

# PEP 517 Workaround for edittable user installs
site.ENABLE_USER_SITE = True

# version = imp.load_source('afterglowpy.version', 'afterglowpy/version.py')

version = {}
with open("afterglowpy/version.py", "r") as f:
    exec(f.read(), version)

with open("README.md", "r") as f:
    long_description = f.read()

inc = [np.get_include()]
libs = []
libdirs = []

jetsources = ["afterglowpy/jetmodule.c", "afterglowpy/offaxis_struct_funcs.c",
              "afterglowpy/integrate.c", "afterglowpy/shockEvolution.c",
              "afterglowpy/interval.c"]
jetdepends = ["afterglowpy/offaxis_struct_funcs.h",
              "afterglowpy/shockEvolution.h", "afterglowpy/interval.h"]

shocksources = ["afterglowpy/shockmodule.c", "afterglowpy/shockEvolution.c"]
shockdepends = ["afterglowpy/shockEvolution.h",
                "afterglowpy/offaxis_struct_funcs.h"]

jetmodule = Extension('afterglowpy.jet', sources=jetsources, include_dirs=inc,
                      depends=jetdepends, extra_compile_args=['-Wall'])
shockmodule = Extension('afterglowpy.shock', sources=shocksources,
                        include_dirs=inc, depends=shockdepends,
                        extra_compile_args=['-Wall'])

setup(
    name='afterglowpy',
    version=version['__version__'],
    author="Geoffrey Ryan",
    author_email="gsryan@umd.edu",
    description='GRB Afterglow Models',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='MIT',
    url='https://github.com/geoffryan/afterglowpy',
    packages=['afterglowpy'],
    ext_modules=[jetmodule, shockmodule],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy"],
    install_requires=['numpy>=1.10', 'scipy>=0.14'],
    extras_require={
        'docs': ['numpydoc']
        },
    project_urls={
        "Source Code": "https://github.com/geoffryan/afterglowpy",
        "Documentation": "https://afterglowpy.readthedocs.io"}
    )
