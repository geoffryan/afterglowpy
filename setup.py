import os
import sys
import subprocess
from setuptools import setup, find_packages, Extension
import numpy.distutils.misc_util
import imp

version = imp.load_source('grbpy.version', 'version.py')

#gsl_dir = os.system("gsl-config --prefix")
gsl_dir = subprocess.check_output(["gsl-config", "--prefix"])[:-1]
if sys.version_info >= (3,0):
    gsl_dir = gsl_dir.decode("utf-8")

inc = numpy.distutils.misc_util.get_numpy_include_dirs()
inc.append(gsl_dir + "/include")
libs = ["gsl"]
libdirs = [gsl_dir + "/lib"]

module = Extension("_grbpy", ["grbpy/_grbpy.c"], include_dirs=inc,
                        libraries=libs, library_dirs=libdirs)

setup(
    name='grbpy',
    version=version.version,
    description='GRBAfterglowModels',
    packages=find_packages(),
    ext_modules = [module],
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy"],
    install_requires=['numpy>=1.10'],
    extras_require={
        'docs': ['numpydoc']
        }
    )



