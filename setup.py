from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
setup(
    name = "swalign",
    description = "cython based swaligner setup to work with python3",
    version = '0.1',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension('swalign.cswalign',
                             ['swalign/cswalign.pyx'],include_dirs=[numpy.get_include()])],
    packages = ["swalign"],
    requires = ['cython (>=0.18)'],
)
                             
                            
