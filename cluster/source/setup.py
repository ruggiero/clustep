from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("rejection", ["source/rejection.pyx"])]

setup(
  name = 'Rejection algorithm',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
