from distutils.core import setup
from Cython.Build import cythonize

setup(name='KeySplitter',
      ext_modules=cythonize("key_split_cy.pyx"))