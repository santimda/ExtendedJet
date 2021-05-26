from setuptools import find_packages
import subprocess
from numpy.distutils.core import setup, Extension

# To generate the .mod although there might be a better way
subprocess.call("make", cwd="./extendedjetpy")

ext1 = Extension(name='extendedjetpy.main',
                 sources=['extendedjetpy/controles.f90','extendedjetpy/dist.f90','extendedjetpy/global.f90','extendedjetpy/radsyn.f90','extendedjetpy/constants.f90'],
                 f2py_options=['--quiet'], include_dirs = ['extendedjetpy'], library_dirs = ['extendedjetpy']
                )

setup(name="ExtendedJetPy",
      version="0.0.1",
      author = "S. del Palacio, E. Gutierrez, FG. Lopez Armengol, L. Combi",
      #package_dir={"jetcode": "src"},
      #license='MIT',
      #url='https://gitlab.com/PI-UoG-RelAstro/kilonova-code',
      #install_requires=['numpy>=1.10', 'scipy>=0.14'],
      packages=['extendedjetpy'], #find_packages(where="src"),
      ext_modules=[ext1])
