from setuptools import find_packages
import subprocess
from numpy.distutils.core import setup, Extension
import glob
import site
## To generate the .mod although there might be a better way
subprocess.call("make", cwd="./extendedjetpy")

site.ENABLE_USER_SITE = True

main = Extension(name='extendedjetpy.main',
        extra_compile_args=['-funroll-loops','-Wno-tabs','-fPIC','-ffree-line-length-none','-fbackslash','-march=native','-fmax-errors=1', '-O3','-fbackslash','-ffast-math'],
        sources = ['./extendedjetpy/controls.f90' ,'./extendedjetpy/constants.f90', './extendedjetpy/system_parameters.f90' ,'./extendedjetpy/globals.f90' ,'./extendedjetpy/initialize.f90' ,'./extendedjetpy/trajectory.f90' ,'./extendedjetpy/evolution.f90', './extendedjetpy/losses.f90','./extendedjetpy/edistribution.f90' ,'./extendedjetpy/synrad.f90' ,'./extendedjetpy/icrad.f90', './extendedjetpy/sscrad.f90', './extendedjetpy/main.f90'],
                 include_dirs = ['extendedjetpy'], library_dirs = ['extendedjetpy']
                )

setup(name="ExtendedJetPy",
      version="0.0.1",
      author = "S. del Palacio, E. Gutierrez, FG. Lopez Armengol, L. Combi",
      #package_dir={"jetcode": "src"},
      #license='MIT',
      #url='https://gitlab.com/PI-UoG-RelAstro/kilonova-code',
      #install_requires=['numpy>=1.10', 'scipy>=0.14'],
      packages=['extendedjetpy'], #find_packages(where="src"),
      ext_modules=[main])
