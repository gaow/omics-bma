#! /usr/bin/env python3
__author__ = "Gao Wang"
__copyright__ = "Copyright 2016, Stephens lab"
__email__ = "gaow@uchicago.edu"
__license__ = "MIT"
__version__ = "0.1.0"
from setuptools import Extension, setup
from distutils import sysconfig
import os, subprocess, platform, glob

# use ccache to speed up build
try:
    if subprocess.call(['ccache'], stderr = open(os.devnull, "w")):
        os.environ['CC'] = 'ccache clang -Qunused-arguments' if platform.system() == 'Darwin' else 'ccache gcc'
except OSError:
    pass

# set c flags
sysconfig.get_config_vars()["CFLAGS"] = \
sysconfig.get_config_vars()["CFLAGS"].replace("-Wstrict-prototypes", "")

# set required files
PACKAGE = 'OmicsBMA'
MODULE = 'pyeqtlbma'
EXTERN_LIB_PATH = os.path.abspath("external")
## lazy setup: simply include everything from external libraries
INCLUDE_DIR = list(set([x[0] for x in os.walk(EXTERN_LIB_PATH + "/eqtlbma/src") if len([y for y in os.path.split(x[0]) if y.startswith('.')]) == 0]))
INCLUDE_DIR.extend([EXTERN_LIB_PATH + "/gsl/include", "pyeqtlbma"])
# LIB_DIR = [EXTERN_LIB_PATH + "/gsl/lib"]
LIB = ["stdc++"]
SOURCE = glob.glob(EXTERN_LIB_PATH + "/eqtlbma/**/*.cpp", recursive=True)
SOURCE += glob.glob(EXTERN_LIB_PATH + "/eqtlbma/**/*.cc", recursive=True)
SOURCE += glob.glob(EXTERN_LIB_PATH + "/eqtlbma/**/*.c", recursive=True)
SOURCE += ['src/pyeqtlbma/pyeqtlbma_wrap.cxx', 'src/pyeqtlbma/libeqtlbma.cpp',
           'src/pyeqtlbma/pyeqtlbma.cpp']
SOURCE = [x for x in SOURCE if not os.path.split(x)[-1] in
          ['eqtlbma_bf.cpp', 'eqtlbma_hm.cpp', 'eqtlbma_avg_bfs.cpp', 'main.c', 'bgzip.c']]
## compiler configs
COMPILE_ARGS = ["-O3", "--std=c++11", "-fPIC", "-fopenmp", "-DVERSION='OmicsBMA'", "-D_FILE_OFFSET_BITS=64", "-DHAVE_SYS_TYPES_H=1", "-DHAVE_SYS_STAT_H=1", "-DHAVE_STDLIB_H=1", "-DHAVE_STRING_H=1", "-DHAVE_MEMORY_H=1", "-DHAVE_STRINGS_H=1", "-DHAVE_INTTYPES_H=1", "-DHAVE_STDINT_H=1", "-DHAVE_UNISTD_H=1", "-DHAVE_DLFCN_H=1", "-DHAVE_LIBZ=1", "-DHAVE_LIBGSLCBLAS=1"]
LINK_ARGS = ["-lgomp", "-lm", "-lz", "-Wl,--no-as-needed"] + [EXTERN_LIB_PATH + x for x in ["/gsl/lib/libgsl.a", "/gsl/lib/libgslcblas.a"]]

# compile
EQTLBMA_MODULE = Extension("{0}.{1}._{1}".format(PACKAGE, 'pyeqtlbma'),
                            SOURCE,
                            include_dirs = INCLUDE_DIR,
                            # library_dirs = LIB_DIR,
                            libraries = LIB,
                            extra_link_args = LINK_ARGS,
                            extra_compile_args = COMPILE_ARGS,
                            swig_opts=['-threads']
                            )

setup(name        = PACKAGE,
      description = "A machinery to interrogate genomics / transcriptomics data using Bayesian Model Averaging",
      author      = "Gao Wang",
      version     = "0.1.0",
      packages    = [PACKAGE, "{}.pyeqtlbma".format(PACKAGE)],
      package_dir = {PACKAGE: "src", "{}.pyeqtlbma".format(PACKAGE): "src/pyeqtlbma"},
      ext_modules = [EQTLBMA_MODULE]
      )
