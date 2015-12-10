#! /usr/bin/env python3
# setup.py
# Gao Wang (c) 2015
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
EXTERN_LIB_PATH = os.path.abspath("../../external")
## lazy setup: simply include everything from external libraries
INCLUDE_DIR = list(set([x[0] for x in os.walk(EXTERN_LIB_PATH + "/eqtlbma/src") if len([y for y in os.path.split(x[0]) if y.startswith('.')]) == 0]))
INCLUDE_DIR.extend([EXTERN_LIB_PATH + "/gsl/include", "./"])
# LIB_DIR = [EXTERN_LIB_PATH + "/gsl/lib"]
LIB = ["stdc++"]
SOURCE = glob.glob(EXTERN_LIB_PATH + "/eqtlbma/**/*.cpp", recursive=True)
SOURCE += glob.glob(EXTERN_LIB_PATH + "/eqtlbma/**/*.cc", recursive=True)
SOURCE += glob.glob(EXTERN_LIB_PATH + "/eqtlbma/**/*.c", recursive=True)
SOURCE += ['pyeqtlbma_wrap.cxx', 'libeqtlbma.cpp', 'pyeqtlbma.cpp']
SOURCE = [x for x in SOURCE if not os.path.split(x)[-1] in ['eqtlbma_bf.cpp', 'eqtlbma_hm.cpp', 'eqtlbma_avg_bfs.cpp', 'main.c', 'bgzip.c']]
## compiler configs
COMPILE_ARGS = ["-O3", "--std=c++11", "-fPIC", "-fopenmp", "-DVERSION='OmicsBMA'", "-D_FILE_OFFSET_BITS=64", "-DHAVE_SYS_TYPES_H=1", "-DHAVE_SYS_STAT_H=1", "-DHAVE_STDLIB_H=1", "-DHAVE_STRING_H=1", "-DHAVE_MEMORY_H=1", "-DHAVE_STRINGS_H=1", "-DHAVE_INTTYPES_H=1", "-DHAVE_STDINT_H=1", "-DHAVE_UNISTD_H=1", "-DHAVE_DLFCN_H=1", "-DHAVE_LIBZ=1", "-DHAVE_LIBGSLCBLAS=1"]
LINK_ARGS = ["-lgomp", "-lm", "-lz", "-Wl,--no-as-needed"] + [EXTERN_LIB_PATH + x for x in ["/gsl/lib/libgsl.a", "/gsl/lib/libgslcblas.a"]]

# compile
cpp_methods_ext = Extension("{}._{}".format(PACKAGE, MODULE),
                            SOURCE,
                            include_dirs = INCLUDE_DIR,
                            # library_dirs = LIB_DIR,
                            libraries = LIB,
                            extra_link_args = LINK_ARGS,
                            extra_compile_args = COMPILE_ARGS,
                            swig_opts=['-threads']
                            )

setup(name        = "{}.{}".format(PACKAGE, MODULE),
      description = "Python interface from eQTLBMA package (the omics-bma branch)",
      author      = "Gao Wang",
      version     = "1.0",
      packages    = [PACKAGE],
      package_dir = {PACKAGE: "."},
      ext_modules = [cpp_methods_ext]
      )
