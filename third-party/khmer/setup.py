#! /usr/bin/env python2
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
""" Setup for khmer project. """

import ez_setup
ez_setup.use_setuptools(version="3.4.1")

import os
import sys
from os import listdir as os_listdir
from os.path import join as path_join

from setuptools import setup
from setuptools import Extension
from setuptools.command.build_ext import build_ext as _build_ext
from distutils.spawn import spawn
from distutils.sysconfig import get_config_vars
from distutils.dist import Distribution
from distutils.errors import DistutilsPlatformError

import versioneer
versioneer.VCS = 'git'
versioneer.versionfile_source = 'khmer/_version.py'
versioneer.versionfile_build = 'khmer/_version.py'
versioneer.tag_prefix = 'v'  # tags are like v1.2.0
versioneer.parentdir_prefix = '.'
CMDCLASS = versioneer.get_cmdclass()

# strip out -Wstrict-prototypes; a hack suggested by
# http://stackoverflow.com/a/9740721
# proper fix coming in http://bugs.python.org/issue1222585
# numpy has a "nicer" fix:
# https://github.com/numpy/numpy/blob/master/numpy/distutils/ccompiler.py
OPT = get_config_vars('OPT')[0]
os.environ['OPT'] = " ".join(
    flag for flag in OPT.split() if flag != '-Wstrict-prototypes'
)

# We bundle tested versions of zlib & bzip2. To use the system zlib and bzip2
# change setup.cfg or use the `--libraries z,bz2` parameter which will make our
# custom build_ext command strip out the bundled versions.

ZLIBDIR = 'third-party/zlib'
BZIP2DIR = 'third-party/bzip2'

BUILD_DEPENDS = []
BUILD_DEPENDS.extend(path_join("lib", bn + ".hh") for bn in [
    "khmer", "khmer_config", "kmer_hash", "hashtable", "counting",
    "hashbits", "labelhash"])

SOURCES = ["khmer/_khmermodule.cc"]
SOURCES.extend(path_join("lib", bn + ".cc") for bn in [
    "khmer_config", "thread_id_map", "trace_logger", "perf_metrics",
    "read_parsers", "kmer_hash", "hashtable", "hashbits", "labelhash",
    "counting", "subset", "read_aligner"])

EXTRA_COMPILE_ARGS = ['-O3']

if sys.platform == 'darwin':
    EXTRA_COMPILE_ARGS.extend(['-arch', 'x86_64'])  # force 64bit only builds

EXTENSION_MOD_DICT = \
    {
        "sources": SOURCES,
        "extra_compile_args": EXTRA_COMPILE_ARGS,
        "depends": BUILD_DEPENDS,
        "language": "c++",
        "define_macros": [("VERSION", versioneer.get_version()), ],
    }

EXTENSION_MOD = Extension("khmer._khmermodule",  # pylint: disable=W0142
                          ** EXTENSION_MOD_DICT)
SCRIPTS = []
SCRIPTS.extend([path_join("scripts", script)
                for script in os_listdir("scripts")
                if script.endswith(".py")])

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: C++",
    "Programming Language :: Python :: 2.7",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
if "-rc" in versioneer.get_version():
    CLASSIFIERS.append("Development Status :: 4 - Beta")
else:
    CLASSIFIERS.append("Development Status :: 5 - Production/Stable")

SETUP_METADATA = \
    {
        "name": "khmer",
        "version": versioneer.get_version(),
        "description": 'khmer k-mer counting library',
        "long_description": open("README.rst").read(),
        "author": 'Michael R. Crusoe, Greg Edvenson, Jordan Fish,'
        ' Adina Howe, Luiz Irber, Eric McDonald, Joshua Nahum, Kaben Nanlohy,'
        ' Humberto Ortiz-Zuazaga, Jason Pell, Jared Simpson, Camille Scott,'
        ' Ramakrishnan Rajaram Srinivasan, Qingpeng Zhang, and C. Titus Brown',
        "author_email": 'khmer-project@idyll.org',
        # "maintainer": 'Michael R. Crusoe', # this overrides the author field
        # "maintainer_email": 'mcrusoe@msu.edu', # so don't include it
        # http://docs.python.org/2/distutils/setupscript.html
        # additiona-meta-data note #3
        "url": 'http://ged.msu.edu/',
        "packages": ['khmer', 'khmer.tests'],
        "package_dir": {'khmer.tests': 'tests'},
        "install_requires": ['screed >= 0.7.1'],
        "extras_require": {':python_version=="2.6"': ['argparse>=1.2.1'],
                           'docs': ['sphinx', 'sphinxcontrib-autoprogram'],
                           'tests': ['nose >= 1.0']},
        "scripts": SCRIPTS,
        "ext_modules": [EXTENSION_MOD, ],
        # "platforms": '', # empty as is conveyed by the classifiers below
        # "license": '', # empty as is conveyed by the classifier below
        "include_package_data": True,
        "zip_safe": False,
        "classifiers": CLASSIFIERS
    }


class KhmerBuildExt(_build_ext):  # pylint: disable=R0904
    """Specialized Python extension builder for khmer project.

    Only run the library setup when needed, not on every invocation.

    Also strips out the bundled zlib and bzip2 libraries if
    `--libraries z,bz2` is specified or the equivalent is in setup.cfg
    """

    def run(self):
        if "%x" % sys.maxsize != '7fffffffffffffff':
            raise DistutilsPlatformError("%s require 64-bit operating system" %
                                         SETUP_METADATA["packages"])

        if "z" not in self.libraries:
            zcmd = ['bash', '-c', 'cd ' + ZLIBDIR + ' && ( test Makefile -nt'
                    ' configure || bash ./configure --static ) && make -f '
                    'Makefile.pic PIC']
            spawn(cmd=zcmd, dry_run=self.dry_run)
            self.extensions[0].extra_objects.extend(
                path_join("third-party", "zlib", bn + ".lo") for bn in [
                    "adler32", "compress", "crc32", "deflate", "gzclose",
                    "gzlib", "gzread", "gzwrite", "infback", "inffast",
                    "inflate", "inftrees", "trees", "uncompr", "zutil"])
        if "bz2" not in self.libraries:
            bz2cmd = ['bash', '-c', 'cd ' + BZIP2DIR + ' && make -f '
                      'Makefile-libbz2_so all']
            spawn(cmd=bz2cmd, dry_run=self.dry_run)
            self.extensions[0].extra_objects.extend(
                path_join("third-party", "bzip2", bn + ".o") for bn in [
                    "blocksort", "huffman", "crctable", "randtable",
                    "compress", "decompress", "bzlib"])
        _build_ext.run(self)

CMDCLASS.update({'build_ext': KhmerBuildExt})

_DISTUTILS_REINIT = Distribution.reinitialize_command


def reinitialize_command(self, command, reinit_subcommands):
    '''
    Monkeypatch distutils.Distribution.reinitialize_command() to match behavior
    of Distribution.get_command_obj()

    This fixes issues with 'pip install -e' and './setup.py nosetests' not
    respecting the setup.cfg configuration directives for the build_ext command
    '''
    cmd_obj = _DISTUTILS_REINIT(self, command, reinit_subcommands)
    options = self.command_options.get(command)
    if options:
        self._set_command_options(  # pylint: disable=protected-access
            cmd_obj, options)
    return cmd_obj
Distribution.reinitialize_command = reinitialize_command

# pylint: disable=W0142
setup(cmdclass=CMDCLASS,
      **SETUP_METADATA)
