r'''
Utility functionality
=====================

.. currentmodule:: micronota.util

This module (:mod:`micronota.util`) provides various utility functionality,
including config config, unit-testing convenience function.

'''

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import shutil
from os import remove
from os.path import exists, isdir, join, abspath, dirname, basename, splitext
from urllib.request import urlopen
from unittest import TestCase
from sqlite3 import connect
from inspect import stack


def _overwrite(fp, overwrite=False, append=False):
    if exists(fp):
        if overwrite:
            if isdir(fp):
                shutil.rmtree(fp)
            else:
                remove(fp)
        elif append:
            return
        else:
            raise FileExistsError('The file path %s exists.' % fp)


def _download(src, dest, **kwargs):
    _overwrite(dest, **kwargs)
    with urlopen(src) as i_f, open(dest, 'wb') as o_f:
        shutil.copyfileobj(i_f, o_f)


def _get_named_data_path(fname):
    # get caller's file path
    caller_fp = abspath(stack()[1][1])
    d = dirname(caller_fp)
    # remove file suffix and prefix of "test_"
    name = splitext(basename(caller_fp))[0][5:]
    return join(d, 'data', name, fname)


class _DBTest(TestCase):
    def _test_eq_db(self, db1, db2):
        '''Test if two database files have the same contents.'''
        with connect(db1) as o, connect(self.exp_db_fp) as e:
            # compare each table
            for table, in o.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"):
                co = o.cursor()
                co.execute('SELECT * from %s' % table)
                ce = e.cursor()
                ce.execute('SELECT * from %s' % table)
                self.assertCountEqual(co.fetchall(), ce.fetchall())
