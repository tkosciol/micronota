# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from tempfile import mkdtemp
from os import getcwd
from shutil import rmtree
from os.path import join
from unittest import TestCase, main
from functools import partial
from skbio.util import get_data_path
# from burrito.util import ApplicationError

from micronota.bfillings.signalp import SignalP, predict_signal


class SignalPTests(TestCase):
    def setUp(self):
        self.temp_dir = mkdtemp()
        self.get_signalp_path = partial(
            get_data_path, subfolder=join('data', 'signalp'))

        # taken from SignalP test files
        self.positive_fps = list(map(self.get_signalp_path,
                                 ['euk10.fsa',
                                  'euk10.fsa',
                                  'euk10.fsa',
                                  'euk10.fsa']))
        self.negative_fps = list(map(get_data_path, [
            # SignalP hangs on whitespace_only
            'whitespace_only',
            'empty']))
        self.positive_params = [
            {'-t': 'euk'},
            {'-t': 'euk', '-f': 'summary'},
            {'-t': 'euk', '-f': 'long', '-m': 'foo'},
            {'-n': 'foo', '-f': 'all'}]
        self.positive_prefix = 'euk10'
        self.positive_suffices = [
            {'out': 'sp_short'},
            {'out': 'sp_summary'},
            {'out': 'sp_long', '-m': 'fasta'},
            {'out': 'sp_all', '-g': 'gff'}]

    def test_base_command(self):
        c = SignalP()
        self.assertEqual(
            c.BaseCommand,
            'cd "%s/"; %s' % (getcwd(), c._command))

    # SignalP accepts any input and does not raise errors
    # def test_predict_signal_wrong_input(self):
    #     for fp in self.negative_fps:
    #         with self.assertRaisesRegex(
    #                 ApplicationError,
    #                 r'Error constructing CommandLineAppResult.'):
    #             predict_signal(fp, self.temp_dir, 'foo')

    def test_predict_signal(self):
        for in_fp, params, suffix in zip(self.positive_fps,
                                         self.positive_params,
                                         self.positive_suffices):
            prefix = self.positive_prefix
            res = predict_signal(in_fp, self.temp_dir, prefix, params=params)
            self.assertEqual(res['ExitStatus'], 0)

            for suff_keys in suffix.keys():
                suff_val = suffix[suff_keys]
                fp = self.get_signalp_path('.'.join([prefix, suff_val]))
                if suff_val.startswith('sp_'):
                    suff_val = 'StdOut'

                with open(fp) as f:
                    self.assertEqual(
                        # skip comment lines
                        [i for i in f.readlines()
                         if not i.startswith('#')],
                        [j for j in res[suff_val].readlines()
                         if not j.startswith('#')])
                    res[suff_val].close()
            res['StdOut'].close()
            res['StdErr'].close()

    def tearDown(self):
        # remove the tempdir and contents
        rmtree(self.temp_dir)

if __name__ == '__main__':
    main()
