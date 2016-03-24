# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

from skbio.util import get_data_path
from skbio.metadata import IntervalMetadata, Feature

from micronota.parsers.gff import (
    _gff_sniffer, _gff_to_metadata, _gff_to_generator)


class GffIOTests(TestCase):
    def setUp(self):
        self.multi_fp = get_data_path('full.gff')
        self.simple_fp = get_data_path('simple.gff')
        self.attributes = {
                            'SEQID': 'gi|15282445|ref|NC_000918.1|',
                            'SOURCE': 'minced:0.2.0',
                            'TYPE': 'CRISPR',
                            'SCORE': 5,
                            'STRAND': '.',
                            'PHASE': '.',
                            'ATTR': ('ID', 'CRISPR1')
                            }
        self.intervals = (156460, 156767)
        self.single_exp_im = IntervalMetadata()
        self.single_exp_im.add(Feature(**self.attributes), self.intervals)


class SnifferTests(TestCase):
    def setUp(self):
        self.positive_fps = list(map(get_data_path, [
            'simple.gff',
            'full.gff']))
        self.negative_fps = list(map(get_data_path, [
            'blank.sam',
            'wrong.gff']))

    def test_positive(self):
        for fp in self.positive_fps:
            self.assertEqual(_gff_sniffer(fp), (True, {}))
        for fp in self.negative_fps:
            self.assertEqual(_gff_sniffer(fp), (False, {}))


class ReaderTests(GffIOTests):
    def test_gff_to_metadata(self):
        obs = _gff_to_metadata(self.multi_fp)
        exp = self.single_exp_im
        self.assertEqual(obs, exp)

    def test_gff_to_generator(self):
        exp = self.single_exp_im
        for obs in _gff_to_generator(self.simple_fp):
            self.assertEqual(obs, exp)

if __name__ == '__main__':
    main()
