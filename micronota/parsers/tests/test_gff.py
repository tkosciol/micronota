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
        self.attribute1 = {
                            'SEQID': 'gi|15282445|ref|NC_000918.1|',
                            'SOURCE': 'minced:0.2.0',
                            'TYPE': 'CRISPR',
                            'SCORE': 5,
                            'STRAND': '.',
                            'PHASE': '.',
                            'ATTR': ('ID', 'CRISPR1')
                            }
        self.attribute2 = {
                            'SEQID': 'gi|15282445|ref|NC_000918.1|',
                            'SOURCE': 'minced:0.2.0',
                            'TYPE': 'CRISPR',
                            'SCORE': 4,
                            'STRAND': '.',
                            'PHASE': '.',
                            'ATTR': ('ID', 'CRISPR2')
                            }
        self.interval1 = (156460, 156767)
        self.interval2 = (244561, 244791)
        self.simple_exp_im = IntervalMetadata(features={
                                              Feature(**self.attribute1):
                                              self.interval1})
        self.full_exp_im = IntervalMetadata(features={
                                            Feature(**self.attribute1):
                                            self.interval1,
                                            Feature(**self.attribute2):
                                            self.interval2})

        # self.sequence_fp = get_data_path('sequence.gff')
        # self.sequence_exp = \
        #     (
        #      'cttctgggcgtacccgattctcggagaacttgccgcaccattccgccttg'
        #      'tgttcattgctgcctgcatgttcattgtctacctcggctacgtgtggcta'
        #      'tctttcctcggtgccctcgtgcacggagtcgagaaaccaaagaacaaaaa'
        #      'aagaaattaaaatatttattttgctgtggtttttgatgtgtgttttttat'
        #      'aatgatttttgatgtgaccaattgtacttttcctttaaatgaaatgtaat'
        #      'cttaaatgtatttccgacgaattcgaggcctgaaaagtgtgacgccattc', {
        #          'SEQID': 'ctg123',
        #          'SOURCE': '.',
        #          'TYPE': 'exon',
        #          'START': 130,
        #          'END': 150,
        #          'SCORE': '.',
        #          'STRAND': '+',
        #          'PHASE': '.',
        #          'ATTR': ('ID', 'exon00001')
        #      })


class SnifferTests(TestCase):
    def setUp(self):
        self.positive_fps = list(map(get_data_path, [
            'simple.gff',
            'full.gff',
            'sequence.gff']))
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
        obs_simple = _gff_to_metadata(self.simple_fp, rec_num=1)
        exp_simple = self.simple_exp_im
        obs_full = _gff_to_metadata(self.multi_fp, rec_num=2)
        exp_full = self.full_exp_im
        self.assertEqual(obs_simple, exp_simple)
        self.assertEqual(obs_full, exp_full)

    def test_gff_to_generator(self):
        exp = self.simple_exp_im
        for obs in _gff_to_generator(self.simple_fp):
            self.assertEqual(obs, exp)

if __name__ == '__main__':
    main()
