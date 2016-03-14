r'''
GFF3 Parser
=====================

GFF is a standard file format for storing genomic features in a text file.
GFF stands for Generic Feature Format. GFF files are plain text, 9 column,
tab-delimited files [#]_.
GFF3 format is a flat tab-delimited file.
The first line of the file is a comment that identifies the format and version.
This is followed by a series of data lines, each one of which corresponds
to an annotation.
The 9 columns of the annotation section are as follows:

    +--------------------------------+----------------------------------+
    | SEQID - ID of the landmark used| (begins each entry; 1 per entry) |
    +--------------------------------+----------------------------------+
    | SOURCE - algorithms used       | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | TYPE - type of the feature     | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | START - start of the feature   | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | END - end of the feature       | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | SCORE - floating point score   | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | STRAND - +/-/./?               | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | PHASE - only for TYPE="CDS"    | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | ATTR - feature attributes (tag)| 1 per entry                      |
    +--------------------------------+----------------------------------+


Column 9 (ATTR) tags have predefined meanings:

ID
Indicates the unique identifier of the feature.
IDs must be unique within the scope of the GFF file.

Name
Display name for the feature. This is the name to be displayed to the user.
Unlike IDs, there is no requirement that the Name be unique within the file.

Alias
A secondary name for the feature.
It is suggested that this tag be used whenever a secondary identifier for
the feature is needed, such as locus names and accession numbers.
Unlike ID, there is no requirement that Alias be unique within the file.

Parent
Indicates the parent of the feature.
A parent ID can be used to group exons into transcripts, transcripts into genes
and so forth. A feature may have multiple parents. Parent can *only* be used
to indicate a partof relationship.

Target
Indicates the target of a nucleotide-to-nucleotide or
protein-to-nucleotide alignment. The format of the value is "target_id start
end [strand]", where strand is optional and may be "+" or "-".
If the target_id contains spaces, they must be escaped as hex escape %20.

Gap
The alignment of the feature to the target if the two are not collinear
(e.g. contain gaps).
The alignment format is taken from the CIGAR format described in the Exonerate
documentation.

Derives_from
Used to disambiguate the relationship between one feature and another when the
relationship is a temporal one rather than a purely structural "part of" one.
This is needed for polycistronic genes.

Note
A free text note.

Dbxref
A database cross reference. See the GFF3 specification for more information.

Ontology_term
A cross reference to an ontology term.

Multiple attributes of the same type are indicated by separating
the values with the comma ","


Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |:mod:`skbio.sequence.Sequence`                                 |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.DNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.RNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.Protein`                                  |
+------+------+---------------------------------------------------------------+
|Yes   |No    |generator of :mod:`skbio.sequence.Sequence` objects            |
+------+------+---------------------------------------------------------------+


Reference
---------
.. [#] http://gmod.org/wiki/GFF3
'''

from skbio.util._misc import merge_dicts
from skbio.io import create_format, FileFormatError
from skbio.sequence import Sequence, DNA, RNA, Protein
from skbio.io.format._base import (
    _line_generator, _get_nth_sequence, _too_many_blanks)


class GFFFormatError(FileFormatError):
    pass

gff = create_format('gff')

# Annotation headers
_ANNOTATION_HEADERS = [
          'SEQID',
          'SOURCE',
          'TYPE',
          'START',
          'END',
          'SCORE',
          'STRAND',
          'PHASE',
          'ATTR'
          ]

# Attribute headers
_ATTR_HEADERS = [
        'ID',
        'Name',
        'Alias',
        'Parent',
        'Target',
        'Gap',
        'Derives_from',
        'Note',
        'Dbxref',
        'Ontology_term'
        ]


@gff.sniffer()
def _gff_sniffer(fh):
    # check the 1st real line is a valid ID line
    if _too_many_blanks(fh, 5):
        return False, {}

    try:
        line = next(_line_generator(fh, skip_blanks=True, strip=False))
    except StopIteration:
        return False, {}

    try:
        assert line.startswith('##')
    except GFFFormatError:
        return False, {}
    return True, {}


def _is_float(input):
    try:
        float(input)
    except ValueError:
        return False
    return True


def parse_optional(s):
    _field, _type, _val = s.split(':')
    if _type == 'i':
        return {_field: int(_val)}
    elif _type == 'f':
        return {_field: float(_val)}
    else:
        return {_field: _val}


def parse_required(s):
    if s.isdigit():
        return int(s)
    elif _is_float(s):
        return float(s)
    else:
        return s


def _construct(record, constructor=None, **kwargs):
    seq, md = record
    if constructor is None:
        constructor = Sequence
    if constructor == RNA:
        return DNA(
            seq, metadata=md, **kwargs).transcribe()
    else:
        return constructor(
            seq, metadata=md, **kwargs)


@gff.reader(None)
def _gff_to_generator(fh, constructor=None, **kwargs):
    for record in _parse_records(fh):
        yield _construct(record, constructor, **kwargs)


@gff.reader(Sequence)
def _gff_to_sequence(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh), seq_num)
    return _construct(record, Protein, **kwargs)


@gff.reader(Protein)
def _gff_to_protein(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh), seq_num)
    return _construct(record, Protein, **kwargs)


@gff.reader(DNA)
def _gff_to_DNA(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh), seq_num)
    return _construct(record, DNA, **kwargs)


@gff.reader(RNA)
def _gff_to_RNA(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh), seq_num)
    return _construct(record, RNA, **kwargs)


def _parse_records(fh, constructor=None, **kwargs):
    metadata = {}
    optional_headers = []
    headers = _ANNOTATION_HEADERS
    for line in _line_generator(fh, skip_blanks=True, strip=True):
        # parse the header (would be nice to abstract this pattern out)
        if line.startswith('@'):
            tabs = line.split('\t')
            key = tabs[0][1:]
            # FIXME:  The vals variable needs to be explicitly tested
            vals = tabs[1:]
            if key == 'CO':
                val = vals[0]
                optional_headers = val.split(',')
                headers = _ANNOTATION_HEADERS + optional_headers
                headers = headers[:9] + headers[10:]
            else:
                vals = tabs[1:]
                if len(vals) > 1:
                    metadata[key] = vals
                else:
                    metadata[key] = vals[0]

        # parse the actual sequences
        else:
            tabs = line.split('\t')

            # extract sequence
            seq = tabs[9]
            tabs = tabs[:9] + tabs[10:]

            req = list(map(parse_required, tabs[:10]))
            opt = list(map(parse_optional, tabs[10:]))
            req = dict(zip(_ANNOTATION_HEADERS, req))

            md = merge_dicts(metadata, req, *opt)
            yield seq, md
