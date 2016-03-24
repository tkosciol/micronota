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
|Yes   |No    |:mod:`skbio.metadata.IntervalMetadata` objects                 |
+------+------+---------------------------------------------------------------+
|Yes   |No    |generator of :mod:`skbio.metadata.IntervalMetadata` objects    |
+------+------+---------------------------------------------------------------+


Reference
---------
.. [#] http://gmod.org/wiki/GFF3
'''

from skbio.util._misc import merge_dicts
from skbio.io import create_format, FileFormatError
from skbio.metadata import IntervalMetadata, Feature
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
        assert line.startswith('##gff-version 3')
    except GFFFormatError:
        return False, {}
    return True, {}


def _is_float(input):
    try:
        float(input)
    except ValueError:
        return False
    return True


def parse_attr(s):
    _field = s.split(';')
    _attr = ()
    for f in _field:
        try:
            _type, _val = f.split('=')
        except ValueError:
            _type, _val = 0, 0

        _attr = _attr + (_type, _val)
    return _attr


def parse_required(s):
    if s.isdigit():
        return int(s)
    elif _is_float(s):
        return float(s)
    else:
        return s


def IntervalMetadata_construct(_attr, _interval):
    im = IntervalMetadata()
    im.add(Feature(**_attr), _interval)
    return im


def _construct(record, constructor=None, **kwargs):
    attr, meta = record
    if constructor is None:
        constructor = IntervalMetadata_construct
    return constructor(attr, meta)


@gff.reader(None)
def _gff_to_generator(fh, constructor=None, **kwargs):
    for record in _parse_records(fh):
        yield _construct(record, constructor, **kwargs)


@gff.reader(IntervalMetadata)
def _gff_to_metadata(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh), seq_num)
    return _construct(record, IntervalMetadata_construct, **kwargs)


# generator passed to gff.reader()
# OUT: yield annotation `annot` list of dicts and `intervals` list of tuples
def _parse_records(fh, constructor=None, **kwargs):
    for line in _line_generator(fh, skip_blanks=True, strip=True):
        if not line.startswith('#'):
            tabs = line.split('\t')

            # extract annotations
            _headers = tabs[:3] + tabs[5:8]

            # extract attributes
            _attributes = tabs[8]

            # extract intervals
            _intervals = tabs[3:5]

            annot = list(map(parse_required, _headers))
            intervals = tuple(map(parse_required, _intervals))
            attr = parse_attr(_attributes)
            annot.append(attr)

            annot = dict(zip(_ANNOTATION_HEADERS, annot))

            yield annot, intervals
