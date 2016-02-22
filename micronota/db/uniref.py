r'''
UniRef
======

.. currentmodule:: micronota.db.uniref

This module create databases from UniRef for micronota usage. Using UniRef
for functional assignment of coding genes, we can get the clustering of
the annotated proteins for free. Using UniRef100 as the reference database,
we can easily collapse the clusters down to the similarity levels
of 90% or 50%. For those UniRef records from UniProKB, we also transfer
the metadata associated with those UniProKB records.

UniProtKB (UniProt Knowledgebase) [#]_ contains two parts, UniProtKB/Swiss-Prot
and UniProtKB/TrEMBL. The first part contains protein records that are
manually annotated and reviewed while the second part is done computationally
and not manually reviewed.

UniParc (UniProt Archive) [#]_ is also part of UniProt project. It is a
is a comprehensive and non-redundant database that contains most of the
publicly available protein sequences in the world.

The UniProt Reference Clusters (UniRef) [#]_ provide clustered sets (UniRef100,
UniRef90 and UniRef50 clusters) of sequences from the UniProtKB
and selected UniParc records, in order to obtain complete coverage of sequence
space at several resolutions (100%, >90% and >50%) while hiding redundant
sequences (but not their descriptions) from view.

If you used this database through micronota, you should also `cite UniProt
<http://www.uniprot.org/help/publications>`_ properly.

Release
-------
:version: 2016_01
:date:    01/20/2016
:link:    ftp://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt

Files
-----
* UniRef files

  UniRef fasta files clustered at the similarity levels of 50, 90 or 100.
  The UniRef100 identifier is generated by placing ``UniRef100_`` prefix
  before the UniProtKB accession number or UniParc identifier of the
  representative UniProtKB or UniParc entry. It is similarly done for
  UniRef50 or UniRef90.

  * uniref50.fasta.gz [#]_

  * uniref90.fasta.gz [#]_

  * uniref100.fasta.gz [#]_

* UniProtKB files

  The .dat.gz files contain the UniProtKB records in a variant format of EMBL.
  The records are divided into taxa groups [#]_.

  * uniprot_sprot.dat.gz

  * uniprot_trembl.dat.gz

Reference
---------
.. [#] http://www.uniprot.org/help/uniparc
.. [#] http://www.uniprot.org/help/uniprotkb
.. [#] http://www.ncbi.nlm.nih.gov/pubmed/17379688
.. [#] ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/README
.. [#] ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/README
.. [#] ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/README
.. [#] ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/README
'''

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join, basename
from sqlite3 import connect
from xml.etree import ElementTree as ET
from itertools import product
import gzip

from skbio import read, Sequence

from ..util import _overwrite, _download


_status = ['Swiss-Prot', 'TrEMBL']
_kingdom = ['Bacteria', 'Archaea', 'Viruses', 'Eukaryota', 'other']


def prepare_db(out_d, downloaded, force=False,
               sprot='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz',
               trembl='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.xml.gz',
               uniref100='ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz'):
    '''Prepare reference database for UniRef.

    Parameters
    ----------
    out_d : str
        The output directory
    force : boolean
        Force overwrite the files
    downloaded : str
        File directory. If the files are already downloaded,
        just use the files in the directory and skip _download.

    Notes
    -----
    This function creates the following files:

    * ``uniref100_Swiss-Prot_Archaea.dmnd``
    * ``uniref100_Swiss-Prot_Bacteria.dmnd``
    * ``uniref100_Swiss-Prot_Viruses.dmnd``
    * ``uniref100_Swiss-Prot_other.dmnd``
    * ``uniref100_Swiss-Prot_Eukaryota.dmnd``
    * ``uniref100_TrEMBL_Archaea.dmnd``
    * ``uniref100_TrEMBL_Bacteria.dmnd``
    * ``uniref100_TrEMBL_Viruses.dmnd``
    * ``uniref100_TrEMBL_other.dmnd``
    * ``uniref100_TrEMBL_Eukaryota.dmnd``
    * ``uniref100__other.dmnd``

    * ``uniprotkb.db``
    '''
    metadata_db = join(out_d, 'uniprotkb.db')

    sprot_raw = join(downloaded, basename(sprot))
    trembl_raw = join(downloaded, basename(trembl))
    uniref100_raw = join(downloaded, basename(uniref100))
    try:
        _download(sprot, sprot_raw, overwrite=force)
        _download(trembl, trembl_raw, overwrite=force)
        _download(uniref100, uniref100_raw, overwrite=force)
    except FileExistsError:
        pass

    prepare_metadata([sprot_raw, trembl_raw], metadata_db)

    sort_uniref(metadata_db, uniref100_raw, out_d, force)


def sort_uniref(db_fp, uniref_fp, out_d, overwrite=False):
    '''Sort UniRef sequences into different partitions.

    This will sort UniRef100 seq into following partitions based on both
    quality and taxon:

    * ``uniref100_Swiss-Prot_Archaea.fasta``
    * ``uniref100_Swiss-Prot_Bacteria.fasta``
    * ``uniref100_Swiss-Prot_Viruses.fasta``
    * ``uniref100_Swiss-Prot_other.fasta``
    * ``uniref100_Swiss-Prot_Eukaryota.fasta``
    * ``uniref100_TrEMBL_Archaea.fasta``
    * ``uniref100_TrEMBL_Bacteria.fasta``
    * ``uniref100_TrEMBL_Viruses.fasta``
    * ``uniref100_TrEMBL_other.fasta``
    * ``uniref100_TrEMBL_Eukaryota.fasta``
    * ``uniref100__other.fasta``

    Parameters
    ----------
    db_fp : str
        The database file created by ``prepare_metadata``.
    uniref_fp : str
        The UniRef100 fasta file. gzipped or not.
    out_d : str
        The output directory.
    '''
    fns = ['%s_%s' % (i, j) for i, j in product(_status, _kingdom)]
    fns.append('_other')
    fps = [join(out_d, 'uniref100_%s.fasta') % f for f in fns]

    for f in fns:
        _overwrite(f, overwrite)

    files = {fn: open(fp, 'w') for fp, fn in zip(fps, fns)}

    with connect(db_fp) as conn:
        cursor = conn.cursor()
        for seq in read(uniref_fp, format='fasta', constructor=Sequence):
            id = seq.metadata['id']
            ac = id.replace('UniRef100_', '')
            group = ['', 'other']
            cursor.execute('''SELECT * FROM metadata
                              WHERE ac = ?''',
                           (ac,))
            for _, s, k in cursor.fetchall():
                group[0] = _status[s]
                group[1] = _kingdom[k]
            seq.write(files['_'.join(group)])

    for f in files:
        files[f].close()


def prepare_metadata(in_fps, db_fp, **kwargs):
    '''
    Parameters
    ----------
    in_fps : list of str
        The gzipped files of either UniProtKB Swiss-Prot or TrEMBLE.
    db_fp : str
        The output database file. See ``Notes``.
    kwargs : dict
        keyword args passed to ``_overwrite``

    Returns
    -------
    int
        The number of records processed.

    Notes
    -----
    The schema of the database file contains one table named `tigrfam` that
    has following columns:

    1. ``ac``. TEXT. UniProtKB primary accession.

    2. ``status``. INT. The index in ``_status``. ``0`` is 'Swiss-Prot'
       and ``1`` is 'TrEMBL'.

    3. ``kingdom``. INT. The index in ``_kingdom``. ``0``, ``1``, ``2``,
       ``3``, and ``4`` represent 'Bacteria', 'Archaea', 'Viruses',
       'Eukaryota', and 'other', respectively.

    The table in the database file will be dropped and re-created if
    the function is re-run.
    '''
    _overwrite(db_fp, **kwargs)
    status_map = {k: i for i, k in enumerate(_status)}
    kingdom_map = {k: i for i, k in enumerate(_kingdom)}
    # this is the namespace for uniprot xml files.
    ns_map = {'xmlns': 'http://uniprot.org/uniprot',
              'xsi': 'http://www.w3.org/2001/XMLSchema-instance'}
    n = 0
    with connect(db_fp) as conn:
        table_name = 'metadata'
        conn.execute('''CREATE TABLE IF NOT EXISTS {t} (
                            ac       TEXT    NOT NULL,
                            status   INT     NOT NULL,
                            kingdom  INT     NOT NULL);'''.format(
                                t=table_name))
        insert = '''INSERT INTO {t} (ac, status, kingdom)
                    VALUES (?,?,?);'''.format(t=table_name)

        for fp in in_fps:
            for i, elem in enumerate(_parse_xml(fp, ns_map), 1):
                ac, status, kingdom = _process_entry(elem, ns_map)
                status = status_map[status]
                kingdom = kingdom_map.get(kingdom, 4)
                conn.execute(insert, (ac, status, kingdom))
            n += i
        # don't forget to index the column to speed up query
        conn.execute('CREATE UNIQUE INDEX IF NOT EXISTS ac ON {t} (ac);'.format(
            t=table_name))
        conn.commit()
    return n


def _parse_xml(in_fp, ns_map):
    def fixtag(ns, tag, nsmap):
        return '{%s}%s' % (nsmap[ns], tag)

    # it is very important to set the events to 'end'; otherwise,
    # elem would be an incomplete record.
    for event, elem in ET.iterparse(gzip.open(in_fp), events=['end']):
        if elem.tag == fixtag('xmlns', 'entry', ns_map):
            yield elem
            # this is necessary for garbage collection
            elem.clear()


def _process_entry(root, ns_map):
    group = root.attrib['dataset']
    try:
        accession = root.find('./xmlns:accession', ns_map).text
    except AttributeError as e:
        for child in root:
            print(child.tag, child.text)
        raise e
    try:
        taxon = root.find('./xmlns:organism/xmlns:lineage/xmlns:taxon',
                          ns_map).text
    except AttributeError:
        taxon = None

    return accession, group, taxon
