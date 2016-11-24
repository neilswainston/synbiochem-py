'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-many-arguments
from xml.etree import ElementTree
import copy
import math
import re
import uuid

from synbiochem.utils import seq_utils


SO_CDS = 'http://purl.obolibrary.org/obo/SO_0000316'
SO_PROM = 'http://purl.obolibrary.org/obo/SO_0000167'

_RDF_NS = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
_NS = {'ns': 'http://sbols.org/v1#',
       'rdf': _RDF_NS}


class Dna(object):
    '''Class to represent a DNA object.'''

    def __init__(self, dna_id=None, name=None, desc=None, typ=None, seq=None,
                 start=float('NaN'), end=float('NaN'), forward=True):
        self.__dict__ = {'id': dna_id if dna_id is not None else _get_uid(),
                         'sequence': seq,
                         'name': name,
                         'description': desc,
                         'type': typ,
                         'start': start if not math.isnan(start) else 1,
                         'end': end if not math.isnan(end) else
                         (len(seq) if seq is not None else end),
                         'forward': forward,
                         'features': []
                         }

    def set_seq(self, seq):
        '''Sets sequence.'''
        self.__dict__['sequence'] = seq

        if math.isnan(self.__dict__['end']):
            self.__dict__['end'] = len(seq)

    def clone(self):
        '''Clones the DNA object (making a deep copy).'''
        clone_dna = Dna()
        clone_dna.__dict__ = copy.deepcopy(self.__dict__)
        return clone_dna

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return self.__dict__['id']


def read(filename):
    '''Parses SBOL v1.'''
    tree = ElementTree.parse(filename)
    root = tree.getroot()

    dna_comp = root.find('ns:DnaComponent', _NS)

    dna = Dna(*_read_dna_comp(dna_comp))

    dna_seq = dna_comp.find('ns:dnaSequence', _NS)
    dna_seq = dna_seq.find('ns:DnaSequence', _NS)
    dna.set_seq(dna_seq.find('ns:nucleotides', _NS).text)

    for annot in dna_comp.findall('ns:annotation', _NS):
        seq_annot = annot.find('ns:SequenceAnnotation', _NS)
        sub_comp = seq_annot.find('ns:subComponent', _NS)
        dna_comp = sub_comp.find('ns:DnaComponent', _NS)

        feat = Dna(*_read_dna_comp(dna_comp),
                   start=int(seq_annot.find('ns:bioStart', _NS).text),
                   end=int(seq_annot.find('ns:bioEnd', _NS).text),
                   forward=seq_annot.find('ns:strand', _NS).text == '+')

        dna.features.append(feat)

    return dna


def write(dna, filename=None):
    '''Writes a Dna object to SBOL v1.'''
    root = ElementTree.Element('ns2:RDF', {'xmlns': 'http://sbols.org/v1#',
                                           'xmlns:ns2': _RDF_NS})

    dna_comp = _write_dna_comp(root, dna)

    dna_seq = _write(dna_comp, 'dnaSequence')
    dna_seq = _write(dna_seq, 'DnaSequence', _get_about())
    _write(dna_seq, 'nucleotides', text=dna.sequence)

    for feature in dna.features:
        annot = _write(dna_comp, 'annotation')
        annot = _write(annot, 'SequenceAnnotation', _get_about())
        _write(annot, 'bioStart', text=str(feature.start))
        _write(annot, 'bioEnd', text=str(feature.end))
        _write(annot, 'strand', text='+' if feature.forward else '-')
        sub_comp = _write(annot, 'subComponent')
        _write_dna_comp(sub_comp, feature)

    sbol = ElementTree.tostring(root)

    with open(filename, 'w') as outfile:
        outfile.write(sbol)

    return sbol


def concat(dnas):
    '''Concatenates a list of DNA objects into a single DNA object.'''
    concat_dna = dnas[0].clone()

    for dna in dnas[1:]:
        concat_dna = _add(concat_dna, dna)

    return concat_dna


def apply_restricts(dna, restricts, circular=False):
    '''Applies restriction site cleavage to forward and reverse strands.'''
    out_dnas = [dna]

    for restrict in restricts:
        out_dnas = _apply_restrict_to_dnas(out_dnas, restrict)

    if circular and len(out_dnas) > 1:
        return [concat([out_dnas[-1], out_dnas[0]])] + out_dnas[1:-1]
    else:
        return out_dnas


def _read_dna_comp(dna_comp):
    '''Read DNAComponent node.'''
    disp_id = _find_text(dna_comp, 'ns:displayId')
    name = _find_text(dna_comp, 'ns:name')
    desc = _find_text(dna_comp, 'ns:description')
    typ_node = dna_comp.find('rdf:type', _NS)
    typ = typ_node.attrib['{' + _RDF_NS + '}resource'] \
        if typ_node is not None else None
    return disp_id, name, desc, typ


def _find_text(parent, field):
    '''Finds text from node.'''
    node = parent.find(field, _NS)
    return None if node is None else node.text


def _write_dna_comp(parent, dna):
    '''Write DNAComponent node.'''
    dna_comp = ElementTree.SubElement(parent, 'DnaComponent',
                                      _get_about(dna.id))

    if dna.type:
        _write(dna_comp, 'ns2:type', {'ns2:resource': dna.type})

    _write(dna_comp, 'displayId', text=dna.id)

    if dna.name:
        _write(dna_comp, 'name', text=dna.name)

    if dna.description:
        _write(dna_comp, 'description', text=dna.description)

    return dna_comp


def _write(parent, name, params=None, text=None):
    '''Writes a node.'''
    if params is None:
        params = {}

    node = ElementTree.SubElement(parent, name, params)
    node.text = text
    return node


def _add(dna1, dna2):
    '''Adds two DNA objects together.'''
    # Add names, etc.
    dna1.id = _concat([dna1.id, dna2.id])
    dna1.name = _concat([dna1.name, dna2.name])
    dna1.description = _concat([dna1.description, dna2.description])

    # Add sequences:
    orig_seq_len = len(dna1.sequence)
    dna1.sequence += dna2.sequence

    # Update features:
    for feature in dna2.features:
        feature = feature.clone()
        feature.start += orig_seq_len
        feature.end += orig_seq_len
        dna1.features.append(feature)

    return dna1


def _concat(strs):
    '''Concatenates non-None strings.'''
    return '_'.join([string for string in strs if string is not None])


def _apply_restrict_to_dnas(dnas, restrict):
    '''Applies restriction site cleavage to forward and reverse strands.'''
    out_dnas = []

    for dna in dnas:
        parent_seq = dna.sequence

        for forw in _apply_restrict_to_seq(parent_seq, restrict):
            for rev in _apply_restrict_to_seq(seq_utils.get_rev_comp(forw[0]),
                                              restrict):
                rev_comp = seq_utils.get_rev_comp(rev[0])
                start = forw[1] + len(forw[0]) - rev[1] - len(rev[0])
                end = start + len(rev_comp)
                out_dnas.append(_get_concat_dna(dna, rev_comp, start, end))

    return out_dnas


def _apply_restrict_to_seq(seq, restrict):
    '''Applies restriction site cleavage to a sequence.'''
    sub_seqs = [(match.group(0), match.start())
                for match in re.finditer(restrict, seq)]
    end = sub_seqs[0][1] if len(sub_seqs) > 0 else len(seq)
    return [(seq[:end], 0)] + sub_seqs


def _get_concat_dna(parent_dna, seq, start, end):
    '''Returns a DNA object from the supplied subsequence from a parent DNA
    object.'''
    frag_str = '_' + str(start) + '_' + str(end) \
        if parent_dna.sequence != seq else ''

    dna = Dna(parent_dna.id + frag_str,
              parent_dna.name + frag_str,
              parent_dna.description + frag_str,
              seq=seq)

    # TODO: This may not work for sub-sequences arriving from circular DNA:
    for feature in parent_dna.features:
        if feature.start >= start and feature.end <= end:
            clone_feature = feature.clone()
            clone_feature.start -= start
            clone_feature.end -= start
            dna.features.append(clone_feature)

    return dna


def _get_about(uid=None):
    '''Gets about attributes.'''
    if uid is None:
        uid = _get_uid()

    return {'ns2:about': 'https://www.synbiochem.co.uk#' + uid}


def _get_uid():
    '''Gets a unique (valid) id.'''
    return '0' + str(uuid.uuid4()).replace('-', '')
