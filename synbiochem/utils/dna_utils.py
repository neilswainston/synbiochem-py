'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-arguments
from xml.etree import ElementTree
import math
import uuid

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

    def add_feature(self, feature):
        '''Adds feature.'''
        self.__dict__['features'].append(feature)

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

        dna.add_feature(feat)

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


def _get_about(uid=None):
    '''Gets about attributes.'''
    if uid is None:
        uid = _get_uid()

    return {'ns2:about': 'https://www.synbiochem.co.uk#' + uid}


def _get_uid():
    '''Gets a unique (valid) id.'''
    return '0' + str(uuid.uuid4()).replace('-', '')

write(read('sbol.xml'), 'out.xml')
