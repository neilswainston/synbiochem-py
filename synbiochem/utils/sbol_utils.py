'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  pablocarbonell / neilswainston / alanwilliams
'''
import uuid

from sbol.sbol import Collection, DNAComponent, DNASequence, Document, \
    SequenceAnnotation


def concatenate(sbol_docs, uri_prefix='http://synbiochem.co.uk#'):
    '''Concatenates a list of Documents into a single Document.'''
    concat = clone(sbol_docs[0], uri_prefix)

    for sbol_doc in sbol_docs[1:]:
        concat = _add(concat, sbol_doc)

    return concat


def clone(orig_doc, uri_prefix='http://synbiochem.co.uk#'):
    '''Clones an sbol Document.'''
    clone_doc = Document()

    for obj in orig_doc.components:
        _clone_comp(clone_doc, obj, uri_prefix)

    for obj in orig_doc.collections:
        coll = Collection(clone_doc, (_get_uri(uri_prefix)
                                      if uri_prefix is not None
                                      else obj.uri))
        coll.description = obj.description
        coll.display_id = obj.display_id
        coll.name = obj.name

    return clone_doc


def _clone_comp(owner_doc, orig_comp, uri_prefix):
    '''Clones a DNAComponent.'''
    for comp in owner_doc.components:
        if comp.uri == orig_comp.uri:
            return comp

    comp = DNAComponent(owner_doc, (_get_uri(uri_prefix)
                                    if uri_prefix is not None
                                    else orig_comp.uri))
    comp.description = orig_comp.description
    comp.display_id = orig_comp.display_id
    comp.name = orig_comp.name
    comp.type = orig_comp.type

    if orig_comp.sequence is not None:
        comp.sequence = DNASequence(owner_doc, (_get_uri(uri_prefix)
                                                if uri_prefix is not None
                                                else orig_comp.sequence.uri))
        comp.sequence.nucleotides = orig_comp.sequence.nucleotides

    for obj in orig_comp.annotations:
        annot = _clone_annotation(owner_doc, obj)
        comp.annotations += annot

    return comp


def _clone_annotation(owner_doc, annot):
    '''Clones a SequenceAnnotation.'''
    clone_annot = SequenceAnnotation(owner_doc, annot.uri)
    clone_annot.start = annot.start
    clone_annot.end = annot.end
    clone_annot.isDownstream = annot.isDownstream
    clone_annot.isUpstream = annot.isUpstream
    clone_annot.strand = annot.strand

    if annot.subcomponent is not None:
        clone_annot.subcomponent = _clone_comp(owner_doc, annot.subcomponent,
                                               None)
    return clone_annot


def _add(sbol_doc1, sbol_doc2):
    '''Adds two sbol Documents together.'''
    # Add names, etc.
    sbol_doc1.components[0].description += ' + ' + \
        sbol_doc2.components[0].description
    sbol_doc1.components[0].display_id += ' + ' + \
        sbol_doc2.components[0].display_id
    sbol_doc1.components[0].name += ' + ' + sbol_doc2.components[0].name

    # Add sequences:
    orig_seq_len = len(sbol_doc1.sequences[0].nucleotides)
    sbol_doc1.sequences[0].nucleotides += sbol_doc2.sequences[0].nucleotides

    # Update SequenceAnnotations:
    for obj in sbol_doc2.annotations:
        annot = obj  # _clone_annotation(sbol_doc1, obj)
        annot.start += orig_seq_len
        annot.end += orig_seq_len

        if annot.subcomponent is not None:
            annot.subcomponent = _clone_comp(sbol_doc1, obj.subcomponent,
                                             None)

        sbol_doc1.components[0].annotations += annot

    return sbol_doc1


def _get_uri(uri_prefix):
    '''Returns a new unique URI.'''
    return uri_prefix + str(uuid.uuid4())
