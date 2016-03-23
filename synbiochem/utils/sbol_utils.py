'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  pablocarbonell / neilswainston / alanwilliams
'''
import re
import uuid

from sbol.sbol import Collection, DNAComponent, DNASequence, Document, \
    SequenceAnnotation
import synbiochem.utils.sequence_utils as seq_utils


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


def apply_restrict(doc, restrict, uri_prefix='http://synbiochem.co.uk#'):
    '''Applies restriction site cleavage to forward and reverse strands.'''
    sbol_docs = []
    parent_seq = doc.sequences[0].nucleotides.upper()

    for forw in _apply_restrict(parent_seq, restrict):
        for rev in _apply_restrict(seq_utils.get_rev_comp(forw[0]), restrict):
            sbol_docs.append(_get_sbol(doc,
                                       seq_utils.get_rev_comp(rev[0]),
                                       forw[1],
                                       uri_prefix))
    return sbol_docs


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

    for annot in orig_comp.annotations:
        clone_annot = _clone_annotation(owner_doc, annot)
        comp.annotations += clone_annot

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


def _apply_restrict(seq, restrict):
    '''Applies restriction site cleavage to a sequence.'''
    sub_seqs = [(match.group(0), match.start())
                for match in re.finditer(restrict, seq)]
    end = sub_seqs[0][1] if len(sub_seqs) > 0 else len(seq)
    return [(seq[:end], 0)] + sub_seqs


def _add(sbol_doc1, sbol_doc2):
    '''Adds two sbol Documents together.'''
    # Add names, etc.
    comp1 = sbol_doc1.components[0]
    comp2 = sbol_doc2.components[0]
    comp1.description = _concat([comp1.description, comp2.description])
    comp1.display_id = _concat([comp1.display_id, comp2.display_id])
    comp1.name = _concat([comp1.name, comp2.name])

    # Add sequences:
    orig_seq_len = len(sbol_doc1.sequences[0].nucleotides)
    sbol_doc1.sequences[0].nucleotides += sbol_doc2.sequences[0].nucleotides

    # Update SequenceAnnotations:
    for annot in sbol_doc2.annotations:
        clone_annot = _clone_annotation(sbol_doc1, annot)
        clone_annot.start += orig_seq_len
        clone_annot.end += orig_seq_len

        if clone_annot.subcomponent is not None:
            clone_annot.subcomponent = _clone_comp(sbol_doc1,
                                                   annot.subcomponent,
                                                   None)

        sbol_doc1.components[0].annotations += clone_annot

    return sbol_doc1


def _get_sbol(parent_doc, seq, start, uri_prefix):
    '''Returns a sbol Document from the supplied subsequence from a parent sbol
    Document.'''
    doc = Document()
    display_id = str(uuid.uuid4())
    comp = DNAComponent(doc, uri_prefix + display_id)
    comp.display_id = display_id
    comp.sequence = DNASequence(doc, _get_uri(uri_prefix))
    comp.sequence.nucleotides = seq.lower()

    end = start + len(seq)

    for annot in parent_doc.annotations:
        if annot.start >= start and annot.end <= end:
            clone_annot = _clone_annotation(doc, annot)
            clone_annot.start -= start
            clone_annot.end -= start
            comp.annotations += clone_annot

    return doc


def _get_uri(uri_prefix):
    '''Returns a new unique URI.'''
    return uri_prefix + str(uuid.uuid4())


def _concat(strs):
    '''Concatenates non-None strings.'''
    return ' + '.join([string for string in strs if string is not None])
