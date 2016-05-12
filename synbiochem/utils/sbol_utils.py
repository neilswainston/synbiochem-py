'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  pablocarbonell / neilswainston / alanwilliams
'''
import re
import uuid

from sbol.sbol import DNAComponent, DNASequence, Document, \
    SequenceAnnotation

import synbiochem.utils.sequence_utils as seq_utils


SO_CDS = 'http://purl.obolibrary.org/obo/SO_0000316'
SO_PROM = 'http://purl.obolibrary.org/obo/SO_0000167'
_DEFAULT_URI_PREFIX = 'http://synbiochem.co.uk/'


def get_seq(sbol_doc):
    '''Returns the sequence from an sbol Document.'''
    return sbol_doc.sequences[0].nucleotides \
        if len(sbol_doc.sequences) > 0 else None


def concat(sbol_docs):
    '''Concatenates a list of Documents into a single Document.'''
    concat_doc = clone(sbol_docs[0])

    for sbol_doc in sbol_docs[1:]:
        concat_doc = _add(concat_doc, sbol_doc)

    return concat_doc


def clone(orig_doc):
    '''Clones an sbol Document.'''
    clone_doc = Document()

    for comp in orig_doc.components:
        _clone_comp(clone_doc, comp)

    return clone_doc


def apply_restricts(doc, restricts, circular=False,
                    uri_prefix=_DEFAULT_URI_PREFIX):
    '''Applies restriction site cleavage to forward and reverse strands.'''
    out_docs = [doc]

    for restrict in restricts:
        out_docs = _apply_restrict_to_docs(out_docs, restrict, uri_prefix)

    if circular and len(out_docs) > 1:
        return [concat([out_docs[-1], out_docs[0]])] + out_docs[1:-1]
    else:
        return out_docs


def apply_pcr(doc, for_primer, rev_primer, circular=False,
              uri_prefix=_DEFAULT_URI_PREFIX):
    '''Applies restriction site cleavage to forward and reverse strands.'''
    sbol_docs = []
    seq = doc.sequences[0].nucleotides
    starts = [m.start() for m in re.finditer(for_primer, seq)]
    ends = sorted([len(seq) - m.start()
                   for m in re.finditer(rev_primer,
                                        seq_utils.get_rev_comp(seq))])

    for start in starts:
        for end in ends:
            if start < end:
                sbol_docs.append(_get_sbol(doc, seq[start:end], start, end,
                                           uri_prefix))
            elif circular:
                sbol_docs.append(_get_sbol(doc, seq[start:] + seq[:end], end,
                                           start, uri_prefix))

    return sbol_docs


def _clone_comp(doc, orig_comp):
    '''Clones a DNAComponent.'''
    for comp in doc.components:
        if comp.uri == orig_comp.uri:
            return comp

    return _build_comp(doc, orig_comp.uri, orig_comp)


def _clone_sub_comp(doc, orig_comp, uri_prefix=_DEFAULT_URI_PREFIX):
    '''Clones a DNAComponent.'''
    uri = orig_comp.uri
    display_id = None

    for comp in doc.components:
        if comp.uri == orig_comp.uri:
            display_id = str(uuid.uuid4())
            uri = uri_prefix + display_id
            break

    return _build_comp(doc, uri, orig_comp, display_id)


def _build_comp(doc, uri, orig_comp, disp_id=None):
    '''Builds (essentially copies) a DNAComponent.'''
    comp = DNAComponent(doc, uri)
    comp.description = orig_comp.description
    comp.display_id = orig_comp.display_id if disp_id is None else disp_id
    comp.name = orig_comp.name
    comp.type = orig_comp.type

    if orig_comp.sequence is not None:
        comp.sequence = DNASequence(doc, orig_comp.sequence.uri)
        comp.sequence.nucleotides = orig_comp.sequence.nucleotides

    for annot in orig_comp.annotations:
        clone_annot = _clone_annotation(doc, annot)
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
        clone_annot.subcomponent = _clone_sub_comp(owner_doc,
                                                   annot.subcomponent)

    return clone_annot


def _apply_restrict_to_docs(docs, restrict, uri_prefix=_DEFAULT_URI_PREFIX):
    '''Applies restriction site cleavage to forward and reverse strands.'''
    out_docs = []

    for doc in docs:
        parent_seq = doc.sequences[0].nucleotides

        for forw in _apply_restrict_to_seq(parent_seq, restrict):
            for rev in _apply_restrict_to_seq(seq_utils.get_rev_comp(forw[0]),
                                              restrict):
                rev_comp = seq_utils.get_rev_comp(rev[0])

                if rev_comp == forw[0]:
                    start = forw[1]
                    end = start + len(rev_comp)
                else:
                    end = forw[1] + len(forw[0])
                    start = end - len(rev_comp)

                out_docs.append(_get_sbol(doc,
                                          rev_comp,
                                          start,
                                          end,
                                          uri_prefix))

    return out_docs


def _apply_restrict_to_seq(seq, restrict):
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
        sbol_doc1.components[0].annotations += clone_annot

    return sbol_doc1


def _get_sbol(parent_doc, seq, start, end, uri_prefix):
    '''Returns a sbol Document from the supplied subsequence from a parent sbol
    Document.'''
    parent_comp = parent_doc.components[0]

    frag_str = ' [' + str(start) + ':' + str(end) + ']'

    doc = Document()
    comp = DNAComponent(doc, uri_prefix + str(uuid.uuid4()))
    comp.display_id = parent_comp.display_id + frag_str
    comp.name = parent_comp.name + frag_str
    comp.description = parent_comp.description + frag_str
    comp.sequence = DNASequence(doc, _get_uri(uri_prefix))
    comp.sequence.nucleotides = seq.lower()

    # TODO: This may not work for sub-sequences arriving from circular DNA:
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
