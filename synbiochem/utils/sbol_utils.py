'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  pablocarbonell / neilswainston / alanwilliams
'''
from sbol.sbol import Collection, DNAComponent, DNASequence, Document, \
    SequenceAnnotation


def concatenate(sbol_docs):
    '''Concatenates a list of Documents into a single Document.'''
    concat = clone(sbol_docs[0])

    for sbol_doc in sbol_docs[1:]:
        concat = _add(concat, sbol_doc)

    return concat


def clone(orig_doc):
    '''Clones an sbol Document.'''
    clone_doc = Document()

    for obj in orig_doc.components:
        _clone_comp(clone_doc, obj)

    for obj in orig_doc.collections:
        coll = Collection(clone_doc, obj.uri)
        coll.description = obj.description
        coll.display_id = obj.display_id
        coll.name = obj.name

    return clone_doc


def _clone_comp(owner_doc, orig_comp):
    '''Clones a DNAComponent.'''
    for comp in owner_doc.components:
        if comp.uri == orig_comp.uri:
            return comp

    comp = DNAComponent(owner_doc, orig_comp.uri)
    comp.description = orig_comp.description
    comp.display_id = orig_comp.display_id
    comp.name = orig_comp.name
    comp.type = orig_comp.type

    if orig_comp.sequence is not None:
        comp.sequence = _clone_seq(owner_doc, orig_comp.sequence)

    for obj in orig_comp.annotations:
        annot = SequenceAnnotation(owner_doc, obj.uri)
        annot.start = obj.start
        annot.end = obj.end
        annot.isDownstream = obj.isDownstream
        annot.isUpstream = obj.isUpstream
        annot.strand = obj.strand

        if obj.subcomponent is not None:
            annot.subcomponent = _clone_comp(owner_doc, obj.subcomponent)

        comp.annotations += annot
    return comp


def _clone_seq(owner_doc, orig_seq):
    '''Clones a DNASequence.'''
    sequence = DNASequence(owner_doc, orig_seq.uri)
    sequence.nucleotides = orig_seq.nucleotides
    return sequence


def _add(sbol_doc1, sbol_doc2):
    '''Adds two sbol Documents together.'''
    return sbol_doc1


def main():
    '''main method.'''
    sbol_doc = Document()
    sbol_doc.read('/Users/neilswainston/Downloads/ahpC L177Q.xml')

    clone_doc = clone(sbol_doc)
    clone_doc.write('/Users/neilswainston/Downloads/ahpC L177Q clone.xml')

    print concatenate([sbol_doc] * 3)

if __name__ == '__main__':
    main()
