'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=ungrouped-imports
from _collections import defaultdict
import re
import sys
import urllib


from libsbml import CVTerm, SBMLDocument, writeSBMLToFile, \
    BIOLOGICAL_QUALIFIER, BQB_IS
from synbiochem import biochem4j
from synbiochem.utils import seq_utils
import libsbml

KCAT = 'kcat'
KCAT_FOR = 'kcat_for'
KCAT_REV = 'kcat_rev'
K_M = 'KM'
BIOCHEM_REACT = 'BIOCHEM_REACT'
VMAX = 'vmax'
CONC = 'CONC'
COMP_INHIB = 'COMP_INHIB'
NON_COMP_INHIB = 'NON_COMP_INHIB'
SIMPLE_CHEM = 'SIMPLE_CHEM'
PROTEIN = 'PROTEIN'
K_I = 'KI'
KEQ = 'KEQ'
COMPARTMENT = 'COMPARTMENT'
ENZYME = 'ENZYME'
VOLUME = 'volume'
ENZYME_CONC = 'Et'

SBO_TERMS = {KCAT: 25,
             KCAT_FOR: 320,
             KCAT_REV: 321,
             K_M: 27,
             BIOCHEM_REACT: 176,
             VMAX: 186,
             CONC: 196,
             COMP_INHIB: 206,
             NON_COMP_INHIB: 207,
             SIMPLE_CHEM: 247,
             PROTEIN: 252,
             K_I: 261,
             KEQ: 281,
             COMPARTMENT: 290,
             ENZYME: 460,
             VOLUME: 468}


UNITS = {SBO_TERMS[KCAT]: 'per_sec',
         SBO_TERMS[KCAT_FOR]: 'per_sec',
         SBO_TERMS[KCAT_REV]: 'per_sec',
         SBO_TERMS[CONC]: 'M',
         SBO_TERMS[K_M]: 'M',
         SBO_TERMS[K_I]: 'M',
         SBO_TERMS[VMAX]: 'M_per_sec',
         SBO_TERMS[KEQ]: 'dimensionless'}


def get_document(params):
    '''Gets a model.'''
    document = SBMLDocument(2, 5)
    model = document.createModel()
    _set_units(model)
    _add_compartment(model, 'c')

    nodes = defaultdict(dict)
    rels = []
    react_to_uniprot = {}

    for param in params:
        ids = param.split(':')
        data = _get_reaction(ids[0])
        _parse(data, nodes, rels)

        if len(ids) > 1:
            react_to_uniprot[_get_id(ids[0])] = ids[1]

    for cid, chemical in nodes['c'].iteritems():
        _add_species(model, cid, chemical)

    for cid, reaction in nodes['r'].iteritems():
        _add_reaction(model, cid, reaction)

    for rel in rels:
        _add_species_ref(model, str(rel[0]), str(rel[2]),
                         rel[1]['stoichiometry'])

    for react_id, uniprot_id in react_to_uniprot.iteritems():
        _add_modifier(model, react_id, uniprot_id)

    return document


def get_annotations(sbase):
    '''Gets list of annotations for a given sbase.'''
    annotations = []

    if sbase.isSetAnnotation():
        cv_terms = sbase.getCVTerms()
        for i in range(cv_terms.getSize()):
            cv_term = cv_terms.get(i)
            qual_type = cv_term.getQualifierType()
            specific_qual_type = cv_term.getModelQualifierType() \
                if qual_type == libsbml.MODEL_QUALIFIER \
                else cv_term.getBiologicalQualifierType()

            for j in range(cv_term.getNumResources()):
                resource = cv_term.getResourceURI(j)
                annotations.append((resource, qual_type, specific_qual_type))

    return annotations


def _set_units(model):
    '''Sets default units on an SBML Model.'''

    # kcat: per second
    kcat_unit_id = 'per_sec'
    unit_def = model.createUnitDefinition()
    unit_def.setId(kcat_unit_id)
    unit_def.setName(kcat_unit_id)
    unit = unit_def.createUnit()
    unit.setExponent(-1)
    unit.setKind(libsbml.UNIT_KIND_SECOND)

    # Concentration, K_M, K_I: M
    conc_unit_id = 'M'
    unit_def = model.createUnitDefinition()
    unit_def.setId(conc_unit_id)
    unit_def.setName(conc_unit_id)
    unit = unit_def.createUnit()
    unit.setKind(libsbml.UNIT_KIND_MOLE)
    unit = unit_def.createUnit()
    unit.setExponent(-1)
    unit.setKind(libsbml.UNIT_KIND_LITRE)

    # vmax: M per second
    v_unit_id = 'M_per_sec'
    unit_def = model.createUnitDefinition()
    unit_def.setId(v_unit_id)
    unit_def.setName(v_unit_id)
    unit = unit_def.createUnit()
    unit.setKind(libsbml.UNIT_KIND_MOLE)
    unit = unit_def.createUnit()
    unit.setExponent(-1)
    unit.setKind(libsbml.UNIT_KIND_LITRE)
    unit = unit_def.createUnit()
    unit.setExponent(-1)
    unit.setKind(libsbml.UNIT_KIND_SECOND)


def _get_reaction(reac_id):
    '''Gets a reaction.'''
    qry = 'MATCH (r:Reaction {id: {reac_id}})-[rel]-(c:Chemical) ' + \
        'RETURN r, rel, c'

    parameters = {'reac_id': reac_id}

    return biochem4j.run_query(qry, parameters)


def _parse(data, nodes, rels):
    '''Parses data.'''
    if len(data['errors']) > 0:
        raise ValueError(str(data['errors']))

    columns = data['results'][0]['columns']

    for datum in data['results'][0]['data']:
        from_node = None
        to_node = None
        rel = None

        for idx, meta_row in enumerate(zip(datum['meta'], datum['row'])):
            if meta_row[0]['type'] == 'node':
                cid = _get_id(meta_row[1]['id'])
                nodes[columns[idx]][cid] = meta_row[1]

                if from_node is None:
                    from_node = cid
                else:
                    to_node = cid

            elif meta_row[0]['type'] == 'relationship':
                rel = meta_row[1]

        if rel is not None:
            rels.append((from_node, rel, to_node))


def _get_id(cid):
    '''Gets an SBML-friendly id.'''
    return '_' + re.sub('\\W', '_', cid)


def _add_compartment(model, cid):
    '''Adds a compartment.'''
    cmpt = model.createCompartment()
    cmpt.setId(cid)
    cmpt.setSize(1)


def _add_species(model, cid, data, sbo=SBO_TERMS[SIMPLE_CHEM], conc=0):
    '''Adds a species.'''
    spec = model.createSpecies()
    _init_sbase(spec, cid, data, sbo)
    spec.setSBOTerm(sbo)

    if 'name' in data:
        spec.setName(str(data['name']))

    spec.setCompartment('c')
    spec.setInitialConcentration(conc)

    return spec


def _add_reaction(model, cid, data, reversible=False):
    '''Adds a reaction.'''
    react = model.createReaction()
    react.setReversible(reversible)

    _init_sbase(react, cid, data, SBO_TERMS[BIOCHEM_REACT])

    if not data.get('balance', False):
        print 'WARNING: Reaction %id unbalanced' % data['id']

    return react


def _add_species_ref(model, react_id, spec_id, stoic):
    '''Adds species reference.'''
    reaction = model.getReaction(react_id)

    if stoic > 0:
        reaction.addProduct(model.getSpecies(spec_id), stoic)

        for ref in reaction.getListOfProducts():
            ref.setSBOTerm(11)
    else:
        reaction.addReactant(model.getSpecies(spec_id), abs(stoic))

        for ref in reaction.getListOfReactants():
            ref.setSBOTerm(10)


def _add_modifier(model, react_id, uniprot_id):
    '''Adds a modifier.'''
    reaction = model.getReaction(react_id)
    cid = _get_id(uniprot_id)
    spec = model.getSpecies(cid)

    if spec is None:
        data = {'uniprot': uniprot_id}

        uniprot_vals = seq_utils.get_uniprot_values([uniprot_id],
                                                    ['protein names'])

        names = uniprot_vals[uniprot_id].get('Protein names', [])

        if len(names) > 0:
            data['name'] = names[0]

        spec = _add_species(model, cid, data, sbo=SBO_TERMS[PROTEIN], conc=1)

    reaction.addModifier(spec)

    for ref in reaction.getListOfModifiers():
        ref.setSBOTerm(SBO_TERMS[ENZYME])


def _init_sbase(sbase, cid, data, sbo):
    '''Initialises an sbase.'''
    sbase.setId(str(cid))
    sbase.setMetaId('_meta' + str(cid))
    sbase.setSBOTerm(sbo)
    _add_identifiers(data, sbase)


def _add_identifiers(properties, sbase):
    '''Gets semantic identifiers.'''
    ec_code = properties.get('ec', None)

    if ec_code:
        properties['ec-code'] = ec_code

    for key, value in properties.iteritems():

        for val in str(value).split(';'):
            url = 'http://identifiers.org/' + key + '/' + val

            if urllib.urlopen(url).getcode() == 200:
                _add_cv_term(url, sbase)


def _add_cv_term(url, sbase):
    '''Adds a CVTerm.'''
    cv_term = CVTerm()
    cv_term.setQualifierType(BIOLOGICAL_QUALIFIER)
    cv_term.setBiologicalQualifierType(BQB_IS)
    cv_term.addResource(str(url))
    sbase.addCVTerm(cv_term)


def main(args):
    '''main method'''
    document = get_document(args[1:])
    document.checkConsistency()
    document.printErrors()
    writeSBMLToFile(document, args[0])


if __name__ == '__main__':
    main(sys.argv[1:])
