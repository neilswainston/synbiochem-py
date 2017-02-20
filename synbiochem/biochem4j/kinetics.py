'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys
import uuid

import libsbml

from synbiochem.biochem4j import model_utils


_SBO_TERM_DEFAULT = {model_utils.SBO_TERMS[model_utils.SIMPLE_CHEM]: 1e-4,
                     model_utils.SBO_TERMS[model_utils.PROTEIN]: 1e-6,
                     model_utils.SBO_TERMS[model_utils.KCAT]: 10,
                     model_utils.SBO_TERMS[model_utils.KCAT_FOR]: 10,
                     model_utils.SBO_TERMS[model_utils.KCAT_REV]: 10,
                     model_utils.SBO_TERMS[model_utils.K_M]: 1e-4}

_SPECIES_TO_IGNORE = ['http://identifiers.org/chebi/CHEBI:15377',
                      'http://identifiers.org/chebi/CHEBI:15378',
                      'http://identifiers.org/chebi/CHEBI:16526',
                      'http://identifiers.org/chebi/CHEBI:24636',
                      'http://identifiers.org/chebi/CHEBI:28938',
                      'http://identifiers.org/chebi/CHEBI:33019',
                      'http://identifiers.org/chebi/CHEBI:456215']

_FORMULA_TO_ID = {}


def add_kinetics(document):
    '''Adds convenience kinetics to all reactions in the document.'''
    model = document.getModel()

    for species in model.getListOfSpecies():
        for annotation in model_utils.get_annotations(species):
            if annotation[0] in _SPECIES_TO_IGNORE:
                species.setConstant(True)
                species.setBoundaryCondition(True)

    for reaction in model.getListOfReactions():
        _add_kinetics_reaction(model, reaction)


def _add_kinetics_reaction(model, reaction):
    '''Adds a convenience kinetic law to a given reaction.'''
    is_reversible = reaction.isSetReversible() and reaction.getReversible()

    enzymes = [modifier for modifier in reaction.getListOfModifiers()
               if modifier.getSBOTerm() ==
               model_utils.SBO_TERMS[model_utils.ENZYME]]

    formula, parameters, num_react, num_prods = \
        _get_formula(model,
                     reaction.getListOfReactants(),
                     reaction.getListOfProducts()
                     if is_reversible else [],
                     model.getSpecies(enzymes[0].getSpecies())
                     if enzymes else None)

    print formula

    func_def = _get_func_def(model,
                             formula,
                             [parameter[0] for parameter in parameters],
                             num_react,
                             num_prods)

    _set_kinetic_law(reaction, func_def, parameters)


def _get_formula(model, reactants, products, enzyme):
    '''Returns formula, react_terms, prod_terms
    for supplied number of reactants and products.'''
    react_terms = _get_terms(model, reactants, 'S')
    prod_terms = _get_terms(model, products, 'P')
    irreversible = len(prod_terms) == 0

    react_numer_terms = model_utils.KCAT_FOR + ' * ' + \
        _get_numer_terms(react_terms)
    react_denom_terms = _get_denom_terms(react_terms)

    prod_numer_terms = '' if irreversible \
        else ' - ( ' + model_utils.KCAT_REV + ' * ' + \
        _get_numer_terms(prod_terms) + ' )'
    prod_denom_terms = '' if irreversible \
        else ' + ( ' + _get_denom_terms(prod_terms) + ' ) - 1'

    numer = '( ' + react_numer_terms + prod_numer_terms + ' )'
    denom = '( ' + react_denom_terms + prod_denom_terms + ' )'

    formula = model_utils.VOLUME + ' * ' + \
        model_utils.ENZYME_CONC + ' * ' + numer + ' / ' + denom

    params = []
    params.append((model_utils.VOLUME,
                   model_utils.SBO_TERMS[model_utils.VOLUME],
                   enzyme.getCompartment(),
                   None))
    params.append((model_utils.ENZYME_CONC,
                   model_utils.SBO_TERMS[model_utils.CONC],
                   enzyme.getId(),
                   None))
    params.append((model_utils.KCAT_FOR,
                   model_utils.SBO_TERMS[model_utils.KCAT_FOR],
                   model_utils.KCAT_FOR,
                   None))

    if not irreversible:
        params.append((model_utils.KCAT_REV,
                       model_utils.SBO_TERMS[model_utils.KCAT_REV],
                       model_utils.KCAT_REV,
                       None))

    params.extend(_get_parameters(react_terms))
    params.extend(_get_parameters(prod_terms))

    return formula, params, len(react_terms), len(prod_terms)


def _get_terms(model, participants, prefix):
    ''''Get list of tuples of (id, stoichiometry, parameter_id,
    sbo_term).'''
    terms = [[participant.getSpecies(),
              participant.getStoichiometry(),
              (prefix + str(i + 1), model_utils.SBO_TERMS[model_utils.CONC]),
              ('KM_' + prefix + str(i + 1),
               model_utils.SBO_TERMS[model_utils.K_M])]
             for i, participant in enumerate(participants)
             if not model.getSpecies(participant.getSpecies()).getConstant()]

    return terms


def _get_numer_terms(terms):
    '''Returns numerator terms in the form S1/KM_S1 * S2/KM_S2.'''
    lst = [['( ' + term[2][0] + ' / ' + term[3][0] + ' )'] * int(term[1])
           for term in terms]
    return ' * '.join([item for sublist in lst for item in sublist])


def _get_denom_terms(terms):
    '''Returns denominator terms in the form
    (((S1/KM_S1)^0 + (S1/KM_S1)^1)) * ((S2/KM_S2)^0 + (S2/KM_S2)^1)).'''
    lst = [' + '.join(['( ( ' + term[2][0] + ' / ' + term[3][0] + ' ) ^ ' +
                       str(x) + ' )'
                       for x in range(int(term[1]) + 1)])
           for term in terms]
    return '( ' + ' ) * ( '.join(lst) + ' )'


def _get_func_def(model, formula, parameters, num_reacts, num_prods):
    '''Gets existing or creates new functionDefinition from given
    parameters.'''
    if formula in _FORMULA_TO_ID:
        formula_id = _FORMULA_TO_ID[formula]
        function_definition = model.getFunctionDefinition(formula_id)
    else:
        function_definition = model.createFunctionDefinition()
        function_definition.setId(_get_unique_id())
        function_definition.setName(_get_func_name(num_reacts, num_prods))
        function_definition.setMath(_get_math(formula, parameters))
        _FORMULA_TO_ID[formula] = function_definition.getId()

    return function_definition


def _get_parameters(terms):
    '''Gets parameters derived from terms.'''
    lst = [[(term[2][0], term[2][1], term[0], None),
            (term[3][0], term[3][1], 'KM_' + term[0], term[0])]
           for term in terms]
    return [item for sublist in lst for item in sublist]


def _set_kinetic_law(reaction, function_definition, parameters):
    '''Sets kineticLaw element to reaction.'''

    mathml = '<?xml version="1.0" encoding="UTF-8"?>'
    mathml += '<math xmlns="http://www.w3.org/1998/Math/MathML">'
    mathml += '<apply>'
    mathml += '<ci>' + function_definition.getId() + '</ci>'

    for parameter in parameters:
        mathml += '<ci>' + parameter[2] + '</ci>'

    mathml += '</apply>'
    mathml += '</math>'

    kinetic_law = reaction.createKineticLaw()
    kinetic_law.setMath(libsbml.readMathMLFromString(mathml))

    for parameter in parameters:
        sbo_term = parameter[1]
        if sbo_term != model_utils.SBO_TERMS[model_utils.CONC] and \
                sbo_term != model_utils.SBO_TERMS[model_utils.VOLUME]:
            param = kinetic_law.createParameter()
            param.setId(parameter[2])
            param.setValue(_SBO_TERM_DEFAULT[sbo_term])
            param.setUnits(model_utils.UNITS[sbo_term])
            param.setSBOTerm(sbo_term)

    return kinetic_law


def _get_func_name(num_reactants, num_products):
    '''Returns a function name based on the number of reactants and
    products.'''
    is_reversible = num_products > 0
    reversible = 'reversible' if is_reversible else 'irreversible'
    name = 'Convenience (' + reversible + '): ' + \
        str(num_reactants) + ' reactants'

    if is_reversible:
        name += ', ' + str(num_products) + ' products'

    return name


def _get_math(formula, parameters):
    '''Returns math element from formula and parameters.'''
    math_elem = '<math xmlns="http://www.w3.org/1998/Math/MathML">'
    param_str = math_elem + '<lambda>'

    for parameter in parameters:
        param_str += '<bvar><ci>' + parameter + '</ci></bvar>'

    mathml = libsbml.writeMathMLToString(libsbml.parseFormula(formula))
    mathml = mathml.replace(math_elem, param_str)
    mathml = mathml.replace('</math>', '</lambda></math>')
    return libsbml.readMathMLFromString(mathml)


def _get_unique_id():
    '''Returns unique and SBML-valid id.'''
    return '_' + uuid.uuid4().get_hex().replace('-', '_')


def main(args):
    '''main method'''
    document = libsbml.readSBMLFromFile(args[0])
    add_kinetics(document)
    libsbml.writeSBMLToFile(document, args[1])


if __name__ == '__main__':
    main(sys.argv[1:])
