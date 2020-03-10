'''
SYNBIOCHEM-DB (c) University of Manchester 2015

SYNBIOCHEM-DB is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import json
import requests


def get_synonyms_by_id(tax_id):
    '''Gets a taxonomy's children.'''
    qry = 'MATCH (o:Organism {taxonomy: {tax_id}}) RETURN o'

    parameters = {'tax_id': tax_id}

    results = _parse(_run_query(qry, parameters))
    return results[0]['names'] if results else ''


def get_children_by_id(tax_id):
    '''Gets a taxonomy's children.'''
    qry = 'MATCH (p:Organism {taxonomy: {tax_id}})<--(c:Organism) ' + \
        'RETURN c'

    parameters = {'tax_id': tax_id}

    return _parse(_run_query(qry, parameters))


def get_parent_by_id(tax_id):
    '''Gets a taxonomy parent.'''
    qry = 'MATCH (c:Organism {taxonomy: {tax_id}})-->(p:Organism) ' + \
        'RETURN p'

    parameters = {'tax_id': tax_id}

    results = _parse(_run_query(qry, parameters))
    return results[0]


def get_parent_by_name(tax_name):
    '''Gets a taxonomy parent.'''

    qry = 'MATCH (c:Organism)-->(p:Organism) WHERE \'' + \
        tax_name + \
        '\' IN c.names RETURN p'

    parameters = {}

    results = _parse(_run_query(qry, parameters))
    return results[0]


def _parse(data):
    '''Parses data.'''
    if data['errors']:
        raise ValueError(str(data['errors']))

    return [datum['row'][0] for datum in data['results'][0]['data']]


def _run_query(query, parameters):
    '''Gets a reaction.'''
    url = 'http://www.biochem4j.org/db/data/transaction/commit'

    data = {'statements': [
        {
            'statement': query,
            'parameters': parameters
        }
    ]}

    headers = {'Accept': 'application/json; charset=UTF-8',
               'Content-Type': 'application/json'}

    resp = requests.post(url, data=json.dumps(data), headers=headers)

    return json.loads(resp.text)
