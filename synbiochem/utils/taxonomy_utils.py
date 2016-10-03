'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import json
import urllib
import urllib2


def get_taxonomy_id(name):
    '''Gets the NCBI Taxonomy id from supplied name.'''
    url = 'http://www.ebi.ac.uk/ols/api/search?ontology=ncbitaxon' + \
        '&exact=true&queryFields=label&q=' + urllib.quote(name)

    response = urllib2.urlopen(url)
    data = json.loads(response.read())

    if data['response']['numFound'] == 1:
        term = data['response']['docs'][0]['obo_id']
        return term[term.find(':') + 1:]

    raise RuntimeError(name + ' not found in NCBI Taxonomy')


def search(term, exact=False):
    '''Gets the NCBI Taxonomy id from supplied name.'''
    url = 'http://www.ebi.ac.uk/ols/api/search?ontology=ncbitaxon' + \
        '&exact=' + str(exact) + '&queryFields=label&q=' + urllib.quote(term)

    response = urllib2.urlopen(url)
    data = json.loads(response.read())

    return [{'id': doc['obo_id'].replace('NCBITaxon:', ''),
             'name': doc['label']}
            for doc in data['response']['docs']]
