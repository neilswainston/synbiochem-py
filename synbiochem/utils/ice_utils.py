'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  pablocarbonell / neilswainston / alanwilliams
'''
# pylint: disable=too-many-arguments
import copy
import json
import tempfile
import urllib

from synbiochem.utils import net_utils as net_utils
from synbiochem.utils import sbol_utils as sbol_utils


class ICEClient(object):
    '''Class representing an ICE session.'''

    def __init__(self, url, username, password, id_prefix='SBC'):
        self.__url = url + '/rest'
        self.__headers = {'Accept': 'application/json',
                          'Content-Type': 'application/json'}
        resp = _read_resp(net_utils.post(self.__url + '/accesstoken',
                                         json.dumps({'email': username,
                                                      'password': password}),
                                         self.__headers))

        self.__sid = resp['sessionId']
        self.__headers[
            'X-ICE-Authentication-SessionId'] = self.__sid
        self.__id_prefix = id_prefix

    def get_meta_data(self, part_id):
        '''Returns an ICE entry metadata.'''
        return _read_resp(net_utils.get(
            self.__url + '/parts/' + self.__get_part_number(part_id),
            self.__headers))

    def get_sequence(self, part_id):
        '''Gets the sequence ICE entry.'''
        url = self.__url + '/file/' + self.__get_part_number(part_id) + \
            '/sequence/sbol?sid=' + self.__sid
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        urllib.urlretrieve(url, temp_file.name)
        return sbol_utils.SbolDocument(temp_file.name).get_dna_sequence()

    def create_part(self, data):
        '''Creates a new ICE part.'''
        url = self.__url + '/parts'
        return _read_resp(
            net_utils.post(url, json.dumps(data), self.__headers))

    def set_sequence(self, part_number, data):
        '''Sets an ICE part sequence.'''
        url = self.__url + '/parts/' + str(part_number) + '/sequence'
        return _read_resp(
            net_utils.post(url, json.dumps(data), self.__headers))

    def update_part(self, part_number, data):
        '''Updates an ICE part.'''
        url = self.__url + '/parts/' + str(part_number)
        return _read_resp(net_utils.put(url, json.dumps(data), self.__headers))

    # def get_sequence(self, part_number):
    #    '''Gets an ICE part sequence.'''
    #    url = self.__url + '/parts/' + str(part_number) + '/sequence'
    #    return _read_resp(net_utils.get(url, self.__headers))

    def upload_sequence(self, part_number, seqfile):
        '''Uploads an annotated sequence.'''
        url = self.__url + '/uploads/' + str(part_number) + '/sequence'
        files = {'file': open(seqfile, 'rb'), 'entryId': str(part_number)}
        return _read_resp(net_utils.post(url, files, self.__headers))

    def __get_part_number(self, part_id):
        '''Maps part number to part id, i.e. from SBC000123 to 123.'''
        return str(int(part_id.replace(self.__id_prefix, '')))


def get_template_newplasmid(name, description, creator, owner, investigator,
                            creator_email=None, owner_email=None,
                            investigator_email=None):
    '''Gets the template for a new ICE plasmid.'''
    return {u'accessPermissions': [],
            u'basePairCount': 0,
            u'bioSafetyLevel': 1,
            u'canEdit': True,
            u'creator': creator,
            u'creatorEmail': creator_email,
            u'featureCount': 0,
            u'fundingSource': u'BBSRC',
            u'hasAttachment': False,
            u'hasOriginalSequence': False,
            u'hasSample': False,
            u'hasSequence': False,
            u'index': 0,
            u'linkedParts': [],
            u'links': [],
            u'name': name,
            u'owner': owner,
            u'ownerEmail': owner_email,
            u'parameters': [],
            u'parents': [],
            u'plasmidData': {u'circular': True},
            u'principalInvestigator': investigator,
            u'principalInvestigatorEmail': investigator_email,
            u'publicRead': False,
            u'selectionMarkers': [],
            u'shortDescription': description,
            u'status': u'Planned',
            u'type': u'PLASMID',
            u'visible': u'OK'
            }


def template_newpart(name, description, creator, owner, investigator,
                     creator_email=None, owner_email=None,
                     investigator_email=None):
    '''Gets the template for a new ICE part.'''
    return {u'accessPermissions': [],
            u'basePairCount': 0,
            u'bioSafetyLevel': 1,
            u'canEdit': True,
            u'creator': creator,
            u'creatorEmail': creator_email,
            u'featureCount': 0,
            u'fundingSource': u'BBSRC',
            u'hasAttachment': False,
            u'hasOriginalSequence': False,
            u'hasSample': False,
            u'hasSequence': False,
            u'index': 0,
            u'linkedParts': [],
            u'links': [],
            u'name': name,
            u'owner': owner,
            u'ownerEmail': owner_email,
            u'parameters': [],
            u'parents': [],
            u'principalInvestigator': investigator,
            u'principalInvestigatorEmail': investigator_email,
            u'publicRead': False,
            u'selectionMarkers': [],
            u'shortDescription': description,
            u'status': u'Complete',
            u'type': u'PART',
            u'visible': u'OK'
            }


def _read_resp(response):
    '''Parses a string response into json.'''
    return json.loads(str(response))


def _join_sequence(seq1, seq2):
    '''Joins two sequences and keeps features.'''
    seq3 = copy.deepcopy(seq1)
    seq3['sequence'] = seq1['sequence'] + seq2['sequence']
    len1 = len(seq1['sequence'])
    for feature in seq2['features']:
        new_feature = copy.deepcopy(feature)
        for location in new_feature['locations']:
            location['genbankStart'] += len1
            location['end'] += len1
        seq3['features'].append(new_feature)
    return seq3


def main():
    '''main method.'''
    url = 'https://testice.synbiochem.co.uk:8443'
    username = 'Administrator'
    password = 'synbiochem'

    ice = ICEClient(url, username, password)
    print ice.get_meta_data('SBC00008')
    print ice.get_sequence('SBC8')


if __name__ == '__main__':
    main()
