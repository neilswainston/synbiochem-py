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
from synbiochem.utils.sbol_utils import SBOLDocument


class ICEEntry(object):
    '''Class to represent an ICE entry.'''

    def __init__(self, uri_prefix=None, sbol_doc=None, typ=None,
                 metadata=None):

        assert typ is not None or metadata is not None

        if sbol_doc is None:
            self.__sbol_doc = SBOLDocument(uri_prefix)
        else:
            self.__sbol_doc = sbol_doc

        if metadata is None:
            self.__metadata = {'type': typ}
        else:
            self.__metadata = metadata

    def get_ice_number(self):
        '''Gets the ICE id.'''
        return self.__metadata['id'] if 'id' in self.__metadata else None

    def get_type(self):
        '''Gets the ICE type.'''
        return self.__metadata['type']

    def get_name(self):
        '''Gets the ICE name.'''
        return self.__metadata['name'] if 'name' in self.__metadata else ''

    def get_metadata(self):
        '''Gets the metadata.'''
        return self.__metadata

    def get_sbol_doc(self):
        '''Gets the SBOLDocument.'''
        return self.__sbol_doc

    def get_sequence(self):
        '''Gets the sequence.'''
        return self.__sbol_doc.get_sequence()

    def set_values(self, new_metadata):
        '''Sets multiple metadata values.'''
        self.__metadata.update(new_metadata)

    def set_value(self, key, value):
        '''Sets a metadata value.'''
        self.__metadata[key] = value

    def __repr__(self):
        return str(self.__metadata) + '\n' + str(self.__sbol_doc)


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

    def get_ice_entry(self, ice_id):
        '''Gets an ICEEntry object from the ICE database.'''
        return ICEEntry(sbol_doc=self.__get_sbol_doc(ice_id),
                        metadata=self.__get_meta_data(ice_id))

    def get_sequence(self, ice_id):
        '''Gets sequence from ICEEntry object from the ICE database.'''
        return self.get_ice_entry(ice_id).get_sequence()

    def get_name(self, ice_id):
        '''Gets name from ICEEntry object from the ICE database.'''
        return self.get_ice_entry(ice_id).get_name()

    def set_ice_entry(self, ice_entry):
        '''Saves an ICEEntry object in the ICE database.'''
        if ice_entry.get_ice_number() is None:
            self.__create_entry(ice_entry)
        else:
            self.__update_entry(ice_entry.get_ice_number(),
                                ice_entry.get_metadata())

        self.__upload_sbol(ice_entry.get_ice_number(),
                           ice_entry.get_type(),
                           ice_entry.get_sbol_doc())

    def __get_meta_data(self, ice_id):
        '''Returns an ICE entry metadata.'''
        return _read_resp(net_utils.get(
            self.__url + '/parts/' + self.__get_ice_number(ice_id),
            self.__headers))

    def __get_sbol_doc(self, ice_id):
        '''Gets the sequence ICE entry.'''
        url = self.__url + '/file/' + self.__get_ice_number(ice_id) + \
            '/sequence/sbol?sid=' + self.__sid
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        urllib.urlretrieve(url, temp_file.name)
        return SBOLDocument(filename=temp_file.name)

    def __create_entry(self, ice_entry):
        '''Creates a new ICE entry in the database.'''
        url = self.__url + '/parts'
        response = _read_resp(
            net_utils.post(url, json.dumps(ice_entry.get_metadata()),
                           self.__headers))

        ice_entry.set_values(response)

    def __update_entry(self, ice_id, metadata):
        '''Updates an ICE entry in the database.'''
        ice_number = self.__get_ice_number(ice_id)
        url = self.__url + '/parts/' + str(ice_number)
        _read_resp(net_utils.put(url, json.dumps(metadata), self.__headers))

    def __upload_sbol(self, ice_number, typ, sbol_doc):
        '''Uploads an SBOLDocument to ICE database.'''
        url = self.__url + '/file/sequence'
        files = {'entryRecordId': str(ice_number),
                 'entryType': typ,
                 'file': str(sbol_doc)}
        headers = self.__headers
        headers['Content-Type'] = 'multipart/form-data'
        return _read_resp(net_utils.post(url, files, headers))

    def __get_ice_number(self, ice_id):
        '''Maps ICE number to ICE id, i.e. from SBC000123 to 123.'''
        return str(int(ice_id.replace(self.__id_prefix, '')))

    def __get_ice_id(self, ice_number):
        '''Maps ICE id to ICE number, i.e. from 123 to SBC000123.'''
        return self.__id_prefix + format(ice_number, '06')


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

    ice = ICEClient(url, username, password, id_prefix='TEST_')
    print ice.get_ice_entry('TEST_000008').get_sbol_doc().get_sequence()

    sbol_doc = SBOLDocument(
        filename='/Users/neilswainston/Downloads/ahpC L177Q.xml')
    ice_entry = ICEEntry(sbol_doc=sbol_doc, typ='PART')
    ice.set_ice_entry(ice_entry)

    ice_entry = ICEEntry(uri_prefix='http://synbiochem.co.uk#', typ='PART')
    print ice_entry.get_sbol_doc().get_sequence()
    print ice_entry.get_sbol_doc()
    ice.set_ice_entry(ice_entry)
    print ice.get_ice_entry(ice_entry.get_ice_number())
    ice_entry.set_value('owner', 'Fabrizio Ravanelli')
    ice.set_ice_entry(ice_entry)
    print ice.get_ice_entry(ice_entry.get_ice_number())

if __name__ == '__main__':
    main()
