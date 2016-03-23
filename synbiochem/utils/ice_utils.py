'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  pablocarbonell / neilswainston / alanwilliams
'''
# pylint: disable=too-many-arguments
import json
import tempfile
import urllib
import requests

import sbol

from synbiochem.utils import net_utils as net_utils

_SESSION_KEY = 'X-ICE-Authentication-SessionId'


class ICEEntry(object):
    '''Class to represent an ICE entry.'''

    def __init__(self, sbol_doc=None, typ=None, metadata=None):

        assert typ is not None or metadata is not None

        if sbol_doc is None:
            self.__sbol_doc = sbol.Document()
        else:
            self.__sbol_doc = sbol_doc

        if metadata is None:
            self.__metadata = {'type': typ}
        else:
            self.__metadata = metadata

    def get_ice_number(self):
        '''Gets the ICE number.'''
        return self.__metadata['id'] if 'id' in self.__metadata else None

    def get_record_id(self):
        '''Gets the ICE record id.'''
        return self.__metadata['recordId'] if 'recordId' in self.__metadata \
            else None

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
        self.__headers[_SESSION_KEY] = self.__sid
        self.__id_prefix = id_prefix

    def get_ice_entry(self, ice_id):
        '''Gets an ICEEntry object from the ICE database.'''
        return ICEEntry(sbol_doc=self.__get_sbol_doc(ice_id),
                        metadata=self.__get_meta_data(ice_id))

    def get_sbol_doc(self, ice_id):
        '''Gets SBOLDocument from ICEEntry object from the ICE database.'''
        return self.get_ice_entry(ice_id).get_sbol_doc()

    def set_ice_entry(self, ice_entry):
        '''Saves an ICEEntry object in the ICE database.'''
        if ice_entry.get_ice_number() is None:
            self.__create_entry(ice_entry)
        else:
            self.__update_entry(ice_entry.get_ice_number(),
                                ice_entry.get_metadata())

        return self.__upload_sbol(ice_entry.get_record_id(),
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
        document = sbol.Document()
        document.read(temp_file.name)
        return document

    def __create_entry(self, ice_entry):
        '''Creates a new ICE entry in the database.'''
        url = self.__url + '/parts'
        response = _read_resp(
            net_utils.post(url, json.dumps(ice_entry.get_metadata()),
                           self.__headers))

        metadata = self.__get_meta_data(self.__get_ice_id(response['id']))
        ice_entry.set_values(metadata)

    def __update_entry(self, ice_id, metadata):
        '''Updates an ICE entry in the database.'''
        ice_number = self.__get_ice_number(ice_id)
        url = self.__url + '/parts/' + str(ice_number)
        _read_resp(net_utils.put(url, json.dumps(metadata), self.__headers))

    def __upload_sbol(self, record_id, typ, sbol_doc):
        '''Uploads an SBOLDocument to ICE database.'''
        sbol_file = tempfile.NamedTemporaryFile()
        sbol_doc.write(sbol_file.name)

        url = self.__url + '/file/sequence'
        response = requests.post(url,
                                 headers={_SESSION_KEY: self.__sid},
                                 files={'file': open(sbol_file.name, 'r'),
                                        'entryType': typ,
                                        'entryRecordId': record_id},
                                 verify=False)

        return _read_resp(response)

    def __get_ice_number(self, ice_id):
        '''Maps ICE number to ICE id, i.e. from SBC000123 to 123.'''
        return str(int(ice_id.replace(self.__id_prefix, '')))

    def __get_ice_id(self, ice_number):
        '''Maps ICE id to ICE number, i.e. from 123 to SBC000123.'''
        return self.__id_prefix + format(ice_number, '06')


def _read_resp(response):
    '''Parses a string response into json.'''
    return json.loads(str(response))


def main():
    '''main method.'''
    url = 'https://testice.synbiochem.co.uk:8443'
    username = 'Administrator'
    password = 'synbiochem'

    ice = ICEClient(url, username, password, id_prefix='TEST_')
    sbol_doc = ice.get_ice_entry('TEST_000008').get_sbol_doc()
    print sbol_doc

    ice_entry = ICEEntry(sbol_doc=sbol_doc, typ='PART')
    ice.set_ice_entry(ice_entry)

    ice_entry = ICEEntry(typ='PART')
    print ice_entry.get_sbol_doc().sequences
    print ice_entry.get_sbol_doc()
    ice.set_ice_entry(ice_entry)
    print ice.get_ice_entry(ice_entry.get_ice_number())
    ice_entry.set_value('owner', 'Fabrizio Ravanelli')
    ice.set_ice_entry(ice_entry)
    print ice.get_ice_entry(ice_entry.get_ice_number())

if __name__ == '__main__':
    main()
