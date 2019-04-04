'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  pablocarbonell / neilswainston / alanwilliams
'''
# pylint: disable=bad-option-value
# pylint: disable=superfluous-parens
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes
import codecs
import copy
import json
import tempfile
import threading
import time
import traceback
from xml.etree.ElementTree import ParseError

from synbiochem.utils import dna_utils, net_utils, sbol_utils


_DEFAULT_ID_PREFIX = 'SBC'
_SESSION_KEY = 'X-ICE-Authentication-SessionId'


class ICEClientFactory():
    '''Factory class for controlling ICE factories.'''

    def __init__(self):
        self.__ice_clients = {}
        self.__running = True
        threading.Thread(target=self.__check_clients).start()

    def close(self):
        '''Close factory.'''
        self.__running = False

    def get_ice_client(self, url, username, psswrd,
                       id_prefix=_DEFAULT_ID_PREFIX, group_names=None):
        '''Get ICEClient.'''
        if group_names is None:
            group_names = []

        params = (url, username, psswrd, id_prefix, tuple(group_names))

        if params in self.__ice_clients:
            return self.__ice_clients[params][0]

        ice_client = ICEClient(url, username, psswrd, id_prefix, group_names)
        self.__ice_clients[params] = [ice_client, time.time()]

        return ice_client

    def __check_clients(self, max_age=60 * 15, sleep_time=60):
        '''Check ICE clients.'''
        while self.__running:
            to_remove = []

            for ice_params, values in self.__ice_clients.items():
                if time.time() - values[1] > max_age:
                    to_remove.append(ice_params)

            for ice_params in to_remove:
                del self.__ice_clients[ice_params]

            print('ICE clients: %s' % len(self.__ice_clients))
            time.sleep(sleep_time)


class ICEEntry():
    '''Class to represent an ICE entry.'''

    def __init__(self, dna=None, typ=None, metadata=None):

        assert typ is not None or metadata is not None

        self.__dna = dna
        self.__dna_updated = dna is not None

        if metadata is None:
            self.__metadata = {'type': typ}
        else:
            if 'type' not in metadata:
                metadata['type'] = typ

            self.__metadata = metadata

    def get_ice_number(self):
        '''Gets the ICE number.'''
        return self.__metadata['id'] if 'id' in self.__metadata else None

    def get_ice_id(self):
        '''Gets the ICE id.'''
        ice_number = self.get_ice_number()
        return get_ice_id(ice_number) if ice_number is not None else None

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

    def set_parameter(self, key, value):
        '''Sets a parameter.'''
        parameters = self.get_parameters()
        parameters[key] = value

        self.set_value('parameters',
                       [{'name': key, 'value': value}
                        for key, value in parameters.items()])

    def get_parameters(self):
        '''Gets parameters.'''
        if 'parameters' in self.__metadata:
            return {parameter['name']: parameter['value']
                    for parameter in self.__metadata['parameters']}

        return {}

    def get_parameter(self, key):
        '''Gets parameter.'''
        parameters = self.get_parameters()
        return parameters.get(key, None)

    def get_metadata(self):
        '''Gets the metadata.'''
        return self.__metadata

    def get_dna(self):
        '''Gets the DNA object.'''
        return self.__dna

    def get_seq(self):
        '''Gets the sequence.'''
        return self.__dna['seq'] if self.__dna else ''

    def get_dna_updated(self):
        '''Gets the DNA object updated flag.'''
        return self.__dna_updated

    def set_values(self, new_metadata):
        '''Sets multiple metadata values.'''
        self.__metadata.update(new_metadata)

    def set_value(self, key, value):
        '''Sets a metadata value.'''
        self.__metadata[key] = value

    def set_dna(self, dna):
        '''Sets the DNA object.'''
        self.__dna_updated = self.__dna is not None \
            or dna is not None
        self.__dna = dna

    def unset_dna_updated(self):
        '''Sets the DNA object updated flag.'''
        self.__dna_updated = False

    def copy(self):
        '''Copy ICEEntry, making a duplicate.'''
        metadata = copy.deepcopy(self.__metadata)
        metadata.pop('id')
        metadata.pop('partId')
        metadata.pop('recordId')
        return ICEEntry(self.__dna, metadata.get('type', None), metadata)

    def __repr__(self):
        return str(self.__metadata) + \
            ('\n' + str(self.__dna)
             if self.__dna is not None
             else '')


class ICEClient():
    '''Class representing an ICE client.'''

    def __init__(self, url, username, psswrd, id_prefix=_DEFAULT_ID_PREFIX,
                 group_names=None):
        self.__url = url[:-1] if url[-1] == '/' else url
        self.__username = username
        self.__psswrd = psswrd
        self.__id_prefix = id_prefix

        if group_names is None:
            group_names = []

        self.__headers = {'Accept': 'application/json',
                          'Content-Type': 'application/json'}

        self.__sid, self.__user, self.__email = self.reconnect()
        self.__headers[_SESSION_KEY] = self.__sid

        self.__group_ids = [group_id
                            for name, group_id
                            in self.get_groups().items()
                            if name in group_names]

    def reconnect(self):
        '''Reconnects to ICE server.'''
        try:
            resp = self.__get_access_token(
                '/accesstoken', self.__username, self.__psswrd)
        except net_utils.NetworkError:
            resp = self.__get_access_token(
                '/accesstokens', self.__username, self.__psswrd)

        return resp['sessionId'], \
            resp['firstName'] + ' ' + resp['lastName'], \
            resp['email']

    def get_ice_entry(self, ice_id):
        '''Gets an ICEEntry object from the ICE database.'''
        metadata = self.__get_meta_data(ice_id)

        try:
            dna = self.__get_dna(ice_id) if metadata['hasSequence'] \
                else None
        except ParseError:
            # Take no action
            dna = None

        return ICEEntry(dna, metadata=metadata)

    def set_ice_entry(self, ice_entry):
        '''Saves an ICEEntry object in the ICE database.'''
        if ice_entry.get_ice_number() is None:
            response = self.__create_entry(ice_entry)
        else:
            response = self.__update_entry(ice_entry.get_ice_number(),
                                           self.__form_metadata(ice_entry))

        metadata = self.__get_meta_data(self.__get_ice_id(response['id']))
        ice_entry.set_values(metadata)

        if ice_entry.get_dna_updated():

            if 'hasSequence' in metadata and metadata['hasSequence']:
                self.__delete_seq(ice_entry.get_ice_number())

            dna = ice_entry.get_dna()

            if dna is not None:
                self.__upload_dna(ice_entry.get_record_id(),
                                  ice_entry.get_type(),
                                  dna)

                ice_entry.unset_dna_updated()

        metadata = self.__get_meta_data(self.__get_ice_id(response['id']))
        ice_entry.set_values(metadata)

        for group_id in self.__group_ids:
            self.add_permission(ice_entry.get_ice_id(), group_id)

        return response['id']

    def do_blast(self, seq, max_num=1024, sort_field='RELEVANCE'):
        '''Performs BLAST search against database.'''
        data = {'blastQuery': {'blastProgram': 'BLAST_N',
                               'sequence': seq.lower()},
                'parameters': {'retrieveCount': max_num,
                               'sortField': sort_field}}
        return _read_resp(net_utils.post(self.__url + '/rest/search',
                                         json.dumps(data), self.__headers))

    def get_ice_entries_by_seq(self, seq):
        '''Returns entries matching supplied sequence.'''
        entries = []

        try:
            response = self.do_blast(seq)

            if 'results' in response:
                for result in response['results']:
                    if result['queryLength'] == int(result['alignment']):
                        entry = self.get_ice_entry(result['entryInfo']['id'])

                        if entry.get_seq().lower() == seq.lower():
                            entries.append(entry)
        except net_utils.NetworkError:
            # Sometime BLAST fails inexplicably (on ICE side).
            # Report error, return empty entries.
            traceback.print_exc()

        return entries

    def add_link(self, parent_id, child_id):
        '''Adds a link between a parent and child.'''
        url = self.__url + '/rest/parts/' + get_ice_number(parent_id) + \
            '/links'
        data = {'id': get_ice_number(child_id)}
        return _read_resp(net_utils.post(url,
                                         data=json.dumps(data),
                                         headers=self.__headers))

    def get_groups(self):
        '''Gets a group name to id map.'''
        url = self.__url + '/rest/groups?offset=0&limit=4096'
        groups = _read_resp(net_utils.get(url, headers=self.__headers))
        return {group['label']: group['id'] for group in groups['data']}

    def search_groups(self, term):
        '''Gets groups from search terms.'''
        url = self.__url + '/rest/groups/autocomplete?token=' + term
        return _read_resp(net_utils.get(url, headers=self.__headers))

    def add_permission(self, ice_id, group_id, read=True):
        '''Adds user permissions to a given ICE entry.'''
        url = self.__url + '/rest/parts/' + self.__get_ice_number(ice_id) + \
            '/permissions'
        data = {'type': 'READ_ENTRY' if read else 'WRITE_ENTRY',
                'article': 'GROUP',
                'articleId': group_id}

        return _read_resp(net_utils.post(url, json.dumps(data),
                                         self.__headers))

    def search(self, term, limit=5):
        '''Searches ICE.'''
        url = self.__url + '/rest/search?offset=0&limit=' + str(limit) + \
            '&sort=relevance&q="' + term + '"'

        return _read_resp(net_utils.get(url, self.__headers))

    def advanced_search(self, term, entry_type, limit=5):
        '''Searches ICE.'''
        data = {
            'queryString': term,
            'entryTypes': [entry_type],
            'parameters': {
                'sortField': 'RELEVANCE',
                'retrieveCount': limit,
            }
        }

        return _read_resp(net_utils.post(self.__url + '/rest/search',
                                         json.dumps(data),
                                         self.__headers))

    def search_name(self, name, entry_type, limit=5):
        '''Search by name and type.'''
        resp = self.advanced_search(name, entry_type, limit)

        return [result for result in resp['results']
                if name in result['entryInfo']['name']]

    def search_design(self, design_number):
        '''Search plasmids by design.'''
        return self.search_name('SBCDE%05d_PL' % (int(design_number)),
                                'PLASMID', limit=128)

    def get_genbank(self, ice_id, out=None):
        '''Get Genbank file.'''
        url = self.__url + '/rest/file/' + self.__get_ice_number(ice_id) + \
            '/sequence/genbank'
        genbank = net_utils.get(url, self.__headers)

        if out:
            with open(out, 'w') as out_file:
                out_file.write(genbank)

        return genbank

    def __get_access_token(self, service, username, psswrd):
        '''Gets access token response.'''
        return _read_resp(net_utils.post(self.__url + '/rest' + service,
                                         json.dumps({'email': username,
                                                     'password': psswrd}),
                                         self.__headers))

    def __form_metadata(self, ice_entry):
        '''Forms metadata dictionary.'''
        metadata = ice_entry.get_metadata()

        if 'creator' not in metadata:
            metadata['creator'] = self.__user

        if 'creatorEmail' not in metadata:
            metadata['creatorEmail'] = self.__email

        return metadata

    def __upload_seq_file(self, record_id, typ, filename):
        '''Uploads a sequence file (not necessarily SBOL).'''
        return _read_resp(net_utils.post_file(self.__url +
                                              '/rest/file/sequence',
                                              {'file': open(filename, 'r'),
                                               'entryType': typ,
                                               'entryRecordId': record_id},
                                              {_SESSION_KEY: self.__sid}))

    def __delete_seq(self, ice_number):
        '''Deletes the sequence associated with supplied record id.'''
        net_utils.delete(self.__url + '/rest/parts/' + str(ice_number) +
                         '/sequence',
                         headers={_SESSION_KEY: self.__sid})

    def __get_meta_data(self, ice_id):
        '''Returns an ICE entry metadata.'''
        return _read_resp(net_utils.get(
            self.__url + '/rest/parts/' + self.__get_ice_number(ice_id),
            self.__headers))

    def __get_dna(self, ice_id):
        '''Gets the sequence ICE entry.'''
        url = self.__url + '/rest/file/' + self.__get_ice_number(ice_id) + \
            '/sequence/sbol1?sid=' + self.__sid
        temp_file = tempfile.NamedTemporaryFile(delete=False)

        with codecs.open(temp_file.name, 'w', 'utf-8') as text_file:
            text_file.write(net_utils.get(url))

        return sbol_utils.read(temp_file.name)

    def __create_entry(self, ice_entry):
        '''Creates a new ICE entry in the database.'''
        url = self.__url + '/rest/parts'
        return _read_resp(
            net_utils.post(url, json.dumps(self.__form_metadata(ice_entry)),
                           self.__headers))

    def __update_entry(self, ice_id, metadata):
        '''Updates an ICE entry in the database.'''
        ice_number = self.__get_ice_number(ice_id)
        url = self.__url + '/rest/parts/' + str(ice_number)

        return _read_resp(net_utils.put(url, json.dumps(metadata),
                                        self.__headers))

    def __upload_dna(self, record_id, typ, dna):
        '''Uploads an SBOLDocument to ICE database.'''
        sbol_file = tempfile.NamedTemporaryFile(suffix='.xml')
        sbol_utils.write(dna, sbol_file.name)
        return self.__upload_seq_file(record_id, typ, sbol_file.name)

    def __get_ice_number(self, ice_identifier):
        '''Maps ICE number to ICE id, i.e. from SBC000123 to 123,
        or if a number is supplied, returns the number.'''
        return get_ice_number(ice_identifier, self.__id_prefix)

    def __get_ice_id(self, ice_number):
        '''Maps ICE id to ICE number, i.e. from 123 to SBC000123.'''
        return get_ice_id(ice_number, self.__id_prefix)


class DNAWriter():
    '''Class for writing DNA objects to ICE.'''

    def __init__(self, ice_client):
        self.__ice_client = ice_client
        self.__cache = {}

    def submit(self, dna):
        '''Submits DNA to ICE (or checks for existing entry.'''
        if dna['seq'] in self.__cache:
            return self.__cache[dna['seq']]

        try:
            ice_entry = self.__ice_client.get_ice_entry(dna['disp_id'])
            ice_id, typ = ice_entry.get_ice_id(), ice_entry.get_type()
        except ValueError:
            # If disp_id is not a valid ICE id, try BLAST:
            ice_entries = self.__ice_client.get_ice_entries_by_seq(dna['seq'])

            if ice_entries:
                ice_entry = ice_entries[0]
                ice_id, typ = ice_entry.get_ice_id(), ice_entry.get_type()
            else:
                ice_id, typ = self.__write(dna)

        self.__cache[dna['seq']] = (ice_id, typ)
        return ice_id, typ

    def __write(self, dna):
        '''Writes DNA document to ICE.'''
        typ = 'PLASMID' if dna.get('typ', None) == dna_utils.SO_PLASMID \
            else 'PART'

        dna['parameters']['Full name'] = dna['name']

        ice_entry = ICEEntry(dna, typ)

        ice_entry.set_value('name', dna['name'][:125])
        ice_entry.set_value('shortDescription', dna['desc'])

        _add_params(ice_entry, dna)
        links = set(dna['links'])

        if dna.get('typ', None):
            links.add(dna['typ'])

        for feature in dna['features']:
            _add_params(ice_entry, feature)
            links |= set(feature['links'])

        ice_entry.set_value('links', list(links))

        entry_id = self.__ice_client.set_ice_entry(ice_entry)

        for child in dna['children']:
            par_ice_entry, _ = self.submit(child)
            self.__ice_client.add_link(entry_id, par_ice_entry)

        return entry_id, ice_entry.get_type()


def get_ice_number(ice_identifier, id_prefix=_DEFAULT_ID_PREFIX):
    '''Maps ICE number to ICE id, i.e. from SBC000123 to 123,
    or if a number is supplied, returns the number.'''
    try:
        ice_number = int(ice_identifier.replace(id_prefix, ''))
    except AttributeError:
        # "Ask forgiveness, not permission" and assume ice_identifier is
        # the ice_number:
        ice_number = ice_identifier

    return str(ice_number)


def get_ice_id(ice_number, id_prefix=_DEFAULT_ID_PREFIX):
    '''Maps ICE id to ICE number, i.e. from 123 to SBC000123.'''
    return id_prefix + format(int(get_ice_number(ice_number)), '06')


def _read_resp(response):
    '''Parses a string response into json.'''
    return json.loads(response, encoding='utf-8') if response else None


def _add_params(ice_entry, dna):
    '''Adds parameter values to ICEENtry.'''
    for key, value in dna['parameters'].items():
        ice_entry.set_parameter(key, (', '.join([str(val) for val in value])
                                      if isinstance(value, list)
                                      else value))
