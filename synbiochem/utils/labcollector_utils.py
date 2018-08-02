'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  alanwilliams
'''
# pylint: disable=missing-docstring
# pylint: disable=too-few-public-methods
# pylint: disable=useless-object-inheritance
import json

import requests


_API_KEY = 'X-LC-APP-Auth'
_CHARSET = 'X-LC-APP-Charset'
_CHARSET_VALUE = 'UTF-8'


class LabCollectorClient(object):
    '''Class representing a LabCollector client.'''

    def __init__(self, url, apikey):
        self.__url = url + '/webservice/v1'
        self.__headers = {'Accept': 'application/json',
                          'Content-Type': 'multipart/form-data',
                          _API_KEY: apikey,
                          _CHARSET: _CHARSET_VALUE}

    @property
    def animals(self):
        response = requests.get(
            self.__url + '/animals', headers=self.__headers)
        return json.loads(response.content)
