'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-locals
# pylint: disable=too-many-arguments

from collections import defaultdict
from functools import partial
from itertools import count
import numpy
import random

from sklearn.metrics import classification_report, confusion_matrix, f1_score
import theanets


class TheanetsBase(object):
    '''Base class for Classifier and Regressor.'''

    def __init__(self, network, data, outputs):
        self._network = network
        self._x_data = _pad(data[0])
        self._y_data = _pad(data[1])
        self._outputs = outputs
        self._exp = None

        # Check lengths of x_data and y_data are equal:
        assert len(self._x_data) == len(self._y_data)

    def train(self, split=0.75, hidden_layers=None, hyperparams=None):
        '''Train the network.'''
        if hidden_layers is None:
            hidden_layers = [1024]

        if hyperparams is None:
            hyperparams = {}

        layers = [len(self._x_data[0])] + hidden_layers + [self._outputs]
        self._exp = theanets.Experiment(self._network, layers=layers)

        # Split data into training and validation:
        ind = int(split * len(self._x_data))
        self._exp.train((self._x_data[:ind], self._y_data[:ind]),
                        (self._x_data[ind:], self._y_data[ind:]),
                        **hyperparams)


class Classifier(TheanetsBase):
    '''Simple classifier in Theanets.'''

    def __init__(self, x_data, y_data):
        y_enum = _enumerate(y_data)
        y_data = numpy.array([y[1] for y in y_enum], dtype=numpy.int32)
        self.__y_map = dict(set(y_enum))

        super(Classifier, self).__init__(theanets.Classifier, (x_data, y_data),
                                         len(self.__y_map))

    def classify(self, x_test, y_test):
        '''Classifies and analyses test data.'''
        y_pred = self._exp.network.classify(_pad(x_test))

        y_test = numpy.array([self.__y_map[y]
                              for y in y_test], dtype=numpy.int32)

        inv_y_map = {v: k for k, v in self.__y_map.items()}

        return [inv_y_map[y] for y in y_pred], inv_y_map, \
            classification_report(y_test, y_pred), \
            confusion_matrix(y_test, y_pred), f1_score(y_test, y_pred,
                                                       average='macro')


class Regressor(TheanetsBase):
    '''Simple regressor in Theanets.'''

    def __init__(self, x_data, y_data):
        super(Regressor, self).__init__(theanets.Regressor, (x_data, y_data),
                                        len(y_data[0]))

    def predict(self, x_test):
        '''Classifies and analyses test data.'''
        return self._exp.network.predict(_pad(x_test))


def randomise_order(x_data, y_data):
    '''Assumes x_data and y_data are paired (such that x_data[i] is paired with
    y_data[i]) and then randomises their orders such that this pairing is
    maintained.'''
    data = zip(x_data, y_data)
    random.shuffle(data)
    return zip(*data)


def _pad(data):
    '''Pad data with average values if sublists of different lengths.'''
    try:
        max_len = max([len(x) for x in data])
        mean_val = numpy.mean([x for sublist in data for x in sublist])
        return numpy.array([x + [mean_val] * (max_len - len(x)) for x in data],
                           dtype=numpy.float32)
    except TypeError:
        # For 1D data, elements cannot be passed to len:
        return data


def _enumerate(lst):
    '''Returns enumeration of supplied list.'''
    label_to_number = defaultdict(partial(next, count()))
    return [(item, label_to_number[item]) for item in lst]
