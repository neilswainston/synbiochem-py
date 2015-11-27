'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-few-public-methods

from collections import defaultdict
from functools import partial
from itertools import count
import numpy
import random

from sklearn.metrics import classification_report, confusion_matrix
import theanets


class TheanetsBase(object):
    '''Base class for Classifier and Regressor.'''

    def __init__(self, network, optimize='sgd', learning_rate=0.01,
                 momentum=0.5):
        self._network = network
        self._optimize = optimize
        self._learning_rate = learning_rate
        self._momentum = momentum
        self._exp = None

    def _train(self, data, outputs, split=0.75, hidden_layers=None):
        '''Train the network.'''
        if hidden_layers is None:
            hidden_layers = [2]

        x_data = _pad(data[0])
        y_data = _pad(data[1])

        # Check lengths of x_data and y_data are equal:
        assert len(x_data) == len(y_data)

        layers = [len(x_data[0])] + hidden_layers + [outputs]
        self._exp = theanets.Experiment(self._network, layers=layers)
        x_data = numpy.array(x_data, dtype=numpy.float32)

        # Split data into training and validation:
        ind = int(split * len(x_data))
        self._exp.train((x_data[:ind], y_data[:ind]),
                        (x_data[ind:], y_data[ind:]),
                        optimize=self._optimize,
                        learning_rate=self._learning_rate,
                        momentum=self._momentum)


class Classifier(TheanetsBase):
    '''Simple classifier in Theanets.'''

    def __init__(self, optimize='sgd', learning_rate=0.01,
                 momentum=0.5):
        super(Classifier, self).__init__(theanets.Classifier, optimize,
                                         learning_rate, momentum)
        self.__y_map = None

    def train(self, x_data, y_data, split=0.75, hidden_layers=None):
        '''Train the network.'''
        y_enum = _enumerate(y_data)
        y_data = numpy.array([y[1] for y in y_enum], dtype=numpy.int32)
        self.__y_map = dict(set(y_enum))
        return super(Classifier, self)._train((x_data, y_data),
                                              len(self.__y_map), split,
                                              hidden_layers)

    def classify(self, x_test, y_test):
        '''Classifies and analyses test data.'''
        y_pred = self._exp.network.classify(_pad(x_test))

        y_test = numpy.array([self.__y_map[y]
                              for y in y_test], dtype=numpy.int32)

        inv_y_map = {v: k for k, v in self.__y_map.items()}

        return [inv_y_map[y] for y in y_pred], inv_y_map, \
            classification_report(y_test, y_pred), \
            confusion_matrix(y_test, y_pred)


class Regressor(TheanetsBase):
    '''Simple regressor in Theanets.'''

    def __init__(self, optimize='sgd', learning_rate=0.01,
                 momentum=0.5):
        super(Regressor, self).__init__(theanets.Regressor, optimize,
                                        learning_rate, momentum)

    def train(self, x_data, y_data, split=0.75, hidden_layers=None):
        '''Train the network.'''
        y_data = numpy.array(y_data, dtype=numpy.float32)
        return super(Regressor, self)._train((x_data, y_data), len(y_data[0]),
                                             split, hidden_layers)

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
        return [x + [mean_val] * (max_len - len(x)) for x in data]
    except TypeError:
        # For 1D data, elements cannot be passed to len:
        return data


def _enumerate(lst):
    '''Returns enumeration of supplied list.'''
    label_to_number = defaultdict(partial(next, count()))
    return [(item, label_to_number[item]) for item in lst]
