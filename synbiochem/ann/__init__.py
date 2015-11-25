'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member

from collections import defaultdict
from functools import partial
from itertools import count
import numpy
import random

from sklearn.metrics import classification_report, confusion_matrix
import theanets


class Classifier(object):
    '''Simple classifier in Theanets.'''

    def __init__(self, optimize='sgd', learning_rate=0.01,
                 momentum=0.5):
        self.__optimize = optimize
        self.__learning_rate = learning_rate
        self.__momentum = momentum
        self.__exp = None
        self.__y_map = None

    def train(self, x_data, y_data, split=0.75, hidden_layers=None):
        '''Train the network. Accepts list of floats as x_data,
        and list of anything as y_data.'''
        if hidden_layers is None:
            hidden_layers = [2]

        # Check lengths of x_data and y_data are equal,
        # assume all tuples in x_data are of the same length.
        assert len(x_data) == len(y_data)

        layers = [len(x_data[0])] + hidden_layers + [len(set(y_data))]
        self.__exp = theanets.Experiment(theanets.Classifier, layers=layers)
        x_data = numpy.array(x_data, dtype=numpy.float32)
        y_enum = _enumerate(y_data)
        y_data = numpy.array([y[1] for y in y_enum], dtype=numpy.int32)
        self.__y_map = dict(set(y_enum))

        # Split data into training and validation:
        ind = int(split * len(x_data))
        self.__exp.train((x_data[:ind], y_data[:ind]),
                         (x_data[ind:], y_data[ind:]),
                         optimize=self.__optimize,
                         learning_rate=self.__learning_rate,
                         momentum=self.__momentum)

    def classify(self, x_test, y_test):
        '''Classifies and analyses test data.'''
        y_test = numpy.array([self.__y_map[y]
                              for y in y_test], dtype=numpy.int32)

        y_pred = self.__exp.network.classify(x_test)

        inv_y_map = {v: k for k, v in self.__y_map.items()}

        return [inv_y_map[y] for y in y_pred], inv_y_map, \
            classification_report(y_test, y_pred), \
            confusion_matrix(y_test, y_pred)


def randomise_order(x_data, y_data):
    '''Assumes x_data and y_data are paired (such that x_data[i] is paired with
    y_data[i]) and then randomises their orders such that this pairing is
    maintained.'''
    data = zip(x_data, y_data)
    random.shuffle(data)
    return zip(*data)


def _enumerate(lst):
    '''Returns enumeration of supplied list.'''
    label_to_number = defaultdict(partial(next, count()))
    return [(item, label_to_number[item]) for item in lst]
