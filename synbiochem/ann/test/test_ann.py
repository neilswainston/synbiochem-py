'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-public-methods
# import climate
import unittest

from sklearn.datasets.samples_generator import make_blobs

import synbiochem.ann


class TestClassifier(unittest.TestCase):
    '''Tests the Classifier class.'''

    def test_classify(self):
        '''Tests the classify method.'''
        # climate.enable_default_logging()

        # Generate some data, convert to list of floats (inputs) and string
        # 'names' for output classifications:
        x_data, y_data = make_blobs(n_samples=1000, centers=5, n_features=3,
                                    cluster_std=1.0, random_state=0)
        x_data = x_data.tolist()
        y_data = [str(unichr(y + ord('A'))) for y in y_data]

        # Split data into training and classifying:
        ind = int(0.8 * len(x_data))

        classifier = synbiochem.ann.Classifier()
        classifier.train(x_data[:ind], y_data[:ind])

        y_test = y_data[ind:]
        y_pred, _, _, _ = classifier.classify(x_data[ind:], y_test)
        self.assertTrue(sum([i == j for i, j in zip(y_test, y_pred)]) /
                        float(len(y_test)) > 0.9)
