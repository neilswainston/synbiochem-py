'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-public-methods
# import climate
import unittest

from sklearn import datasets, linear_model
from sklearn.datasets.samples_generator import make_blobs

import matplotlib.pyplot as plt
import numpy as np
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


class TestRegressor(unittest.TestCase):
    '''Tests the Regressor class.'''

    def test_regressor(self):
        # Load the diabetes dataset
        diabetes = datasets.load_diabetes()

        # Use only one feature
        diabetes_X = diabetes.data[:, np.newaxis, 2]

        # Split the data into training/testing sets
        diabetes_X_train = diabetes_X[:-20]
        diabetes_X_test = diabetes_X[-20:]

        # Split the targets into training/testing sets
        diabetes_y_train = diabetes.target[:-20]
        diabetes_y_test = diabetes.target[-20:]

        # Create linear regression object
        regr = linear_model.LinearRegression()

        # Train the model using the training sets
        regr.fit(diabetes_X_train, diabetes_y_train)

        # The coefficients
        print('Coefficients: \n', regr.coef_)
        # The mean square error
        print("Resid sum of squares: %.2f"
              % np.mean((regr.predict(diabetes_X_test) - diabetes_y_test) **
                        2))
        # Explained variance score: 1 is perfect prediction
        print('Variance score: %.2f' %
              regr.score(diabetes_X_test, diabetes_y_test))

        regressor = synbiochem.ann.Regressor()
        y_train = [[y] for y in diabetes_y_train]
        regressor.train(diabetes_X_train, y_train)
        y_pred = regressor.predict(diabetes_X_test)

        # Plot outputs
        plt.scatter(diabetes_X_test, diabetes_y_test,  color='black')
        plt.scatter(
            diabetes_X_test, regr.predict(diabetes_X_test), color='blue')
        plt.scatter(diabetes_X_test, y_pred, color='red')

        plt.xticks(())
        plt.yticks(())

        plt.show()
