'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import math
import random
from synbiochem.utils.job import EventHandler


class SimulatedAnnealer(object):
    '''Class to perform simulated annealing method.'''

    def __init__(self, verbose=False):
        self.__verbose = verbose
        self.__ev_hand = EventHandler()

    def add_listener(self, listener):
        '''Adds an event listener.'''
        self.__ev_hand.add_listener(listener)

    def remove_listener(self, listener):
        '''Removes an event listener.'''
        self.__ev_hand.remove_listener(listener)

    def optimise(self, solution, acceptance=0.25, max_iter=10000):
        '''Optmises a solution with simulated annealing.'''

        # Initialization:
        counter = 0
        accepts = 0
        rejects = 0
        r_temp = 0.6
        cooling_rate = 1e-3

        energy = solution.get_energy()

        if self.__verbose:
            print str(counter) + '\t' + str(energy) + '\t' + \
                str(solution)

        while energy > acceptance and counter < max_iter:
            counter += 1
            energy_new = solution.mutate()

            if energy_new < energy:
                # Accept move immediately:
                energy = energy_new
                self.__accept(solution, counter, max_iter, energy)
            elif energy == energy_new:
                # Take no action:
                continue
            elif math.exp((energy - energy_new) / r_temp) > random.random():
                # Accept move based on conditional probability:
                energy = energy_new
                self.__accept(solution, counter, max_iter, energy)
                accepts += 1
            else:
                # Reject move:
                rejects += 1

            # Simulated annealing control:
            if accepts + rejects > 50:
                if float(accepts) / float(accepts + rejects) > 0.2:
                    # Too many accepts, reduce r_temp:
                    r_temp /= 2.0
                    accepts = 0
                    rejects = 0
                elif float(accepts) / float(accepts + rejects) < 0.01:
                    # Too many rejects, increase r_temp:
                    r_temp *= 2.0
                    accepts = 0
                    rejects = 0

            r_temp *= 1 - cooling_rate

        if counter == max_iter:
            raise ValueError('Unable to optimise: ' + str(solution))

        return (solution, counter)

    def __accept(self, solution, counter, max_iter, energy):
        '''Accept the current solution.'''
        solution.accept()

        event = dict({'progress': float(counter) / max_iter}.items() +
                     solution.get_json().items())
        self.__ev_hand.fire_event(event)

        if self.__verbose:
            print str(counter) + '\t' + str(energy) + '\t' + \
                str(solution)
