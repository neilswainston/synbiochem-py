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

    def __init__(self, solution, acceptance=0.25, max_iter=10000,
                 verbose=False):
        self.__solution = solution
        self.__acceptance = acceptance
        self.__max_iter = max_iter
        self.__verbose = verbose
        self.__cancelled = False
        self.__ev_hand = EventHandler()

    def cancel(self):
        '''Cancel job.'''
        self.__cancelled = True

    def add_listener(self, listener):
        '''Adds an event listener.'''
        self.__ev_hand.add_listener(listener)

    def remove_listener(self, listener):
        '''Removes an event listener.'''
        self.__ev_hand.remove_listener(listener)

    def optimise(self):
        '''Optmises a solution with simulated annealing.'''

        # Initialization:
        iteration = 0
        accepts = 0
        rejects = 0
        r_temp = 0.6
        cooling_rate = 1e-3

        energy = self.__solution.get_energy()

        if self.__verbose:
            print str(iteration) + '\t' + str(energy) + '\t' + \
                str(self.__solution)

        while not self.__cancelled \
                and energy > self.__acceptance \
                and iteration < self.__max_iter:
            iteration += 1
            energy_new = self.__solution.mutate()

            if energy_new < energy:
                # Accept move immediately:
                energy = energy_new
                self.__accept(iteration, energy)
            elif energy == energy_new:
                # Take no action:
                continue
            elif math.exp((energy - energy_new) / r_temp) > random.random():
                # Accept move based on conditional probability:
                energy = energy_new
                self.__accept(iteration, energy)
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

        if iteration == self.__max_iter:
            raise ValueError('Unable to optimise in ' + self.__max_iter +
                             ' iterations')

        self.__fire_event(iteration, 100, result=True)

    def __accept(self, iteration, energy):
        '''Accept the current solution.'''
        self.__solution.accept()
        self.__fire_event(iteration, float(iteration) / self.__max_iter)
        if self.__verbose:
            print str(iteration) + '\t' + str(energy) + '\t' + \
                str(self.__solution)

    def __fire_event(self, iteration, progress, result=False):
        '''Fires an event.'''
        event = dict({'progress': progress,
                      'iteration': iteration,
                      'max_iter': self.__max_iter}.items() +
                     self.__solution.get_status().items() +
                     (self.__solution.get_result() if result else {}).items())

        self.__ev_hand.fire_event(event)
