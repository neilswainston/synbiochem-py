'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-arguments
import math
import random

from synbiochem.utils.job import EventHandler


class SimulatedAnnealer(object):
    '''Class to perform simulated annealing method.'''

    def __init__(self, solution, acceptance=0.1, max_iter=10000,
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
        r_temp = 0.0025
        cooling_rate = 1e-6

        energy = self.__solution.get_energy()

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
                print '\t'.join([str(energy - energy_new),
                                 str(math.exp((energy - energy_new) / r_temp)),
                                 str(r_temp)])
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
            message = 'Unable to optimise in ' + str(self.__max_iter) + \
                ' iterations'
            self.__fire_event('error', 100, iteration, message=message)
        elif self.__cancelled:
            self.__fire_event('cancelled', 100, iteration,
                              message='Job cancelled')
        else:
            self.__fire_event('finished', 100, iteration,
                              message='Job completed', result=True)

    def __accept(self, iteration, energy):
        '''Accept the current solution.'''
        self.__solution.accept()
        self.__fire_event('running', float(iteration) / self.__max_iter * 100,
                          iteration)
        if self.__verbose:
            print str(iteration) + '\t' + str(energy) + '\t' + \
                str(self.__solution)

    def __fire_event(self, status, progress, iteration, message='',
                     result=False):
        '''Fires an event.'''
        event = dict({'status': status,
                      'message': message,
                      'progress': progress,
                      'iteration': iteration,
                      'max_iter': self.__max_iter}.items() +
                     self.__solution.get_query().items() +
                     self.__solution.get_update().items() +
                     (self.__solution.get_result() if result else {}).items())

        self.__ev_hand.fire_event(event)
