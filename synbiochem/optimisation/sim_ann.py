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
        r_temp = 0.0075
        cooling_rate = 1e-3

        energy = self.__solution.get_energy()

        while not self.__cancelled \
                and energy > self.__acceptance \
                and iteration < self.__max_iter:
            iteration += 1
            energy_new = self.__solution.mutate()

            if energy_new < energy:
                # Accept move immediately:
                energy = energy_new
                self.__accept(iteration)
                print '\t'.join(['T', str(iteration), str(energy_new),
                                 str(self.__solution)])
            elif energy_new == energy:
                # Reject move:
                self.__solution.reject()
                rejects += 1
                print '\t'.join([' ', str(iteration), str(energy_new),
                                 str(self.__solution)])
            elif math.exp((energy - energy_new) / r_temp) > random.random():
                # Accept move based on conditional probability:
                energy = energy_new
                self.__accept(iteration)
                accepts += 1
                print '\t'.join(['*', str(iteration), str(energy_new),
                                 str(self.__solution)])
            else:
                # Reject move:
                self.__solution.reject()
                rejects += 1
                print '\t'.join([' ', str(iteration), str(energy_new),
                                 str(self.__solution)])

            # Heartbeat:
            if float(iteration) % 10 == 0:
                self.__fire_event('running',
                                  float(iteration) / self.__max_iter * 100,
                                  iteration,
                                  'Running...')

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
                              message='Job completed')

    def __accept(self, iteration):
        '''Accept the current solution.'''
        self.__solution.accept()

        self.__fire_event('running', float(iteration) / self.__max_iter * 100,
                          iteration, 'Running...')

    def __fire_event(self, status, progress, iteration, message=''):
        '''Fires an event.'''
        event = {'update': {'status': status,
                            'message': message,
                            'progress': progress,
                            'iteration': iteration,
                            'max_iter': self.__max_iter,
                            'values': self.__solution.get_values()},
                 'query': self.__solution.get_query()
                 }

        if status == 'finished':
            event['result'] = self.__solution.get_result()

        self.__ev_hand.fire_event(event)
