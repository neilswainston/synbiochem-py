'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=broad-except
# pylint: disable=too-many-arguments
import math
import random
import traceback

from synbiochem.utils.job import JobThread


class SimulatedAnnealer(JobThread):
    '''Class to perform simulated annealing method.'''

    def __init__(self, solution,
                 acceptance=0.01, max_iter=10000,
                 r_temp=0.025, cooling_rate=2.5e-4,
                 heartbeat=1, verbose=False):
        self.__solution = solution
        self.__acceptance = acceptance
        self.__max_iter = max_iter
        self.__r_temp = r_temp
        self.__cooling_rate = cooling_rate
        self.__verbose = verbose
        self.__heartbeat = heartbeat
        JobThread.__init__(self)

    def run(self):
        '''Optimises a solution with simulated annealing.'''
        try:
            if self.__init():
                # Initialization:
                iteration = 0
                accepts = 0
                rejects = 0
                r_temp = self.__r_temp

                energy = self.__solution.get_energy()

                while not self._cancelled \
                        and energy > self.__acceptance \
                        and iteration < self.__max_iter:
                    iteration += 1
                    energy_new = self.__solution.mutate()

                    if energy_new < energy:
                        # Accept move immediately:
                        self.__accept(iteration)
                        self.__log(iteration, r_temp, energy, energy_new, 'T')
                        energy = energy_new
                    elif energy_new == energy:
                        # Reject move:
                        self.__solution.reject()
                        rejects += 1
                    elif _get_exp((energy - energy_new) / r_temp) > \
                            random.random():
                        # Accept move based on conditional probability:
                        self.__accept(iteration)
                        self.__log(iteration, r_temp, energy, energy_new, '*')
                        energy = energy_new
                        accepts += 1
                    else:
                        # Reject move:
                        self.__solution.reject()
                        rejects += 1

                    # Heartbeat:
                    self.__check_heartbeat(iteration)

                    # Simulated annealing control:
                    accepts, rejects, r_temp = \
                        _control_r_temp(accepts, rejects, r_temp)

                    r_temp *= 1 - self.__cooling_rate

                self.__log(iteration, r_temp, energy, energy_new, 'E')
                self.__exit(iteration)
        except Exception:
            self.__fire_event('error', 100, 0, message=traceback.format_exc())

    def __init(self):
        '''Initialise.'''
        self.__fire_event('running', 0, 0, message='Job initialising...')

        try:
            self.__solution.init()
        except ValueError:
            self.__fire_event('error', 100, 0, message=traceback.format_exc())
            return False

        self.__fire_event('running', 0, 0, message='Job initialised')
        return True

    def __log(self, iteration, r_temp, energy, energy_new, symbol):
        '''Log data.'''
        if self.__verbose:
            print('\t'.join([symbol,
                             str(iteration),
                             str(energy),
                             str(energy_new),
                             '%.2f' % r_temp,
                             '%.2f' % _get_exp((energy - energy_new) / r_temp),
                             str(self.__solution)]))

    def __check_heartbeat(self, iteration):
        '''Heartbeat.'''
        if float(iteration) % self.__heartbeat == 0:
            self.__fire_event('running',
                              float(iteration) / self.__max_iter * 100,
                              iteration,
                              'Running...')

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

        self._fire_event(event)

    def __exit(self, iteration):
        '''Exit.'''
        if iteration == self.__max_iter:
            message = 'Unable to optimise in ' + str(self.__max_iter) + \
                ' iterations'
            self.__fire_event('error', 100, iteration, message=message)
        elif self._cancelled:
            self.__fire_event('cancelled', 100, iteration,
                              message='Job cancelled')
        else:
            self.__fire_event('finished', 100, iteration,
                              message='Job completed')


def _control_r_temp(accepts, rejects, r_temp):
    '''Control temperature.'''
    if accepts + rejects > 50:
        if float(accepts) / float(accepts + rejects) > 0.2:
            # Too many accepts, reduce r_temp:
            r_temp /= 2.0
            accepts = 0
            rejects = 0
        elif float(accepts) / float(accepts + rejects) < 0.025:
            # Too many rejects, increase r_temp:
            r_temp *= 2.0
            accepts = 0
            rejects = 0

    return accepts, rejects, r_temp


def _get_exp(value):
    '''Get exp safely.'''
    try:
        ans = math.exp(value)
    except OverflowError:
        ans = float('inf')

    return ans
