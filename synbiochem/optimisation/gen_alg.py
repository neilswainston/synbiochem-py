'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
import abc
import random

import numpy


class Chromosome(object):
    '''Class to represent a chromosome.'''

    def __init__(self, chromosome, len_chromosome, mutation_rate=0.01):
        self.__chromosome = chromosome
        self.__len_chromosome = len_chromosome
        self.__mask = 2 ** (self.__len_chromosome) - 1
        self.__mutation_rate = mutation_rate

    def get_chromosome(self):
        '''Gets the chromosome.'''
        return self.__chromosome

    def mutate(self):
        '''Mutates the chromosome.'''
        for i in range(self.__len_chromosome):
            if random.random() < self.__mutation_rate:
                # Mutate:
                mask = 1 << i
                self.__chromosome = self.__chromosome ^ mask

    def breed(self, partner):
        '''Breeds chromosome with a partner Chromosome.'''
        i = int(random.random() * self.__len_chromosome)
        end = 2 ** i - 1
        start = self.__mask - end
        return (partner.get_chromosome() & start) + (self.__chromosome & end)

    def __repr__(self):
        return format(self.__chromosome, '0' + str(self.__len_chromosome) +
                      'b') + '\t' + str(self.__chromosome)


class GeneticAlgorithm(object):
    '''Base class to run a genetic algorithm.'''
    __metaclass__ = abc.ABCMeta

    def __init__(self, pop_size, args, retain=0.2, random_select=0.05,
                 mutate=0.01, verbose=False):
        self.__pop_size = pop_size
        self.__args = args
        self.__retain = retain
        self.__random_select = random_select
        self.__mutate = mutate
        self._verbose = verbose
        self.__population = []

        while len(list(numpy.unique(numpy.array(self.__population)))) < \
                pop_size:
            self.__population.append(self.__get_individual())

    def run(self, max_iter=1024, max_tries=1024):
        '''Runs the genetic algorithm.'''
        for _ in range(max_iter):
            result = self.__evolve(max_tries)

            if result is not None:
                return result

        raise ValueError('Unable to optimise in ' + str(max_iter) +
                         ' iterations.')

    @abc.abstractmethod
    def _fitness(self, individual):
        '''Determine the fitness of an individual.'''
        return

    def __get_individual(self):
        '''Create a member of the population.'''
        return {k: _get_arg(args) for k, args in self.__args.iteritems()}

    def __evolve(self, max_tries):
        '''Performs one round of evolution.'''
        graded = sorted([(self._fitness(x), x) for x in self.__population])

        if graded[0][0] == 0:
            return graded[0][1]

        if self._verbose:
            print graded[0]

        graded = [x[1] for x in graded]
        retain_length = int(self.__pop_size * self.__retain)

        # Retain best and randomly add other individuals to promote genetic
        # diversity:
        self.__population = graded[:retain_length] + \
            [ind for ind in graded[retain_length:]
             if self.__random_select > random.random()]

        # Mutate some individuals:
        for individual in self.__population:
            if self.__mutate > random.random():
                key = random.choice(individual.keys())

                if key in self.__args:
                    individual[key] = _get_arg(self.__args[key])

        # Ensure uniqueness in population:
        self.__population = list(numpy.unique(numpy.array(self.__population)))

        self.__breed(max_tries)

    def __breed(self, max_tries):
        '''Breeds parents to create children.'''
        new_population = []
        tries = 0

        while len(new_population) < self.__pop_size:
            male = random.choice(self.__population)
            female = random.choice(self.__population)

            if male != female:
                pos = random.randint(0, len(male))

                male_parts = {k: male[k]
                              for i, k in enumerate(male.keys())
                              if i < pos}
                female_parts = {k: female[k]
                                for i, k in enumerate(female.keys())
                                if i >= pos}

                child = dict(male_parts.items() +
                             female_parts.items())

                new_population.append(child)
                new_population = \
                    list(numpy.unique(numpy.array(new_population)))

            tries += 1

            if tries == max_tries:
                raise ValueError('Unable to generate unique population.')

        self.__population = new_population


def _get_arg(args):
    '''Gets a random argument.'''
    return random.randint(args[0], args[1]) if isinstance(args, tuple) \
        else random.choice(args)
