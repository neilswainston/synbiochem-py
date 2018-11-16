'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-self-use
# pylint: disable=protected-access
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
import random

import numpy


class Chromosome(object):
    '''Class to represent a chromosome.'''

    def __init__(self, chromosome, mutation_rate=0.1):
        self._chromosome = chromosome
        self.__mutation_rate = mutation_rate

    def get_chromosome(self):
        '''Gets the chromosome.'''
        return self._chromosome

    def mutate(self):
        '''Mutate.'''
        pass

    def breed(self, partner):
        '''Breed.'''
        pass

    def fitness(self):
        '''Get fitness.'''
        pass

    def __repr__(self):
        return str(self._chromosome)


class GeneticAlgorithm(object):
    '''Base class to run a genetic algorithm.'''

    def __init__(self, population, retain=0.2, random_select=0.05,
                 mutate=0.1, verbose=False):
        self.__population = population
        self.__retain = retain
        self.__random_select = random_select
        self.__mutate = mutate
        self._verbose = verbose

    def run(self, max_iter=1024):
        '''Run.'''
        for _ in range(max_iter):
            result = self.__evolve()

            if result is not None:
                return result

        raise ValueError('Unable to optimise in %d iterations.' % max_iter)

    def __evolve(self):
        graded = sorted([(chrm.fitness(), chrm)
                         for chrm in self.__population])

        if graded[0][0] == 0:
            return graded[0][1]

        if self._verbose:
            print graded[0]

        graded = [x[1] for x in graded]
        retain_length = int(len(self.__population) * self.__retain)

        # Retain best and randomly add other individuals to promote genetic
        # diversity:
        new_population = graded[:retain_length] + \
            [ind for ind in graded[retain_length:]
             if self.__random_select > random.random()]

        while len(new_population) < len(self.__population):
            male = random.choice(self.__population)
            female = random.choice(self.__population)

            if male != female:
                child = male.breed(female)

                new_population.append(child)

                new_population = list(
                    numpy.unique(numpy.array(new_population)))

        # Mutate some individuals:
        for individual in new_population:
            if self.__mutate > random.random():
                individual.mutate()

        self.__population = new_population

        # print numpy.mean([self._fitness(ind) for ind in self.__population])

        return None
