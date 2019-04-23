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


class Chromosome():
    '''Class to represent a chromosome.'''

    def __init__(self, chromosome, mutation_rate=0.1):
        self._chromosome = chromosome
        self.__mutation_rate = mutation_rate

    def get_chromosome(self):
        '''Gets the chromosome.'''
        return self._chromosome

    def mutate(self):
        '''Mutate.'''
        return None

    def breed(self, partner):
        '''Breed.'''
        pass

    def fitness(self):
        '''Get fitness.'''
        pass

    def __repr__(self):
        return str(self._chromosome)


class GeneticAlgorithm():
    '''Base class to run a genetic algorithm.'''

    def __init__(self, population, retain=0.05, random_select=0.2,
                 mutate=0.25, verbose=False):
        self.__population = population
        self.__retain = retain
        self.__random_select = random_select
        self.__mutate = mutate
        self._verbose = verbose

    def run(self, max_iter=1024):
        '''Run.'''
        for _ in range(max_iter):
            results = self.__evolve()

            if results is not None:
                return results

        raise ValueError('Unable to optimise in %d iterations.' % max_iter)

    def get_results(self):
        '''Get results.'''
        return sorted([(chrm.fitness(), chrm) for chrm in self.__population])

    def __evolve(self):
        results = self.get_results()

        if results[0][0] == 0:
            return results

        if self._verbose:
            print(results[0])

        results = [x[1] for x in results]
        retain_length = int(len(self.__population) * self.__retain)

        # Retain best and randomly add other individuals to promote genetic
        # diversity:
        new_population = results[:retain_length]

        # Mutate some individuals:
        for individual in new_population:
            if self.__mutate > random.random():
                mutated = individual.mutate()

                if mutated:
                    new_population.append(mutated)

        new_population.extend([ind for ind in results[retain_length:]
                               if self.__random_select > random.random()])

        while len(new_population) < len(self.__population):
            male = random.choice(self.__population)
            female = random.choice(self.__population)

            if male != female:
                child = male.breed(female)

                if child:
                    new_population.append(child)

                new_population = list(
                    numpy.unique(numpy.array(new_population)))

        self.__population = new_population

        # print numpy.mean([self._fitness(ind) for ind in self.__population])

        return None
