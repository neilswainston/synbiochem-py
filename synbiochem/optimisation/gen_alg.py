'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import abc
import numpy
import random


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
        for i in range(self.__len_chromosome):
            if random.random() < self.__mutation_rate:
                # Mutate:
                mask = 1 << i
                self.__chromosome = self.__chromosome ^ mask

    def breed(self, partner):
        i = int(random.random() * self.__len_chromosome)
        end = 2 ** i - 1
        start = self.__mask - end
        return (partner.__chromosome & start) + (self.__chromosome & end)

    def __repr__(self):
        return format(self.__chromosome, '0' + str(self.__len_chromosome) +
                      'b') + '\t' + str(self.__chromosome)


class GeneticAlgorithm(object):
    '''Base class to run a genetic algorithm.'''
    __metaclass__ = abc.ABCMeta

    def __init__(self, pop_size, args):
        self.__pop_size = pop_size
        self.__args = args
        self.__population = [self.__get_individual()
                             for _ in xrange(pop_size)]

    def run(self, max_iter=1024):
        for _ in range(max_iter):
            result = self.__evolve()

            if result is not None:
                return result

        raise ValueError('Unable to optimise in ' + str(max_iter) +
                         ' iterations.')

    @abc.abstractmethod
    def _fitness(self, individual):
        '''Determine the fitness of an individual.'''
        return

    def __get_individual(self):
        'Create a member of the population.'
        return {k: self.__get_arg(args) for k, args in self.__args.iteritems()}

    def __get_arg(self, args):
        return random.randint(args[0], args[1]) if isinstance(args, tuple) \
            else random.choice(args)

    def __evolve(self, retain=0.2, random_select=0.05, mutate=0.01):
        # Ensure uniqueness in population:
        self.__population = list(numpy.unique(numpy.array(self.__population)))

        graded = sorted([(self._fitness(x), x) for x in self.__population])

        if graded[0][0] == 0:
            return graded[0][1]

        graded = [x[1] for x in graded]
        retain_length = int(self.__pop_size * retain)

        # Retain best and randomly add other individuals to promote genetic
        # diversity:
        self.__population = graded[:retain_length] + \
            [ind for ind in graded[retain_length:]
             if random_select > random.random()]

        # Mutate some individuals:
        for individual in self.__population:
            if mutate > random.random():
                pos = random.randint(0, len(individual) - 1)
                individual[pos] = self.__get_arg(self.__args[pos])

        # Breed parents to create children:
        children = []

        while len(children) < self.__pop_size - len(self.__population):
            male = random.choice(self.__population)
            female = random.choice(self.__population)

            if male != female:
                pos = random.randint(0, len(male))

                child = dict({k: v for k, v in male.iteritems()
                              if k < pos}.items() +
                             {k: v for k, v in female.iteritems()
                              if k >= pos}.items())

                children.append(child)

        self.__population.extend(children)

        # print numpy.mean([self._fitness(ind) for ind in self.__population])

        return None
