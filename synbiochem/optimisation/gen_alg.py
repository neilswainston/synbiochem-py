'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
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

    def breed(self, other):
        i = int(random.random() * self.__len_chromosome)
        temp = self.__chromosome
        end = 2 ** i - 1
        start = self.__mask - end
        self.__chromosome = (other.__chromosome & start) + (temp & end)
        other.__chromosome = (temp & start) + (other.__chromosome & end)

    def __repr__(self):
        return format(self.__chromosome, '0' + str(self.__len_chromosome) +
                      'b') + '\t' + str(self.__chromosome)
