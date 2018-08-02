'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=broad-except
# pylint: disable=superfluous-parens
# pylint: disable=useless-object-inheritance
from threading import Thread


try:
    import Queue as queue
except ImportError:
    pass


class Worker(Thread):
    '''Thread executing tasks from a given tasks queue.'''

    def __init__(self, tasks):
        Thread.__init__(self)
        self.__tasks = tasks
        self.daemon = True
        self.start()

    def run(self):
        while True:
            func, args, kargs = self.__tasks.get()
            try:
                func(*args, **kargs)
            except Exception as err:
                print(err)
            finally:
                # Mark this task as done, whether an exception happened or not
                self.__tasks.task_done()


class ThreadPool(object):
    '''Pool of threads consuming tasks from a queue.'''

    def __init__(self, num_threads):
        self.__tasks = queue.Queue(num_threads)

        for _ in range(num_threads):
            Worker(self.__tasks)

    def add_task(self, func, *args, **kargs):
        '''Add a task to the queue.'''
        self.__tasks.put((func, args, kargs))

    def map(self, func, args_list):
        '''Add a list of tasks to the queue.'''
        for args in args_list:
            self.add_task(func, args)

    def wait_completion(self):
        '''Wait for completion of all the tasks in the queue.'''
        self.__tasks.join()
