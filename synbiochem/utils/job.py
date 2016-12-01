'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from threading import Thread
import uuid


class JobThread(Thread):
    '''Wraps a job into a thread, and fires events.'''

    def __init__(self):
        Thread.__init__(self)
        self._cancelled = False
        self.__job_id = str(uuid.uuid4())
        self.__listeners = set()

    def get_job_id(self):
        '''Gets thread job id.'''
        return self.__job_id

    def cancel(self):
        '''Cancels the current job.'''
        self._cancelled = True

    def add_listener(self, listener):
        '''Adds an event listener.'''
        self.__listeners.add(listener)

    def remove_listener(self, listener):
        '''Removes an event listener.'''
        self.__listeners.remove(listener)

    def _fire_event(self, event):
        '''Event listener, passes events on to registered listeners.'''
        event.update({'job_id': self.__job_id})

        for listener in self.__listeners:
            listener.event_fired(event)
