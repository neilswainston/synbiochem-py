'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from threading import Thread


class JobThread(Thread):
    '''Wraps a job into a thread, and fires events.'''

    def __init__(self, job_id):
        Thread.__init__(self)
        self.__job_id = job_id
        self.__cancelled = False
        self._event_handler = EventHandler()

    def cancel(self):
        '''Cancels the current job.'''
        self.__cancelled = True

    def add_listener(self, listener):
        '''Adds an event listener.'''
        self._event_handler.add_listener(listener)

    def remove_listener(self, listener):
        '''Removes an event listener.'''
        self._event_handler.remove_listener(listener)

    def event_fired(self, event):
        '''Event listener, passes events on to registered listeners.'''
        event = dict({'job_id': self.__job_id}.items() + event.items())
        for listener in self._event_handler.get_listeners():
            listener.event_fired(event)


class EventHandler(object):
    '''Simple event handler class.'''

    def __init__(self):
        self.__listeners = set()

    def add_listener(self, listener):
        '''Adds an event listener.'''
        self.__listeners.add(listener)

    def remove_listener(self, listener):
        '''Removes an event listener.'''
        self.__listeners.remove(listener)

    def get_listeners(self):
        '''Returns a copied list of the event listeners.'''
        return list(self.__listeners)

    def fire_event(self, event):
        '''Fires an event to event listeners.'''
        for listener in self.__listeners:
            listener.event_fired(event)
