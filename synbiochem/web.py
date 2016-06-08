'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import json
import time
import uuid

from flask import Response


class FlaskManager(object):
    '''Class to run Flask application.'''

    def __init__(self):
        self.__status = {}
        self.__threads = {}

    def submit(self, req):
        '''Responds to submission.'''
        query = json.loads(req.data)
        self.process_query(query)

        # Do job in new thread, return result when completed:
        job_id = str(uuid.uuid4())
        self.__status[job_id] = {'job_id': job_id, 'progress': 0}
        thread = self.get_thread(job_id, query)
        thread.add_listener(self)
        self.__threads[job_id] = thread
        thread.start()

        return json.dumps({'job_id': job_id})

    def get_progress(self, job_id):
        '''Returns progress of job.'''
        def _check_progress(job_id):
            '''Checks job progress.'''
            while self.__status[job_id]['progress'] < 100:
                time.sleep(1)
                yield 'data:' + self.__get_response(job_id) + '\n\n'

            yield 'data:' + self.__get_response(job_id) + '\n\n'

        if job_id in self.__status:
            return Response(_check_progress(job_id),
                            mimetype='text/event-stream')
        else:
            return 'Job ' + job_id + ' unknown or finished.'

    def cancel(self, job_id):
        '''Cancels job.'''
        self.__threads[job_id].cancel()
        return 'Cancelled: ' + job_id

    def event_fired(self, event):
        '''Responds to event being fired.'''
        self.__status[event['job_id']] = event

    def process_query(self, query):
        '''Perform application-specific pre-processing of query.'''
        # Take no action: needs to be overridden.

    def get_thread(self, job_id, query):
        '''Gets a thread to run query.'''
        # Needs to be overridden.
        return None

    def __get_response(self, job_id):
        '''Returns current progress for job id.'''
        return json.dumps(self.__status[job_id])
