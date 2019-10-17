'''
Created on Nov 22, 2018

@author: mluksza
'''

import time

class Log(object):
    '''
    classdocs
    '''

    VERBOSE=3
    TSTART = time.time()

    @staticmethod
    def log(msg, level, warn=False):
        """
        Print log message msg to stdout.

        Parameters
        -----------

         msg : str
            String to print on the screen

         level : int
            Log-level. Only the messages with a level higher than the
            current verbose level will be shown.

         warn : bool
            Warning flag. If True, the message will be displayed
            regardless of its log-level.

        """
        if level<Log.VERBOSE or (warn and level<=Log.VERBOSE):
            msg=str(msg)
            dt = time.time() - Log.TSTART
            outstr = '\n' if level<2 else ''
            outstr += format(dt, '4.2f')+'\t'
            outstr += level*'-'
            outstr += msg
            print outstr#, file=sys.stdout)

