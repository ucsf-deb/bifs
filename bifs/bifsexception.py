"""
This module defines Exception classes to indicate errors or warnings and functions to
install global handlers for those exceptions.  Global handlers are important since
exceptions may be raised from code executing inside the Qt event loop.

If used from a Qt program pops up a message box with exception info.  
Perhaps too much info (full traceback)

Exceptions are also logged to sys.stdout.

Exceptions
----------
All classes are subclasses of BifsException.  But note the code may raise other kinds of
exceptions.   BifsFatal is a base class for exceptions that should terminate the program, and they
always include an error code with which to terminate.  Other exceptions will typically function
as notifications.

Logging
-------

Goes to a logger named "Bifs".

Useage
------

In your main program
====================
import bifs.bifsexception
# will activate the error handlers
# if you need to use the exceptions or the log, use the patterns below instead.

Other modules
=============
from bifs.bifsexception import *
# and then raise any of the exceptions
# if you need to log something
import logging
log = logging.getLogger("Bifs")

Implementation Note
-------------------
Python's import machinery should assure the code here is executed only once even if multiple
modules import this module.  For that to work you must use `import` rather than other strategies to access
the code.

Future Directions
-----------------
* ability to suppress traceback in dialog box
* separate the general and GUI specific parts
* more use of logging

Credits
-------
Based on code from Tim Lehr https://timlehr.com/python-exception-hooks-with-qt-message-box/
He indicated no explicit license but did write, just below it
"Just paste it somewhere in one of your modules"

Modifications by Ross Boylan <ross.boylan@ucsf.edu> under GPL-3 or BSD-ish license for the whole package.

"""

__all__ = ("BifsException", "BifsBadParameter", "BifsBadIndex", "BifsLookupFailure", "BifsBadInputs")


import sys
import traceback
import logging
from PyQt5 import QtCore, QtWidgets

# basic logger functionality
log = logging.getLogger("Bifs")
handler = logging.StreamHandler(stream=sys.stdout)
log.addHandler(handler)

class BifsException(Exception):
    pass

class BifsBadParameter(BifsException):
    pass

class BifsBadIndex(BifsException):
    pass

class BifsLookupFailure(BifsException):
    pass

class BifsBadInputs(BifsException):
    pass

def show_exception_box(log_msg):
    """Checks if a QApplication instance is available and shows a messagebox with the exception message. 
    If unavailable (non-console application), log an additional notice.
    """
    if QtWidgets.QApplication.instance() is not None:
            errorbox = QtWidgets.QMessageBox()
            errorbox.setText("Oops. An unexpected error occured:\n{0}".format(log_msg))
            errorbox.exec_()
    else:
        log.debug("No QApplication instance available.")
 
class UncaughtHook(QtCore.QObject):
    _exception_caught = QtCore.Signal(object)
 
    def __init__(self, *args, **kwargs):
        super(UncaughtHook, self).__init__(*args, **kwargs)

        # this registers the exception_hook() function as hook with the Python interpreter
        sys.excepthook = self.exception_hook

        # connect signal to execute the message box function always on main thread
        self._exception_caught.connect(show_exception_box)
 
    def exception_hook(self, exc_type, exc_value, exc_traceback):
        """Function handling uncaught exceptions.
        It is triggered each time an uncaught exception occurs. 
        """
        if issubclass(exc_type, KeyboardInterrupt):
            # ignore keyboard interrupt to support console applications
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
        else:
            exc_info = (exc_type, exc_value, exc_traceback)
            log_msg = '\n'.join([''.join(traceback.format_tb(exc_traceback)),
                                 '{0}: {1}'.format(exc_type.__name__, exc_value)])
            if issubclass(exc_type, BifsException):
                log_msg = str(exc_value)+"\nProceed with caution.\n"+log_msg
                log.warn(log_msg)
            else:
                log.critical("Uncaught exception:\n {0}".format(log_msg), exc_info=exc_info)

            # trigger message box show
            self._exception_caught.emit(log_msg)
 
# create a global instance of our class to register the hook
qt_exception_hook = UncaughtHook()
