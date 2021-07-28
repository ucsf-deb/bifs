# Added this file on the theory the import machinery requires it
# to access stuff further down
print("__name__ is ", __name__)
from .bifs import *
from . import pset_dialogs
