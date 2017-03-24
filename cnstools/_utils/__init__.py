from .callcommands import *
from .headerprint import *
from .multiprinter import *
from .multitracker import *
from .safeprint import *

def make_dir(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path): raise