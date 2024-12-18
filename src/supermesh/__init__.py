import ctypes
import os
from ctypes.util import find_library
from pathlib import Path

_libsupermesh_dir = Path(__file__).parent

def get_library():
    """Retrieve the shared library as a string
    """
    return str(next(_libsupermesh_dir.joinpath("lib").glob("*supermesh*")))

def get_include():
    """Retrieve the path to the headers as a string
    """
    return str(_libsupermesh_dir.joinpath("include"))

def load():
    """Load the shared library and return a handle
    """
    path = get_library()
    return ctypes.cdll.LoadLibrary(path)

# libsupermesh = load()
# access funcions with
# libsupermesh.libsupermesh_print_backtrace
