# load inputs
import src.main.java.sc_sizing.sizing as sz
import jpype
import jpype.imports
from jpype.types import *

file_name = 'test_input.json'
instrument_lists = sz.get_instrument_lists(file_name)
orbit_list = sz.get_orbit_lists(file_name)

# print inputs
print(instrument_lists)
print(orbit_list)

jpype.startJVM()
jpype.addClassPath("../../")

x = 1