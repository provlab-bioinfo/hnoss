import sys # added!
sys.path.append("..")
from analysis.functions import *

def test_sigfig():
    assert sigfig(3,0) == "3"
    assert sigfig(3,1) == "3.0"
    assert sigfig(3,2) == "3.00"
    assert sigfig(3,3) == "3.000"
    assert sigfig(3.1,3) == "3.100"
    assert sigfig(3.14,3) == "3.140"
    assert sigfig(3.141,3) == "3.141"
    assert sigfig(3.1415,3) == "3.142"
    assert sigfig(3.14159,3) == "3.142"
    assert sigfig(3.14159,10) == "3.1415900000"
