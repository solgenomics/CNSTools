import time
from safeprint import safe_print

#a variety fancy headers!
def header_print(header,h_type=0):
    wall = lambda char: char*70
    if h_type==0:
        safe_print("\n%s\n%s \n%s" % (wall(" "),header,wall("~")))
    if h_type==1:
        safe_print("\n%s\n%s \n%s" % (wall("~"),header,wall("~")))
    if h_type==2:
        safe_print("\n%s\n%s \n%s" % (wall("="),header,wall("~")))
    if h_type==3:
        safe_print("\n%s\n%s \n%s" % (wall("="),header,wall("=")))