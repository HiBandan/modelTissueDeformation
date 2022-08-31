import os
import sys
import glob
import numpy
import shutil
import warnings
warnings.filterwarnings("ignore")

# suppress-warning
sys.dont_write_bytecode=True
warnings.filterwarnings("ignore")

def recreateDirectory(path):        
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)
    
def decomment(infile):
    for row in infile:
        raw = row.split('#')[0].strip()
        if raw: yield raw
        
def listdir_nohidden(path):
    return glob.glob(os.path.join(path,'*'))

# intersection-point-of-two-line-segments
def intersectionLineSegments(v1,v2,u1,u2):
    encounter = False
    u = [a_i - b_i for a_i, b_i in zip(u2,u1)] 
    v = [a_i - b_i for a_i, b_i in zip(v2,v1)] 
    w = [a_i - b_i for a_i, b_i in zip(u1,v1)] 
    sval1 = (v[1]*w[0]-v[0]*w[1])/(v[0]*u[1]-v[1]*u[0]) 
    w = [a_i - b_i for a_i, b_i in zip(v1,u1)]
    sval2 = (u[1]*w[0]-u[0]*w[1])/(u[0]*v[1]-u[1]*v[0])
    interSecPoint = None
    if(((sval1 >= 0.0) and (sval1 <= 1.0)) and ((sval2 >= 0.0) and (sval2 <= 1.0))):
        encounter = True
        interSecPoint = [u1[0]+sval1*(u2[0]-u1[0]),u1[1]+sval1*(u2[1]-u1[1])]  
    return(encounter,numpy.array(interSecPoint))
        
