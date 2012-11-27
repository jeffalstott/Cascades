from os import listdir
from Helix import biowulf

n_rewires = 2
n_iters = 5

data_directory = '/data/alstottjd/Sini/'

swarm = biowulf.Swarm(memory_requirement=16)

files = listdir(data_directory)

for filename in files:
    if filename.startswith(str(n_rewires)+','+str(n_iters)) and filename.endswith(".mat") and 'tl_'+filename not in files:
        print filename
        jobstring = """from numpy import shape, arange, ceil, floor
import sys
import os
sys.path.append('/home/alstottjd/Code/Cascade/')

from scipy.io import loadmat
from Timelines import Timelines, Timeline
data_directory = %r
filename = %r
mat = loadmat(data_directory+filename)
n_nets = shape(mat['pnets'])[1]
n_runs = shape(mat['pnets'][0,0])[1]
n_controls = shape(mat['pnets_spr_out'][0,0])[0]
step_size = 5
n_samples = ceil(n_runs/step_size)

T_out = Timelines()
T_in = Timelines()
for i in range(n_nets):
    tl = Timeline()
    CT_out = Timelines(with_control=False)
    CT_in = Timelines(with_control=False)    
    for j in arange(0,n_runs,step_size):
        if floor(j%%100)==0:
            print "Generation = %%i"%%j

        tl.add_gen(mat['pnets'][0,i][0,j].toarray().tolist())
        
        if j==0:
                c_out=[]
                c_in=[]
        for k in range(n_controls):
            #print "k=%%i"%%k
            if j==0:
                c_out.append(Timeline())
                c_in.append(Timeline())
         
            control_out = mat['pnets_spr_out'][i,j,k].toarray().tolist()
            c_out[k].add_gen(control_out)
            control_in = mat['pnets_spr_in'][i,j,k].toarray().tolist()
            c_in[k].add_gen(control_in)
            if j==n_runs-step_size:
                c_out[k].close_gens()
                c_in[k].close_gens()
        
        CT_out.add_timelines(c_out)
        CT_in.add_timelines(c_in)
    
    tl.close_gens()
    T_out.add_timeline(tl, CT_out) 
    T_in.add_timeline(tl, CT_in)

import pickle
data = {"T_out": T_out, "T_in": T_in}
pickle.dump(data, open( data_directory+"tl_"+filename[:-4]+".p", "wb" ) )
        """%(data_directory, filename)
        swarm.add_job(jobstring)

swarm.submit()
