from os import listdir
from Helix import biowulf

n_rewires = 10
n_iters = 5

data_directory = '/data/alstottjd/Cascade/Controlled_Data/'
output_directory = '/data/alstottjd/Cascade/Timelines/'
swarm = biowulf.Swarm(memory_requirement=16)

files = listdir(data_directory)

for filename in files:
    if filename.startswith(str(n_rewires)+','+str(n_iters))\
        and filename.endswith(".mat") and 'tl_'+filename not in files\
        and ('mlast' in filename or 'mo-1' in filename\
        or 'mo-2' in filename or 'mall' in filename):
        print filename
        jobstring = """from numpy import shape, arange, ceil, floor
import sys
import os
sys.path.append('/home/alstottjd/Code/Cascade/')

from scipy.io import loadmat
from Timelines import Timelines, Timeline
data_directory = %r
output_directory = %r
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
    tl = Timeline(mat['pnets'][0,i][0,::step_size])
    tl.calculate()

    CT_out = Timelines(with_control=False)
    CT_in = Timelines(with_control=False)
    for k in range(n_controls):
        co = Timeline(mat['pnets_spr_out'][i,::step_size,k])
        co.calculate()
        CT_out.add_timeline(co)

        ci = Timeline(mat['pnets_spr_in'][i,::step_size,k])
        ci.calculate()
        CT_in.add_timeline(ci)

    T_out.add_timeline(tl, CT_out) 
    T_in.add_timeline(tl, CT_in)

import pickle
data = {"T_out": T_out, "T_in": T_in}
pickle.dump(data, open( output_directory+"tl_"+filename[:-4]+".p", "wb" ) )
        """%(data_directory, output_directory, filename)
        swarm.add_job(jobstring)

swarm.submit()
