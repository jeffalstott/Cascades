from os import listdir
from Helix import biowulf

n_rewires = 10
randomizations = ['in', 'out']
n_iters = 5

data_directory = '/data/alstottjd/Sini/'

swarm = biowulf.Swarm(memory_requirement=72)

files = listdir(data_directory)

for filename in files:
    if filename.startswith('NL_m') and '_60_' in filename and str(n_rewires)+','+str(n_iters)+'_'+filename not in files:
        print filename
        jobstring = """data_directory = %r
filename = %r
print filename
n_rewires = %i
randomizations = ['in', 'out']
n_iters = %i
from numpy import shape, empty
import richclub
from scipy.io import loadmat, savemat
from igraph import Graph
from scipy.sparse import csc, csc_matrix
mat = loadmat(data_directory+filename)
n_nets = shape(mat['pnets'])[1]
n_runs = shape(mat['pnets'][0,0])[1]
for randomization in randomizations:
    print randomization
    mat['pnets_spr_'+randomization] = empty([n_nets, n_runs, n_iters], dtype=csc.csc_matrix)
    for i in range(n_nets):
        for j in range(n_runs):
            if j%%100==0:
                print i, j
            g = Graph.Weighted_Adjacency(mat['pnets'][0,i][0,j].toarray().tolist())
            for k in range(n_iters):
                random_graph = richclub.directed_spr(g, n_rewires=n_rewires, preserve=randomization)
                mat['pnets_spr_'+randomization][i, j, k] = csc_matrix(random_graph.get_adjacency(attribute='weight').data)
savemat(data_directory+str(n_rewires)+','+str(n_iters)+'_'+filename, mat)
        """%(data_directory, filename, n_rewires, n_iters)
        swarm.add_job(jobstring)

swarm.submit()
