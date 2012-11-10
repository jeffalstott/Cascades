# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

data_directory = '/data/alstottjd/Sini/'

# <codecell>

from scipy.io import loadmat, savemat
from igraph import Graph
from scipy.sparse import csc, csc_matrix

# <codecell>

def directed_spr(G, n=1000, weighted='out'):
    #des = []
    #tes = []
    g = G.copy()
    nes = len(g.es)
    for i in range(n):
        e1 = randint(nes)
        e2 = randint(nes)
        #In case we select the same edge twice, roll again.
        #If either of the to-be-newly-wired connections already exist, roll again.
        #This prevents multiple edges going in the same direction between two nodes.
        while e1==e2
            or : 
            e2 = randint(nes)
            
        s1 = g.es[e1].source
        t1 = g.es[e1].target
        a1 = g.es[e1].attributes()
        s2 = g.es[e2].source
        t2 = g.es[e2].target
        a2 = g.es[e2].attributes()
        
        n = len(g.es)
        g.delete_edges([e1, e2])
        m = len(g.es)
        #print(m-n)
        if weighted=='out': #Rewire the outgoing connections
            g.add_edge(s1, t2, **a1)
            g.add_edge(s2, t1, **a2)
        elif weighted=='in': #Rewire the incoming connections
            #Only difference is in the weight assignments
            g.add_edge(s1, t2, **a2)
            g.add_edge(s2, t1, **a1)
        #Note that in the unweighted condition these two methods are equivalent, so either option
        #for 'weighted' produces the same (correct) randomization, while preserving in degree
        # and out degree for each node
            
            
        #des.append(len(g.es)-n)
        #tes.append(len(g.es))
    return g

# <codecell>

filename = 'NL_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_ws_-0.0333333-4_1_5_ncorr.mat'

# <codecell>

print data_directory
print filename

# <codecell>

mat = loadmat(data_directory+filename)#, variable_names = ['pnets',])

# <codecell>

normalized = 5000
randomizations = ['in', 'out']
n_iters = 10

# <codecell>

n_nets = shape(mat['pnets'])[1]
n_runs = shape(mat['pnets'][0,0])[1]

# <codecell>

for randomization in randomizations:
    print randomization
    mat['pnets_spr_'+randomization] = empty([n_nets, n_runs, n_iters], dtype=csc.csc_matrix)           
    for i in range(n_nets):
        for j in range(n_runs):
            if j%100==0:
                print i, j
            g = Graph.Weighted_Adjacency(mat['pnets'][0,i][0,j].toarray().tolist())
            for k in range(n_iters):
                random_graph = directed_spr(g, n=normalized, weighted=randomization)
                print random_graph.ecount()
                print csc_matrix(random_graph.get_adjacency(attribute='weight').data)
                mat['pnets_spr_'+randomization][i, j, k] = csc_matrix(random_graph.get_adjacency(attribute='weight').data)
                imshow(mat['pnets_spr_in'][0,0,0].toarray())
                break
            break
        break
    break

# <codecell>

x = asarray(ceil(random_graph.get_adjacency(attribute='weight').data))-ceil(asarray(random_graph.get_adjacency().data))
print sum(x)
imshow(x)

# <codecell>

savemat(data_directory+'test_'+filename, mat)

