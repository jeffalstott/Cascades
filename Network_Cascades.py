# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from scipy.io import loadmat
from igraph import Graph, summary
import networkx as nx
from scipy.signal import correlate

# <codecell>

def directed_spr(G, n=1000, weighted='out'):
    #des = []
    #tes = []
    g = G.copy()
    nes = len(g.es)
    for i in range(n):
        e1 = randint(nes)
        e2 = randint(nes)
        while e1==e2: #In case we select the same edge twice, roll again.
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

def rich_nodes(graph, fraction=0.1, highest=True, scores=None):
    """Extracts the "rich club" of the given graph, i.e. the subgraph spanned
    between vertices having the top X% of some score.
 
    Scores are given by the vertex degrees by default.
 
    @param graph:    the graph to work on
    @param fraction: the fraction of vertices to extract; must be between 0 and 1.
    @param highest:  whether to extract the subgraph spanned by the highest or
                     lowest scores.
    @param scores:   the scores themsestrength_sprlves. C{None} uses the vertex degrees.
    """
 
    if scores is None:
        scores = graph.degree()
 
    indices = range(graph.vcount())
    indices.sort(key=scores.__getitem__)
 
    n = int(round(graph.vcount() * fraction))
    if highest:
        indices = indices[-n:]
    else:
        indices = indices[:n]
 
    return indices

# <codecell>

def rich_club_coefficient(graph, fraction=None, highest=True, scores_name=None, rewire=1000, average=1, control=None):
    if type(fraction)==float:
        fraction = [fraction]
        
    if fraction is None:
        from numpy import arange
        fraction = arange(.9,0, -.1)

    from numpy import zeros
    
    if scores_name is None or scores_name=='degree':
        scores = graph.degree()
        randomization = 'out'
    elif scores_name=='out_strength':
        scores = graph.strength(graph.vs, mode=2,weights=graph.es["weight"])
        randomization = 'out'
    elif scores_name=='in_strength':
        scores = graph.strength(graph.vs, mode=1,weights=graph.es["weight"])
        randomization = 'in'
    
    rc_coefficient = zeros(len(fraction))
    
    for i in range(len(fraction)):
    
        node_indices = rich_nodes(graph, fraction=fraction[i], highest=highest, scores=scores)
        
        if scores_name is None or scores_name=='degree':
            numerator = len(g.es.select(_within=node_indices))
            denominator = sum(g.es.select(_target_in=node_indices)["weight"])
        elif scores_name=='out_strength':
            numerator = sum(g.es.select(_within=node_indices)["weight"])
            denominator = sum(g.es.select(_source_in=node_indices)["weight"])
        elif scores_name=='in_strength':    
            numerator = sum(g.es.select(_within=node_indices)["weight"])
            denominator = sum(g.es.select(_target_in=node_indices)["weight"])
            
            
        rc_coefficient[i] = numerator/denominator
    
    if control!=None:
        from numpy import ndarray
        if ~(type(control)==list or type(control)==ndarray):
            control = list(control)
            
        control_rc_coefficient = zeros(len(fraction))
        for i in range(len(control)):
            random_graph = Graph.Weighted_Adjacency(control[i].toarray().tolist())
            control_rc_coefficient = control_rc_coefficient +\
                rich_club_coefficient(random_graph, fraction=fraction, highest=highest, \
                    scores_name=scores_name, rewire=False)
        control_rc_coefficient = control_rc_coefficient/len(control)
        
        return rc_coefficient/control_rc_coefficient
    elif rewire:
        control_rc_coefficient = zeros(len(fraction))
        for i in range(average):
            random_graph = directed_spr(graph, n=rewire, weighted=randomization)
            control_rc_coefficient = control_rc_coefficient +\
                rich_club_coefficient(random_graph, fraction=fraction, highest=highest, \
                    scores_name=scores_name, rewire=False)
        control_rc_coefficient = control_rc_coefficient/average
            
        return rc_coefficient/control
    else:
        return rc_coefficient

# <codecell>

cd /data/alstottjd/Sini

# <codecell>

mat = loadmat('test_NL_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_oho_-6_1_5_ncorr.mat')

# <codecell>

#q = mat['pnets'][0,0][0,0]
g = Graph.Weighted_Adjacency(mat['pnets'][0,0][0,0].toarray().tolist())

# <codecell>

n_nets = shape(mat['pnets'])[1]
n_runs = shape(mat['pnets'][0,0])[1]

q_betweeness = zeros([n_nets, n_runs])
n_betweeness = zeros([n_nets, n_runs])
q_walktrap = zeros([n_nets, n_runs])
n_walktrap = zeros([n_nets, n_runs])
q_infomap = zeros([n_nets, n_runs])
n_infomap = zeros([n_nets, n_runs])
codelength = zeros([n_nets, n_runs])
rc_out = zeros([n_nets, n_runs, 9])
rc_in = zeros([n_nets, n_runs, 9])
rc_out_int = zeros([n_nets, n_runs])
rc_in_int = zeros([n_nets, n_runs])
close('all')
node_size = 100
alpha = .9
width = .4

for i in range(n_nets):
    for j in arange(0,n_runs):#,10):
        if floor(j%50)=0:
            print j

        g = Graph.Weighted_Adjacency(mat['pnets'][0,i][0,j].toarray().tolist())
        
        
        #b = g.community_edge_betweenness(directed=True, weights=g.es["weight"])
        #n_betweeness.append(b.optimal_count)
        #q_betweeness.append(b.as_clustering().q)
        
        #w = g.community_walktrap(weights=g.es["weight"], steps=100)
        #n_walktrap.append(w.optimal_count)
        #q_walktrap.append(w.as_clustering().q)
        
        infomap = g.community_infomap(edge_weights=g.es["weight"])
        n_infomap[i,j] = infomap.cluster_graph().vcount()
        q_infomap[i,j] = infomap.q
        codelength[i, j] = infomap.codelength
        
        rc_out[i,j,:] = rich_club_coefficient(g, scores_name='out_strength', control = mat['pnets_spr_out'][i,j,:])
        rc_in[i,j,:] = rich_club_coefficient(g, scores_name='in_strength', control = mat['pnets_spr_in'][i,j,:])
        rc_out_int[i,j] = sum(rc_out[i,j,:]-1)
        rc_in_int[i,j] = sum(rc_in[i,j,:]-1)
        
        
            
        if j in [0, 180, 480]:
            savetxt('inrichclub_frame%i.txt'%j, rc_in[i,j,:])
            savetxt('outrichclub_frame%i.txt'%j, rc_out[i,j,:])
            
            figure()
            plot(rc_in[i,j,:], 'b')
            plot(rc_out[i,j,:], 'g')
            #show()
            savefig('richclub_frame%i.pdf'%j)
            
        #net = nx.DiGraph(mat['pnets'][0,i][0,j])
        #pos=nx.spring_layout(net)
        #figure()
        #title(str(i)+', '+str(z.as_clustering().q))
        #nx.draw(net,pos,node_size=node_size,alpha=alpha, width=width, with_labels=False)
        #show()

# <codecell>

increment = 1.0/10.0
x_vals = arange(.1, 1, .1)
f = figure()

ax1 = f.add_subplot(131)
handles = {}
handles['In Strength'] = ax1.plot(x_vals, rc_in[0,0, :], label='In Strength Rich Club', color='b')
handles['In Strength'] = ax1.plot(x_vals, rc_out[0,0,:], label='Out Strength Rich Club', color='g')
plt.setp(ax1.get_xticklabels(), visible=False)
ylabel("Rich Club Coefficient")
text(.5, .9, '0 Cascades', transform = ax1.transAxes, horizontalalignment='center')
#xlabel('Strength Decile')
plt.xticks(x_vals[::2], x_vals[::2])
#handles, labels = ax1.get_legend_handles_labels()
#ax1.legend(handles, labels, loc=6)

ax2 = f.add_subplot(132, sharey=ax1)
handles = {}
handles['In Strength'] = ax2.plot(x_vals, rc_in[0,180,:], label='In Strength Rich Club', color='b')
handles['In Strength'] = ax2.plot(x_vals, rc_out[0,180,:], label='Out Strength Rich Club', color='g')
plt.setp(ax2.get_yticklabels(), visible=False)
xlabel('Strength Decile')
plt.xticks(x_vals[::2], (x_vals[::2]*100).astype(int))
text(.5, .9, '36x10$^{5}$ Cascades', transform = ax2.transAxes, horizontalalignment='center')
#handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles, labels, loc=1)

ax3 = f.add_subplot(133, sharey=ax1)
handles = {}
handles['In Strength'] = ax3.plot(x_vals, rc_in[0,480,:], label='In Strength', color='b')
handles['In Strength'] = ax3.plot(x_vals, rc_out[0,480,:], label='Out Strength', color='g')
plt.setp(ax3.get_yticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
text(.5, .9, '96x10$^{5}$ Cascades', transform = ax3.transAxes, horizontalalignment='center')
#xlabel('Strength Decile')
plt.xticks(x_vals[::2], x_vals[::2])

handles, labels = ax3.get_legend_handles_labels()
ax3.legend(handles, labels, loc=5)

suptitle('Rich Club Growth and Death')
savefig('RichClubSamples.pdf')

# <codecell>

close('all')
increment = 10.0/500.0
x_vals = arange(0, 10+increment, increment)

# <codecell>

#Modularity Metrics
#plot(n_betweeness)
figure()
#plot(q_betweeness)
y_vals = mean(n_infomap, axis=0)
error = std(n_infomap, axis=0)
plot(x_vals, y_vals)
fill_between(x_vals, y_vals-error, y_vals+error, alpha=.5)
title("Modularity with Learning")
xlabel("Cascades (n x 10$^{6}$)")
ylabel("Number of Modules")
inset_ax = axes([0.58,0.58,0.3,0.3])
y_vals = mean(q_infomap, axis=0)
error = std(q_infomap, axis=0)
inset_ax.plot(x_vals, y_vals)
inset_ax.fill_between(x_vals, y_vals-error, y_vals+error, alpha=.5)
ylabel("Modularity Index")

savetxt('modules.txt', mean(n_infomap, axis=0))
savetxt('modularity.txt', mean(q_infomap, axis=0))
savefig('Modularity.pdf')

# <codecell>

figure()
y_vals = mean(codelength, axis=0)
error = std(codelength, axis=0)
plot(x_vals, y_vals)
fill_between(x_vals, y_vals-error, y_vals+error, alpha=.5)
title("Network Compression with Learning")
xlabel("Cascades (n x 10$^{6}$)")
ylabel("Bits to Represent Network")

savetxt('compression.txt', mean(codelength, axis=0))
savefig('Compression.pdf')

# <codecell>

f = figure()
ax = f.add_subplot(111)
handles = {}
y_vals = mean(rc_in_int, axis=0)
error = std(rc_in_int, axis=0)
handles['In Strength'] = ax.plot(x_vals, y_vals, label='In Strength Rich Club', color='b')
fill_between(x_vals, y_vals-error, y_vals+error, color='b', alpha=.5)

y_vals = mean(rc_out_int, axis=0)
error = std(rc_out_int, axis=0)
handles['Out Strength'] = ax.plot(x_vals, y_vals, label='Out Strength Rich Club', color='g')
fill_between(x_vals, y_vals-error, y_vals+error, color='g', alpha=.5)

title("Rich Club with Learning")
xlabel("Cascades (n x 10$^{6}$)")
ylabel("Integrated Rich Club Index")

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=1)

savetxt('inrichclubint.txt', rc_in_int)
savetxt('outrichclubint.txt', rc_out_int)
savefig('RichClubInt.pdf')

# <codecell>

x_vals = arange(-1*(shape(rc_in_int)[1]-1), shape(rc_in_int)[1], 1)

# <codecell>

##Different Rich Club Metrics' Correlation with Each Other
f = figure()
ax = f.add_subplot(111)
handles = {}

corr_out_in = zeros([shape(rc_in_int)[0], shape(rc_in_int)[1]*2-1])

for i in range(shape(rc_in_int)[0]):
    corr_out_in[i] = correlate(array(rc_in_int[i]), array(rc_out_int[i]))
##"When does in strength look like out strength?"
#ax.plot(x, mean(corr_out_in, axis=0))
y_vals = mean(corr_out_in, axis=0)
error = std(corr_out_in, axis=0)
plot(x_vals, y_vals)
fill_between(x_vals, y_vals-error, y_vals+error, alpha=.5)

max_corr_out_in = x_vals[argmax(mean(corr_out_in, axis=0))]
ax.plot((max_corr_out_in, max_corr_out_in), ylim(), 'b--')
title("In Strength Rich Club and Out Strength Rich Club Cross Correlation")
xlabel("Cascades (n x 10$^{6}$)")
ylabel("Cross-Correlation")
handles = {}

# <codecell>

def plot_line(x, y, ax, label, color='b'):

    y_vals = mean(y, axis=0)
    error = std(y, axis=0)
    ax.plot(x_vals, y_vals, color, label=label)
    fill_between(x_vals, y_vals-error, y_vals+error, alpha=.5, color = color)

# <codecell>

def correlation_figure(a, b, ax, color='b', label=None):
    
    
    corr = zeros([shape(a)[0], shape(a)[1]*2-1])
    for i in range(shape(a)[0]):
        corr[i] = correlate(array(a[i]), array(b[i]))
    
    x_vals = arange(-1*(shape(a)[1]-1), shape(a)[1], 1)
    
    
    plot_line(x_vals, corr, ax, color = color, label=label)
    
    max_corr = x_vals[argmax(mean(corr, axis=0))]
    
    ax.plot((max_corr, max_corr), ylim(), color+'--')
    xlabel("Cascades (n x 10$^{6}$)")
    ylabel("Cross-Correlation")
    
    

# <codecell>

##Rich Club correlation with Modularity
f = figure()
ax = f.add_subplot(111)
correlation_figure(rc_in_int, q_infomap, ax, 'b', 'In Strength')
correlation_figure(rc_out_int, q_infomap, ax, 'g', 'Out Strength')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=1)
title("Rich Club and Modularity Cross Correlation")

# <codecell>

##Rich Club correlation with Modularity
f = figure()
ax = f.add_subplot(111)
correlation_figure(rc_in_int, n_infomap, ax, 'b', 'In Strength')
correlation_figure(rc_out_int, n_infomap, ax, 'g', 'Out Strength')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=1)
title("Rich Club and Number of Modules Cross Correlation")

# <codecell>

##Rich Club correlation with Modularity
f = figure()
ax = f.add_subplot(111)
correlation_figure(rc_in_int, 1/codelength, ax, 'b', 'In Strength')
correlation_figure(rc_out_int, 1/codelength, ax, 'g', 'Out Strength')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=1)
title("Rich Club and Network Compression Cross Correlation")

# <markdowncell>

# In Strength Rich Clubs trail Out Strength Rich Club
# 
# Both Rich Clubs trail Modularity

# <codecell>

from rpy2.robjects import r
r("source('/home/alstottjd/alluvial/code_for_jeff_20Aug2012.R')")
change_in_community_structure=asarray(r('OUTPUT1'))

# <codecell>

print change_in_community_structure
print diff(rc_in_int[:40])
scatter(change_in_community_structure, diff(rc_in_int[:40]))

# <codecell>

import Aaron

# <codecell>

plot(change_in_community_structure)

# <codecell>

f = figure()
ax = f.add_subplot(111)

ax.scatter(change_in_community_structure, abs(diff(rc_in_int[:40])), c='b', label='In Strength')
ax.scatter(change_in_community_structure, abs(diff(rc_out_int[:40])), c='g', label='Out Strength')

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=1)

