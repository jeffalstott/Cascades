# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from matplotlib import rc
rc('text', usetex=False)
rc('font', family='serif')

# <codecell>

from scipy.io import loadmat
from igraph import Graph
from scipy.signal import correlate
from Timelines import Timelines, Timeline
import richclub

# <codecell>

data_dir = "/data/alstottjd/Sini/"

# <codecell>

mat = loadmat(data_dir+'test_NL_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_ws_-0.0333333-4_1_5_ncorr.mat')
#mat = loadmat(data_dir+'test_NL_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_oho_-6_1_5_ncorr.mat')

# <codecell>

n_nets = shape(mat['pnets'])[1]
n_runs = shape(mat['pnets'][0,0])[1]
n_controls = shape(mat['pnets_spr_out'][0,0])[0]
step_size = 50

T_out = Timelines()
T_in = Timelines()
for i in range(n_nets):
    tl = Timeline()
    CT_out = Timelines(with_control=False)
    CT_in = Timelines(with_control=False)    
    for j in arange(0,n_runs,step_size):
        if floor(j%100)==0:
            print "Generation = %i"%j

        tl.add_gen(mat['pnets'][0,i][0,j].toarray().tolist())
        
        if j==0:
                c_out=[]
                c_in=[]
        for k in range(n_controls):
            #print "k=%i"%k
            if j==0:
                c_out.append(Timeline())
                c_in.append(Timeline())
         
            control_out = mat['pnets_spr_out'][i,j,k].toarray().tolist()
            c_out[k].add_gen(control_out)
            control_in = mat['pnets_spr_in'][i,j,k].toarray().tolist()
            c_in[k].add_gen(control_in)
        
        CT_out.add_timelines(c_out)
        CT_in.add_timelines(c_in)
        
    T_out.add_timeline(tl, CT_out) 
    T_in.add_timeline(tl, CT_in)
        

# <rawcell>

# n_nets = shape(mat['pnets'])[1]
# n_runs = shape(mat['pnets'][0,0])[1]
# 
# q_betweeness = zeros([n_nets, n_runs])
# n_betweeness = zeros([n_nets, n_runs])
# q_walktrap = zeros([n_nets, n_runs])
# n_walktrap = zeros([n_nets, n_runs])
# q_infomap = zeros([n_nets, n_runs])
# n_infomap = zeros([n_nets, n_runs])
# codelength = zeros([n_nets, n_runs])
# mean_path_length = zeros([n_nets, n_runs])
# mean_clustering = zeros([n_nets, n_runs])
# betweeness_change = zeros([n_nets, n_runs])
# 
# rc_out = zeros([n_nets, n_runs, 9])
# rc_in = zeros([n_nets, n_runs, 9])
# rc_out_int = zeros([n_nets, n_runs])
# rc_in_int = zeros([n_nets, n_runs])
# close('all')
# node_size = 100
# alpha = .9
# width = .4
# 
# for i in range(n_nets):
#     last_betweeness = 0
#     for j in arange(0,n_runs):#,10):
#         if floor(j%50)==0:
#             print j
# 
#         g = Graph.Weighted_Adjacency(mat['pnets'][0,i][0,j].toarray().tolist())
#         
#         
#         #b = g.community_edge_betweenness(directed=True, weights=g.es["weight"])
#         #n_betweeness.append(b.optimal_count)
#         #q_betweeness.append(b.as_clustering().q)
#         
#         #w = g.community_walktrap(weights=g.es["weight"], steps=100)
#         #n_walktrap.append(w.optimal_count)
#         #q_walktrap.append(w.as_clustering().q)
#         
#         infomap = g.community_infomap(edge_weights=g.es["weight"])
#         n_infomap[i,j] = infomap.cluster_graph().vcount()
#         q_infomap[i,j] = infomap.q
#         codelength[i, j] = infomap.codelength
# 
#         mean_path_length[i,j] = mean(g.shortest_paths(weights='weight'))
#         mean_clustering[i,j] = mean(g.transitivity_local_undirected(weights='weight'))
#         
#         betweeness_sequence = g.edge_betweenness(weights='weight')
#         if last_betweeness==0:
#             last_betweeness = betweeness_sequence
#         else:
#             betweeness_change[i,j] = kendalltau(last_betweeness, betweeness_sequence)[0]
#         
#         rc_out[i,j,:] = rich_club_coefficient(g, scores_name='out_strength', control = mat['pnets_spr_out'][i,j,:])
#         rc_in[i,j,:] = rich_club_coefficient(g, scores_name='in_strength', control = mat['pnets_spr_in'][i,j,:])
#         rc_out_int[i,j] = sum(rc_out[i,j,:]-1)
#         rc_in_int[i,j] = sum(rc_in[i,j,:]-1)
#         
#         
#             
#         if j in [0, 180, 480]:
#             savetxt('inrichclub_frame%i.txt'%j, rc_in[i,j,:])
#             savetxt('outrichclub_frame%i.txt'%j, rc_out[i,j,:])
#             
#             figure()
#             plot(rc_in[i,j,:], 'b')
#             plot(rc_out[i,j,:], 'g')
#             #show()
#             savefig('richclub_frame%i.pdf'%j)
#             
#         #net = nx.DiGraph(mat['pnets'][0,i][0,j])
#         #pos=nx.spring_layout(net)
#         #figure()
#         #title(str(i)+', '+str(z.as_clustering().q))
#         #nx.draw(net,pos,node_size=node_size,alpha=alpha, width=width, with_labels=False)
#         #show()

# <codecell>

increment = 1.0/10.0
x_vals = arange(.1, 1, .1)
w= 3.375
h = w/1.6180
#f = figure(figsize=(w,h))
f = figure()

ax1 = f.add_subplot(131)
handles = {}
handles['In Strength'] = ax1.plot(x_vals, T_in.rc_in[0]rc_in[0,0, :], label='In Strength Rich Club', color='b')
handles['In Strength'] = ax1.plot(x_vals, rc_out[0,0,:], label='Out Strength Rich Club', color='g')
plt.setp(ax1.get_xticklabels(), visible=False)
ylabel("Rich Club Coefficient")
text(.5, .9, '0 Cascades', transform = ax1.transAxes, horizontalalignment='center', fontsize=10)
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
text(.5, .9, '36x10$^{5}$ Cascades', transform = ax2.transAxes, horizontalalignment='center', fontsize=10)
#handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles, labels, loc=1)

ax3 = f.add_subplot(133, sharey=ax1)
handles = {}
handles['In Strength'] = ax3.plot(x_vals, rc_in[0,480,:], label='In Strength', color='b')
handles['In Strength'] = ax3.plot(x_vals, rc_out[0,480,:], label='Out Strength', color='g')
plt.setp(ax3.get_yticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
text(.5, .9, '96x10$^{5}$ Cascades', transform = ax3.transAxes, horizontalalignment='center', fontsize=10)
#xlabel('Strength Decile')
plt.xticks(x_vals[::2], x_vals[::2])

handles, labels = ax3.get_legend_handles_labels()
ax3.legend(handles, labels, loc=5, fontsize=8)


savefig(data_dir+'RichClubSamples.pdf', bbox_inches='tight')
#suptitle('Rich Club Growth and Death')

# <codecell>

close('all')
x_vals = arange(0, n_runs, step_size)
x_vals = 10.0*x_vals/n_runs

# <codecell>

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#Modularity Metrics
#plot(n_betweeness)
figure()
#plot(q_betweeness)
y_vals, error  = T_out.__getattr__("n_infomap", version='raw')
plot(x_vals, y_vals)
fill_between(x_vals, y_vals-error, y_vals+error, alpha=.5)

xlabel("Cascades (n x 10$^{6}$)")
ylabel("Number of Modules")
inset_ax = axes([0.58,0.58,0.3,0.3])
y_vals, error  = T_out.__getattr__("q_infomap", version='raw')
inset_ax.plot(x_vals, y_vals)
inset_ax.fill_between(x_vals, y_vals-error, y_vals+error, alpha=.5)
ylabel("Modularity Index", fontsize=10)
inset_ax.yaxis.set_major_locator(MultipleLocator(.1))


#savetxt('modules.txt', mean(n_infomap, axis=0))
#savetxt('modularity.txt', mean(q_infomap, axis=0))
savefig(data_dir+'Modularity.pdf', bbox_inches='tight')
#title("Modularity with Learning")

# <codecell>

figure()
y_vals, error = T_out.__getattr__("codelength", version='normalized')
plot(x_vals, y_vals)
fill_between(x_vals, y_vals-error, y_vals+error, alpha=.5)

xlabel("Cascades (n x 10$^{6}$)")
ylabel("Bits to Represent Network")

#savetxt('compression.txt', mean(codelength, axis=0))
savefig(data_dir+'Compression.pdf', bbox_inches='tight')
#title("Network Compression with Learning")

# <codecell>

f = figure()
ax = f.add_subplot(111)
handles = {}
y_vals, error = T_in.__getattr__("rc_in_int", version='raw')
handles['In Strength'] = ax.plot(x_vals, y_vals, label='In Strength', color='b')
fill_between(x_vals, y_vals-error, y_vals+error, color='b', alpha=.5)

y_vals, error = T_out.__getattr__("rc_out_int", version='raw')
handles['Out Strength'] = ax.plot(x_vals, y_vals, label='Out Strength', color='g')
fill_between(x_vals, y_vals-error, y_vals+error, color='g', alpha=.5)


xlabel("Cascades (n x 10$^{6}$)")
ylabel("Integrated Rich Club Index")

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=1)

#savetxt('inrichclubint.txt', rc_in_int)
#savetxt('outrichclubint.txt', rc_out_int)
savefig(data_dir+'RichClubInt.pdf', bbox_inches='tight')
#title("Rich Club with Learning")

# <codecell>

f = figure()
ax = f.add_subplot(111)
handles = {}
y_vals, error = T_out.mean_path_length
handles['In Strength'] = ax.plot(x_vals, y_vals, label='Path Length', color='b')
fill_between(x_vals, y_vals-error, y_vals+error, color='b', alpha=.5)

y_vals, error = T_out.mean_clustering
handles['Out Strength'] = ax.plot(x_vals, y_vals, label='Clustering', color='g')
fill_between(x_vals, y_vals-error, y_vals+error, color='g', alpha=.5)


xlabel("Cascades (n x 10$^{6}$)")
#ylabel("Integrated Rich Club Index")

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=1)

#savetxt('inrichclubint.txt', rc_in_int)
#savetxt('outrichclubint.txt', rc_out_int)
#savefig(data_dir+'RichClubInt.pdf', bbox_inches='tight')
#title("Rich Club with Learning")

# <codecell>

figure()
y_vals, error = T_out.__getattr__("betweeness_change_kendall", version='raw')
plot(x_vals[1:], y_vals)
fill_between(x_vals[1:], y_vals-error, y_vals+error, alpha=.5)

xlabel("Cascades (n x 10$^{6}$)")
ylabel("betweeness_change")

#savetxt('compression.txt', mean(codelength, axis=0))
#savefig(data_dir+'Compression.pdf', bbox_inches='tight')
#title("Network Compression with Learning")

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

xlabel("Cascades (n x 10$^{6}$)")
ylabel("Cross-Correlation")
title("In Strength Rich Club and Out Strength Rich Club Cross Correlation")

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

