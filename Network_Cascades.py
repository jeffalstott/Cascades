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
from scipy.stats import kendalltau
from numpy import mean, std, asarray
import richclub

# <codecell>

data_dir = "/data/alstottjd/Sini/"

# <codecell>

#mat = loadmat(data_dir+'test_NL_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_ws_-0.0333333-4_1_5_ncorr.mat')
#data_name='WSlast'
mat = loadmat(data_dir+'test_NL_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_oho_-6_1_5_ncorr.mat')
data_name='OHOlast'
#mat = loadmat(data_dir+'2,5_NL_mo-3_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_oho_-6_1_5_ncorr.mat')
#data_name = "2,5OHOF3"
#mat = loadmat(data_dir+'2,5_NL_mall_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_ws_-0.0333333-4_1_5_ncorr.mat')
#data_name = "2,5WSall"

# <codecell>

from Timelines import Timelines, Timeline
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
            if j==n_runs-step_size:
                c_out[k].close_gens()
                c_in[k].close_gens()
        
        CT_out.add_timelines(c_out)
        CT_in.add_timelines(c_in)
    
    tl.close_gens()
    T_out.add_timeline(tl, CT_out) 
    T_in.add_timeline(tl, CT_in)

# <rawcell>

# n_nets = shape(mat['pnets'])[1]
# n_runs = shape(mat['pnets'][0,0])[1]
# n_controls = shape(mat['pnets_spr_out'][0,0])[0]
# step_size = 5.0
# 
# n_samples = ceil(n_runs/step_size)
# 
# q_betweeness = zeros([n_nets, n_samples])
# n_betweeness = zeros([n_nets, n_samples])
# q_walktrap = zeros([n_nets, n_samples])
# n_walktrap = zeros([n_nets, n_samples])
# q_infomap = zeros([n_nets, n_samples])
# n_infomap = zeros([n_nets, n_samples])
# codelength = zeros([n_nets, n_samples])
# mean_path_length = zeros([n_nets, n_samples])
# mean_clustering = zeros([n_nets, n_samples])
# betweeness_change = zeros([n_nets, n_samples])
# 
# rc_out = zeros([n_nets, n_samples, 9])
# rc_in = zeros([n_nets, n_samples, 9])
# rc_out_int = zeros([n_nets, n_samples])
# rc_in_int = zeros([n_nets, n_samples])
# 
# node_size = 100
# alpha = .9
# width = .4
# 
# for i in range(n_nets):
#     last_betweeness = 0
#     for j in arange(0,n_runs, step_size):
#         if floor(j%100)==0:
#             print j
# 
#         g = Graph.Weighted_Adjacency(mat['pnets'][0,i][0,j].toarray().tolist())
#         
#         sample = ceil(j/step_size)
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
#         n_infomap[i,sample] = infomap.cluster_graph().vcount()
#         q_infomap[i,sample] = infomap.q
#         codelength[i, sample] = infomap.codelength
# 
#         mean_path_length[i,sample] = mean(g.shortest_paths(weights='weight'))
#         mean_clustering[i,sample] = mean(g.transitivity_local_undirected(weights='weight'))
#         
#         betweeness_sequence = g.edge_betweenness(weights='weight')
#         if last_betweeness==0:
#             last_betweeness = betweeness_sequence
#         else:
#             betweeness_change[i,sample] = kendalltau(last_betweeness, betweeness_sequence)[0]
#         
#         rc_out[i,sample,:] = richclub.rich_club_coefficient(g, scores_name='out_strength', control = mat['pnets_spr_out'][i,j,:])
#         rc_in[i,sample,:] = richclub.rich_club_coefficient(g, scores_name='in_strength', control = mat['pnets_spr_in'][i,j,:])
#         rc_out_int[i,sample] = sum(rc_out[i,sample,:]-1)
#         rc_in_int[i,sample] = sum(rc_in[i,sample,:]-1)
#         
#         
#             
#         if j in [0, 180, 480]:
#             hist(g.es["weight"])
#             savetxt('inrichclub_frame%i.txt'%j, rc_in[i,sample,:])
#             savetxt('outrichclub_frame%i.txt'%j, rc_out[i,sample,:])
#             
#             figure()
#             plot(rc_in[i,sample,:], 'b')
#             plot(rc_out[i,sample,:], 'g')
#             #show()
#             savefig('richclub_frame%i.pdf'%j)
#             
#         #net = nx.DiGraph(mat['pnets'][0,i][0,j])
#         #pos=nx.spring_layout(net)
#         #figure()
#         #title(str(i)+', '+str(z.as_clustering().q))
#         #nx.draw(net,pos,node_size=node_size,alpha=alpha, width=width, with_labels=False)
#         #show()
#         

# <codecell>

increment = 1.0/10.0
x_vals = arange(.1, 1, .1)
w= 3.375
h = w/1.6180
#f = figure(figsize=(w,h))
f = figure()

n_cascades = 10.0
cascades_per_run = n_cascades/n_runs

c1 = 0
c2 = 3.6
c3 = 9.6

r1 = round(c1/cascades_per_run)
r2 = round(c2/cascades_per_run)
r3 = round(c3/cascades_per_run)

#c1 = r1*cascades_per_step/10**6
#c2 = r2*cascades_per_step/10**6
#c3 = r3*cascades_per_step/10**6

s1 = floor((r1/n_runs)*n_samples)
s2 = floor((r2/n_runs)*n_samples)
s3 = floor((r3/n_runs)*n_samples)

rc_in = T_out.raw_normalized("rc_in")
rc_out = T_out.raw_normalized("rc_out")


alpha = .1

xlims = (.1,.9)



ax1 = f.add_subplot(131)

y_vals = mean(rc_in[:,s1,:], axis=0)
error = std(rc_in[:,s1, :], axis=0)
ax1.plot(x_vals, y_vals, label='In Strength', color='b')
ax1.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color='b')

y_vals = mean(rc_out[:,s1,:], axis=0)
error = std(rc_out[:,s1, :], axis=0)
ax1.plot(x_vals, y_vals, label='Out Strength', color='g')
ax1.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color='g')

ax1.set_xlim(xlims)

ax1.plot(xlims, (1,1), 'k--')
plt.setp(ax1.get_xticklabels(), visible=False)
ylabel("Rich Club Coefficient")
text(.5, .9, '0 Cascades', transform = ax1.transAxes, horizontalalignment='center', fontsize=10)
#xlabel('Strength Decile')
plt.xticks(x_vals[::2], x_vals[::2])
#handles, labels = ax1.get_legend_handles_labels()
#ax1.legend(handles, labels, loc=6)

ax2 = f.add_subplot(132, sharey=ax1)
y_vals = mean(rc_in[:,s2,:], axis=0)
error = std(rc_in[:,s2, :], axis=0)
ax2.plot(x_vals, y_vals, label='In Strength', color='b')
ax2.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color='b')

y_vals = mean(rc_out[:,s2,:], axis=0)
error = std(rc_out[:,s2, :], axis=0)
ax2.plot(x_vals, y_vals, label='Out Strength', color='g')
ax2.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color='g')

ax2.set_xlim(xlims)
ax2.plot(xlims, (1,1), 'k--')
plt.setp(ax2.get_yticklabels(), visible=False)
xlabel('Strength Decile')
plt.xticks(x_vals[::2], (x_vals[::2]*100).astype(int))
text(.5, .9, '$%.1f*10^{6}$ Cascades'%c2, transform = ax2.transAxes, horizontalalignment='center', fontsize=10)
#handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles, labels, loc=1)

ax3 = f.add_subplot(133, sharey=ax1)
y_vals = mean(rc_in[:,s3,:], axis=0)
error = std(rc_in[:,s3, :], axis=0)
ax3.plot(x_vals, y_vals, label='In Strength', color='b')
ax3.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color='b')

y_vals = mean(rc_out[:,s3,:], axis=0)
error = std(rc_out[:,s3, :], axis=0)
ax3.plot(x_vals, y_vals, label='Out Strength', color='g')
ax3.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color='g')

ax3.set_xlim(xlims)
ax3.plot(xlims, (1,1), 'k--')
plt.setp(ax3.get_yticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
text(.5, .9, '$%.1f*10^{6}$ Cascades'%c3, transform = ax3.transAxes, horizontalalignment='center', fontsize=10)
#xlabel('Strength Decile')
plt.xticks(x_vals[::2], x_vals[::2])

handles, labels = ax3.get_legend_handles_labels()
ax3.legend(handles, labels, loc=5, fontsize=8)

ax3.set_ylim(0,6)


savefig(data_dir+data_name+'_'+'RichClubSamples.pdf', bbox_inches='tight')
#suptitle('Rich Club Growth and Death')

# <codecell>

close('all')
x_vals = arange(0, n_runs, step_size)
x_vals = 10.0*x_vals/n_runs

# <codecell>

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#Modularity Metrics
figure()
y_vals, error = T_out.data("n_infomap")
#mean(n_infomap, axis=0)
#error = std(n_infomap, axis=0)
plot(x_vals, y_vals)
fill_between(x_vals, y_vals-error, y_vals+error, alpha=.5)

xlabel("Cascades (n x 10$^{6}$)")
ylabel("Number of Modules")
inset_ax = axes([0.58,0.58,0.3,0.3])
y_vals, error = T_out.data("q_infomap")
#y_vals = mean(q_infomap, axis=0)
#error = std(q_infomap, axis=0)
inset_ax.plot(x_vals, y_vals)
inset_ax.fill_between(x_vals, y_vals-error, y_vals+error, alpha=.5)
ylabel("Modularity Index", fontsize=10)
inset_ax.yaxis.set_major_locator(MultipleLocator(.1))



#savetxt(data_dir+data_name+'_'+'modules.txt', mean(n_infomap, axis=0))
#savetxt(data_dir+data_name+'_'+'modularity.txt', mean(q_infomap, axis=0))
savefig(data_dir+data_name+'_'+'Modularity.pdf', bbox_inches='tight')
#title("Modularity with Learning")

# <codecell>

figure()
y_vals, error = T_out.data("codelength")
#y_vals = mean(codelength, axis=0)
#error = std(codelength, axis=0)
plot(x_vals, y_vals)
fill_between(x_vals, y_vals-error, y_vals+error, alpha=.5)

xlabel("Cascades (n x 10$^{6}$)")
ylabel("Bits to Represent Network")

#savetxt(data_dir+data_name+'_'+'compression.txt', mean(codelength, axis=0))
savefig(data_dir+data_name+'_'+'Compression.pdf', bbox_inches='tight')
#title("Network Compression with Learning")

# <codecell>

f = figure()
ax = f.add_subplot(111)


data = sum(T_in.raw_normalized("rc_in")-1, axis=-1)
y_vals = mean(data, axis=0)
error = std(data, axis=0)
ax.plot(x_vals, y_vals, label='In Strength', color='b')
fill_between(x_vals, y_vals-error, y_vals+error, color='b', alpha=.5)

data = sum(T_out.raw_normalized("rc_out")-1, axis=-1)
y_vals = mean(data, axis=0)
error = std(data, axis=0)
ax.plot(x_vals, y_vals, label='Out Strength', color='g')
fill_between(x_vals, y_vals-error, y_vals+error, color='g', alpha=.5)

plot(xlim(), (0,0), 'k--')


xlabel("Cascades (n x 10$^{6}$)")
ylabel("Integrated Rich Club Index")

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=1)

#savetxt(data_dir+data_name+'_'+'inrichclubint.txt', rc_in_int)
#savetxt(data_dir+data_name+'_'+'outrichclubint.txt', rc_out_int)
savefig(data_dir+data_name+'_'+'RichClubInt.pdf', bbox_inches='tight')
#title("Rich Club with Learning")

# <codecell>

f = figure()
ax = f.add_subplot(111)

y_vals, error = T_out.normalized("mean_path_length")
#error = std(mean_path_length, axis=0)
ax.plot(x_vals, y_vals, label='Path Length', color='b')
fill_between(x_vals, y_vals-error, y_vals+error, color='b', alpha=.5)

y_vals, error = T_out.normalized("mean_clustering")
#y_vals = mean(mean_clustering, axis=0)
#error = std(mean_clustering, axis=0)
ax.plot(x_vals, y_vals, label='Clustering', color='g')
fill_between(x_vals, y_vals-error, y_vals+error, color='g', alpha=.5)


y_vals = T_out.normalized("mean_clustering")[0]/T_out.normalized("mean_path_length")[0]
ax.plot(x_vals, y_vals, label='Small World', color='k', linewidth=2)


xlabel("Cascades (n x 10$^{6}$)")
#ylabel("Integrated Rich Club Index")

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=1)

savefig(data_dir+data_name+'_'+'PathLengthClustering.pdf', bbox_inches='tight')

# <codecell>

figure()
y_vals, error = T_out.data("betweeness_change_kendall")
#y_vals = mean(betweeness_change, axis=0)[1:]
#error = std(betweeness_change, axis=0)[1:]
plot(x_vals[1:], y_vals, label='Original Networks', color='b')
fill_between(x_vals[1:], y_vals-error, y_vals+error, alpha=.5)

y_vals, error = T_out.control("betweeness_change_kendall")
#y_vals = mean(betweeness_change, axis=0)[1:]
#error = std(betweeness_change, axis=0)[1:]
plot(x_vals[1:], y_vals, label='Random Controls', color='g')
fill_between(x_vals[1:], y_vals-error, y_vals+error, alpha=.5)

xlabel("Cascades (n x 10$^{6}$)")
ylabel("Correlation in Nodes' Betweeness Centrality Rank Order")

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=1)

#savetxt(data_dir+data_name+'_'+'Rerouting.txt', mean(y_vals, axis=0))
savefig(data_dir+data_name+'_'+'Rerouting.pdf', bbox_inches='tight')

# <codecell>

def plot_line(x, y, ax, label, color='b'):

    y_vals = mean(y, axis=0)
    error = std(y, axis=0)
    ax.plot(x, y_vals, color, label=label)
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

##Different Rich Club Metrics' Correlation with Each Other
##Rich Club correlation with Modularity
f = figure()
ax = f.add_subplot(111)
correlation_figure(rc_in_int, rc_out_int, ax)

xlabel("Cascades (n x 10$^{6}$)")
ylabel("Cross-Correlation")
title("In Strength Rich Club and Out Strength Rich Club Cross Correlation")

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

