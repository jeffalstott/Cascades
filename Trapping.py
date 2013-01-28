# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import matplotlib.gridspec as gridspec
from matplotlib import rc
rc('text', usetex=False)
rc('font', family='serif')

# <codecell>

#data_directory = '/data/alstottjd/Cascade/Controlled_Data/'
#data_file = '10,5_NL_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_oho_-6_1_5_ncorr.mat'
data_directory = '/data/alstottjd/Sini/Original/'
if network_type in ('OHO', 'oho'):
    nets_file = 'NL_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_oho_-6_1_5_ncorr.mat'
    summary_file = 'NL_summary_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_oho_-6_1_5_ncorr.mat'
    activity_file = 'NL_activity_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_oho_-6_1_5_ncorr.mat'
elif network_type in ('WS', 'ws', 'wn', 'WN'):
    nets_file = 'NL_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_ws_-0.0333333-4_1_5_ncorr.mat'
    summary_file = 'NL_summary_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_ws_-0.0333333-4_1_5_ncorr.mat'
    activity_file = 'NL_activity_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_ws_-0.0333333-4_1_5_ncorr.mat'

from scipy.io import loadmat
nets = loadmat(data_directory+nets_file)
summary = loadmat(data_directory+summary_file)

# <codecell>

from igraph import Graph

# <codecell>

n_samples = 500
n_runs_per_sample = 20000
step_size = 1
x_vals = arange(0, n_samples+1, step_size)
x_vals = 10.0*x_vals/n_samples
x_vals = x_vals[1:]
alpha = .1

# <codecell>

activity = loadmat(data_directory+activity_file)

# <codecell>

summary['outstrend']

# <codecell>

node_out_strength

# <codecell>

node_os

# <codecell>

#####
replicate = 4
initial_graph = Graph.Weighted_Adjacency(nets['pnets'][0,replicate][0,0].toarray().tolist())
n_nodes = len(initial_graph.vs)
degree_sequence = initial_graph.degree(mode=1)

node_dynamical_degree = zeros([n_samples, n_nodes])
node_dynamical_out_strength = zeros([n_samples, n_nodes])
node_out_strength = zeros([n_samples, n_nodes])

system_dynamical_degree = zeros([n_samples, n_nodes])
system_dynamical_out_strength = zeros([n_samples, n_nodes])
system_out_strength = zeros([n_samples, n_nodes])
system_dynamical_degree_norm = zeros([n_samples, n_nodes])
system_dynamical_out_strength_norm = zeros([n_samples, n_nodes])

for sample in range(n_samples):
    if sample%100==0:
        print sample
    deg = zeros(n_runs_per_sample)
    os = zeros(n_runs_per_sample)
    dyndeg = zeros(n_runs_per_sample)
    dynos = zeros(n_runs_per_sample)
    dyndeg_norm = zeros(n_runs_per_sample)
    dynos_norm = zeros(n_runs_per_sample)

    node_dyndeg = [[] for i in range(n_nodes)]
    node_dyndeg_norm = [[] for i in range(n_nodes)]
    node_dynos = [[] for i in range(n_nodes)]
    node_dynos_norm = [[] for i in range(n_nodes)]
    node_os = [[] for i in range(n_nodes)]
    
    this_graph = Graph.Weighted_Adjacency(nets['pnets'][0,replicate][0,sample].toarray().tolist())
    strength_sequence = this_graph.strength(mode=1, weights='weight')
    for run in range(n_runs_per_sample):
        endnodes = atleast_1d(activity['endnodes'][sample*n_runs_per_sample+run][0].squeeze())-1
        activenodes = atleast_1d(activity['activenodes'][sample*n_runs_per_sample+run][0].squeeze())-1
        for endnode in endnodes:
            deg[run] = degree_sequence[endnode]
            dyndeg[run] = deg[run]
            os[run] = strength_sequence[endnode]
            dynos[run] = os[run]
            neighbors = this_graph.neighbors(endnode, mode=1)
            
            for i in neighbors:
                if i in activenodes:
                    dyndeg[run] -= 1
                    dynos[run] -= this_graph.es[this_graph.get_eid(endnode, i)]['weight']
            node_dyndeg[endnode].append(dyndeg[run])
            node_dynos[endnode].append(dynos[run])
            node_os[endnode].append(os[run])

            dyndeg_norm[run] = dyndeg[run]/deg[run]
            dynos_norm[run] = dynos[run]/os[run]
  
    node_dynamical_degree[sample, :] = [mean(i) for i in node_dyndeg]
    node_dynamical_out_strength[sample, :] = [mean(i) for i in node_dynos]
    node_out_strength[sample,:] = [mean(i) for i in node_os]

    system_dynamical_degree[sample] = mean(dyndeg)
    system_dynamical_out_strength[sample] = mean(dynos)

    system_dynamical_degree_norm[sample] = mean(dyndeg_norm)
    system_dynamical_out_strength_norm[sample] = mean(dynos_norm)

# <codecell>

output_directory = '/data/alstottjd/Cascade/'
figures = []
fn = 0

# <codecell>

#####
close('all')
f = figure(figsize=(11,8))
gs = gridspec.GridSpec(4,16)
gs.update(hspace=0.25, wspace=0.3)
alpha = 1
s = 20
linewidths=.15

TA = summary['totactivity'].astype('float')/(10**7/500.0)

clustering_sequence = initial_graph.transitivity_local_undirected()
from matplotlib.cm import get_cmap
cmap = get_cmap('spectral')
colors = [cmap(i) for i in clustering_sequence]
coloring_sequence = array(clustering_sequence)
from matplotlib.collections import LineCollection
from matplotlib.colors import NoNorm
norm=None #NoNorm()


####
termination_plot = f.add_subplot(gs[0,:7])
termination_plot.annotate("A", (0,0.95), xycoords=(termination_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)

segs = list(zip(x_vals, TA[:,i]) for i in range(n_nodes))
line_segments = LineCollection(segs, cmap=cmap, norm=norm, linewidths=linewidths)
line_segments.set_array(coloring_sequence)
termination_plot.add_collection(line_segments)

termination_plot.set_ylim(0, .1)
#termination_plot.set_ylim(200, 400)
termination_plot.set_ylabel(r'$AF$, Activation Freq.')

for i in termination_plot.get_xticklabels():
    i.set_visible(False)

####
strength_plot = f.add_subplot(gs[1,:7], sharex=termination_plot)
strength_plot.annotate("C", (0,0.95), xycoords=(strength_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)

segs = list(zip(x_vals, summary['inpstrend'][:,i]) for i in range(n_nodes))

line_segments = LineCollection(segs, cmap=cmap, norm=norm, linewidths=linewidths)
line_segments.set_array(coloring_sequence)
strength_plot.add_collection(line_segments)

strength_plot.set_ylim(.4, 2)
strength_plot.set_ylabel(r'$IS$, In Strength')
for i in strength_plot.get_xticklabels():
    i.set_visible(False)

####
termination_strength_corr_plot = f.add_subplot(gs[2,:7], sharex=termination_plot)
termination_strength_corr_plot.annotate("E", (0,0.95), xycoords=(termination_strength_corr_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)


termination_strength_corr_plot.plot(x_vals, [pearsonr(summary['inpstrend'][t,:], TA[t,:])[0] for t in range(len(x_vals))])

termination_strength_corr_plot.set_ylim(-1, 1)
termination_strength_corr_plot.set_ylabel(r'$R_{AF-IS}$')
#for i in strength_plot.get_xticklabels():
#    i.set_visible(False)
termination_strength_corr_plot.set_xlabel("Cascades (n x 10$^{6}$)")
termination_strength_corr_plot.set_xlim(x_vals[0], x_vals[-1])
termination_strength_corr_plot.set_xticks(arange(x_vals[-1]+1))

####
from mpl_toolkits.axes_grid.inset_locator import inset_axes
#termination_strength_corr_plot = inset_axes(termination_strength_corr_plot,
#                        width="30%", # width = 30% of parent_bbox
#                        height="50%",
#                        loc=1)

termination_strength_corr_plot = f.add_subplot(gs[3,:7])
termination_strength_corr_plot.annotate("G", (0,0.95), xycoords=(termination_strength_corr_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)
                                            
from matplotlib import colors, cm
norm = colors.normalize(coloring_sequence.min(), coloring_sequence.max())
cmap = cm.spectral

for i in range(n_nodes):
    termination_strength_corr_plot.scatter(coloring_sequence[i], pearsonr(summary['inpstrend'][:,i], TA[:,i])[0], s=s, color=cmap(norm(coloring_sequence[i])))

termination_strength_corr_plot.set_ylim(-1, 1)
termination_strength_corr_plot.set_ylabel(r'$R_{AF-IS}$ by node')
#for i in strength_plot.get_xticklabels():
#    i.set_visible(False)
termination_strength_corr_plot.set_xlabel("Clustering Coefficient")

#####
TF = summary['endactive'].astype('float')/summary['totactivity'].astype('float')
norm=None #NoNorm()


####
termination_plot = f.add_subplot(gs[0,8:15])
termination_plot.annotate("B", (0,0.95), xycoords=(termination_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)

segs = list(zip(x_vals, TF[:,i]) for i in range(n_nodes))
line_segments = LineCollection(segs, cmap=cmap, norm=norm, linewidths=linewidths)
line_segments.set_array(coloring_sequence)
termination_plot.add_collection(line_segments)

termination_plot.set_ylim(.1, .4)
termination_plot.set_xlim(x_vals[0], x_vals[-1])
#termination_plot.set_ylim(200, 400)
termination_plot.set_ylabel(r'$TF$, Termination Freq.')

for i in termination_plot.get_xticklabels():
    i.set_visible(False)

####
strength_plot = f.add_subplot(gs[1,8:15], sharex=termination_plot)
strength_plot.annotate("D", (0,0.95), xycoords=(strength_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)

segs = list(zip(x_vals, summary['outstrend'][:,i]) for i in range(n_nodes))

line_segments = LineCollection(segs, cmap=cmap, norm=norm, linewidths=linewidths)
line_segments.set_array(coloring_sequence)
strength_plot.add_collection(line_segments)

strength_plot.set_ylim(.4, 2)
strength_plot.set_xlim(x_vals[0], x_vals[-1])
strength_plot.set_ylabel(r'$OS$, Out Strength')
for i in strength_plot.get_xticklabels():
    i.set_visible(False)

####
termination_strength_corr_plot = f.add_subplot(gs[2,8:15], sharex=termination_plot)
termination_strength_corr_plot.annotate("F", (0,0.95), xycoords=(termination_strength_corr_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)


termination_strength_corr_plot.plot(x_vals, [pearsonr(summary['outstrend'][t,:], TF[t,:])[0] for t in range(len(x_vals))])

termination_strength_corr_plot.set_ylim(-1, 1)
termination_strength_corr_plot.set_ylabel(r'$R_{TF-OS}$')
#for i in strength_plot.get_xticklabels():
#    i.set_visible(False)
termination_strength_corr_plot.set_xlabel("Cascades (n x 10$^{6}$)")
termination_strength_corr_plot.set_xlim(x_vals[0], x_vals[-1])
termination_strength_corr_plot.set_xticks(arange(x_vals[-1]+1))

####
from mpl_toolkits.axes_grid.inset_locator import inset_axes
#termination_strength_corr_plot = inset_axes(termination_strength_corr_plot,
#                        width="30%", # width = 30% of parent_bbox
#                        height="50%",
#                        loc=1)

termination_strength_corr_plot = f.add_subplot(gs[3,8:15])
termination_strength_corr_plot.annotate("H", (0,0.95), xycoords=(termination_strength_corr_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)
                                            
from matplotlib import colors, cm
norm = colors.normalize(coloring_sequence.min(), coloring_sequence.max())
cmap = cm.spectral

for i in range(n_nodes):
    termination_strength_corr_plot.scatter(coloring_sequence[i], pearsonr(summary['outstrend'][:,i], fTN[:,i])[0], s=s, color=cmap(norm(coloring_sequence[i])))

termination_strength_corr_plot.set_ylim(-1, 1)
termination_strength_corr_plot.set_ylabel(r'$R_{TF-OS}$ by node')
#for i in strength_plot.get_xticklabels():
#    i.set_visible(False)
termination_strength_corr_plot.set_xlabel("Clustering Coefficient")

 
#savefig(output_directory+'DynamicDegreeStrength'+".pdf", bbox_inches='tight')
cax = f.add_subplot(gs[:,-1])
cb = f.colorbar(line_segments, cax=cax)
cb.set_label('Clustering Coefficient')

fn+=1
f.suptitle("Figure %i"%fn)
figures.append(f)

# <codecell>

dynamical_degree_dip = zeros(n_nodes)
dynamical_out_strength_dip = zeros(n_nodes)
for i in range(n_nodes):
    normdd = node_dynamical_degree[:,i]/float(degree_sequence[i])
    dynamical_degree_dip[i]= sum(normdd-normdd[0])
    normos = node_dynamical_out_strength[:,i]/node_out_strength[:,i]
    dynamical_out_strength_dip[i] = sum(normos-normos[0])
    
close('all')
f = figure(figsize=(8,8))
ax = f.add_subplot(121)
ax.scatter(clustering_sequence, dynamical_degree_dip)
print pearsonr(clustering_sequence, dynamical_degree_dip)
ax.scatter(clustering_sequence, dynamical_out_strength_dip, color='r')
print pearsonr(clustering_sequence, dynamical_out_strength_dip)
ax = f.add_subplot(122)
ax.scatter(degree_sequence, dynamical_degree_dip)
print pearsonr(degree_sequence, dynamical_degree_dip)
ax.scatter(degree_sequence, dynamical_out_strength_dip, color='r')
print pearsonr(degree_sequence, dynamical_out_strength_dip)

# <codecell>

close('all')
f = figure(figsize=(11,8))
gs = gridspec.GridSpec(3, 3)
gs.update(hspace=0.1, wspace=0.3)
alpha = .1

#####
replicate = 4
initial_graph = Graph.Weighted_Adjacency(nets['pnets'][0,replicate][0,0].toarray().tolist())
n_nodes = len(initial_graph.vs)
d = asarray(initial_graph.degree(mode=1))
cc = asarray(initial_graph.transitivity_local_undirected())
fTN = summary['endactive'].astype('float')/summary['totactivity'].astype('float')

#####
from scipy.stats import pearsonr
Rtn_d = zeros(500)
Rtn_cc = zeros(500)
for i in range(500):
    Rtn_d[i] = pearsonr(fTN[i], d)[0]
    Rtn_cc[i] = pearsonr(fTN[i], cc)[0]
    
termination_correlation_plot = f.add_subplot(gs[0:2])

termination_correlation_plot.plot(x_vals, Rtn_d, color='k', label=r'$R_{TN-d}$')
termination_correlation_plot.plot(x_vals, Rtn_cc, color='r', label=r'$R_{TN-CC}$')

handles, labels = termination_correlation_plot.get_legend_handles_labels()
leg = termination_correlation_plot.legend(handles, labels, loc=1)
    
termination_correlation_plot.set_ylim(-1,1)
termination_correlation_plot.set_ylabel(r'$R_{TN-d}$, $R_{TN-CC}$', color='k')
#termination_correlation_plot.set_xlabel("Cascades (n x 10$^{6}$)")
#plt.setp(termination_correlation_plot.get_xticklabels(), visible=False)
for i in termination_correlation_plot.get_xticklabels():
    i.set_visible(False)
    
termination_correlation_plot.annotate("A", (0,0.95), xycoords=(termination_correlation_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)

#####
cc_d_corr_plot = f.add_subplot(gs[2])
cc_d_corr_plot.scatter(d,cc, color='k', marker='.', s=2)
cc_d_corr_plot.set_ylim(0,1)
cc_d_corr_plot.set_ylabel('CC')
#cc_d_corr_plot.set_xlabel('Node Degree')
#plt.setp(cc_d_corr_plot.get_xticklabels(), visible=False)
for i in cc_d_corr_plot.get_xticklabels():
    i.set_visible(False)

# <codecell>

#####
d_target = 8
cc_range = (.45, .55)

n_replicates = len(nets['pnets'][0,:])

d = zeros([n_replicates, n_nodes])
cc = zeros([n_replicates, n_nodes])
fixed_d_nodes = []
fixed_cc_nodes = []

for replicate in range(n_replicates):
    initial_graph = Graph.Weighted_Adjacency(nets['pnets'][0,replicate][0,0].toarray().tolist())
    d[replicate,:] = asarray(initial_graph.degree(mode=1))
    cc[replicate,:] = asarray(initial_graph.transitivity_local_undirected())
    fixed_d_nodes.append(where(d[replicate,:]==d_target)[0])
    fixed_cc_nodes.append(where((cc_range[0]<=cc[replicate,:])&(cc[replicate,:]<cc_range[1]))[0])

#####
from scipy.stats import pearsonr
Rs_d_mean = zeros(500)
Rs_cc_mean = zeros(500)
Rs_d_sd = zeros(500)
Rs_cc_sd = zeros(500)
for i in range(500):
    Rs_d = []
    Rs_cc = []
    for replicate in range(n_replicates):
        this_graph = Graph.Weighted_Adjacency(nets['pnets'][0,replicate][0,i].toarray().tolist())
        Rs_d.append(pearsonr(this_graph.strength(weights="weight", mode=1), d[replicate])[0])
        Rs_cc.append(pearsonr(this_graph.strength(weights="weight", mode=1), cc[replicate])[0])
    
    Rs_d_mean[i] = mean(Rs_d)
    Rs_cc_mean[i] = mean(Rs_cc)
    Rs_d_sd[i] = std(Rs_d)
    Rs_cc_sd[i] = std(Rs_cc)
    
strength_correlation_plot = f.add_subplot(gs[3:5], sharex=termination_correlation_plot, sharey=termination_correlation_plot)

strength_correlation_plot.plot(x_vals, Rs_d_mean, color='k', label=r'$R_{s-d}$')
strength_correlation_plot.fill_between(x_vals, Rs_d_mean-Rs_d_sd, Rs_d_mean+Rs_d_sd, alpha=alpha, color='k')
strength_correlation_plot.plot(x_vals, Rs_cc_mean, color='r', label=r'$R_{s-CC}$')
strength_correlation_plot.fill_between(x_vals, Rs_cc_mean-Rs_cc_sd, Rs_cc_mean+Rs_cc_sd, alpha=alpha, color='r')

#Rtn_d = zeros(500)
#Rtn_cc = zeros(500)
#for i in range(500):
#    Rtn_d[i] = pearsonr(fTN[i][fixed_cc_nodes], d[fixed_cc_nodes])[0]
#    Rtn_cc[i] = pearsonr(fTN[i][fixed_d_nodes], cc[fixed_d_nodes])[0]
#strength_correlation_plot.plot(x_vals, Rtn_d, color='k')
#strength_correlation_plot.plot(x_vals, Rtn_cc, color='r')

handles, labels = strength_correlation_plot.get_legend_handles_labels()
leg = strength_correlation_plot.legend(handles, labels, loc=1)

strength_correlation_plot.set_ylim(-1,1)
strength_correlation_plot.set_ylabel(r'$R_{s-d}$, $R_{s-CC}$')
strength_correlation_plot.set_xlabel("Cascades (n x 10$^{6}$)")
strength_correlation_plot.set_xticks(range(11))

strength_correlation_plot.annotate("B", (0,0.95), xycoords=(strength_correlation_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)

#####
cc_d_corr_plot_sample = f.add_subplot(gs[5], sharex=cc_d_corr_plot, sharey=cc_d_corr_plot)

for replicate in range(n_replicates):
    cc_d_corr_plot_sample.scatter(d[replicate], cc[replicate], color='k', marker='.', s=2)

cc_d_corr_plot_sample.set_ylim(0,1)
cc_d_corr_plot_sample.set_ylabel('CC')
cc_d_corr_plot_sample.set_yticks([0, .5, 1.0])
#cc_d_corr_plot_sample.set_xlabel('d')
for i in cc_d_corr_plot_sample.get_xticklabels():
    i.set_visible(False)


#####
f.show()
figures.append(f)
#savefig(output_directory+'TerminationCorrelation'+".pdf", bbox_inches='tight')

# <codecell>

#####
if network_type in ('OHO', 'oho'):
    d_target = 8
    cc_range = (.45, .55)
elif network_type in ('WS', 'ws', 'wn', 'WN'):
    d_target = 9
    cc_range = (.30, .45)
n_replicates = len(nets['pnets'][0,:])

d = zeros([n_replicates, n_nodes])
cc = zeros([n_replicates, n_nodes])
fixed_d_nodes = []
fixed_cc_nodes = []

for replicate in range(n_replicates):
    initial_graph = Graph.Weighted_Adjacency(nets['pnets'][0,replicate][0,0].toarray().tolist())
    d[replicate,:] = asarray(initial_graph.degree(mode=1))
    cc[replicate,:] = asarray(initial_graph.transitivity_local_undirected())
    fixed_d_nodes.append(where(d[replicate,:]==d_target)[0])
    fixed_cc_nodes.append(where((cc_range[0]<=cc[replicate,:])&(cc[replicate,:]<cc_range[1]))[0])

#####
from scipy.stats import pearsonr
Rs_d_mean = zeros(500)
Rs_cc_mean = zeros(500)
Rs_d_sd = zeros(500)
Rs_cc_sd = zeros(500)
for i in range(500):
    Rs_d = []
    Rs_cc = []
    for replicate in range(n_replicates):
        this_graph = Graph.Weighted_Adjacency(nets['pnets'][0,replicate][0,i].toarray().tolist())
        Rs_d.append(pearsonr(this_graph.strength(fixed_cc_nodes[replicate], weights="weight", mode=1), d[replicate,fixed_cc_nodes[replicate]])[0])
        Rs_cc.append(pearsonr(this_graph.strength(fixed_d_nodes[replicate], weights="weight", mode=1), cc[replicate, fixed_d_nodes[replicate]])[0])
    
    Rs_d_mean[i] = mean(Rs_d)
    Rs_cc_mean[i] = mean(Rs_cc)
    Rs_d_sd[i] = std(Rs_d)
    Rs_cc_sd[i] = std(Rs_cc)
    
strength_correlation_plot = f.add_subplot(gs[6:8], sharex=termination_correlation_plot, sharey=termination_correlation_plot)

strength_correlation_plot.plot(x_vals, Rs_d_mean, color='k', label=r'$R_{s-d}$ fixed $CC$')
strength_correlation_plot.fill_between(x_vals, Rs_d_mean-Rs_d_sd, Rs_d_mean+Rs_d_sd, alpha=alpha, color='k')
strength_correlation_plot.plot(x_vals, Rs_cc_mean, color='r', label=r'$R_{s-CC}$ fixed $d$')
strength_correlation_plot.fill_between(x_vals, Rs_cc_mean-Rs_cc_sd, Rs_cc_mean+Rs_cc_sd, alpha=alpha, color='r')

#Rtn_d = zeros(500)
#Rtn_cc = zeros(500)
#for i in range(500):
#    Rtn_d[i] = pearsonr(fTN[i][fixed_cc_nodes], d[fixed_cc_nodes])[0]
#    Rtn_cc[i] = pearsonr(fTN[i][fixed_d_nodes], cc[fixed_d_nodes])[0]
#strength_correlation_plot.plot(x_vals, Rtn_d, color='k')
#strength_correlation_plot.plot(x_vals, Rtn_cc, color='r')

handles, labels = strength_correlation_plot.get_legend_handles_labels()
leg = strength_correlation_plot.legend(handles, labels, loc=1)

strength_correlation_plot.set_ylim(-1,1)
strength_correlation_plot.set_ylabel(r'$R_{s-d}$, $R_{s-CC}$')
strength_correlation_plot.set_xlabel("Cascades (n x 10$^{6}$)")
strength_correlation_plot.set_xticks(range(11))

strength_correlation_plot.annotate("C", (0,0.95), xycoords=(strength_correlation_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)

#####
cc_d_corr_plot_sample = f.add_subplot(gs[8], sharex=cc_d_corr_plot, sharey=cc_d_corr_plot)

for replicate in range(n_replicates):
    cc_d_corr_plot_sample.scatter(d[replicate, fixed_d_nodes[replicate]], cc[replicate, fixed_d_nodes[replicate]], color='r')
    cc_d_corr_plot_sample.scatter(d[replicate, fixed_cc_nodes[replicate]], cc[replicate, fixed_cc_nodes[replicate]], color='k')

cc_d_corr_plot_sample.set_ylim(0,1)
cc_d_corr_plot_sample.set_ylabel('CC')
cc_d_corr_plot_sample.set_xlabel('d')
cc_d_corr_plot_sample.set_yticks([0, .5, 1.0])

#####
fn+=1
f.suptitle("Figure %i"%fn)
f.show()
figures.append(f)
#savefig(output_directory+'TerminationCorrelation'+".pdf", bbox_inches='tight')

# <codecell>

#####
close('all')
f = figure(figsize=(11,8))
gs = gridspec.GridSpec(2, 16)
gs.update(hspace=0.05, wspace=0.25)
alpha = 1

clustering_sequence = initial_graph.transitivity_local_undirected()
from matplotlib.cm import get_cmap
cmap = get_cmap('spectral')
colors = [cmap(i) for i in clustering_sequence]
coloring_sequence = array(clustering_sequence)
from matplotlib.collections import LineCollection
from matplotlib.colors import NoNorm
norm=None #NoNorm()
linewidths=.15


dynamical_degree_plot = f.add_subplot(gs[0,:7])
dynamical_degree_plot.annotate("A", (0,0.95), xycoords=(dynamical_degree_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)

segs = list(zip(x_vals, node_dynamical_degree[:,i]) for i in range(n_nodes))
line_segments = LineCollection(segs, cmap=cmap, norm=norm, linewidths=linewidths)
line_segments.set_array(coloring_sequence)
dynamical_degree_plot.add_collection(line_segments)

#for i in range(n_nodes):
#    dynamical_degree_plot.plot(x_vals, node_dynamical_degree[:,i], alpha=alpha, color=cmap(clustering_sequence[i]))
dynamical_degree_plot.plot(x_vals, system_dynamical_degree, color='k', linewidth=3)

dynamical_degree_plot.set_ylabel('Dynamic Degree')
for i in dynamical_degree_plot.get_xticklabels():
    i.set_visible(False)
    
dynamical_degree_norm_plot = f.add_subplot(gs[0,8:15])
dynamical_degree_norm_plot.annotate("B", (0,0.95), xycoords=(dynamical_degree_norm_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)
segs = list(zip(x_vals, node_dynamical_degree[:,i]/float(degree_sequence[i])) for i in range(n_nodes))
line_segments = LineCollection(segs, cmap=cmap, norm=norm, linewidths=linewidths)
line_segments.set_array(coloring_sequence)
dynamical_degree_norm_plot.add_collection(line_segments)

#for i in range(n_nodes):
#    dynamical_degree_norm_plot.plot(x_vals, node_dynamical_degree[:,i]/float(degree_sequence[i]), alpha=alpha, color=cmap(clustering_sequence[i]))
dynamical_degree_norm_plot.plot(x_vals, system_dynamical_degree_norm, color='k', linewidth=3)

dynamical_degree_norm_plot.set_ylabel('Dynamic Degree, Normalized')
for i in dynamical_degree_norm_plot.get_xticklabels():
    i.set_visible(False)
    
    

dynamical_strength_plot = f.add_subplot(gs[1,:7], sharex=dynamical_degree_plot)
dynamical_strength_plot.annotate("C", (0,0.95), xycoords=(dynamical_strength_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)
segs = list(zip(x_vals, node_dynamical_out_strength[:,i]) for i in range(n_nodes))
line_segments = LineCollection(segs, cmap=cmap, norm=norm, linewidths=linewidths)
line_segments.set_array(coloring_sequence)
dynamical_strength_plot.add_collection(line_segments)

#for i in range(n_nodes):
#    dynamical_strength_plot.plot(x_vals, node_dynamical_out_strength[:,i], alpha=alpha, color=cmap(clustering_sequence[i]))
dynamical_strength_plot.plot(x_vals, system_dynamical_out_strength, color='k', linewidth=3)
dynamical_strength_plot.set_ylabel('Dynamic Out Strength')
dynamical_strength_plot.set_xlabel("Cascades (n x 10$^{6}$)")
dynamical_strength_plot.set_xlim(x_vals[0], x_vals[-1])
dynamical_strength_plot.set_xticks(arange(x_vals[-1]+1))


dynamical_strength_norm_plot = f.add_subplot(gs[1,8:15], sharex=dynamical_degree_norm_plot)
dynamical_strength_norm_plot.annotate("D", (0,0.95), xycoords=(dynamical_strength_norm_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)
segs = list(zip(x_vals, node_dynamical_out_strength[:,i]/node_out_strength[:,i]) for i in range(n_nodes))
line_segments = LineCollection(segs, cmap=cmap, norm=norm, linewidths=linewidths)
line_segments.set_array(coloring_sequence)
dynamical_strength_norm_plot.add_collection(line_segments)
#for i in range(n_nodes):
#    dynamical_strength_norm_plot.plot(x_vals, node_dynamical_out_strength[:,i]/node_out_strength[:,i], alpha=alpha, color=cmap(clustering_sequence[i]))
dynamical_strength_norm_plot.plot(x_vals, system_dynamical_out_strength_norm, color='k', linewidth=3)
dynamical_strength_norm_plot.set_ylabel('Dynamic Out Strength, Normalized')
dynamical_strength_norm_plot.set_xlabel("Cascades (n x 10$^{6}$)")
dynamical_strength_norm_plot.set_xlim(x_vals[0], x_vals[-1])
dynamical_strength_norm_plot.set_xticks(arange(x_vals[-1]+1))

#savefig(output_directory+'DynamicDegreeStrength'+".pdf", bbox_inches='tight')
cax = f.add_subplot(gs[:,15])
cb = f.colorbar(line_segments, cax=cax)
cb.set_label('Clustering Coefficient')

fn+=1
f.suptitle("Figure %i"%fn)
figures.append(f)

# <codecell>

close('all')
f = figure(figsize=(11,8))
gs = gridspec.GridSpec(2, 16)
gs.update(hspace=0.2, wspace=0.6)
coloring_sequence = array(clustering_sequence)

from matplotlib import colors, cm
norm = colors.normalize(coloring_sequence.min(), coloring_sequence.max())
cmap = cm.spectral


####
s=.05
degree_strength_plot = f.add_subplot(gs[0,:7])
degree_strength_plot.annotate("A", (0,0.95), xycoords=(dynamical_strength_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)

for i in range(n_nodes):
    degree_strength_plot.scatter(node_dynamical_degree[:,i], node_dynamical_out_strength[:,i], s=s, color=cmap(norm(coloring_sequence[i])))


degree_strength_plot.set_ylabel('Dynamic Out Strength')
degree_strength_plot.set_xlabel("Dynamic Degree")
degree_strength_plot.annotate("A", (0,0.95), xycoords=(dynamical_strength_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)

####
degree_strength_norm_plot = f.add_subplot(gs[0,8:15])
degree_strength_norm_plot.annotate("B", (0,0.95), xycoords=(degree_strength_norm_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)
for i in range(n_nodes):
    degree_strength_norm_plot.scatter(node_dynamical_degree[:,i]/float(degree_sequence[i]), node_dynamical_out_strength[:,i]/node_out_strength[:,i], s=s, color=cmap(norm(coloring_sequence[i])))

    degree_strength_norm_plot.set_ylabel('Dynamic Out Strength, Normalized')
degree_strength_norm_plot.set_xlabel("Dynamic Degree, Normalized")

########
####
from scipy.stats import pearsonr
s=10
degree_strength_plot = f.add_subplot(gs[1,:7])
degree_strength_plot.annotate("C", (0,0.95), xycoords=(dynamical_strength_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)

for i in range(n_nodes):
    degree_strength_plot.scatter(coloring_sequence[i], pearsonr(node_dynamical_degree[:,i], node_dynamical_out_strength[:,i])[0], s=s, color=cmap(norm(coloring_sequence[i])))

degree_strength_plot.set_xlim(0,1)
degree_strength_plot.set_ylabel('Dynamic Degree/Out Strength Correlation')
degree_strength_plot.set_xlabel("Clustering Coefficient")


####
degree_strength_norm_plot = f.add_subplot(gs[1,8:15])
degree_strength_norm_plot.annotate("D", (0,0.95), xycoords=(degree_strength_norm_plot.get_yaxis().get_label(), "axes fraction"), fontsize=14)

for i in range(n_nodes):
    degree_strength_norm_plot.scatter(coloring_sequence[i], pearsonr(node_dynamical_degree[:,i]/float(degree_sequence[i]), node_dynamical_out_strength[:,i]/node_out_strength[:,i])[0], s=s, color=cmap(norm(coloring_sequence[i])))

degree_strength_norm_plot.set_xlim(0,1)
degree_strength_norm_plot.set_ylabel('Dynamic Degree/Out Strength\nCorrelation, Normalized')
degree_strength_norm_plot.set_xlabel("Clustering Coefficient")


####
from matplotlib import colorbar
cax = f.add_subplot(gs[:,15])
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
cb.set_label('Clustering Coefficient')

fn+=1
f.suptitle("Figure %i"%fn)
figures.append(f)

# <codecell>

from matplotlib.backends.backend_pdf import PdfPages
plots = PdfPages('/data/alstottjd/Cascade/%s_Figures.pdf'%network_type)
for i in figures:
    plots.savefig(i)
plots.close()

# <codecell>

#Significance testing

# <codecell>


