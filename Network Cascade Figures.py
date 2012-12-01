# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from matplotlib import rc
rc('text', usetex=False)
rc('font', family='serif')

# <codecell>

from scipy.signal import correlate
import cPickle as pickle
from os import listdir

# <codecell>

data_directory = '/data/alstottjd/Sini/'
files = listdir(data_directory)

# <codecell>

n_rewires = 10
n_iters = 5
network_type = 'oho'

# <codecell>

datanames = [#'mlast',
 #'mL1',
 #'mL2',
 #'mL3',
 #'mL4',
 #'mL5',
 'mo-1',
 'mall',
 'mo-2',
 'mlast'
 #'mo-2',
 #'mo-3',
 #'mo-4',
 #'mo-5',
]

# <codecell>

data = {
#'mlast': {'name': 'Last', 'color': 'k', 'line': '-'},
'mlast': {'name': 'Last', 'color': 'k', 'line': '-'},
'mL1': {'name': 'L1', 'color': 'r', 'line': '-'},
'mL2': {'name': 'L2', 'color': 'g', 'line': '-'},
'mL3': {'name': 'L3', 'color': 'b', 'line': '-'},
'mL4': {'name': 'L4', 'color': 'c', 'line': '-'},
'mL5': {'name': 'L5', 'color': 'm', 'line': '-'},
#'mall': {'name': 'All', 'color': 'k', 'line': '--'},
'mall': {'name': 'All', 'color': 'b', 'line': '-'},
#'mo-1': {'name': 'F1', 'color': 'r', 'line': '--'},
'mo-1': {'name': 'Random', 'color': 'r', 'line': '-'},
#'mo-2': {'name': 'F2', 'color': 'g', 'line': '--'},
'mo-2': {'name': 'Second', 'color': 'g', 'line': '-'},
'mo-3' : {'name': 'F3', 'color': 'b', 'line': '--'},
'mo-4': {'name': 'F4', 'color': 'c', 'line': '--'},
'mo-5': {'name': 'F5', 'color': 'm', 'line': '--'}
}

# <codecell>

for filename in files:
    if filename.startswith('tl_'+str(n_rewires)+','+str(n_iters)):
        dataname = filename.split('_')[3]
        if network_type in filename and dataname in datanames:
            print dataname
            data[dataname]['data'] = pickle.load( open( data_directory+filename, "rb" ))

# <codecell>

n_runs = 500
step_size = 5
close('all')
x_vals = arange(0, n_runs+1, step_size)
x_vals = 10.0*x_vals/n_runs
alpha = .1

# <codecell>

f = figure()
ax = f.add_subplot(111)

for i in datanames:
    y_vals, error = data[i]['data']["T_out"].data("codelength")
    ax.plot(x_vals, y_vals, 
        label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])
    fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])

handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles, labels, loc=4)
leg.get_frame().set_alpha(0) 


xlabel("Cascades (n x 10$^{6}$)")
ylabel("Bits to Represent Network")

f.tight_layout()
savefig(data_directory+'Compression.pdf', bbox_inches='tight')
#title("Network Compression with Learning")

# <codecell>

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#Modularity Metrics
f = figure()
ax = f.add_subplot(111)

for i in datanames:
    y_vals, error = data[i]['data']["T_out"].data("q_infomap")
    ax.plot(x_vals, y_vals, 
        label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])
    fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])

ax.set_ylim(0,.701)
ax.set_xlabel("Cascades (n x 10$^{6}$)")
ax.set_ylabel("Modularity Index")



handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles, labels, loc=1, ncol=2)
leg.get_frame().set_alpha(0) 

inset_ax = axes([0.62,0.22,0.3,0.3])
for i in datanames:
    y_vals, error = data[i]['data']["T_out"].data("n_infomap")
    inset_ax.plot(x_vals, y_vals, 
        label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])
    inset_ax.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])
ylabel("Number of Modules", fontsize=10)
#inset_ax.yaxis.set_major_locator(MultipleLocator(.1))

f.tight_layout()

savefig(data_directory+'Modularity.pdf', bbox_inches='tight')
#title("Modularity with Learning")

# <codecell>

f = figure()
ax = f.add_subplot(111)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax.set_xlabel("Cascades (n x 10$^{6}$)")
ax.set_ylabel("Integrated Rich Club Index")

ax1 = f.add_subplot(211)
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
ax1.yaxis.set_major_locator(MultipleLocator(10))
ax2 = f.add_subplot(212, sharey=ax1, sharex=ax1)
plt.setp(ax1.get_xticklabels(), visible=False)

ax1.text(.5, .8, 'In Strength',
 transform = ax1.transAxes, horizontalalignment='center', fontsize=10)
ax2.text(.5, .8, 'Out Strength',
 transform = ax2.transAxes, horizontalalignment='center', fontsize=10)

ax1.set_xlim(x_vals[0], x_vals[-1])
#ax1.set_ylim(-5, 40)
ax1.plot(ax1.get_xlim(), (0,0), 'k')
ax2.plot(ax2.get_xlim(), (0,0), 'k')

for i in datanames:
    y = sum(data[i]['data']["T_in"].raw_normalized("rc_in")-1, axis=-1)
    y_vals = mean(y, axis=0)
    error = std(y, axis=0)
    ax1.plot(x_vals, y_vals, 
        label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])
    ax1.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])
    
    y = sum(data[i]['data']["T_out"].raw_normalized("rc_out")-1, axis=-1)
    y_vals = mean(y, axis=0)
    error = std(y, axis=0)
    ax2.plot(x_vals, y_vals, 
        label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])
    ax2.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])

ax1.plot(ax1.get_xlim(), (0,0), 'k')
ax2.plot(ax2.get_xlim(), (0,0), 'k')



handles, labels = ax2.get_legend_handles_labels()
leg = ax2.legend(handles, labels, loc=1)
leg.get_frame().set_alpha(0) 

f.tight_layout()

savefig(data_directory+'RichClubInt.pdf', bbox_inches='tight')
#title("Rich Club with Learning")

# <codecell>

f = figure()
ax = f.add_subplot(111)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax.set_xlabel("Cascades (n x 10$^{6}$)")
ax.set_ylabel("Integrated Rich Club Index")

ax1 = f.add_subplot(121)
ax2 = f.add_subplot(122, sharey=ax1, sharex=ax1)
plt.setp(ax2.get_yticklabels(), visible=False)

#ax1.text(.5, .8, 'In Strength',
# transform = ax1.transAxes, horizontalalignment='center', fontsize=10)
#ax2.text(.5, .8, 'Out Strength',
# transform = ax2.transAxes, horizontalalignment='center', fontsize=10)
ax1.set_title("In Strength")
ax2.set_title("Out Strength")

ax1.set_xlim(x_vals[0], x_vals[-1])
ax1.plot(ax1.get_xlim(), (0,0), 'k')
ax2.plot(ax2.get_xlim(), (0,0), 'k')

for i in datanames:
    y = sum(data[i]['data']["T_in"].raw_normalized("rc_in")-1, axis=-1)
    y_vals = mean(y, axis=0)
    error = std(y, axis=0)
    ax1.plot(x_vals, y_vals, 
        label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])
    ax1.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])
    
    y = sum(data[i]['data']["T_out"].raw_normalized("rc_out")-1, axis=-1)
    y_vals = mean(y, axis=0)
    error = std(y, axis=0)
    ax2.plot(x_vals, y_vals, 
        label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])
    ax2.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])

ax1.plot(ax1.get_xlim(), (0,0), 'k')
ax2.plot(ax2.get_xlim(), (0,0), 'k')

handles, labels = ax2.get_legend_handles_labels()
leg = ax2.legend(handles, labels, loc=1)
leg.get_frame().set_alpha(0) 

f.tight_layout()

savefig(data_directory+'RichClubInt.pdf', bbox_inches='tight')
#title("Rich Club with Learning")

# <codecell>

f = figure()
ax = f.add_subplot(111)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax.set_xlabel("Cascades (n x 10$^{6}$)")

ax1 = f.add_subplot(131)
ax1.set_ylabel('Mean Path Length')
ax2 = f.add_subplot(132, sharey=ax1, sharex=ax1)
ax2.set_ylabel('Mean Clustering')
plt.setp(ax2.get_yticklabels(), visible=False)
ax3 = f.add_subplot(133, sharey=ax1, sharex=ax1)
ax3.set_ylabel('Small World Index')
plt.setp(ax3.get_yticklabels(), visible=False)

ax1.set_xlim(x_vals[0], x_vals[-1])#len(x_vals)*2/3])

for i in datanames:
    y_vals, error = data[i]['data']["T_out"].normalized("mean_path_length")
    ax1.plot(x_vals, y_vals, 
        label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])
    ax1.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])
    
    y_vals, error = data[i]['data']["T_out"].normalized("mean_clustering")
    ax2.plot(x_vals, y_vals, 
        label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])
    ax2.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])
    
    y_vals = data[i]['data']["T_out"].normalized("mean_clustering")[0]/data[i]['data']["T_out"].normalized("mean_path_length")[0]
    ax3.plot(x_vals, y_vals, 
        label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])

handles, labels = ax3.get_legend_handles_labels()
leg = ax3.legend(handles, labels, loc=1)
leg.get_frame().set_alpha(0) 

f.tight_layout()

savefig(data_directory+'PathLengthClustering.pdf', bbox_inches='tight')

# <codecell>

f = figure()
ax = f.add_subplot(111)

for i in datanames:
    y_vals, error = data[i]['data']["T_out"].data("betweeness_change_kendall")
    ax.plot(x_vals, y_vals, 
        label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])
    ax.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])

    y_vals, error = data[i]['data']["T_out"].control("betweeness_change_kendall")
    ax.plot(x_vals, y_vals, color=data[i]['color'], linestyle=data[i]['line']+'.')
#        label=data[i]['name'], 
    #ax.fill_between(x_vals[1:], y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])

ax.set_ylim(ax.get_ylim()[0], 1.01)
handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles, labels, loc=3, ncol=4, mode='expand')
leg.get_frame().set_alpha(0) 

ax.set_xlabel("Cascades (n x 10$^{6}$)")
ax.set_ylabel("Correlation in Betweeness Centrality Rank Order")

f.tight_layout()

savefig(data_directory+'Rerouting.pdf', bbox_inches='tight')
#title("Network Compression with Learning")

# <codecell>


