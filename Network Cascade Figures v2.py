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

data_directory = '/data/alstottjd/Cascade/Timelines/'
output_directory = '/data/alstottjd/Cascade/Figures/'
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

def plot_intrc(measure, normalized=True, offset=1):
            
    f = figure()
    ax = f.add_subplot(111)
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    ax.set_xlabel("Cascades (n x 10$^{6}$)")
    ax.set_ylabel("Integrated Rich Club Index, "+measure)
    
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
        if normalized:
            y = sum(data[i]['data']["T_in"].raw_normalized("rc_"+measure+"_in")-offset, axis=-1)
        else:
            y = sum(data[i]['data']["T_in"].raw_data("rc_"+measure+"_in")-offset, axis=-1)
        y_vals = mean(y, axis=0)
        error = std(y, axis=0)
        ax1.plot(x_vals, y_vals, 
            label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])
        ax1.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])
        
        if normalized:
            y = sum(data[i]['data']["T_out"].raw_normalized("rc_"+measure+"_out")-offset, axis=-1)
        else:
            y = sum(data[i]['data']["T_out"].raw_data("rc_"+measure+"_out")-offset, axis=-1)
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
    
    savefig(output_directory+'RichClub'+measure+".pdf", bbox_inches='tight')
    #title("Rich Club with Learning")

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
savefig(output_directory+'Compression.pdf', bbox_inches='tight')
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

savefig(output_directory+'Modularity.pdf', bbox_inches='tight')
#title("Modularity with Learning")

# <codecell>

plot_intrc('iNl')

# <codecell>

plot_intrc('iNpwml')

# <codecell>

plot_intrc('q_infomap', normalized=False, offset=0)

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

savefig(output_directory+'PathLengthClustering.pdf', bbox_inches='tight')

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

savefig(output_directory+'BetweenessDerivative.pdf', bbox_inches='tight')
#title("Network Compression with Learning")

# <codecell>

f = figure()
ax = f.add_subplot(111)

for i in datanames:
    y_vals, error = data[i]['data']["T_out"].data("betweeness_change_from_first")
    ax.plot(x_vals, y_vals, 
        label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])
    ax.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])

    y_vals, error = data[i]['data']["T_out"].control("betweeness_change_from_first")
    ax.plot(x_vals, y_vals, color=data[i]['color'], linestyle=data[i]['line']+'.')
#        label=data[i]['name'], 
    #ax.fill_between(x_vals[1:], y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])

ax.set_ylim(ax.get_ylim()[0], 1.01)
handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles, labels, loc=3, ncol=4, mode='expand')
leg.get_frame().set_alpha(0) 

ax.set_xlabel("Cascades (n x 10$^{6}$)")
ax.set_ylabel("Correlation in Betweeness Centrality Rank Order\n From Original Network")

f.tight_layout()

savefig(output_directory+'BetweenessCorrtoFirst.pdf', bbox_inches='tight')
#title("Network Compression with Learning")

# <codecell>

f = figure()
ax = f.add_subplot(111)

for i in datanames:
    y_vals, error = data[i]['data']["T_out"].data("mincut_value")
    ax.plot(x_vals, y_vals, 
        label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])
    ax.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])

    y_vals, error = data[i]['data']["T_out"].control("mincut_value")
    ax.plot(x_vals, y_vals, color=data[i]['color'], linestyle=data[i]['line']+'.')
#        label=data[i]['name'], 
    #ax.fill_between(x_vals[1:], y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])

handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles, labels, loc=3, ncol=4, mode='expand')
leg.get_frame().set_alpha(0) 

ax.set_xlabel("Cascades (n x 10$^{6}$)")
ax.set_ylabel("Mincut Value")

f.tight_layout()

savefig(output_directory+'Mincut.pdf', bbox_inches='tight')
#title("Network Compression with Learning")

# <codecell>

f = figure()
ax = f.add_subplot(111)

for i in datanames:
    y_vals, error = data[i]['data']["T_out"].data("diameter")
    ax.plot(x_vals, y_vals, 
        label=data[i]['name'], color=data[i]['color'], linestyle=data[i]['line'])
    ax.fill_between(x_vals, y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])

    y_vals, error = data[i]['data']["T_out"].control("diameter")
    ax.plot(x_vals, y_vals, color=data[i]['color'], linestyle=data[i]['line']+'.')
#        label=data[i]['name'], 
    #ax.fill_between(x_vals[1:], y_vals-error, y_vals+error, alpha=alpha, color=data[i]['color'])

handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles, labels, loc=3, ncol=4, mode='expand')
leg.get_frame().set_alpha(0) 

ax.set_xlabel("Cascades (n x 10$^{6}$)")
ax.set_ylabel("Diameter")

f.tight_layout()

savefig(output_directory+'Diameter.pdf', bbox_inches='tight')
#title("Network Compression with Learning")

# <codecell>

increment = 1.0/10.0
x = arange(.1, 1, .1)
w= 3.375
h = w/1.6180
#f = figure(figsize=(w,h))
f = figure()

n_cascades = 10.0
cascades_per_run = n_cascades/n_runs
step_size = 5
n_samples = ceil(n_runs/step_size)

c1 = 0
c2 = 2
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

rc_in = data['mlast']['data']["T_in"].raw_normalized("rc_itl_in")
rc_out = data['mlast']['data']["T_out"].raw_normalized("rc_itl_out")


alpha = .1

xlims = (.1,.9)



ax1 = f.add_subplot(131)

y_vals = mean(rc_in[:,s1,:], axis=0)
error = std(rc_in[:,s1, :], axis=0)
ax1.plot(x, y_vals, label='In Strength', color='b')
ax1.fill_between(x, y_vals-error, y_vals+error, alpha=alpha, color='b')

y_vals = mean(rc_out[:,s1,:], axis=0)
error = std(rc_out[:,s1, :], axis=0)
ax1.plot(x, y_vals, label='Out Strength', color='g')
ax1.fill_between(x, y_vals-error, y_vals+error, alpha=alpha, color='g')

ax1.set_xlim(xlims)

ax1.plot(xlims, (1,1), 'k--')
plt.setp(ax1.get_xticklabels(), visible=False)
ylabel("Normalized Rich Club Coefficient")
text(.5, .9, '0 Cascades', transform = ax1.transAxes, horizontalalignment='center', fontsize=10)
#xlabel('Strength Decile')
plt.xticks(x[::2], x[::2])
#handles, labels = ax1.get_legend_handles_labels()
#ax1.legend(handles, labels, loc=6)

ax2 = f.add_subplot(132, sharey=ax1)
y_vals = mean(rc_in[:,s2,:], axis=0)
error = std(rc_in[:,s2, :], axis=0)
ax2.plot(x, y_vals, label='In Strength', color='b')
ax2.fill_between(x, y_vals-error, y_vals+error, alpha=alpha, color='b')

y_vals = mean(rc_out[:,s2,:], axis=0)
error = std(rc_out[:,s2, :], axis=0)
ax2.plot(x, y_vals, label='Out Strength', color='g')
ax2.fill_between(x, y_vals-error, y_vals+error, alpha=alpha, color='g')

ax2.set_xlim(xlims)
ax2.plot(xlims, (1,1), 'k--')
plt.setp(ax2.get_yticklabels(), visible=False)
xlabel('Strength Percentile')
plt.xticks(x[::2], (x[::2]*100).astype(int))
text(.5, .9, '$%.1f*10^{6}$ Cascades'%c2, transform = ax2.transAxes, horizontalalignment='center', fontsize=10)
#handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles, labels, loc=1)

ax3 = f.add_subplot(133, sharey=ax1)
y_vals = mean(rc_in[:,s3,:], axis=0)
error = std(rc_in[:,s3, :], axis=0)
ax3.plot(x, y_vals, label='In Strength', color='b')
ax3.fill_between(x, y_vals-error, y_vals+error, alpha=alpha, color='b')

y_vals = mean(rc_out[:,s3,:], axis=0)
error = std(rc_out[:,s3, :], axis=0)
ax3.plot(x, y_vals, label='Out Strength', color='g')
ax3.fill_between(x, y_vals-error, y_vals+error, alpha=alpha, color='g')

ax3.set_xlim(xlims)
ax3.plot(xlims, (1,1), 'k--')
plt.setp(ax3.get_yticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
text(.5, .9, '$%.1f*10^{6}$ Cascades'%c3, transform = ax3.transAxes, horizontalalignment='center', fontsize=10)
#xlabel('Strength Decile')
plt.xticks(x[::2], x[::2])

handles, labels = ax3.get_legend_handles_labels()
ax3.legend(handles, labels, loc=5, fontsize=8)

ax3.set_ylim(0,6)

savefig(output_directory+'RichClubSamples.pdf', bbox_inches='tight')
#suptitle('Rich Club Growth and Death')

# <codecell>


# <codecell>


