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

n_rewires = 2
n_iters = 5

# <codecell>

data = {}
for filename in files:
    if filename.startswith('tl_'+str(n_rewires)+','+str(n_iters)):
        dataname = filename.split('_')[3]
        data[dataname] = pickle.load( open( data_directory+filename, "rb" ))

# <codecell>

n_runs = 500
step_size = 5
close('all')
x_vals = arange(0, n_runs, step_size)
x_vals = 10.0*x_vals/n_runs

# <codecell>

figure()

for i in data:
    y_vals, error = data[i]["T_out"].data("codelength")
    plot(x_vals, y_vals, label=i)
#    fill_between(x_vals, y_vals-error, y_vals+error, alpha=.5)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=1)

xlabel("Cascades (n x 10$^{6}$)")
ylabel("Bits to Represent Network")

savefig(data_dir+'Compression.pdf', bbox_inches='tight')
#title("Network Compression with Learning")

# <codecell>

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#Modularity Metrics
figure()
for i in data:
    y_vals, error = data[i]["T_out"].data("q_infomap")
    plot(x_vals, y_vals, label=i)
#    fill_between(x_vals, y_vals-error, y_vals+error, alpha=.5)

xlabel("Cascades (n x 10$^{6}$)")
ylabel("Number of Modules")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=2)

inset_ax = axes([0.58,0.58,0.3,0.3])
for i in data:
    y_vals, error = data[i]["T_out"].data("n_infomap")
    inset_axplot(x_vals, y_vals, label=i)
    #inset_ax.fill_between(x_vals, y_vals-error, y_vals+error, alpha=.5)
ylabel("Modularity Index", fontsize=10)
inset_ax.yaxis.set_major_locator(MultipleLocator(.1))


savefig(data_dir+'Modularity.pdf', bbox_inches='tight')
#title("Modularity with Learning")

