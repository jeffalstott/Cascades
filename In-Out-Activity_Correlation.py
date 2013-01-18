# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from scipy.io import loadmat
from igraph import Graph

# <codecell>

data_directory = "/data/alstottjd/Sini/"

# <codecell>

#mat = loadmat(data_dir+'test_NL_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_ws_-0.0333333-4_1_5_ncorr.mat')
#data_name='WSlast'
mat = loadmat(data_directory+'test_NL_mlast_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_oho_-6_1_5_ncorr.mat')
data_name='OHOlast'
#mat = loadmat(data_dir+'2,5_NL_mo-3_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_oho_-6_1_5_ncorr.mat')
#data_name = "2,5OHOF3"
#mat = loadmat(data_dir+'2,5_NL_mall_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_ws_-0.0333333-4_1_5_ncorr.mat')
#data_name = "2,5WSall"

# <codecell>

mat = loadmat(data_directory+'NL_summary_mo-2_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_ws_-0.0333333-4_1_5_ncorr.mat')
f = figure()
ax = f.add_subplot(111)

from scipy.stats import pearsonr
n_runs = 500
act_in = ones(n_runs-1)
act_out = ones(n_runs-1)
in_act = ones(n_runs-1)
out_act = ones(n_runs-1)
out_in = ones(n_runs-1)
in_out = ones(n_runs-1)
in_in = ones(n_runs-1)
out_out = ones(n_runs-1)
for i in range(1,n_runs):
    in_act[i-1]=pearsonr(mat['totactivity'][i:].flatten(),mat['inpstrend'][:-i].flatten())[0]
    out_act[i-1]=pearsonr(mat['totactivity'][i:].flatten(),mat['outstrend'][:-i].flatten())[0]
    act_in[i-1]=pearsonr(mat['totactivity'][:-i].flatten(),mat['inpstrend'][i:].flatten())[0]
    act_out[i-1]=pearsonr(mat['totactivity'][:-i].flatten(),mat['outstrend'][i:].flatten())[0]
    in_out[i-1]=pearsonr(mat['outstrend'][i:].flatten(),mat['inpstrend'][:-i].flatten())[0]
    out_in[i-1]=pearsonr(mat['inpstrend'][i:].flatten(),mat['outstrend'][:-i].flatten())[0]
    in_in[i-1]=pearsonr(mat['inpstrend'][i:].flatten(),mat['inpstrend'][:-i].flatten())[0]
    out_out[i-1]=pearsonr(mat['outstrend'][i:].flatten(),mat['outstrend'][:-i].flatten())[0]
ax.plot(act_in, label="Activity Predicting In")
ax.plot(act_out, label="Activity Predicting Out")
ax.plot(act_in, label="In Predicting Activity")
ax.plot(act_out, label="Out Predicting Activity")
ax.plot(in_out, label="In Predicting Out")
ax.plot(out_in, label="Out predicting In")
ax.plot(in_in, label="In Strength Autocorrelation")
ax.plot(out_out, label="Out Strength Autocorrelation")
ax.set_xlabel("Lag (Cascades)")
ax.set_ylim(ax.get_ylim()[0],1)
ax.legend(bbox_to_anchor=(1.1, 1.05),loc=2)
print "Collapsing across cascades, impact of present structure on future activity"
"or structure decreases the further in the future you look."

# <codecell>

f = figure()
ax = f.add_subplot(111)

from scipy.stats import pearsonr
n_runs = 400
act_in = zeros(n_runs)
act_out = zeros(n_runs)
out_in = zeros(n_runs)
for j in [0,20]:#range(0,50,10):
    for i in range(n_runs):
        act_in[i]=pearsonr(mat['totactivity'][i],mat['inpstrend'][i])[0]
        act_out[i]=pearsonr(mat['totactivity'][i],mat['outstrend'][i])[0]
        out_in[i]=pearsonr(mat['inpstrend'][i],mat['outstrend'][i+j])[0]
    #ax.plot(act_in)
    #ax.plot(act_out)
    ax.plot(out_in, label="Lag: "+str(j))
ax.legend(loc=3)

ax.set_ylim(0,1)

# <codecell>

from scipy.signal import correlate
f = figure()
ax = f.add_subplot(111)
for i in range(len(mat['inpstrend'][0])):
    ax.plot(correlate(mat['inpstrend'][:,i], mat['outstrend'][:,i])[400:600], label=str(i))
ylabel("Cross-Correlation")
title("In Strength and Out Strength Cross Correlation")

# <codecell>


