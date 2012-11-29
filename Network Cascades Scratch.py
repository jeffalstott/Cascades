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


