# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from scipy.io import loadmat
data_directory = "/data/alstottjd/Sini/"
import Timelines

# <codecell>

mat = loadmat(data_directory+'2,5_NL_mall_nruns70_pinc0.001_addinc0_mnwt7_mxwt5_cf99_R1000000_netprms:_60_oho_-6_1_5_ncorr.mat')
data_name='OHOlast'

# <codecell>

step_size=5

# <codecell>

tl = Timelines.Timeline(mat['pnets'][0,0][0,::step_size])
tlr = Timelines.Timeline(mat['pnets_spr_out'][0,::step_size,0])

# <codecell>


