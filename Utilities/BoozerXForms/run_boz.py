#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 17:50:56 2022

This file initialises a booz_xForms instance
and runs it before writing all outputs to 
an hdf5 file.

@author: S.Guinchard
"""

############
# PACKAGES #
############

from init_from_spec import init_from_spec
import matplotlib.pyplot as plt
import booz_xform as bx
import numpy as np
import h5py as h5

##############
# INIT + RUN #
##############

file     = 'Name'
filename  = file + '.sp.h5'
bozout   = file +'.boz.h5'
b        = init_from_spec('path_to_file'+filename)
b.verbose = 2
print("Selected surfaces:", b.compute_surfs)
b.run()

##########################
# CREATE OUTPUT DATASETS #
##########################

f   = h5.File('path_to_Boozer_out_file'+ bozout, 'w')
print(f.filename)
grp = f.create_group('Booz_xForms')
grp.create_group('Inputs')
inputs = grp['Inputs']
inputs['verbose'] = b.verbose
inputs['asym']   = b.asym
inputs['nfp']    = b.nfp
inputs['ntor'] = b.ntor
inputs['mpol'] = b.mpol
inputs['mnmax'] = b.mnmax
inputs['mnmax_nyq'] = b.mnmax_nyq
inputs['mpol_nyq'] = b.mpol_nyq
inputs['ntor_nyq'] = b.ntor_nyq
inputs['xn'] = b.xn
inputs['xm'] = b.xm
inputs['xn_nyq'] = b.xn_nyq
inputs['xm_nyq'] = b.xm_nyq
inputs['ns_in'] = b.ns_in
inputs['iota'] = b.iota
inputs['rmnc'] = b.rmnc
inputs['rmns'] = b.rmns
inputs['zmnc'] = b.zmnc
inputs['zmns'] = b.zmns
inputs['bmn'] = b.bmnc
inputs['bsubumnc'] = b.bsubumnc
inputs['bsubvmnc'] = b.bsubvmnc
inputs['compute_surfs'] = b.compute_surfs
inputs['aspect'] = b.aspect
inputs['toroidal_flux'] = b.toroidal_flux
inputs['mboz'] = b.mboz
inputs['nboz'] = b.nboz
inputs['lmns'] = b.lmns
inputs['lmnc'] = b.lmnc

grp.create_group('Outputs')
outputs = grp['Outputs']
outputs.create_dataset('ns_b', data = b.ns_b)
outputs.create_dataset('s_b', data = b.s_b)
outputs.create_dataset('mnboz', data = b.mnboz)
outputs.create_dataset('xm_b', data = b.xm_b)
outputs.create_dataset('xn_b', data = b.xn_b)
outputs.create_dataset('bmnc_b', data = b.bmnc_b)
outputs.create_dataset('bmns_b', data = b.bmns_b)
outputs.create_dataset('gmnc_b', data = b.gmnc_b)
outputs.create_dataset('gmns_b', data = b.gmns_b)
outputs.create_dataset('rmnc_b', data = b.rmnc_b)
outputs.create_dataset('zmnc_b', data = b.zmnc_b)
outputs.create_dataset('rmns_b', data = b.rmns_b)
outputs.create_dataset('zmns_b', data = b.zmns_b)
outputs.create_dataset('numnc_b', data = b.numnc_b)
outputs.create_dataset('numns_b', data = b.numns_b)
outputs.create_dataset('boozer_g', data = b.Boozer_G)
outputs.create_dataset('boozer_g_all', data = b.Boozer_G_all)
outputs.create_dataset('boozer_i', data = b.Boozer_I)
outputs.create_dataset('boozer_i_all', data = b.Boozer_I_all)

    
f.close()
    


###########
# FIGURES #
###########

plt.figure
bx.surfplot(b, js=0, fill=False, cmap=plt.cm.jet, ntheta=50, nphi=90, ncontours=25)
plt.savefig('Contour.eps')
plt.show()


plt.figure
bx.surfplot(b, js=0, cmap=plt.cm.jet, shading = 'gouraud')
plt.savefig('Filled.eps')
plt.show()

# ANOTHER EXAMPLE
# plt.figure
# bx.surfplot(b, js=0, fill=False, cmap=plt.cm.jet, levels=np.arange(0.8, 1.3, 0.05))
# plt.show()


