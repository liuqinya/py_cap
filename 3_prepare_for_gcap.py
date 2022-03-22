#!/usr/bin/env python3

import sys, os, glob, subprocess
import data_utils
# use absolute path for gcap_dir, event_dir 

# run directory for gcap code, data dir is imbedded as a link
gcap_dir='/home/lqy/gcap_test_20220303/gcap'

# event data and greens function directory
event_dir='/home/lqy/gcap_test_20220303/20080418093700'
gcap_command_file=event_dir+'/cap_auto.bash'
green_dir=event_dir+'/' # no need to add model for the cap.pl command
data_dir=event_dir+'/data'
dist_list=data_dir+'/'+'dist.list'

##  ****** adjust the following parameters ******* 
model='cus'
#write depths into a loop, should be consistent with the depths list in process_data_to_sac.py
depths=[10, 15, 20] # avoid depth on model interfaces
deltat=0.05
mw=5.0 # event test magnitude
include_iso_clvd=True  # also search for iso and clvd componenents in the inversion (requires fk already computed grn.[abc])
niter=1  # number of iterations

# filter parameters: tweak these to see fits
pnl_tmin=2.5; pnl_tmax=12.5 # sec
surf_tmin=10; surf_tmax=20 # sec
pnl_fmin=1./pnl_tmax; pnl_fmax=1./pnl_tmin
surf_fmin=1./surf_tmax; surf_fmax=1./surf_tmin
dist_scale='-D0.5/1/0.5' # pnl weight (2)/scaling (1), surf scaling (0.5)
plot_scale='-P0.2/45' # use bigger number to increase plot amplitude
body_surf_shift='-S2/5/0'; window='-T30/70'

# minimum distance to start separating Pnl and surface waves
# check screen_event_data.py output to allow at least 2 cycles between tp and ts
min_dist_pnl=85 # km
# check screen_event_data.py output for the SNR deterioration with distance
max_dist=300 # km

# inversion parameters:
if include_iso_clvd:
    iso_clvd='-J-0.9/0.05/-0.45/0.05'
else:
    iso_clvd=''
## ************************************************

bad_fit_station_file=data_dir+'/bad-fit-stations.txt'
if os.path.isfile(bad_fit_station_file):
    check_bad_fit=True
else:
    check_bad_fit=False

for depth in depths:
    model_dep=model+'_'+str(depth)
    model_dir=green_dir+'/'+model+'/'+model_dep

    weight_file='weight.dat.'+str(depth) # no need to add data_dir/
    print('write '+data_dir+'/'+weight_file+' for depth '+str(depth)+' km and data in dir '+data_dir+' and dist_list ' + dist_list+' ...\n')
    
    data_utils.write_weight_file(data_dir+'/'+weight_file,dist_list,model_dir,check_bad_fit,bad_fit_station_file,min_dist_pnl,max_dist)
 
# write command
# use my modified cap_plt.pl file
print('Write gcap command file '+gcap_command_file+' ...')

# link data dir for cap command in cap.pl: /cap data vmn_5
fw=open(gcap_command_file,'w')
depths_str = " ".join(list(map(str,depths))) # convert the depths list to str
fw.write('#!/bin/bash\ncd '+gcap_dir+'\n')
# the gcap code is finicky on how long the data_dir can be, so always link data/
fw.write('rm -f data\nln -s '+data_dir+' data\n')
fw.write('for h in '+depths_str+'; do \n ./cap.pl -M'+model+'_$h'+'/'+str(mw)+ \
         '        -H'+str(deltat)+ \
         '        -C'+str(pnl_fmin)+'/'+str(pnl_fmax)+'/'
         +str(surf_fmin)+'/'+str(surf_fmax)+ \
         ' -W1 -X10 '+\
         dist_scale+' '+body_surf_shift+' '+plot_scale+' '+window+ \
         ' -Z'+'weight.dat.$h'+' '+iso_clvd+' -N'+str(niter) +\
         ' -G'+green_dir+'  data'+\
         '\ndone\n')
fw.write('grep -h Event data/'+model+'_*.out > data/junk.out\n')
fw.write('./depth.pl data/junk.out data > data/'+model+'.ps\n')
fw.write('ps2pdf data/'+model+'.ps data/'+model+'.pdf\n')
fw.close()
