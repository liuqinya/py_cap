#!/usr/bin/env python

import sys, os, glob, commands

green_dir='/data2/gcap-inv/induced-earthquakes/'
gcap_dir='/data2/gcap-inv/gcap'
data_dir='/data2/gcap-inv/induced-earthquakes/data/'
data_link_dir_in_gcap='data' #important to strip out the '/'
data_link=gcap_dir+'/'+data_link_dir_in_gcap
model='vmn'
depth=5  # write this into a loop?
deltat=0.05

gcap_command_file=gcap_dir+'/cap_auto.bash'
dist_list='dist.list'

model_dir=model+'_'+str(depth)
weight_file='weight.dat.'+str(depth)
print('write '+weight_file+' for depth '+str(depth)+' km and data in '+data_dir+'...\n')
fw=open(weight_file,'w')
for line in open(dist_list,'r'):
    [file,dist]=line.split()
    file1=file[:-2] # distance
    #dist=os.system('grep '+file+' '+dist_list+' | awk \'{print $2}\'')
    dist_km=int(round(float(dist)))
    #print(file1, dist_km)
    sac_grn=green_dir+'/'+model+'/'+model_dir+'/'+dist+'.grn.0'
    if not os.path.isfile(sac_grn):
        sys.exit('No such greens function file '+sac_grn)
    tp=commands.getoutput('saclst t1 t2 f '+sac_grn).split()[1]
    #print(tp)
    # note the cap c code is very finicky about the format of the weight input
    # need to read the source code to understand why
    fw.write("%-10s %5d %.0f %.0f %.0f %.0f %.0f %.1f %.0f\n" %
             (file1, dist_km, 1, 1, 1, 1, 1, float(tp), 0))
fw.close()

# write command
# use my modified cap_plt.pl file
print('Write gcap command file '+gcap_command_file+' ...')
mw=3.8 # event test magnitude
# filter parameters: tweak these to see fits
pnl_fmin=0.08; pnl_fmax=0.4
surf_fmin=0.05; surf_fmax=0.1
dist_scale='-D2/1/0.5'; plot_scale='-P0.3/45'
body_surf_weight='-S2/5/0'; window='-T30/70'

# link data dir for cap command in cap.pl: /cap data vmn_5
fw=open(gcap_command_file,'w')
fw.write('#!/bin/bash\ncd '+gcap_dir+'\n')
if os.path.isfile(data_link):
    sys.exit('data link dir '+data_dir+' already exists as a file')
if os.path.islink(data_link):
    os.unlink(data_link)
os.symlink(data_dir, data_link)
fw.write('./cap.pl -M'+model_dir+'/'+str(mw)+ \
         '        -H'+str(deltat)+ \
         '        -C'+str(pnl_fmin)+'/'+str(pnl_fmax)+'/'
         +str(surf_fmin)+'/'+str(surf_fmax)+ \
         ' -W1 -X10 '+\
         dist_scale+' '+body_surf_weight+' '+plot_scale+' '+window+ \
         ' -Z'+weight_file+' ' +\
         ' -G'+green_dir+' '+\
         data_link_dir_in_gcap+'\n')
fw.close()

# why does it fail
