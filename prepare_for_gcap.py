#!/usr/bin/env python

import sys, os, glob, commands
# use absolute path for gcap_dir, event_dir 

# run directory for gcap code, data dir is imbedded as a link
gcap_dir='/data2/gcap-inv/gcap'

# event data and greens function directory
event_dir='/data2/gcap-inv/test'
gcap_command_file=event_dir+'/cap_auto.bash'
model='vmn'
green_dir=event_dir+'/' # no need to add model for the cap.pl command
data_dir=event_dir+'/data'
dist_list=data_dir+'/'+'dist.list'
if not os.path.isfile(dist_list):
    sys.exit('No file: '+dist_list)

depths=['1','3','5','7','9'] #write depths into a loop, should be consistent with the depths list in process_data_to_sac.py
deltat=0.05

model_dep=model+'_'+str(depth)
model_dir=green_dir+'/'+model+'/'+model_dep
if not os.path.isdir(model_dir):
    sys.exit('No greens function dir: '+model_dir)

weight_file='weight.dat.'+str(depth) # no need to add data_dir/
print('write '+data_dir+'/'+weight_file+' for depth '+str(depth)+' km and data in '+data_dir+'...\n')
fw=open(data_dir+'/'+weight_file,'w')
for line in open(dist_list,'r'):
    [file,dist]=line.split()
    file1=file[:-2] # distance
    #dist=os.system('grep '+file+' '+dist_list+' | awk \'{print $2}\'')
    dist_km=int(round(float(dist)))
    #print(file1, dist_km)
    sac_grn=model_dir+'/'+dist+'.grn.0'
    if not os.path.isfile(sac_grn):
        sys.exit('No such greens function file '+sac_grn)
    tp=commands.getoutput('saclst t1 t2 f '+sac_grn).split()[1]
    #print(tp)
    # note the gcap c code is very finicky about the format of the weight input
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
depths_str = " ".join(depths) # convert the depths list to str
fw.write('#!/bin/bash\ncd '+gcap_dir+'\n')
fw.write('for h in '+depths_str+'; do ./cap.pl -M'+model+'_$h'+'/'+str(mw)+ \ 
         '        -H'+str(deltat)+ \
         '        -C'+str(pnl_fmin)+'/'+str(pnl_fmax)+'/'
         +str(surf_fmin)+'/'+str(surf_fmax)+ \
         ' -W1 -X10 '+\
         dist_scale+' '+body_surf_weight+' '+plot_scale+' '+window+ \
         ' -Z'+weight_file+' ' +\
         ' -G'+green_dir+' '+\
         data_dir+'; done')
fw.close()

