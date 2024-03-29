#!/usr/bin/env python3
"""
Written in March-April 2018
updated August 2018
"""

from obspy import Catalog, Stream
import matplotlib.pyplot as plt
from obspy.geodetics import gps2dist_azimuth
import sac_utils
import obspy.taup as taup
import sys,math,glob, os
import pickle

# this code can be only run in the data directory where data.pkl sits.
data_dir=os.getcwd() 
parent_dir=os.path.dirname(os.path.normpath(data_dir))
data_pkl=data_dir+'/data.pkl'
print('Unpacking '+data_pkl+' file ...')
try:
    f=open(data_pkl,'rb')
except FileNotFoundError:
    raise Exception('Error opening file '+data_pkl)
ev=pickle.load(f)
inv=pickle.load(f)
stream=pickle.load(f)
f.close()

# anti-aliasing filter before interpolation (delta=0.05)
T_min=0.25 # sec
T_max=100 # sec
f_min=1./T_max # Hz
f_max=1./T_min # Hz
filter=True
print('Processing data (remove instrument response, apply anti-aliasing filter) ...')
# water-level of 100 is clearly too big, the default 60 is probably alright.
# use plot=True to see the effect of water-level. water level 10 is also ok.
stream.remove_response(water_level=50,output='VEL')
if filter:
    stream.detrend("linear")
    stream.taper(max_percentage=0.05, type="hann")
    stream.filter('bandpass',freqmin=f_min,freqmax=f_max,corners=2,zerophase=True)
print('Write data to sac files ...')
sac_utils.stream_add_stats(stream,inv,ev,write_sac=True)

#------------
fk_dir='/directory-for-fk-source-code' # absolute path to fk run code
if not os.path.isdir(fk_dir):
    raise Exception('Error finding directory '+fk_dir)

model='cus'
green_dir=parent_dir+'/'+model # absolute path to greens function output dir
model_file=green_dir+'/'+model # absolute path to the model file 
if not os.path.isfile(model_file):
    raise Exception('No such file '+model_file)

depths=[10, 15, 20] # avoid depth on model interfaces
inv_iso_clvd=False  # also compute greens functions .[abc] for iso+clvd inversions
run_bash=False
deltat=0.05
syn_rec_length=160 # in seconds
nt= int(math.pow(2,round(math.log(syn_rec_length/deltat,2))))

pre_origin_length=50
record_length = 180 # 3 min

print('\nWriting fk_bash.cmd for processing in sac ...')
# rotate sac (STA.NET.LOC.CMP.sac) in SAC
print('  rotate seismograms ...\n')
fp=open('fk_bash.cmd','w')
fp.write('#!/bin/bash\ncd '+data_dir+'\nsac <<EOF\n') #set -x to debug bash code
fp.write('  cut  o '+str(-pre_origin_length+1)+' '+str(record_length-1)+'\n')

for efile in glob.glob(data_dir+'/*.*.*.??[1E].sac'):
    nfile=efile.replace('1.sac','2.sac').replace('E.sac','N.sac')
    rfile=efile.replace('1.sac','R.sac').replace('E.sac','R.sac')
    tfile=nfile.replace('2.sac','T.sac').replace('N.sac','T.sac')
#    print(efile, nfile, rfile, tfile)
    if not os.path.isfile(nfile):
        sys.exit('Error finding matching '+nfile+' to '+efile)
#   use the cut command in sac to ensure the same length of e,n or 1,2 files
    fp.write('  r '+efile+' '+nfile+'\n')
#    fp.write('  cut off\n')
    fp.write('  rotate\n')
    fp.write('  w '+rfile+' '+tfile+'\n')

# convert from m to cm, interpolate to the deltat used in fk
print('  convert units to cm and interpolate with '+str(deltat)+'...\n')
fp.write('r *.*.*.??[RTZ].sac\nmul 100\ninterp delta '+str(deltat)+'\nw over\n')
fp.write('quit\nEOF\n')

fp.write('mkdir -p EN12; mv -f *[EN12].sac EN12/\n');
# rename sac files  NET_STA.[RTZ]
print('  rename sac files ...\n')
fp.write("""for file in *.??[RTZ].sac; do
  sta=`echo $file | cut -d . -f 1`
  net=`echo $file | cut -d . -f 2`
  cmp=`echo $file | cut -d . -f 4`
  cmp1=`echo ${cmp:2:1} | tr [:upper:] [:lower:]`
  mv -f $file ${net}_${sta}.$cmp1
done\n""")

# write distance list
dist_list='dist.list'
print('  get distance list for fk ...\n')
fp.write('saclst dist f *.z | sort -k 2 -g | awk \'{printf "%s  %.0f \\n", $1, $2}\'> '+dist_list+'\n')
fp.write('dist=`awk \'{print $2}\' '+dist_list+' | sort -g -u | awk  \'{printf " %.0f", $1}\'`\n')

# write fk command
print('  write fk command ...\n')
fp.write('cp '+model_file+' '+fk_dir+'\n')
for depth in depths:
    fk_command='./fk.pl -M'+model+'/'+str(depth)+'/k -N'+str(nt)+'/'+str(deltat)+' $dist\n'
    fp.write('\ncd '+fk_dir+'\necho *****'+fk_command+fk_command)
    if inv_iso_clvd:
        fk_command='./fk.pl -M'+model+'/'+str(depth)+'/k -N'+str(nt)+'/'+str(deltat)+' -S0 $dist\n'
        fp.write('echo *****'+fk_command+fk_command)
    fp.write('rm -rf '+model_file+'_'+str(depth)+'\nmv '+model+'_'+str(depth)+' '+green_dir+'\n')
fp.close()

# Question: any need to pre-pad greens functions?
if run_bash:
    print('Run fk command ...\n')
    os.system('bash fk_bash.cmd')

#====================== more processing tips==============================#
# could choose to decimate in obspy instead of sac later
#tr.decimate(factor=4, strict_length=False)
# pick first P and S arrivals: the arrivals can be picked through taup package by defining a taup model based on tvel ASCII file (depth pVel sVel Density) ... skip here ...
#model='vmn'
#build taupy model ... write a subroutine
#taupy_model=TauPyModel(model='../models/vmn.npz')
#taupy_model.get_travel_times(source_depth_in_km=5,distance_in_km=)

# stream_add_stats() adds proper headers for plotting record section (trace.stats.distance/latitude/longitude; get ev_coord from event catalog
# plot record section
# stream.select(component='Z').plot(type='section')
