#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written in Aug 2018
Use a bad station screening file to eliminate the bad stations from stream
Right now removes all three components (maybe refined to remove certain
components in the future)
"""

import os,sys,subprocess
from obspy.core.stream import Stream

# add stats headers including trace.stats.distance/latitude/longitude
# figure out channel code number before select stream ...
def remove_station_from_stream(data_stream,inv,bad_station_file='bad-quality-station.txt',verbose=True):
    if not os.path.isfile(bad_station_file):
        sys.exit('No bad_quality_station_file: '+bad_station_file)
    lines=open(bad_station_file,'r').readlines()
    for i in range(len(lines)):
        net_sta=lines[i].split()[0]
        sta_info=net_sta.split('.')
        if len(sta_info) == 2:
            (net,sta)=sta_info[0:2]
        elif len(sta_info) == 4:
            (net,sta,loc,chan1)=sta_info[0:4]
        else:
            sys.exit('Incorrect format in bad_quality_station_file: net.sta or net.sta.chan.cmp')
        if verbose:
            print(net, sta)
# confirm this net.sta[.loc.chan1] indeed exist
        found=False
        for nt in inv:
            for st in nt:
                for ch in st:
                    if len(sta_info) == 2:
                        if nt.code == net and st.code==sta:
                            found=True; break
                    else:
                        if nt.code == net and st.code==sta and ch.location_code == loc and ch.code[0:2] == chan1:
                            found=True
                            break
        if not found:
            sys.exit('Check if '+net_sta+' exists in the inventory')

# screen out the bad traces
        if len(sta_info) == 2:
            str_bad=data_stream.select(network=net,station=sta)
        elif len(sta_info) == 4:
            str_bad=data_stream.select(network=net,station=sta,location=loc,channel=chan1+'?')
        for trace in str_bad:
            data_stream.remove(trace)
            if verbose: 
                print('Remove trace ',trace)
    return

# --------------------------------------------------------
def write_weight_file(weight_file,dist_list,model_dir,check_bad_fit=False,bad_fit_station_file='bad-fit-stations.txt',min_dist_pnl=0., max_dist=1000.):
    
#    print('Write weight file '+ weight_file+'...')
    if not os.path.isfile(dist_list):
        sys.exit('No file: '+dist_list)
    if not os.path.isdir(model_dir):
        sys.exit('No greens function dir: '+model_dir)
        
    if check_bad_fit:
        print('Eliminate bad-fitting components based on '+bad_fit_station_file+'...')
        if not os.path.isfile(bad_fit_station_file):
            sys.exit('Check bad_fit_station_file '+bad_fit_station_file)
        # step 1: read in bad_fit station/componentlist
        tmp=open(bad_fit_station_file).readlines()
        bad_fit_station_list=list(map(str.rstrip,tmp))
        bad_fit={}
        for fit in bad_fit_station_list:
            [net_sta,cmps]=fit.split()
            bad_fit[net_sta]=cmps

    # step 2: write weight file based on dist.list and bad sta/comp file
    fw=open(weight_file,'w')
    for line in open(dist_list,'r'):
        [file,dist]=line.split()  #e.g., GS_OK035.z 22
        file1=file[:-2]  #e.g., GS_OK035
        dist_km=int(round(float(dist)))
        # check if greens function exist and obtain tp arrival
        sac_grn=model_dir+'/'+dist+'.grn.0'
        if not os.path.isfile(sac_grn):
            sys.exit('No such greens function file '+sac_grn)
        tp=subprocess.getoutput('saclst t1 t2 f '+sac_grn).split()[1]
        #print(tp)
        use_cmps={'P':1,'Z':1,'R':1,'T':1}

        # eliminate bad fit station (only used in a second run)
        if check_bad_fit:
            if file1 in bad_fit.keys():
                cmps=bad_fit[file1]
                for cmp in list(cmps):
                    if cmp not in ['P','Z','R','T']:
                        sys.exit('Error component names in '+bad_fit_station_file)
                    use_cmps[cmp]=0
        # adjust based on distance
        if dist_km < min_dist_pnl:
            use_cmps['P'] = 0
        if dist_km > max_dist:
            use_cmps={'P':0,'Z':0,'R':0,'T':0}

# note the gcap c code is very finicky about the format of the weight input
# need to read the source code to understand why 
        fw.write("%-10s %5d %.0f %.0f %.0f %.0f %.0f %.1f %.0f\n" %
                 (file1, dist_km, use_cmps['P'], use_cmps['P'], use_cmps['Z'], use_cmps['R'], use_cmps['T'], float(tp), 0))

    fw.close()
