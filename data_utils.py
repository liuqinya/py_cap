#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written in Aug 2018
Use a bad station screening file to eliminate the bad stations from stream
Right now removes all three components (maybe refined to remove certain
components in the future)
"""

import os,sys
from obspy.core.stream import Stream


# add stats headers including trace.stats.distance/latitude/longitude
# figure out channel code number before select stream ...
def remove_station_from_stream(data_stream,inv,bad_station_file='bad-station.txt',verbose=True):
    if not os.path.isfile(bad_station_file):
        sys.exit('No bad_station_file: '+bad_station_file)
    lines=open(bad_station_file,'r').readlines()
    for i in range(len(lines)):
        net_sta=lines[i].split()[0]
        (net,sta)=net_sta.split('.')[0:2]
        if verbose:
            print(net, sta)
   
        found=False
        for nt in inv:
            for st in nt:
                if nt.code == net and st.code==sta:
                    found=True
                    break
        if not found:
            sys.exit('Check if '+net_sta+' exists in the inventory')

        str_bad=data_stream.select(network=net,station=sta)
        for trace in str_bad:
            data_stream.remove(trace)
            if verbose: 
                print('Remove trace ',trace)
    return
