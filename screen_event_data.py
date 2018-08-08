#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A data screening script that works with get_event_data.py
written by Q. Liu
August 2018
"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy import Catalog,Stream
import matplotlib.pyplot as plt
from obspy.geodetics import gps2dist_azimuth
import pickle
import data_utils

data_pkl='./data.pkl'
print('Unpacking '+data_pkl+' file ...')
f=open(data_pkl,'rb')
ev=pickle.load(f)
inv=pickle.load(f)
stream=pickle.load(f)
f.close()

#print('Remove pre-screened bad stations from stream ...')
#station_screen_file='bad-stations.txt'
#data_utils.remove_station_from_stream(stream,inv,bad_station_file=station_screen_file)

# prefilter to facilitate interpolation later

filter=True
f_min=0.01 # Hz
f_max=4.0 # Hz
print('Processing data (remove instrument response, anti-aliasing filter) ...')
stream.remove_response(water_level=50,output='VEL')
if filter:
    stream.detrend("linear")
    stream.taper(max_percentage=0.05, type="hann")
    stream.filter('bandpass',freqmin=f_min,freqmax=f_max,corners=2,zerophase=True)

for i in range(int(len(stream)/3)):
  str=stream[i*3:i*3+3]
  str.plot()
