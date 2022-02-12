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
import matplotlib.pyplot as plt 
from matplotlib.dates import date2num
from obspy.taup import TauPyModel
import sac_utils

model = TauPyModel(model="ak135")

data_pkl='./data.pkl'
print('Unpacking '+data_pkl+' file ...')
f=open(data_pkl,'rb')
ev=pickle.load(f)
inv=pickle.load(f)
stream=pickle.load(f)
f.close()

# print('Remove pre-screened bad stations from stream ...')
use_existing_screen_file=False
if use_existing_screen_file:
    station_screen_file='bad-quality-stations.txt'
    data_utils.remove_station_from_stream(stream,inv,bad_station_file=station_screen_file)

print('Processing data (e.g., remove instrument response, anti-aliasing filter) ...')
# to screen raw data for clipping, set remove_response and filter to be both False
# remove instrument response
remove_response=True
if remove_response: 
    stream.remove_response(water_level=50,output='VEL')

# apply filter to screen specific type of data: surface waves
bp_filter=True
if bp_filter:
    # use [2.5, 20] sec to cover the bands for pnl and surf
    T_min=2.5  # 0.5
    T_max=20. # 40
    f_min=1./T_max # Hz
    f_max=1./T_min # Hz

    stream.detrend("linear")
    stream.taper(max_percentage=0.05, type="hann")
    stream.filter('bandpass',freqmin=f_min,freqmax=f_max,corners=2,zerophase=True)

#screen_by_snr=True; pnl_snr_min=1.5; surf_snr_min=2.5
surf_length=60.
# plot three-component data with event origin time indicated.
print('Plotting three-component traces: event origin time, P and S arrival based on ak135 model ...')
for i in range(int(len(stream)/3)):
  str=stream[i*3:i*3+3]
  fig = plt.figure()
  str.plot(fig=fig,show=False)
  Ptime=0; Stime=0
  for j in range(3):
      # add origin time
      ax=fig.axes[j]
      ax.axvline(date2num(ev.origins[0].time.datetime),lw=2)
#      print('depth',ev.origins[0].depth/1000,'dist',str[j].stats.distance/1000./111.195)
      if j == 0: 
          # calculate P and S arrival times based on ak135
          arrivals = model.get_travel_times(source_depth_in_km=ev.origins[0].depth/1000.,distance_in_degree=str[j].stats.distance/1000./111.195)
          for k in range(len(arrivals)):
              if (arrivals[k].name == 'P' or arrivals[k].name == 'p') and Ptime == 0:
                  Ptime=arrivals[k].time
              elif (arrivals[k].name == 'S' or arrivals[k].name == 's') and Stime == 0:
                  Stime=arrivals[k].time
          print('%s.%s Dist %4.1f   Ptime %4.1f   Stimes %4.1f' % (str[j].stats.network, str[j].stats.station, str[j].stats.distance/1000., Ptime,Stime))
          Parr=ev.origins[0].time+Ptime
          Sarr=ev.origins[0].time+Stime
      ax.axvline(date2num(Parr.datetime),lw=1,ls='dashed',c='g')
      ax.axvline(date2num(Sarr.datetime),lw=1,ls='dashed',c='m')

      ## display snr
      ts=str[j].stats.starttime; t0=ev.origins[0].time
      t1=Parr; t2=Sarr; t3=Sarr+surf_length
      dt=str[j].stats.delta
      n0=int((t0-ts)/dt); n1=int((t1-ts)/dt)
      n2=int((t2-ts)/dt); n3=int((t3-ts)/dt)
      if j==0:
          pnl_snr=abs(str[j].data[n1:n2]).mean()/abs(str[j].data[n0:n1]).mean()
          print('Pnl SNR    = %4.1f' %  (pnl_snr))
      
      surf_snr=abs(str[j].data[n2:n3]).mean()/abs(str[j].data[n0:n1]).mean()
      print('Surf SNR %2d = %4.1f' % (j, surf_snr))

  plt.show()
  
  # use the output info to determine at what distance Pnl starts to separate
  # out from surface waves (e.g., dist_pnl_from_surf=85 km)

  # to do list:
  # * add picked P and S arrival to trace.stats.
  # * rotation in obspy
  # * screen data automatically by SNR
  # * indicate a single component of a station to be kept
  # * translate SNR for Pnl and surf into weight.dat file
  

