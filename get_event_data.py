#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written in March-April 2018
rewrite from yiruzhou's original script
"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy import Catalog,Stream
import matplotlib.pyplot as plt
from obspy.geodetics import gps2dist_azimuth
import pickle

# ============= get event metadata =====================
# data server
client = Client("IRIS")

# events info
start_time = ["2014-08-05"]
end_time = ["2018-01-30"]
min_mag= [3.3]
# define search region
ev_min_lat=49.0; ev_max_lat=55.0; ev_min_lon=-120.0; ev_max_lon=-100.0 

# search for all events
catalog=Catalog() # catalog is a list of events
for i in range(len(min_mag)):
    print('Getting events ...')
    catalog+=client.get_events(
        minmagnitude=min_mag[i],
    	starttime=UTCDateTime(start_time[i]),
        endtime=UTCDateTime(end_time[i]),
        minlatitude = ev_min_lat, maxlatitude = ev_max_lat, 
        minlongitude = ev_min_lon, maxlongitude =ev_max_lon)
    
# get event interested in
ev = catalog[-1]
#print(event)
ev_time=ev.origins[0].time
ev_lat=ev.origins[0].latitude
ev_lon=ev.origins[0].longitude

# ============= get station metadata =====================
sta_max_radius=2 # degrees
print('Getting stations ....')
inv=client.get_stations(
    level='response',starttime=ev_time, endtime=ev_time + 24*60*60,
    latitude=ev_lat, longitude=ev_lon, maxradius=sta_max_radius)
#print(inv)
#inv.plot(projection='local')
# for net in inv:
#     for sta in net:
#         (dist,az,baz)=gps2dist_azimuth(ev_lat,ev_lon,sta.latitude,sta.longitude)
#         print(net.code,sta.code,dist/(100*1000))

# ###===========get the waveform================
bulk=[];
pre_origin_length=50
record_length = 180 # 3 min
exclude_sta=['TD06A']

for net in inv:
    for sta in net:
        # exclude the station with missing data (even after merge)
        if sta.code not in exclude_sta:
            bulk.append((net.code,sta.code,'*','?H?',ev_time-pre_origin_length,ev_time+record_length))
#print(bulk)
print('Requesting waveforms ...')
stream=client.get_waveforms_bulk(bulk, attach_response=True)
stream.merge(fill_value='interpolate')
stream.sort(keys=['network','station','location','channel'])
#print(stream)
#for i in range(len(stream)/3):
#    str=stream[i*3:i*3+3]
#    str.plot()

data_pkl='data.pkl'
print('Dump event, inventory and stream to '+data_pkl+'...')
f=open(data_pkl,'wb')
pickle.dump(ev,f)
pickle.dump(inv,f)
pickle.dump(stream,f)
f.close()
