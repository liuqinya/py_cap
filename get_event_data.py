#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written in March-April 2018 by Q. Liu
rewrite from yiruzhou's original script
"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy import Catalog,Stream
import matplotlib.pyplot as plt
from obspy.geodetics import gps2dist_azimuth
import pickle
import data_utils
import os

# ============= get event metadata =====================
# data server
client = Client("IRIS")

# events info
start_time = ["2016-01-01"]
end_time = ["2016-12-31"]
min_mag= [4.0]
# define search region (OK induced seismicity zone --> a square box)
ev_min_lat=33.70; ev_max_lat=37.02; ev_min_lon=-98.96; ev_max_lon=-96.11
# event info for Montney
# events info
#start_time = ["2014-08-05"]
#end_time = ["2018-01-30"]
#min_mag= [3.3]
# define search region
#ev_min_lat=49.0; ev_max_lat=55.0; ev_min_lon=-120.0; ev_max_lon=-100.0 


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
    
# list all the events to get the event interested in
#for i in range(len(catalog)):
#    print(i, catalog[i].origins[0].time, catalog[i].magnitudes[0].mag)

event_number=-1
ev = catalog[event_number]
#print(event)
ev_time=ev.origins[0].time
ev_lat=ev.origins[0].latitude
ev_lon=ev.origins[0].longitude
ev_dep=ev.origins[0].depth/1000.
print('event info:', ev_time, ev_lat, ev_lon, ev_dep)
# event info: 2016-01-01T11:39:39.800000Z 35.6688 -97.4065 5.825
# ============= get station metadata =====================
sta_max_radius=3 # degrees
print('Getting broadband stations ....')
inv=client.get_stations(
    level='response',starttime=ev_time, endtime=ev_time + 24*60*60,
    latitude=ev_lat, longitude=ev_lon, maxradius=sta_max_radius,
    channel="BH?,HH?")
# uncomment this part for station verbose output
# print(inv)
# inv.plot(projection='local')
# nsta=0
# for net in inv:
#     nsta+=len(net)
#     for sta in net:
#         (dist,az,baz)=gps2dist_azimuth(ev_lat,ev_lon,sta.latitude,sta.longitude)
#         print(net.code,sta.code,dist/1000)
# print('Total number of available broadband stations: '+str(nsta))
        
# # # ###===========get the waveform================
bulk=[];
pre_origin_length=50
record_length = 180 # 3 min
# examine waveforms to exclude bad stations
exclude_sta=[] 

for net in inv:
    for sta in net:
        # exclude the station with missing data (even after merge)
        if sta.code not in exclude_sta:
            bulk.append((net.code,sta.code,'*','BH?,HH?',ev_time-pre_origin_length,ev_time+record_length))
#print(bulk)
print('Requesting waveforms ...')
stream=client.get_waveforms_bulk(bulk, attach_response=True)
stream.merge(fill_value='interpolate')
stream.sort(keys=['network','station','location','channel'])
print(stream)
print('Total number of stations acquired: ',int(len(stream)/3))

# # first run: skip this part and save data.pkl first,
# # then use screen_event_data.py
# # to quickly screen three-component data
# # at this point you need to screen the waveform, either to
# # go back and fill the exclude_sta list, or assemble the bad stations
# # in bad-quality-stations.txt and call the following

station_screen_file='bad-quality-stations.txt'
if os.path.isfile(station_screen_file):
    print('Remove pre-screened bad stations from stream ...')
    data_utils.remove_station_from_stream(stream,inv,bad_station_file=station_screen_file,verbose=False)
    nstation=int(len(stream)/3)
    print('Total number of stations left after screening: ',nstation)
else:
    print('Try run screen_event_data.py to screen data quality')

# double check station after elimination
# for i in range(nstation):
#    str=stream[i*3:i*3+3]
#    str.plot()

data_pkl='data.pkl'
print('Dump event, inventory and stream to '+data_pkl+'...')
f=open(data_pkl,'wb')
pickle.dump(ev,f)
pickle.dump(inv,f)
pickle.dump(stream,f)
f.close()

