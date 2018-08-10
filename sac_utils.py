#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written in March-April 2018
rewrite from yiruzhou's original script
"""

import numpy as np
import os,sys,glob,subprocess,copy, pickle

from obspy.clients.fdsn import Client

from obspy.core.stream import Stream
from obspy.core.util.attribdict import AttribDict

from obspy import UTCDateTime
from obspy.geodetics.base import locations2degrees
from obspy.geodetics.base import gps2dist_azimuth
from obspy.taup import TauPyModel
from obspy.signal.rotate import rotate2zne

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
#from sets import Set

# add stats headers including trace.stats.distance/latitude/longitude
# figure out channel code number before select stream ...
def stream_add_stats(data_stream,inv,evt,write_sac=False,rotate_in_obspy=False):
    for net in inv:
        for sta in net:
            str1=data_stream.select(network=net.code,station=sta.code)
#            print(str(net.code),str(sta.code),len(str1))
            if len(str1) == 0:
                continue
            # update in future to deal with multiple channel (total_number_of channels)
            if len(str1) % 3 !=0:
                print('Problem: missing components', str1); exit()
                
            for tr in str1:
                for chan in sta:
                    if tr.stats.channel == chan.code and tr.stats.location == chan.location_code:
                        break
                else:
                    print('Problem finding channel in inventory',tr); exit()
                tr.stats.coordinates={'latitude':chan.latitude,'longitude':chan.longitude}
                (tr.stats.distance,tr.stats.azimuth,tr.stats.back_azimuth)=gps2dist_azimuth(
                    chan.latitude, chan.longitude, evt.origins[0].latitude, evt.origins[0].longitude)
                if write_sac==True:
                    sac= AttribDict()
                    sac.kstnm=str(sta.code);
                    sac.knetwk=str(net.code);
                    sac.kcmpnm=str(chan.code)
                    sac.khole=str(chan.location_code)
                    sac.stla=chan.latitude; sac.stlo=chan.longitude; sac.stel=chan.elevation
                    sac.evla=evt.origins[0].latitude; sac.evlo=evt.origins[0].longitude;
                    sac.evdp=evt.origins[0].depth/1000. # in km
                    sac.mag=evt.magnitudes[0].mag; time=evt.origins[0].time
    
                    sac.nzyear,  sac.nzjday,  sac.nzhour,  sac.nzmin,  sac.nzsec,  sac.nzmsec=time.year, time.julday, time.hour, time.minute, time.second,  time.microsecond/1000
                    sac.o=0.
                    sac.b=tr.stats.starttime-time # this is very important!!
                    sac.kevnm=str(time)
                    sac.cmpaz=chan.azimuth
                    # dip is from horizontal downward; inc is from vertical downward
                    sac.cmpinc=chan.dip+90
                    sac.gcarc = locations2degrees(evt.origins[0].latitude, evt.origins[0].longitude, chan.latitude, chan.longitude)
                    sac.dist,sac.az,sac.baz= tr.stats.distance/1000,tr.stats.azimuth,tr.stats.back_azimuth
                    tr.stats.sac=sac
                    tr_name=sta.code+'.'+net.code+'.'+chan.location_code+'.'+chan.code+'.sac'
                    tr.write(tr_name,format='SAC')
                    



# def sac_header_from_event_station(evt,  net,  sta):
    
#     # note   sac.b=start time - evt.origins[0].time and sac.e needs to be set when the trace has been acquired
#     sac=AttribDict()
#     # station location: obspy missing station burial information!!
#     sac.kstnm=str(sta.code) 
#     sac.knetwk=str(net.code)
#     sac.stla=sta.latitude
#     sac.stlo=sta.longitude
#     sac.stel=sta.elevation

#     # component info set in sac_header_from_channel: sac.kcmpnm=channel.code ; cmpaz, cmpinc,
#     # distance/ azimuth
#     gcarc = locations2degrees(evt.origins[0].latitude, evt.origins[0].longitude, sta.latitude, sta.longitude)
#     azi_baz = gps2dist_azimuth(evt.origins[0].latitude, evt.origins[0].longitude, sta.latitude, sta.longitude)
#     sac.dist, sac.az,  sac.baz,  sac.gcarc=  azi_baz[0]/1000.0, azi_baz[1], azi_baz[2], gcarc 
#     # event location
#     sac.evla=evt.origins[0].latitude; sac.evlo=evt.origins[0].longitude; sac.evdp=evt.origins[0].depth/1000.
#     sac.mag=evt.magnitudes[0].mag; time=evt.origins[0].time
#     sac.kevnm=str(time)
#     sac.nzyear,  sac.nzjday,  sac.nzhour,  sac.nzmin,  sac.nzsec,  sac.nzmsec=time.year, time.julday, time.hour, time.minute, time.second,  time.microsecond/1000
#     sac.o=0
#     return sac

# def sac_header_from_channel(chan, sac):
#     sac.kcmpnm=str(chan.code)# convert unicode
#     #sac.kcmpnm = str(chan)
#     sac.cmpaz=chan.azimuth
#     sac.cmpinc=chan.dip+90 # dip is from horizontal downward; inc is from vertical downward
#     sac.khole=str(chan.location_code)
