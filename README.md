# py_cap
Python utilities for running GCAP (a source mechanism inversion code by Zhu L.)

Currently a manual page for using FK and GCAP has been setup at:
https://www.overleaf.com/read/dtxzwrqntztd

These utility scripts need to modified and run in the following order:

get_event_data.py -- select events based on start/end time, region and magnitude; get station meta data in the region within given epicentral distance; get bulk waveform for all stations and pickle event catalog, station inventory and stream objects.

process_data_to_sac.py -- load in waveform stream, remove instrument response, filter, add stats and write into sac files. Write `fk_bash.cmd' bash script to rotate horizontal components, convert displacement to cm, rename sac files to gcap expected format, write distance list, and write (and possibly execute) fk commands to generate the corresponding fk Greens function synthetics for the distance list.

prepare_for_gcap.py -- prepare the weight file (weight.dat.depth) and create bash script `cap_auto.bash' in the gcap directory to write the `cap.pl' command line based on choices of inversion parameters (filters, distance scale, plot scale, weights, window selection, etc). Create all the necessary data links.

