# py_cap
This is a collection of Python utilities our group members wrote for running GCAP (a source mechanism inversion code by Zhu L. http://www.eas.slu.edu/People/LZhu/home.html) efficiently. Note that FK and GCAP3D codes from Dr. Lupei Zhu (downloaded on 2022/02/12)'s website are included for reference. We will take these codes off if requested. Here is Prof. Zhu's short-course on the FK code: http://home.ustc.edu.cn/~zhx0906/AdvancedGeophysics/_document/note.pdf

For a complete manual page on how to setup and use FK, GCAP, and this utility tool, please consult to the group research handbook (Chapter 2). Here is a short description.

For your own source inversion, these utility scripts need to be slightly modified and run in the following order:

get_event_data.py -- select events based on start/end time, region and magnitude; identify one event; get station meta data (i.e., inventory) in the region within given epicentral distance; get bulk waveform for all stations; (in the second round eliminate traces from stream based on bad-quality-stations.txt); pickle event catalog, station inventory and stream objects.

select_event_data.py -- read pickled data file, plot three-component station data with P and S arrivals (based on ak135 model). Open a separate text file named 'bad-quality-stations.txt', and every row has the format 'NT.STA.loc.[BH]H  Z12  Noisy' to indicate the bad stations (right now copmonents Z12 are not used, the entire station is eliminated). Remove instrument response and filter are optional. Try run this multiple times (e.g., without removing response, or with proper filter) to get rid of stations that are clipped, extremely noisy or with glitches.

Now you can rerun get_event_data.py --- eliminate the bad quality station, and repickle data.

Now move data.pkl into data/ sub-directory

process_data_to_sac.py (run under data/ as ../process_data_to_sac.py)-- load in waveform stream, remove instrument response, filter, add stats and write into sac files. Write 'fk_bash.cmd' bash script to rotate horizontal components, convert displacement to cm, rename sac files to gcap expected format, write distance list, and write (and possibly execute) fk commands to generate the corresponding fk Greens function synthetics for the distance list. Parameters to be adjusted: T_min, Tmax, fk_dir, model, depths list

prepare_for_gcap.py -- prepare the weight file (weight.dat.depth) and create bash script 'cap_auto.bash' in the gcap directory to write the 'cap.pl' command line based on choices of inversion parameters (filters, distance scale, plot scale, weights, window selection, etc). Create all the necessary data links. Note: a line of system("ps2pdf $outps $outpdf"); has been added to cap_plt.pl to generate output plots in pdf format which is much more friendly to view. depth.pl has been updated to get the mechanism misfit plot vs depth. Parameters to be adjusted: green_dir, event_dir, model, depths list, test mw, pnl/surf_tmin/max, pnl weight, dist scale, plot_scale, body_surf_shift, min_dist_pnl, max_dist. This can be also done for several iterations, everytime update the 'data/bad-fit-stations.txt', every row has the format of 'NT_STA  [PZRT]', and eliminate those stations/components with bad fits.

Before running prepare_for_gcap.py, update cap.pl in the gcap directory to reflect the correct path to cap_plt.pl: 
require "$HOME/gcap/cap_plt.pl
If no non-double couple inversion needed, ignore warnings like this:
Unable to open /home/xyz/gcap/OK_event/Ok/Ok_3/33.grn.b
Warning: non-DC components are ignored



