
#  Usage: dist_list_from_sacfiles.py sac-files

import sys
from obspy import read

# read input arguments
if len(sys.argv) <= 1:
    raise BaseException('Usage: dist_list_from_sacfiles.py sac-files')

dist_list=[]
az_list=[]
# loop over events
for sacfile in sys.argv[1:]:
    st=read(sacfile)  # sac header goes into st[0].stats.sac
    tr=st[0]
    if hasattr(tr.stats.sac, 'dist') and hasattr(tr.stats.sac,'az'):
        dist_list.append(float("{:.2f}".format(tr.stats.sac.dist)))
        az_list.append(float("{:.2f}".format(tr.stats.sac.az)))
    else:
        raise BaseException('Problem inquiring dist and az header of file'+sacfile)
    #print(("%s   %.3f   %.3f") %(sacfile, dist, az))

dist_list = sorted(set(dist_list))
print(' '.join(map(str,dist_list)))



    
