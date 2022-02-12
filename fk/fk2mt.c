/******************************************************************
   F-K: @(#) fk2mt.c             1.0 11/24/2013
 
   Copyright (c) 2013 by L. Zhu
   See README file for copying and redistribution conditions.

 Convert FK Green's functions to MT Green's functions

 Revision History:
   11/24/2013	Lupei Zhu	Initial coding.
******************************************************************/
#include <stdio.h> 
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sac.h"
#include "radiats.h"

int main(int argc, char **argv) {
  int	i,npt;
  char	*ccc,nam[128],outnm[128];
  float	dt, tp, ts, az, sa, ca, sa2, ca2, *fk[12], *mt[18];
  SACHEAD	hd;

  if(argc != 4) {
    fprintf(stderr,"Usage: %s station_azimuth_in_deg FirstComponentOfFKGF FirstComponentOfMTGF\n\
 The input GFs are 9- or 12-component outputs from FK\n\
 The output GFs are 18-component moment-tensor GFs, in the order of\n\
 Mxx[3=Up,Radial (away positive),Transverse (CW positive)], 2Mxy[3], 2Mxz[3], Myy[3], 2Myz[3], and Mzz[3] (x=North, y=East, z=Down)\n",argv[0]);
    return -1;
  }

  sscanf(argv[1],"%f",&az);
  strcpy(nam, argv[2]);
  strcpy(outnm, argv[3]);

  ccc = nam + strlen(nam) - 1;
  for(i=0; i<12; i++) {
        if ( (fk[i]=read_sac(nam,&hd)) == NULL ) break;
	if ( i==0 ) {
           npt = hd.npts;
	   dt = hd.delta;
	   tp = hd.t1;
	   ts = hd.t2;
	}
	else if (hd.npts != npt) {
	   fprintf(stderr,"number points in %s not agree with %d\n",nam,npt);
	   return -1;
	}
        (*ccc)++;
	if (*ccc == '9') (*ccc)+=40;	/* explosion components start at 'a' instead of '0' */
  }
  if (i<9) return -1;
  for(;i<12;i++) fk[i] = (float *) calloc(npt, sizeof(float));

  for(i=0;i<18;i++) mt[i] = (float *) malloc(npt*sizeof(float));
  fk2mtg(az, npt, fk, mt);

  /* output */
  ccc = outnm + strlen(outnm) - 1;
  hd.t1 = tp; hd.t2 = ts;
  *ccc = 'a';
  for(i=0;i<18;i++) {
     hd.cmpaz = 0.;
     hd.cmpinc = 90.;
     if ( i%3 == 0 ) hd.cmpinc = 0.;
     if ( i%3 == 1 ) hd.cmpaz = az;
     if ( i%3 == 2 ) hd.cmpaz = az+90.; if (hd.cmpaz>360.) hd.cmpaz -= 360.;
     write_sac(outnm,hd,mt[i]);
     (*ccc)++;
  }

  return 0;

}
