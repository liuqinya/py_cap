#include "cap3D.h"

void    taper(float *aa, int n) {
  int	i, m;
  float	tt, pi1;
  m = rint(0.3*n);
  pi1 = 3.1415926/m;
  for (i=0; i<m; i++) {
    tt = 0.5*(1.-cos(i*pi1));
    aa[i] *= tt;
    aa[n-i-1] *= tt;
  }
}


int discard_bad_data(int nda, DATA *obs, SOLN sol, float sig, float rms_cut[]) {
   int i, j, n;
   float minCC = 50.;
   COMP	*spt;
   n = 0;
   for(i=0;i<nda;i++,obs++) {
     spt = obs->com;
     for(j=0; j<NCP; j++,spt++) {
        if (sol.cfg[i][j]<minCC && sol.error[i][j]/sig > rms_cut[j]) {
	   spt->on_off = 0;
	   n++;
	}
     }
   }
   return(n);
}


float *cutTrace(float *trace, int npts, int offset, int n) {
   int m;
   float *cut;
   cut = (float *) calloc(n, sizeof(float));
   if (cut == NULL) return cut;
   m = n+offset;
   if (offset<0) {
      if (m>npts) m = npts;
      if (m>0) memcpy(cut-offset, trace, m*sizeof(float));
   } else {
      if (m>npts) n = npts-offset;
      if (n>0) memcpy(cut, trace+offset, n*sizeof(float));
   }
   taper(cut, n);
   return cut;
}


int check_first_motion(float mt[3][3], FM *fm, int n, float fm_thr) {
  int	i;
  FM	*pt;
  for(pt=fm,i=0;i<n;i++,pt++) {
    if (pt->type*radpmt(mt, pt->alpha, pt->az, abs(pt->type))/abs(pt->type)<fm_thr)
      return -1;
  }
  return 0;

}
