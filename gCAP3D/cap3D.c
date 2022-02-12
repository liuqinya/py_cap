/****************************************************************
				cap3D.c
  Generalized Cut-and-Paste (gCAP3D) program. The code uses two windows
  (P and S) of 3-component waveforms to determined the moment tensor
    M_ij = M0 * [ sqrt(1-iso*iso)*D_ij  + sqrt(2/3)*iso*I_ij ],
  where I_ij is the unit tensor and D_ij is a deviatoric tensor:
    D_ij = sqrt(1-clvd*clvd)*DC_ij + clvd*CLVD_ij,
  and
    DC_ij   = n_i v_j + n_j v_i,
    CLVD_ij = (2 N_i N_j - v_i v_j - n_i n_j)/sqrt(3),
    (n is the fault normal, v is the slip vector, N=nXv)
    iso = tr(M)/M0/sqrt(6),    -1<= iso <=1,
    clvd = sqrt(3/2)*m2,    -0.5 <= clvd <=0.5,
  where m2 is the intermediate eigenvalue of the normalized deviatoric tensor.

  The solution is given in terms of Mw, strike, dip,
  rake, iso, and clvd. The full moment tensor is also given in the
  output file in the formate of
    M0_in_dyncm Mxx Mxy Mxz Myy Myz Mzz
  where x=North, y=East, z=Down.

  For reference, see
  Zhu and Helmberger, Advancements in source estimation techniques using
      broadband regional seismograms. BSSA, 86, 1634-1641, 1996.
  Zhu and Ben-Zion, Parameterization of general seismic potency and moment
      tensors for source inversion of seismic waveform data. GJI, 2013.
  Zhu and Zhou, Moment tensor solutions of the 2013 Lushan earthquake and its
      aftershocks using 3D crustal-upper mantle velocity models.
      Physics and Chemistry of the Earth. doi:10.1016/j.pce.2016.01.002, 2016.

  Requirements:
	Green's functions -- has P and S arrival time set (t1 and t2 in SAC header)
  Optional:
	Data		 -- P pick (A in SAC header) --> align with t1
			 -- Pnl window (t1 and t2)
			 -- P-SV, SH window (t3 and t4)

  Modify history:
  Nov 24, 2013	Lupei Zhu	modified from cap.c.
  Nov 23, 2014	Lupei Zhu	revised to filter data/GFs before windowing.
  May 05, 2015	Lupei Zhu	delete search for mw. use m0 = sqrt(u2/s2).
  Sept. 18, 2015  LZ		remove using NR subroutines.

  Known bugs:

****************************************************************/
#include "cap3D.h"

int main (int argc, char **argv) {
  int 	i,j,k,k1,l,m,nda,npt,plot,nfm,useDisp,dof,tele,indx,ngf;
  int	ns, mltp, nup, up[3], total_n, n_shft, nqP, nqS;
  int	n1,n2,mm[2],n[NCP],max_shft[NCP],npts[NRC];
  int	repeat, bootstrap, failed;
  char	tmp[128],glib[128],dep[32],loc[32],dst[16],eve[32],*c_pt;
  float	x,x1,x2,y,y1,amp,dt,rad[MAXS],arad[4][3],fm_thr,tie,mtensor[3][3],rec2=0.;
  float	rms_cut[NCP], t0[NCP], tb[NRC], t1, t2, t3, t4, srcDelay;
  float	con_shft[STN], s_shft, shft0[STN][NCP];
  float	tstarP, tstarS, attnP[NFFT], attnS[NFFT];
  float *data[NRC], *green[MAXS][NRC];
  float	bs_body,bs_surf,bs[NCP],weight,nof_per_samp;
  float	w_pnl[NCP];
  float	distance,dmin=100.,vp,vs1,vs2,depSqr=25,vpvs,mu,mw;
  float	*syn,*f_pt,*f_pt0,*f_pt1,*tmpData;
  GRID	grid;
  MTPAR mt[2];	// iso and clvd
  COMP	*spt;
  DATA	*obs, *obs0;
  FM	*fm, *fm0;
  SOLN	sol;
  SACHEAD hd[NRC];
  FILE 	*f_out;
  float tau0, riseTime, *src;
  char type[2] = {'B','P'}, proto[2] = {'B','U'};
  double f1_pnl, f2_pnl, f1_sw, f2_sw;
  float pnl_sn[30], pnl_sd[30], sw_sn[30], sw_sd[30];
  long int order=4, nsects;
#ifdef DIRECTIVITY
  int ns_pnl, ns_sw;
  float *src_pnl, *src_sw;
  float tau, faultStr, faultDip, rupDir, rupV, pVel, sVel, temp;
  scanf("%f%f%f%f%f",&pVel,&sVel,&riseTime,&tau0,&rupDir);
  rupDir *= DEG2RAD;
  rupV = 0.8*sVel;
#endif
  
  if (argc<3) {
     fprintf(stderr,"Usage: %s eventID depthID [elocID] < other_parameters\n",argv[0]);
     return -1;
  } else {
     strcpy(eve,argv[1]);
     strcpy(dep,argv[2]);
  }

  /****** input control parameters *************/
  scanf("%f%f%f%f%d%d%f%f",&x1,&y1,&x,&y,&repeat,&bootstrap,&fm_thr,&tie);
  if (repeat) for(j=0;j<NCP;j++) scanf("%f",rms_cut+4-j);
  scanf("%f%f",&vpvs,&mu);
  scanf("%f%f%f",&vp,&vs1,&vs2);
  scanf("%f%f%f%f",&bs_body,&bs_surf,&x2,&nof_per_samp);
  scanf("%d",&plot);
  scanf("%d%d",&useDisp,&mltp);
  scanf("%s",glib);
  scanf("%d",&ngf);
  fprintf(stderr,"You are using %d fundamental sources\n",ngf);
  for(i=0;i<NRC;i++) ng[i]=ngf;
  if (ngf<6) ng[0]= 2;	// Transverse only has non-zero DD and DS.
  else strcpy(loc,argv[3]);

  /*** input source functions and filters for pnl and sw ***/
  scanf("%f",&dt);
  if (dt>0.) {
     scanf("%f%f",&tau0,&riseTime);
     if ((src = trap(tau0, riseTime, dt, &ns)) == NULL) {
        fprintf(stderr,"fail to make a trapzoid stf\n");
	return -1;
     }
     srcDelay = 0.;
  } else {
     scanf("%s",tmp); scanf("%f",&riseTime);
     if ((src = read_sac(tmp,hd)) == NULL) {
        fprintf(stderr,"fail to read in source time: %s\n",tmp);
        return -1;
     }
     dt = hd->delta;
     ns = hd->npts;
     srcDelay = -hd->b;
  }
  scanf("%lf%lf%lf%lf",&f1_pnl,&f2_pnl,&f1_sw,&f2_sw);
  if (f1_pnl>0.) design(order, type, proto, 1., 1., f1_pnl, f2_pnl, (double) dt, pnl_sn, pnl_sd, &nsects);
  if (f1_sw>0.)  design(order, type, proto, 1., 1., f1_sw, f2_sw, (double) dt, sw_sn, sw_sd, &nsects);

  /** max. window length, shift, and weight for Pnl portion **/
  mm[0]=rint(x1/dt);
  max_shft[3]=max_shft[4]=2*rint(x/dt);
  w_pnl[3]=w_pnl[4]=x2;
  /** max. window length, shift, and weight for P-SV, SH **/
  mm[1]=rint(y1/dt);
  max_shft[0]=max_shft[1]=max_shft[2]=2*rint(y/dt);
  w_pnl[0]=w_pnl[1]=w_pnl[2]=1.;
  /** and tie of time shifts between SH and P-SV **/

  /** input grid-search range **/
  scanf("%f%f",&(mt[0].par),&(mt[0].dd)); mt[0].min = -1.;  mt[0].max = 1.;	// iso
  scanf("%f%f",&(mt[1].par),&(mt[1].dd)); mt[1].min = -0.5; mt[1].max = 0.5;	// clvd
  for(j=0;j<3;j++) {
    scanf("%f%f%f",&x1,&x2,&grid.step[j]);
    grid.n[j] = rint((x2-x1)/grid.step[j]) + 1;
    grid.x0[j] = x1;
  }
  grid.err = (float *) malloc(grid.n[0]*grid.n[1]*grid.n[2]*sizeof(float));
  if (grid.err == NULL ) {
     fprintf(stderr,"fail to allocate memory for storing misfit errors\n");
     return -1;
  }

#ifdef DIRECTIVITY
  faultStr = grid.x0[0]*DEG2RAD;
  faultDip = grid.x0[1]*DEG2RAD;
#endif

  /** input number of stations **/
  scanf("%d",&nda);
  if (nda > STN) {
     fprintf(stderr,"number of station, %d, exceeds max., some stations are discarded\n",nda);
     nda = STN;
  }
  obs = obs0 = (DATA *) malloc(nda*sizeof(DATA));
  fm = fm0 = (FM *) malloc(3*nda*sizeof(FM));
  if (obs == NULL || fm == NULL) {
     fprintf(stderr,"fail to allocate memory for data\n");
     return -1;
  }
  
  /**** loop over stations *****/
  total_n = 0;
  n_shft = 0;
  nfm = 0;
  for(i=0;i<nda;i++) {

    /***** input station name and weighting factor ******/
    failed = 0;
    scanf("%s%s",tmp,dst);
    nup = sscanf(tmp,"%[^/]/%d/%d/%d",obs->stn,&up[0],&up[1],&up[2]);
    if ( fm_thr > 1 ) nup = 1;
    for(k=0,j=0;j<NCP;j++) {
       scanf("%d",&obs->com[4-j].on_off);
       k += obs->com[4-j].on_off;
    }
    scanf("%f%f",&x1,&s_shft);
    tele = 0;
    bs[0] = bs[1] = bs[2] = bs_surf;
    bs[3] = bs[4] = bs_body;
    if (obs->com[3].on_off<0) {
       tele = 1;
       tstarS = obs->com[1].on_off;
       tstarP = obs->com[2].on_off;
       obs->com[1].on_off = obs->com[2].on_off = obs->com[3].on_off = 0;
       k = obs->com[0].on_off + obs->com[4].on_off;
       bs[0] = bs[1] = bs[2] = bs_body;
       j = NFFT;
       if (tstarP>0.) fttq_(&dt, &tstarP, &j, &nqP, attnP);
       if (tstarS>0.) fttq_(&dt, &tstarS, &j, &nqS, attnS);
    }
    if (k==0) {failed++;goto SKIP;}

    /**************input waveforms************/
    strcat(strcat(strcat(strcpy(tmp,eve),"/"),obs->stn), ".t");
    c_pt = strrchr(tmp,(int) 't');
    for(npt=0,j=0;j<NRC;j++){
      *c_pt = cm[j];
      if ((data[j] = read_sac(tmp,&hd[j])) == NULL) {failed++;goto SKIP;}
      tb[j] = hd[j].b-hd[j].o;
      npts[j] = hd[j].npts;
      if (npts[j] > npt) npt=npts[j];
    }
    obs->az = hd->az;
    obs->dist = distance = hd->dist;
    obs->tele = tele;
    if (x1<=0.) x1 = hd[2].a;
    x1 -= hd[2].o;
    if (tele && s_shft>0.) s_shft -= hd[0].o;
    t1 = hd[2].t1-hd[2].o;
    t2 = hd[2].t2-hd[2].o;
    t3 = hd[0].t3-hd[0].o;
    t4 = hd[0].t4-hd[0].o;

    /**************compute source time function***********/
#ifdef DIRECTIVITY
    temp = hd->az*DEG2RAD-faultStr;
    temp = rupV*cos(temp)*cos(rupDir)-sin(temp)*sin(rupDir)*cos(faultDip);
    tau = tau0*(1-temp/pVel);
    src_pnl = trap(tau, riseTime, dt, &ns_pnl);
    tau = tau0*(1-temp/sVel);
    src_sw  = trap(tau, riseTime, dt, &ns_sw);
    if (src_pnl == NULL || src_sw == NULL) {
       fprintf(stderr, "failed to make src for pnl or sw\n");
       return -1;
    }
    fprintf(stderr,"station %s %5.1f tauS %5.1f\n",obs->stn,hd->az,tau);
#endif
    
    /************input green's functions***********/
    float tp = -12345., ts = -12345.;
    if (ngf == 6) sprintf(tmp,"%s/%s/%s/%s.mt.a",glib,dep,loc,obs->stn);
    else          sprintf(tmp,"%s/%s/%s.grn.0",glib,dep,dst);
    c_pt = tmp + strlen(tmp) - 1;
    for(j=0;j<NRC;j++) {
       (*c_pt) = '2'-j; if (ngf==6) (*c_pt) += 49;
       for(k=0;k<ng[j];k++) {
         if (*c_pt == '2') (*c_pt) += 3;
         if ((green[k][j] = read_sac(tmp,&hd[j])) == NULL) {failed++;goto SKIP;}
         if (hd[j].t1>0.) tp=hd[j].t1;
         if (hd[j].t2>tp) ts=hd[j].t2;
	 if (hd[j].npts>npt) npt=hd[j].npts;
	 if (*c_pt=='0' || *c_pt=='5' || *c_pt=='r') obs->alpha = hd[j].user1;
         conv(src, ns, green[k][j], hd[j].npts);
         if (tele) {
            if (tstarP>0. && j==2) conv(attnP, nqP, green[k][j], hd[j].npts);
            if (tstarS>0. && j==0) conv(attnS, nqS, green[k][j], hd[j].npts);
         }
         if (*c_pt == '6' || *c_pt == '7') (*c_pt) += 40;
         (*c_pt) += 3;
       }
    }
    if (tp<0. || ts<0.) {fprintf(stderr,"No P and S arrival times in GFs %s\n", tmp); failed++; goto SKIP;}

    /* generate first-motion polarity data */
    if (nup>1 && (hd[2].user1<0. || hd[0].user2<0.)) {
      fprintf(stderr,"No P/S take-off angle in Greens' function %s\n",tmp);
    } else {
      for(j=1;j<nup;j++) {
        fm->type = up[j-1];
        fm->az = obs->az;
        if (abs(fm->type)==1)	fm->alpha = hd[2].user1;
	else 			fm->alpha = hd[0].user2;
        nfm++;
        fm++;
      }
    }

    /*** calculate time shift needed to align data and GF approximately ****/
    /* positive shift means GF is earlier */
    con_shft[i] = -srcDelay;
    if ( x1 > 0.) {			/* if first-arrival is specified */
       con_shft[i] += x1 - tp;		/* use it to align with greens' fn*/
    }
    if (tele && s_shft > x1 ) {
       s_shft -= ts+con_shft[i];	/* align teleseismic S */
    }

    /** calculate time windows for Pnl and Surface wave portions **/

    /* for Pnl portion */
    if (t1 < 0 || t2 < 0 ) {	/* no time window in the data trace. use default time window in GF */
      if (!tele && vp>0.)
	 t1 = sqrt(distance*distance+depSqr)/vp;	/* use vp to compute t1 */
      else
	 t1 = tp;					/* use tp as t1 */
      t1 = t1 - 0.3*mm[0]*dt + con_shft[i];
      //t2 = ts + 0.2*mm[0]*dt + con_shft[i];		/* ts plus some delay */
      t2 = ts + con_shft[i];	/* ts */
      if (t2-t1<2/f2_pnl+0.3*mm[0]*dt) t2=t1+2/f2_pnl+0.3*mm[0]*dt;	/* if ts-tp is smaller than 2 periods */
    }

    /* do the same for the s/surface wave portion */
    if (t3 < 0 || t4 < 0 ) {
      if (!tele && vs1>0. && vs2> 0.) {
	 t3 = sqrt(distance*distance+depSqr)/vs1 - 0.3*mm[1]*dt;
	 t4 = sqrt(distance*distance+depSqr)/vs2 + 0.7*mm[1]*dt;
      }
      else {
         t3 = ts - 0.3*mm[1]*dt;
         t4 = t3+mm[1]*dt;
      }
      t3 += con_shft[i] + s_shft;
      t4 += con_shft[i] + s_shft;
    }

    /*calculate the time windows */
    n1 = rint((t2 - t1)/dt);	/*Pnl*/
    n2 = rint((t4 - t3)/dt);	/*PSV/SH*/
    if (n1>mm[0]) n1=mm[0];
    if (n2>mm[1]) n2=mm[1];

    /***window data+Greens, do correlation and L2 norms **/
    t0[0]=t3;			/* love wave */
    t0[1]=t0[2]=t4-n2*dt;	/* rayleigh wave */
    t0[3]=t0[4]=t1;		/* Pnl */
    n[0]=n[1]=n[2]=n2;	n[3]=n[4]=n1;
    shft0[i][0] = shft0[i][1] = shft0[i][2] = s_shft;
    shft0[i][3] = shft0[i][4] = 0.;
    if (obs->com[0].on_off>0) n_shft++;
    if (obs->com[1].on_off>0 || obs->com[2].on_off>0) n_shft++;
    if (obs->com[3].on_off>0 || obs->com[4].on_off>0) n_shft++;
    tmpData = (float *) malloc(npt*sizeof(float));
    for(spt=obs->com,j=0;j<NCP;j++,spt++) {
      //if (spt->on_off == 0) continue;
      indx  = kd[j];
      spt->npt = npt = n[j];
      spt->b = t0[j];
      weight = w_pnl[j]*pow(distance/dmin,bs[j]);
      memcpy(tmpData, data[indx], sizeof(float)*npts[indx]);
      if (j<3) {if (f1_sw>0.)  apply(tmpData,(long int) npts[indx],0,sw_sn,sw_sd,nsects);}
      else     {if (f1_pnl>0.) apply(tmpData,(long int) npts[indx],0,pnl_sn,pnl_sd,nsects);}
      if (useDisp==1) cumsum(tmpData, npts[indx], dt); /*use displacement data*/
      f_pt = cutTrace(tmpData, npts[indx], rint((t0[j]-tb[indx])/dt), npt);
      if ( f_pt == NULL ) {
         fprintf(stderr, "fail to window the data\n");
	 failed++; goto SKIP;
      }
      total_n += npt;
      spt->rec = f_pt;
      for(x2=0.,l=0;l<npt;l++,f_pt++) {
	*f_pt *= weight;
	x2+=(*f_pt)*(*f_pt);
      }
      spt->rec2 = x2;
      rec2 += x2;
      for(m=0,k=0;k<ng[indx];k++) {
	memcpy(tmpData, green[k][indx], sizeof(float)*hd[indx].npts);
	if (j<3) {
#ifdef DIRECTIVITY
		conv(src_sw, ns_sw, tmpData, hd[indx].npts);
#endif
		if (f1_sw>0.)  apply(tmpData,(long int) hd[indx].npts,0,sw_sn,sw_sd,nsects);
	} else {
#ifdef DIRECTIVITY
		conv(src_pnl, ns_pnl, tmpData, hd[indx].npts);
#endif
		if (f1_pnl>0.) apply(tmpData,(long int) hd[indx].npts,0,pnl_sn,pnl_sd,nsects);
	}
	if (useDisp) cumsum(tmpData, hd[indx].npts, dt);
	f_pt = cutTrace(tmpData, hd[indx].npts, rint((t0[j]-con_shft[i]-shft0[i][j]-hd[indx].b)/dt), npt);
	if ( f_pt == NULL ) {
	   fprintf(stderr, "fail to window the Greens functions\n");
	   failed++; goto SKIP;
	}
	spt->grn[k] = f_pt;
	for(l=0;l<npt;l++) f_pt[l] *= weight;
	spt->crl[k] = crscrl(npt,spt->rec,f_pt,max_shft[j]);
	for(x=1.,k1=k;k1>=0;k1--,x=2.) {
	  f_pt0=spt->grn[k];
	  f_pt1=spt->grn[k1];
	  for(x2=0.,l=0;l<npt;l++) x2+=(*f_pt0++)*(*f_pt1++);
	  spt->grn2[m++] = x*x2;
	}
      }
      //fprintf(stderr, "%s %d %e\n",obs->stn, j, spt->rec2);
    }

    SKIP:
    if (failed) {/* skip this station */
       fprintf(stderr,"skip this statin ... %s %d\n", obs->stn, failed);
       nda--; i--;
    } else {
       obs++;
    }
    if (!failed || failed>2) for(j=0;j<NRC;j++) free(data[j]);
    if (!failed || failed>3) {free(tmpData); for(j=0;j<ngf;j++) for(k=0;k<NRC;k++) free(green[j][k]);}

  }	/*********end of loop over stations ********/

  if (nda < 1) {
    fprintf(stderr,"No station available for inversion\n");
    return -1;
  }

  /************grid-search for full moment tensor***********/
  INVERSION:
  sol = error(2,ngf,nda,obs0,nfm,fm0,fm_thr,max_shft,tie,mt,grid,0,bootstrap);
  amp = sol.m0;
  mw = (log(amp)/log(10.)+20-16.1)/1.5;

  dof = nof_per_samp*total_n - n_shft;
  x2 = sol.err/dof;		/* data variance */
  /* repeat grid search if needed */
  if ( repeat && discard_bad_data(nda,obs0,sol,x2,rms_cut) ) {
    repeat--;
    goto INVERSION;
  }

  /**************output the results***********************/
  if (sol.flag) fprintf(stderr,"Warning: flag=%d => the minimum %5.1f/%4.1f/%5.1f is at boundary\n",sol.flag,sol.meca.stk,sol.meca.dip,sol.meca.rak);
  for(i=0; i<3; i++) rad[i] = sqrt(2*x2/fabs(sol.dev[i]));
  if (sol.meca.dip>90.) {
     fprintf(stderr,"Warning: dip corrected by %f\n",sol.meca.dip-90);
     sol.meca.dip = 90.;
  }
  //grid.n[0]=grid.n[1]=grid.n[2]=1;
  //grid.x0[0]=sol.meca.stk; grid.x0[1]=sol.meca.dip; grid.x0[2]=sol.meca.rak;
  //sol = error(nda,obs0,nfm,fm0,max_shft,m0,grid,fm_thr,tie);
  strcat(strcat(strcat(strcpy(tmp,eve),"/"),dep),".out");
  f_out=fopen(tmp,"w");
  fprintf(f_out,"Event %s Model %s FM %3d %2d %3d Mw %4.2f E %9.3e %5d ERR %3d %3d %3d ISO %3.2f %3.2f CLVD %3.2f %3.2f\n",eve,dep,
	(int) rint(sol.meca.stk), (int) rint(sol.meca.dip), (int) rint(sol.meca.rak),
	mw, sol.err, dof,
	(int) rint(rad[0]), (int) rint(rad[1]), (int) rint(rad[2]),
        mt[0].par, sqrt(mt[0].sigma*x2),mt[1].par, sqrt(mt[1].sigma*x2));
  fprintf(f_out,"# Variance reduction %4.1f\n",100*(1.-sol.err/rec2));
  nmtensor(mt[0].par,mt[1].par,sol.meca.stk,sol.meca.dip,sol.meca.rak,mtensor);
  fprintf(f_out,"# MomentTensor = %8.3e %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",amp*1.0e13,mtensor[0][0],mtensor[0][1],mtensor[0][2],mtensor[1][1],mtensor[1][2],mtensor[2][2]);
  if (vpvs>0.) {
     x1 = 1.5*vpvs*vpvs-2.;	// eta = (1+nu)/(1-2*nu) = 3/2*vpvs^2 - 2
     x1 = x1*x1;
     y1 = mt[0].par/sqrt(x1+(1-x1)*mt[0].par*mt[0].par);
     fprintf(f_out,"# PotencyTensor ISO = %5.2f",y1);
     if (mu>0.) fprintf(f_out,"   P0 = %8.3e m3\n",amp*1.0e13/mu/sqrt(1-(1-x1)*y1*y1));
     else fprintf(f_out,"\n");
  }
  for(i=1;i<sol.ms;i++) {
     j = sol.others[i];
     if (grid.err[j]-grid.err[sol.others[0]]<mltp*x2) {
	k = j/(grid.n[0]*grid.n[1]);
	k1 = (j-k*grid.n[0]*grid.n[1])/grid.n[0];
	fprintf(f_out,"# %3d %2d %3d %4.2f %9.3e %3.1f\n",
		(int) rint(grid.x0[0]+(j-k1*grid.n[0]-k*grid.n[0]*grid.n[1])*grid.step[0]),
		(int) rint(grid.x0[1]+k1*grid.step[1]),
		(int) rint(grid.x0[2]+k*grid.step[2]),
		mw,grid.err[j],(grid.err[j]-grid.err[sol.others[0]])/x2);
     }
  } 
  for(obs=obs0,i=0;i<nda;i++,obs++) {
    fprintf(f_out,"%-9s %5.1f/%-5.2f",obs->stn, obs->dist, con_shft[i]);
    for(j=0;j<NCP;j++) {
      k = NCP - 1 - j;
      sol.shft[i][k]=sol.shft[i][k] - max_shft[k]/2;
      l = rint(100.*sol.cfg[i][k]); if (l<0) l = 0;
      fprintf(f_out," %1d %8.2e %2d %5.2f",obs->com[k].on_off,sol.error[i][k],l,shft0[i][k]+dt*sol.shft[i][k]);
    }
    fprintf(f_out,"\n");
  }
  fclose(f_out);

  if ( ! plot ) return 0;

  /***** output waveforms for both data and synthetics ****/
  i = mm[1]; if(mm[0]>i) i=mm[0];
  syn = (float *) malloc(i*sizeof(float));
  if (syn == NULL) {
     fprintf(stderr,"fail to allocate memory for output waveforms\n");
     return -1;
  }
  if (ngf==6) for(k=0,i=0;i<3;i++) for(j=i;j<3;j++,k++) rad[k]=mtensor[i][j];
  f_pt = rad;
  for(obs=obs0,i=0;i<nda;i++,obs++){
    if (ngf<6) {
       mt_radiat(obs->az, mtensor, arad);
       for(k=0;k<4;k++) rad[k]=arad[k][0];	// P-SV
       for(   ;k<6;k++) rad[k]=arad[k-3][2];	// only DS and SS have no-zero SH
       f_pt = rad+4;
    }
    strcat(strcat(strcat(strcat(strcat(strcpy(tmp,eve),"/"),dep), "_"),obs->stn),".0");
    c_pt = strrchr(tmp,(int) '0');
    for(spt=obs->com,j=0;j<NCP;j++,spt++,f_pt=rad) {
      //if (spt->on_off == 0) {(*c_pt)+=2; continue;}
      npt=spt->npt;
      hd[0] = sachdr(dt, npt, spt->b);
      hd->dist = obs->dist; hd->az = obs->az; hd->user1 = obs->alpha;
      hd->a = hd->b;
      for(l=0;l<npt;l++) syn[l] = spt->rec[l]/amp;
      write_sac(tmp,hd[0],syn);
      (*c_pt)++;
      for(l=0;l<npt;l++) {
	for(x2=0.,k=0;k<ng[kd[j]];k++) x2 += f_pt[k]*spt->grn[k][l];
	syn[l] = sol.scl[i][j]*x2;
      }
      hd->b -= (shft0[i][j]+con_shft[i]);
      hd->a = hd->b-sol.shft[i][j]*dt;
      write_sac(tmp,hd[0],syn);
      (*c_pt)++;
    }
  }
  return 0;
}


// grid-search for the full moment tensor
SOLN	error(	int		npar,	// 2=iso; 1=clvd; 0=strike/dip/rake
		int		ngf,
		int		nda,
		DATA		*obs0,
		int		nfm,
		FM		*fm,
		float		fm_thr,
		const int	*max_shft,
		float		tie,
		MTPAR		*mt,
		GRID		grid,
		int		interp,
		int		bootstrap
	) {
  int	i, j, k, l, m, k1, z0, z1, z2, debug=0;
  int	i_stk, i_dip, i_rak;
  float	amp, rad[MAXS], arad[4][3], x, x1, x2, y, y1, y2, s3d[9],  u2 ,s2;
  float	*fpt, *frad, **f_pt, **r_pt, **z_pt, *grd_err;
  float dx, mtensor[3][3];
  DATA	*obs;
  COMP	*spt;
  SOLN	sol, sol1, sol2, best_sol;

  if ( npar ) {	// search for mw, iso, and clvd ================================

  npar--;
  dx = mt[npar].dd;
  mt[npar].sigma = 0.;
  if (bootstrap) {	// do full-grid search for bootstrapping
     if (dx<0.001) {mt[npar].min=mt[npar].max=mt[npar].par; dx=1.;}
     sol.err = FLT_MAX;
     for(mt[npar].par = mt[npar].min; mt[npar].par<=mt[npar].max; mt[npar].par+=dx) {
        sol1 = error(npar,ngf,nda,obs0,nfm,fm,fm_thr,max_shft,tie,mt,grid,0,bootstrap);
	if (sol1.err<sol.err) {
	   sol = sol1;
	   x = mt[npar].par;
	}
     }
     mt[npar].par = x;
  } else {		// do line search for efficiency
     i = 1; if (dx>0.001) i = 0;
     sol = error(npar,ngf,nda,obs0,nfm,fm,fm_thr,max_shft,tie,mt,grid,i,bootstrap);
     if (dx>0.001) {
        mt[npar].par += dx;
        sol2 = error(npar,ngf,nda,obs0,nfm,fm,fm_thr,max_shft,tie,mt,grid,0,bootstrap);
        if (sol2.err > sol.err) {	/* this is the wrong direction, turn around */
           dx = -dx;
           sol1 = sol2; sol2 = sol; sol  = sol1; /*swap sol, sol2 */
           mt[npar].par += dx;
        }
        while(sol2.err < sol.err) {	/* keep going until passing by the mininum */
           sol1 = sol;
           sol = sol2;
           mt[npar].par += dx;
           if (mt[npar].par>mt[npar].max || mt[npar].par<mt[npar].min) sol2.err = sol1.err;
           else sol2 = error(npar,ngf,nda,obs0,nfm,fm,fm_thr,max_shft,tie,mt,grid,0,bootstrap);
        }
        mt[npar].sigma = 2*dx*dx/(sol2.err+sol1.err-2*sol.err);
        mt[npar].par -= dx+0.5*dx*(sol2.err-sol1.err)/(sol2.err+sol1.err-2*sol.err);
        sol = error(npar,ngf,nda,obs0,nfm,fm,fm_thr,max_shft,tie,mt,grid,1,bootstrap);
     }
  }
  return(sol);

  } else {	// the base case: grid-search for strike, dip, and rake ========

  best_sol.err = FLT_MAX;
  grd_err = grid.err;
  for(i_rak=0; i_rak<grid.n[2]; i_rak++) {
     sol.meca.rak=grid.x0[2]+i_rak*grid.step[2];
     for(i_dip=0; i_dip<grid.n[1]; i_dip++) {
       sol.meca.dip=grid.x0[1]+i_dip*grid.step[1];
       for(i_stk=0; i_stk<grid.n[0]; i_stk++) {
          sol.meca.stk=grid.x0[0]+i_stk*grid.step[1];
          nmtensor(mt[0].par,mt[1].par,sol.meca.stk,sol.meca.dip,sol.meca.rak,mtensor);
	  if (check_first_motion(mtensor,fm,nfm,fm_thr)<0) {
		*grd_err++ = sol.err = FLT_MAX;
		continue;
	  }
          if (bootstrap) fprintf(stderr,"BOOTSTRAPPING grid %5.2f %5.2f %5.1f %5.1f %5.1f\n", mt[0].par, mt[1].par, sol.meca.stk, sol.meca.dip, sol.meca.rak);
          if (ngf==6) for(k=0,i=0;i<3;i++) for(j=i;j<3;j++,k++) rad[k]=mtensor[i][j];
	  frad = rad;
	  for(obs=obs0,u2=s2=0.,i=0;i<nda;i++,obs++){
	    
	    if (ngf<6) {
		mt_radiat(obs->az,mtensor,arad);
		for(k=0;k<4;k++) rad[k]=arad[k][0];
		for(   ;k<6;k++) rad[k]=arad[k-3][2];
		frad = rad+4;	// for SH 
	    }

	    /*******find the time shift*************/
	    /**SH surface wave**/
	    spt = obs->com;
	    f_pt = spt->crl;
	    z0 = spt->on_off>0?1:0;
	    /**PSV surface wave**/
	    spt++;
	    r_pt = spt->crl;
	    z1 = spt->on_off>0?1:0;
	    spt++;
	    z_pt = spt->crl;
	    z2 = spt->on_off>0?1:0;
	    for(y1=y2=-FLT_MAX,l=0;l<=max_shft[1];l++) {
	      x =0.; if (z0) for(k=0;k<ng[0];k++) x+=frad[k]*f_pt[k][l];
	      x1=0.; if (z1) for(k=0;k<ngf;k++) x1 += rad[k]*r_pt[k][l];
	      x2=0.; if (z2) for(k=0;k<ngf;k++) x2 += rad[k]*z_pt[k][l];
	      y = (1-tie)*z0*x + tie*(z1*x1+z2*x2);
	      if (y>y2) {y2=y;sol.cfg[i][0]=x;sol.shft[i][0]=l;}
	      y = tie*z0*x + (1-tie)*(z1*x1+z2*x2);
	      if (y>y1) {y1=y;sol.cfg[i][1]=x1;sol.cfg[i][2]=x2;m=l;}
	    }
	    sol.shft[i][1]=sol.shft[i][2]=m;
	    /**Pnl*/
	    spt++;
	    r_pt = spt->crl;
	    z1 = spt->on_off>0?1:0;
	    spt++;
	    z_pt = spt->crl;
	    z2 = spt->on_off>0?1:0;
	    for(y1=-FLT_MAX,l=0;l<=max_shft[3];l++) {
	      x1=0.; if (z1) for(k=0;k<ngf;k++) x1 += rad[k]*r_pt[k][l];
	      x2=0.; if (z2) for(k=0;k<ngf;k++) x2 += rad[k]*z_pt[k][l];
	      y = z1*x1 + z2*x2;
	      if (y>y1) {y1=y;sol.cfg[i][3]=x1;sol.cfg[i][4]=x2;m=l;}
	    }
	    sol.shft[i][3]=sol.shft[i][4]=m;
	    
	    /* compute the L2 norm of synthetics */
	    for(spt=obs->com,j=0;j<NCP;j++,spt++,frad=rad) {
	      for(x2=0.,fpt=spt->grn2,k=0;k<ng[kd[j]];k++) {
		for(k1=k;k1>=0;k1--)
		  x2+=frad[k]*frad[k1]*(*fpt++);
	      }
	      spt->syn2 = x2;
	      if (spt->on_off) {u2 += spt->rec2; s2 += x2;}
	    }
          }
	  sol.m0 = amp = sqrt(u2/s2);

	  /***error calculation*****/
	  for(obs=obs0,sol.err=0.,i=0;i<nda;i++,obs++){
	    if (ngf<6) {
		mt_radiat(obs->az,mtensor,arad);
		for(k=0;k<4;k++) rad[k]=arad[k][0];
		for(   ;k<6;k++) rad[k]=arad[k-3][2];
		frad = rad+4;	// for SH 
	    }
	    spt = obs->com;
	    for(j=0;j<NCP;j++,spt++,frad=rad) {
	      y1 = 1.;
	      /* find out the scaling factor for teleseismic distances */
	      if (obs->tele) {
	         if (sol.cfg[i][j]>0.) y1 = sol.cfg[i][j]/spt->syn2;
		 else y1 = 0.;
	      }
	      sol.scl[i][j] = y1;
	      y1 *= amp;

	      x1 = spt->rec2+spt->syn2*y1*y1-2.*sol.cfg[i][j]*y1;
	      sol.error[i][j] = x1;	/*L2 error for this com.*/
	      sol.cfg[i][j] /= sqrt(spt->rec2*spt->syn2);
	      sol.err += spt->on_off*x1;
	      //if (bootstrap) fprintf(stderr,"BOOTSTRAPPING %-10s %d %d %9.3e\n", obs->stn, j, spt->on_off, x1);

	    }

	  } /*-------------------------end of all stations*/
	  *grd_err++ = sol.err;		/*error for this solution*/
	  if (bootstrap) fprintf(stderr,"BOOTSTRAPPING chi2 %9.3e\n", sol.err);
	  if (best_sol.err>sol.err) best_sol = sol;

       }
    }
  }
  if (debug) fprintf(stderr, "iso=%5.2f clvd=%5.2f misfit = %9.3e\n", mt[0].par, mt[1].par, best_sol.err);
  for(i=0;i<6;i++) best_sol.dev[i]  = 1.;
  if (bootstrap || interp == 0) return(best_sol);
  /* do interpolation */
  best_sol.err = grid3d(grid.err,&(grid.n[0]),s3d,&(best_sol.flag),&(best_sol.ms),best_sol.others);
  if (debug) fprintf(stderr, " interpolation  misfit = %9.3e\n", best_sol.err);
  best_sol.meca.stk = grid.x0[0]+s3d[0]*grid.step[0];
  best_sol.meca.dip = grid.x0[1]+s3d[1]*grid.step[1];
  best_sol.meca.rak = grid.x0[2]+s3d[2]*grid.step[2];
  for(i=0;i<3;i++) best_sol.dev[i]  = s3d[3+i]/(grid.step[i]*grid.step[i]);
  best_sol.dev[3] = s3d[6]/(grid.step[0]*grid.step[1]);
  best_sol.dev[4] = s3d[7]/(grid.step[0]*grid.step[2]);
  best_sol.dev[5] = s3d[8]/(grid.step[1]*grid.step[2]);
  return(best_sol);

  }
}
