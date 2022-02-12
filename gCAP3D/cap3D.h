/***************************************************************
	cap3D.h		Head file for cap3D.c
***************************************************************/


#ifndef __CAP_HEAD__
  #define __CAP_HEAD__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "sac.h"
#include "Complex.h"
#define NRANSI
#include "inversion.h"
#include "radiats.h"

/***********************Constants********************************/

#define STN	200	// max number of stations

#define MAXS	6	// max number of fundamental sources
			// for 1D FK GFs, they are DD, DS, SS, EX
			// for 3D MT GFs, they are xx, xy, xz, yy, yz, zz

#define NRC	3	// 3 components of records
static char cm[NRC]={'t','r','z'};
int  ng[NRC];		// number of fundamental sources for each component.

#define NCP	5	// 5 segments, 2*Pnl 2*PSV 1*SH
static int kd[NCP]={0,1,2,1,2};		//in the order of SH, SVr, SVz, Pnlr, Pnlz

#define NFFT	2048	// nfft for Q operator in teleseismic cases.

/*********************** Data Structure***************************/

/* focal mechanism (strike, dip, rake) data structure */
typedef struct {
	float	stk;	// strkie
	float	dip;	// dip
	float	rak;	// rake
} MECA;

/* a portion of observed waveform and corresponding
Green's functions, their cross-correlations, and L2 norms, etc */
typedef struct {
	int	on_off;		// on or off
	int	npt;		// number of points
	float	b;		// beginning time of the data
	float	*rec;		// observation
	float	*grn[MAXS];	// Green's functions
	float	rec2, syn2;	// L2 norm of the data (sum rec*rec) and syn.
	float	grn2[21];	// xc of GFs, c*sum grn[i]*grn[j], c=1 if i=j, 2 if j>i
	float	*crl[MAXS];	// xc of obs and GFs, sum rec*grn[i]
} COMP;

/* data */
typedef struct {
	char	stn[10];	// station name
	float	az;		// azimuth in deg
	float	dist;		// epicentral distance in km
	float	alpha;		// take-off angle
	int	tele;		// 1 if the station is at teleseismic distances
	COMP	com[NCP];	// 5 segments
} DATA;

/* solution */
typedef struct {
	MECA	meca;
	float	m0;			/* scalar moment */
	float	dev[6];			/* uncertainty ellipsis */
	float	err;			/* L2-norm of waveform misfits for this solution */
	float	cfg[STN][NCP];		/* correlation of data and syn for each comp. */
	int	shft[STN][NCP];		/* time shift for each comp. */
	float	error[STN][NCP];	/* L2-norm of waveform misfits for each component */
	float   scl[STN][NCP];		/* amplifications to GF for each component */
	int	ms;			/* number of local minimums < 10 */
	int	others[10];		/* top 10 best solutions */
	int	flag;			/* =1 if the best soln is at boundary */
} SOLN;

/* Moment tensor parameter mw, iso, or clvd */
typedef struct {
	float	par;		// can be mw, iso, or cvld
	float	sigma;		// variance
	float	min, max;	// range
	float	dd;		// search step
} MTPAR;

/* first-motion data */
typedef struct {
	float	az;	/* azimuth */
	float	alpha;	/* take-off angle */
	int	type;	/* 1=P; 2=SV; 3=SH; positive=up; negative=down */
} FM;

/* function declaration */
SOLN	error(int,int,int,DATA *,int,FM *,float,const int *,float,MTPAR *,GRID,int,int);
void    taper(float *aa, int n);
float	*trap(float, float, float, int *);
float	*cutTrace(float *, int, int, int);
int	discard_bad_data(int,DATA *,SOLN,float,float *);
int	check_first_motion(float mt[3][3], FM *fm, int n, float fm_thr);
void	fttq_(float *,float *,int *,int *,float *);


#endif
