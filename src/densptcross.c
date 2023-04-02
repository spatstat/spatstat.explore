#include <Rmath.h>
#include <R_ext/Utils.h>

#include "chunkloop.h"
#include "crossloop.h"
#include "constants.h"
/*

  densptcross.c

  $Revision: 1.6 $     $Date: 2023/04/02 00:18:44 $

  Assumes point patterns are sorted in increasing order of x coordinate

  *crdenspt     Density estimate at points
  *crsmoopt     Smoothed mark values at points

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#define TWOPI M_2PI

double sqrt(double x);
double exp(double x);

#define STD_DECLARATIONS				\
  int i, j, n1, n2, maxchunk, jleft;                    \
  double x1i, y1i, xleft, dx, dy, d2, rmax, r2max;      \
  double *x1, *y1, *x2, *y2;

#define STD_INITIALISE				\
  n1 = *nquery;					\
  x1 = xq; y1 = yq;                             \
  n2 = *ndata;					\
  x2 = xd; y2 = yd;                             \
  rmax = *rmaxi;				\
  r2max = rmax * rmax


/* ----------------- density estimation -------------------- */

void crdenspt(
     /* inputs */
     int *nquery,             /* number of locations to be interrogated */
     double *xq, double *yq,  /* (x,y) coordinates to be interrogated */
     int *ndata,              /* number of data points */
     double *xd, double *yd,  /* (x,y) coordinates of data */
     double *rmaxi,           /* maximum distance at which points contribute */
     double *sig,             /* Gaussian sd */
     int *squared,            /* whether to use the squared kernel */
     /* output */
     double *result           /* vector of computed density values */
) {
  STD_DECLARATIONS;
  double resulti, coef;	
  double sigma, a; 
  STD_INITIALISE;

  if(n1 == 0 || n2 == 0) 
    return;

  sigma = *sig;				      
  a = 1.0/(2.0 * sigma * sigma);	
  coef = 1.0/(TWOPI * sigma * sigma);  

  if(*squared) {
    coef *= coef;
    a *= 2.0;
  }
  
  CROSSLOOP( { resulti = 0.0; },
            { resulti += exp(-d2 * a); } ,
	    { result[i] = coef * resulti; })

}


void wtcrdenspt(
     /* inputs */
     int *nquery,             /* number of locations to be interrogated */
     double *xq, double *yq,  /* (x,y) coordinates to be interrogated */
     int *ndata,              /* number of data points */
     double *xd, double *yd,  /* (x,y) coordinates of data */
     double *wd,              /* weights of data points */
     double *rmaxi,           /* maximum distance at which points contribute */
     double *sig,             /* Gaussian sd */
     int *squared,            /* whether to use the squared kernel */
     /* output */
     double *result           /* vector of computed density values */
) {
  STD_DECLARATIONS;
  double resulti, coef;	
  double sigma, a; 
  STD_INITIALISE;

  if(n1 == 0 || n2 == 0) 
    return;

  sigma = *sig;				      
  a = 1.0/(2.0 * sigma * sigma);	
  coef = 1.0/(TWOPI * sigma * sigma);  

  if(*squared) {
    coef *= coef;
    a *= 2.0;
  }
  
  CROSSLOOP( { resulti = 0.0; },
	    { resulti += wd[j] * exp(-d2 * a); },
	    { result[i] = coef * resulti; } )

 }

/* ------------- anisotropic versions -------------------- */

void acrdenspt(
     /* inputs */
     int *nquery,              /* number of locations to be interrogated */
     double *xq, double *yq,   /* (x,y) coordinates to be interrogated */
     int *ndata,               /* number of data points */
     double *xd, double *yd,   /* (x,y) coordinates of data */
     double *rmaxi,            /* maximum distance at which points contribute */
     double *detsigma,         /* determinant of variance matrix */
     double *sinv,             /* inverse variance matrix (2x2, flattened) */
     int *squared,            /* whether to use the squared kernel */
     /* output */
     double *result            /* vector of computed density values */
) {
  STD_DECLARATIONS;
  double resulti, coef;	
  double detsig, s11, s12, s21, s22;
  STD_INITIALISE;

  if(n1 == 0 || n2 == 0) 
    return;

  detsig = *detsigma;
  coef = 1.0/(TWOPI * sqrt(detsig));
  
  if(*squared) {
    coef *= coef;
    s11 = sinv[0];
    s12 = sinv[1];
    s21 = sinv[2];
    s22 = sinv[3];
  } else {
    s11 = sinv[0]/2.0;
    s12 = sinv[1]/2.0;
    s21 = sinv[2]/2.0;
    s22 = sinv[3]/2.0;
  }

  CROSSLOOP( { resulti = 0.0; },
	    { resulti += exp(- dx * (dx * s11 + dy * s12) \
			     + dy * (dx * s21 + dy * s22)); },
	    { result[i] = coef * resulti; })
}


void awtcrdenspt(
  /* inputs */
  int *nquery,             /* number of locations to be interrogated */
  double *xq, double *yq,  /* (x,y) coordinates to be interrogated */
  int *ndata,              /* number of data points */
  double *xd, double *yd,  /* (x,y) coordinates of data */
  double *wd,              /* weights of data points */
  double *rmaxi,           /* maximum distance at which points contribute */
  double *detsigma,        /* determinant of variance matrix */
  double *sinv,            /* inverse variance matrix (2x2, flattened) */
  int *squared,            /* whether to use the squared kernel */
  /* output */
  double *result           /* vector of weighted density values */
) {
  STD_DECLARATIONS;
  double resulti, coef;	
  double detsig, s11, s12, s21, s22;
  STD_INITIALISE;

  if(n1 == 0 || n2 == 0) 
    return;

  detsig = *detsigma;
  coef = 1.0/(TWOPI * sqrt(detsig));

  if(*squared) {
    coef *= coef;
    s11 = sinv[0];
    s12 = sinv[1];
    s21 = sinv[2];
    s22 = sinv[3];
  } else {
    s11 = sinv[0]/2.0;
    s12 = sinv[1]/2.0;
    s21 = sinv[2]/2.0;
    s22 = sinv[3]/2.0;
  }

  CROSSLOOP( { resulti = 0.0; },
	    { resulti += wd[j] * \
		exp(- dx * (dx * s11 + dy * s12)			\
		    + dy * (dx * s21 + dy * s22)); },
	    { result[i] = coef * resulti; })
 }


/* --------------- smoothing --------------------------- */

void crsmoopt(
  /* inputs */
  int *nquery,             /* number of locations to be interrogated */
  double *xq, double *yq,  /* (x,y) coordinates to be interrogated */
  int *ndata,              /* number of data points */
  double *xd, double *yd,  /* (x,y) coordinates of data */
  double *vd,              /* mark values at data points */
  double *rmaxi,           /* maximum distance at which points contribute */
  double *sig,             /* Gaussian sd */
  /* output */
  double *result           /* vector of computed smoothed values */
) {
  STD_DECLARATIONS;
  double sigma, twosig2;
  double numer, denom, wij; 

  STD_INITIALISE;
  sigma = *sig;
  twosig2 = 2.0 * sigma * sigma;

  if(n1 == 0 || n2 == 0) 
    return;

  CROSSLOOP({ numer = denom = 0.0; },
	   { \
	     wij = exp(-d2/twosig2);		\
	     denom += wij;			\
	     numer += wij * vd[j];		\
	   },					
	   {					\
	     result[i] = numer/denom;		\
	   })
 }


void wtcrsmoopt(
  /* inputs */
  int *nquery,              /* number of locations to be interrogated */
  double *xq, double *yq,   /* (x,y) coordinates to be interrogated */
  int *ndata,               /* number of data points */
  double *xd, double *yd,   /* (x,y) coordinates of data */
  double *vd,               /* mark values at data points */
  double *wd,               /* weights of data points */
  double *rmaxi,            /* maximum distance */
  double *sig,              /* Gaussian sd */
  /* output */
  double *result    /* vector of computed smoothed values */
) {
  STD_DECLARATIONS;
  double sigma, twosig2;
  double numer, denom, wij; 

  STD_INITIALISE;
  sigma = *sig;
  twosig2 = 2.0 * sigma * sigma;

  if(n1 == 0 || n2 == 0) 
    return;

  CROSSLOOP({ numer = denom = 0.0; },
	   {						\
	     wij = wd[j] * exp(-d2/twosig2);	\
	     denom += wij;				\
	     numer += wij * vd[j];			\
	   },						
	   {						\
	     result[i] = numer/denom;			\
	   })
}

/* ------------- anisotropic versions -------------------- */

void acrsmoopt(
  /* inputs */
  int *nquery,             /* number of locations to be interrogated */
  double *xq, double *yq,  /* (x,y) coordinates to be interrogated */
  int *ndata,              /* number of data points */
  double *xd, double *yd,  /* (x,y) coordinates of data */
  double *vd,              /* mark values at data points */
  double *rmaxi,           /* maximum distance at which points contribute */
  double *sinv,            /* inverse variance matrix (2x2, flattened) */
  /* output */
  double *result           /* vector of smoothed values */
) {
  STD_DECLARATIONS;
  double s11, s12, s21, s22;
  double numer, denom, wij; 

  STD_INITIALISE;
  s11 = sinv[0];
  s12 = sinv[1];
  s21 = sinv[2];
  s22 = sinv[3];

  if(n1 == 0 || n2 == 0) 
    return;

  CROSSLOOP({ numer = denom = 0.0; },
	   {							\
	     wij = exp(-(dx * (dx * s11 + dy * s12)		\
			 + dy * (dx * s21 + dy * s22))/2.0);	\
	     denom += wij;					\
	     numer += wij * vd[j];				\
	   },
	   {					\
	     result[i] = numer/denom;		\
	   })
}


void awtcrsmoopt(
  /* inputs */
  int *nquery,              /* number of locations to be interrogated */
  double *xq, double *yq,   /* (x,y) coordinates to be interrogated */
  int *ndata,               /* number of data points */
  double *xd, double *yd,   /* (x,y) coordinates of data */
  double *vd,               /* mark values at data points */
  double *wd,               /* weights of data points */
  double *rmaxi,            /* maximum distance at which points contribute */
  double *sinv,             /* inverse variance matrix (2x2, flattened) */
  /* output */
  double *result            /* vector of smoothed values */
) {
  STD_DECLARATIONS;
  double s11, s12, s21, s22;
  double numer, denom, wij; 

  STD_INITIALISE;

  s11 = sinv[0];
  s12 = sinv[1];
  s21 = sinv[2];
  s22 = sinv[3];

  if(n1 == 0 || n2 == 0) 
    return;

  CROSSLOOP({ numer = denom = 0.0; },
	   {								\
	     wij = wd[j] * exp(-(dx * (dx * s11 + dy * s12)		\
				     + dy * (dx * s21 + dy * s22))/2.0); \
	     denom += wij;						\
	     numer += wij * vd[j];					\
	   },
	   {					\
	     result[i] = numer/denom;		\
	   })
}

