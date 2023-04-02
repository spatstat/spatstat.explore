#include <Rmath.h>
#include <R_ext/Utils.h>

#include "chunkloop.h"
#include "pairloop.h"
#include "constants.h"
/*

  denspt.c

  Calculation of density estimate at data points

  $Revision: 1.26 $     $Date: 2023/04/02 00:17:21 $

  Assumes point pattern is sorted in increasing order of x coordinate

  *denspt*     Density estimate at points
  *smoopt*     Smoothed mark values at points

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#define TWOPI M_2PI

double sqrt(double x);
double exp(double x);

#define STD_DECLARATIONS				\
  int n, i, j, maxchunk;				\
  double xi, yi, rmax, r2max, dx, dy, dx2, d2	

#define STD_INITIALISE				\
  n = *nxy;					\
  rmax = *rmaxi;				\
  r2max = rmax * rmax


/* ----------------- density estimation -------------------- */

void denspt(
  /* inputs */
  int    *nxy,             /* number of (x,y) points */
  double *x, double *y,    /* (x,y) coordinates */
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

  if(n == 0) 
    return;

  sigma = *sig;				      
  a = 1.0/(2.0 * sigma * sigma);	
  coef = 1.0/(TWOPI * sigma * sigma);  

  if(*squared) {
    coef *= coef;
    a *= 2.0;
  }
  
  PAIRLOOP( { resulti = 0.0; },
            { resulti += exp(-d2 * a); } ,
	    { result[i] = coef * resulti; })

}


void wtdenspt(
     /* inputs */
     int *nxy,                /* number of (x,y) points */
     double *x, double *y,    /* (x,y) coordinates */
     double *rmaxi,           /* maximum distance */
     double *sig,             /* Gaussian sd */
     double *weight,          /* vector of weights */
     int *squared,            /* whether to use the squared kernel */
     /* output */
     double *result           /* vector of weighted density values */
) {
  STD_DECLARATIONS;
  double resulti, coef;	
  double sigma, a; 
  STD_INITIALISE;

  if(n == 0) 
    return;

  sigma = *sig;				      
  a = 1.0/(2.0 * sigma * sigma);	
  coef = 1.0/(TWOPI * sigma * sigma);  

  if(*squared) {
    coef *= coef;
    a *= 2.0;
  }
  
  PAIRLOOP( { resulti = 0.0; },
	    { resulti += weight[j] * exp(-d2 * a); },
	    { result[i] = coef * resulti; } )

 }

/* ------------- anisotropic versions -------------------- */

void adenspt(
     /* inputs */
     int *nxy,                /* number of (x,y) points */
     double *x, double *y,    /* (x,y) coordinates */
     double *rmaxi,           /* maximum distance at which points contribute */
     double *detsigma,        /* determinant of variance matrix */
     double *sinv,            /* inverse variance matrix (2x2, flattened) */
     int *squared,            /* whether to use the squared kernel */
     /* output */
     double *result           /* vector of density values */
) {
  STD_DECLARATIONS;
  double resulti, coef;	
  double detsig, s11, s12, s21, s22;
  STD_INITIALISE;

  if(n == 0)
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

  PAIRLOOP( { resulti = 0.0; },
	    { resulti += exp(-dx * (dx * s11 + dy * s12) \
			     + dy * (dx * s21 + dy * s22)); },
	    { result[i] = coef * resulti; })
}


void awtdenspt(
     /* inputs */
     int *nxy,          /* number of (x,y) points */
     double *x,
     double *y,         /* (x,y) coordinates */
     double *rmaxi,     /* maximum distance at which points contribute */
     double *detsigma,  /* determinant of variance matrix */
     double *sinv,      /* inverse variance matrix (2x2, flattened) */
     double *weight,    /* vector of weights */
     int *squared,      /* whether to use the squared kernel */
     /* output */
     double *result    /* vector of weighted density values */
) {
  STD_DECLARATIONS;
  double resulti, coef;	
  double detsig, s11, s12, s21, s22;
  STD_INITIALISE;
  detsig = *detsigma;
  coef = 1.0/(TWOPI * sqrt(detsig));

  if(n == 0) 
    return;

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

  PAIRLOOP( { resulti = 0.0; },
	    { resulti += weight[j] * \
		exp(-(dx * (dx * s11 + dy * s12)			\
		      + dy * (dx * s21 + dy * s22))); },
	    { result[i] = coef * resulti; })
}


/* --------------- smoothing --------------------------- */

void smoopt(
     /* inputs */
     int *nxy,                /* number of (x,y) points */
     double *x, double *y,    /* (x,y) coordinates */
     double *v,               /* vector of mark values to be smoothed */
     int *self,               /* 0 if leave-one-out */
     double *rmaxi,           /* maximum distance at which points contribute */
     double *sig,             /* Gaussian sd */
     /* output */
     double *result           /* vector of computed smoothed values */
) {
  STD_DECLARATIONS;
  int countself;
  double sigma, twosig2;
  double numer, denom, wij; 

  STD_INITIALISE;
  sigma = *sig;
  countself = *self;
  twosig2 = 2.0 * sigma * sigma;

  if(n == 0) 
    return;

  if(countself != 0) {
  PAIRLOOP({ numer = denom = 0.0; },
	   { \
	     wij = exp(-d2/twosig2);		\
	     denom += wij;			\
	     numer += wij * v[j];		\
	   },					
	   {					\
	     denom += 1;			\
	     numer += v[i];			\
	     result[i] = numer/denom;		\
	   })
    } else {
  PAIRLOOP({ numer = denom = 0.0; },
	   { \
	     wij = exp(-d2/twosig2);		\
	     denom += wij;			\
	     numer += wij * v[j];		\
	   },					
	   {					\
	     result[i] = numer/denom;		\
	   })
    }
 }


void wtsmoopt(
     /* inputs */
     int *nxy,                /* number of (x,y) points */
     double *x, double *y,    /* (x,y) coordinates */
     double *v,               /* vector of mark values to be smoothed */
     int *self,               /* 0 if leave-one-out */
     double *rmaxi,           /* maximum distance */
     double *sig,             /* Gaussian sd */
     double *weight,          /* vector of weights */
     /* output */
     double *result    /* vector of computed smoothed values */
) {
  STD_DECLARATIONS;
  int countself;
  double sigma, twosig2;
  double numer, denom, wij; 

  STD_INITIALISE;
  sigma = *sig;
  countself = *self;
  twosig2 = 2.0 * sigma * sigma;

  if(n == 0) 
    return;

  if(countself != 0) {
  PAIRLOOP({ numer = denom = 0.0; },
	   {						\
	     wij = weight[j] * exp(-d2/twosig2);	\
	     denom += wij;				\
	     numer += wij * v[j];			\
	   },						
	   {						\
	     denom += weight[i];			\
	     numer += weight[i] * v[i];		\
	     result[i] = numer/denom;			\
	   })
  } else {
  PAIRLOOP({ numer = denom = 0.0; },
	   {						\
	     wij = weight[j] * exp(-d2/twosig2);	\
	     denom += wij;				\
	     numer += wij * v[j];			\
	   },						
	   {						\
	     result[i] = numer/denom;			\
	   })
    }
}

/* ------------- anisotropic versions -------------------- */

void asmoopt(
     /* inputs */
     int *nxy,                /* number of (x,y) points */
     double *x, double *y,    /* (x,y) coordinates */
     double *v,               /* vector of mark values to be smoothed */
     int *self,               /* 0 if leave-one-out */
     double *rmaxi,           /* maximum distance at which points contribute */
     double *sinv,            /* inverse variance matrix (2x2, flattened) */
     /* output */
     double *result           /* vector of smoothed values */
) {
  STD_DECLARATIONS;
  int countself;
  double s11, s12, s21, s22;
  double numer, denom, wij; 

  STD_INITIALISE;
  countself = *self;
  s11 = sinv[0];
  s12 = sinv[1];
  s21 = sinv[2];
  s22 = sinv[3];

  if(n == 0) 
    return;

  if(countself != 0) {
  PAIRLOOP({ numer = denom = 0.0; },
	   {							\
	     wij = exp(-(dx * (dx * s11 + dy * s12)		\
			 + dy * (dx * s21 + dy * s22))/2.0);	\
	     denom += wij;					\
	     numer += wij * v[j];				\
	   },
	   {					\
	     denom += 1;			\
	     numer += v[i];			\
	     result[i] = numer/denom;		\
	   })
    } else {
  PAIRLOOP({ numer = denom = 0.0; },
	   {							\
	     wij = exp(-(dx * (dx * s11 + dy * s12)		\
			 + dy * (dx * s21 + dy * s22))/2.0);	\
	     denom += wij;					\
	     numer += wij * v[j];				\
	   },
	   {					\
	     result[i] = numer/denom;		\
	   })
    }
}


void awtsmoopt(
     /* inputs */
     int *nxy,                /* number of (x,y) points */
     double *x, double *y,    /* (x,y) coordinates */
     double *v,               /* vector of mark values to be smoothed */
     int *self,               /* 0 if leave-one-out */
     double *rmaxi,           /* maximum distance at which points contribute */
     double *sinv,            /* inverse variance matrix (2x2, flattened) */
     double *weight,          /* vector of weights */
     /* output */
     double *result           /* vector of smoothed values */
) {
  STD_DECLARATIONS;
  int countself;
  double s11, s12, s21, s22;
  double numer, denom, wij; 

  STD_INITIALISE;
  countself = *self;

  s11 = sinv[0];
  s12 = sinv[1];
  s21 = sinv[2];
  s22 = sinv[3];

  if(n == 0) 
    return;

  if(countself != 0) {
  PAIRLOOP({ numer = denom = 0.0; },
	   {								\
	     wij = weight[j] * exp(-(dx * (dx * s11 + dy * s12)		\
				     + dy * (dx * s21 + dy * s22))/2.0); \
	     denom += wij;						\
	     numer += wij * v[j];					\
	   },
	   {					\
	     denom += weight[i];		\
	     numer += weight[i] * v[i];	\
	     result[i] = numer/denom;		\
	   })
    } else {
  PAIRLOOP({ numer = denom = 0.0; },
	   {								\
	     wij = weight[j] * exp(-(dx * (dx * s11 + dy * s12)		\
				     + dy * (dx * s21 + dy * s22))/2.0); \
	     denom += wij;						\
	     numer += wij * v[j];					\
	   },
	   {					\
	     result[i] = numer/denom;		\
	   })
    }
}

/* ----------------- transformed coordinates -------------------- */
/*

   The following functions assume that x, y have been transformed
   by the inverse of the variance matrix,
   and subsequently scaled by 1/sqrt(2) so that
   the Gaussian density is proportional to exp(-(x^2+y^2)). 

   Constant factor in density is omitted.
*/
   
void Gdenspt(
     /* inputs */
     int *nxy,                /* number of (x,y) points */
     double *x, double *y,    /* (x,y) coordinates */
     double *rmaxi,           /* maximum distance at which points contribute */
     /* output */
     double *result           /* vector of computed density values */
) {
  STD_DECLARATIONS;
  double resulti;
  STD_INITIALISE;

  if(n == 0) 
    return;

  PAIRLOOP( { resulti = 0.0; },
            { resulti += exp(-d2); } ,
	    { result[i] = resulti; })
}

void Gwtdenspt(
     /* inputs */
     int *nxy,                /* number of (x,y) points */
     double *x, double *y,    /* (x,y) coordinates */
     double *rmaxi,           /* maximum distance */
     double *weight,          /* vector of weights */
     /* output */
     double *result           /* vector of weighted density values */
) {
  STD_DECLARATIONS;
  double resulti;	
  STD_INITIALISE;

  if(n == 0) 
    return;

  PAIRLOOP( { resulti = 0.0; },
	    { resulti += weight[j] * exp(-d2); },
	    { result[i] = resulti; } )
 }

void Gsmoopt(
     /* inputs */
     int *nxy,                /* number of (x,y) points */
     double *x, double *y,    /* (x,y) coordinates */
     double *v,               /* vector of mark values to be smoothed */
     int *self,               /* 0 if leave-one-out */
     double *rmaxi,           /* maximum distance at which points contribute */
     /* output */
     double *result           /* vector of computed smoothed values */
) {
  STD_DECLARATIONS;
  int countself;
  double numer, denom, wij; 

  STD_INITIALISE;
  countself = *self;

  if(n == 0) 
    return;

  if(countself != 0) {
  PAIRLOOP({ numer = denom = 0.0; },
	   { \
	     wij = exp(-d2);		\
	     denom += wij;			\
	     numer += wij * v[j];		\
	   },					
	   {					\
	     denom += 1;			\
	     numer += v[i];			\
	     result[i] = numer/denom;		\
	   })
    } else {
  PAIRLOOP({ numer = denom = 0.0; },
	   { \
	     wij = exp(-d2);		\
	     denom += wij;			\
	     numer += wij * v[j];		\
	   },					
	   {					\
	     result[i] = numer/denom;		\
	   })
    }
 }


void Gwtsmoopt(
     /* inputs */
     int *nxy,                /* number of (x,y) points */
     double *x, double *y,    /* (x,y) coordinates */
     double *v,               /* vector of mark values to be smoothed */
     int *self,               /* 0 if leave-one-out */
     double *rmaxi,           /* maximum distance */
     double *weight,          /* vector of weights */
     /* output */
     double *result           /* vector of computed smoothed values */
) {
  STD_DECLARATIONS;
  int countself;
  double numer, denom, wij; 

  STD_INITIALISE;
  countself = *self;

  if(n == 0) 
    return;

  if(countself != 0) {
  PAIRLOOP({ numer = denom = 0.0; },
	   {						\
	     wij = weight[j] * exp(-d2);	\
	     denom += wij;				\
	     numer += wij * v[j];			\
	   },						
	   {						\
	     denom += weight[i];			\
	     numer += weight[i] * v[i];		\
	     result[i] = numer/denom;			\
	   })
  } else {
  PAIRLOOP({ numer = denom = 0.0; },
	   {						\
	     wij = weight[j] * exp(-d2);	\
	     denom += wij;				\
	     numer += wij * v[j];			\
	   },						
	   {						\
	     result[i] = numer/denom;			\
	   })
    }
}
