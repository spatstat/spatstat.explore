#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>

/*
  segdens.c

  Convolution of segments with Gaussian kernel

  Entry points:

  segdens    unweighted
  segwdens   weighted

  Copyright (c) Adrian Baddeley, 02 dec 2016 (modified 01 mar 2021)
  Licence: GPL >= 2.0

  $Revision: 1.8 $ $Date: 2022/10/23 05:56:57 $
*/


#define DNORM(X, SIG) dnorm((X), (double) 0.0, (SIG), FALSE)

#define PNORM(X, SIG) pnorm((X), (double) 0.0, (SIG), TRUE, FALSE)

void segdens(
  double *sigma,   /* bandwidth */
  /* input line segments */
  int *ns,         /* number of line segments */
  double *xs,
  double *ys,
  double *alps,
  double *lens,    /* first endpoint, angle, length */
  /* query locations */
  int *np,         /* number of pixels or test locations */
  double *xp,
  double *yp,      /* coordinates of pixels or test locations */
  /* output */
  double *z       /* result, assumed initially 0 */
) {
  int i, j, Ns, Np;
  double Sigma;
  double xsi, ysi, angi, leni, cosi, sini;
  double dx, dy, u1, u2;

  Ns = *ns;
  Np = *np;
  Sigma = *sigma;

  for(i = 0; i < Ns; i++) {
    R_CheckUserInterrupt();
    xsi = xs[i];
    ysi = ys[i];
    angi = alps[i];
    leni = lens[i];
    cosi = cos(angi);
    sini = sin(angi);
    for(j = 0; j < Np; j++) {
      dx = xp[j] - xsi;
      dy = yp[j] - ysi;
      u1 = dx * cosi + dy * sini;
      u2 = -dx * sini + dy * cosi;
      z[j] += DNORM(u2, Sigma) * (PNORM(u1, Sigma) - PNORM(u1-leni, Sigma));
    }
  }
}

void segwdens(
  double *sigma,   /* bandwidth */
  /* input line segments */
  int *ns,         /* number of line segments */
  double *xs,
  double *ys,
  double *alps,
  double *lens,    /* first endpoint, angle, length */
  double *ws,  /* segment weights */
  /* query locations */
  int *np,         /* number of pixels or test locations */
  double *xp,
  double *yp,      /* coordinates of pixels or test locations */
  /* output */
  double *z       /* result, assumed initially 0 */
) {
  int i, j, Ns, Np;
  double Sigma;
  double xsi, ysi, angi, leni, cosi, sini, wi;
  double dx, dy, u1, u2;

  Ns = *ns;
  Np = *np;
  Sigma = *sigma;

  for(i = 0; i < Ns; i++) {
    R_CheckUserInterrupt();
    xsi = xs[i];
    ysi = ys[i];
    angi = alps[i];
    leni = lens[i];
    wi   = ws[i];
    cosi = cos(angi);
    sini = sin(angi);
    for(j = 0; j < Np; j++) {
      dx = xp[j] - xsi;
      dy = yp[j] - ysi;
      u1 = dx * cosi + dy * sini;
      u2 = -dx * sini + dy * cosi;
      z[j] += wi * DNORM(u2, Sigma) * (PNORM(u1, Sigma) - PNORM(u1-leni, Sigma));
    }
  }
}
