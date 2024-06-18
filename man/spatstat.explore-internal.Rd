\name{spatstat.explore-internal} 
\title{Internal spatstat.explore functions}
\alias{spatstat.explore-internal} %DoNotExport
\alias{[.localpcfmatrix}
\alias{[.rat}
\alias{adjust.ratfv}
\alias{ang2rad}
\alias{as.data.frame.bw.optim}
\alias{as.data.frame.fv}
\alias{assemble.plot.objects}       
\alias{bandwidth.is.infinite}
\alias{BartCalc}
\alias{bermantestCalc}
\alias{bermantestEngine}
\alias{bind.ratfv}
\alias{bw.optim}
\alias{calc.DR}
\alias{calc.NNIR}
\alias{calc.SAVE}
\alias{calc.SIR}
\alias{calc.TSE}
\alias{censtimeCDFest}
\alias{check.testfun}
\alias{circticks}
\alias{clarkevansCalc}
\alias{compatible.rat}
\alias{conform.ratfv}
\alias{CressieReadStatistic}
\alias{CressieReadSymbol}
\alias{CressieReadName}
\alias{cutoff2Dkernel}
\alias{CVforPCF}
\alias{Deviation}
\alias{densitycrossEngine}
\alias{densitypointsEngine}
\alias{digestCovariates}
\alias{digital.volume}
\alias{dim.fasp}               
\alias{dimnames.fasp}               
\alias{dimnames<-.fasp}
\alias{distributecbind}
\alias{ensure.listarg}
\alias{envelopeEngine}
\alias{envelopeProgressData}
\alias{envelopeTest}
\alias{envelope.hasenvelope}
\alias{envelope.matrix}
\alias{evaluateCovariate}
\alias{evaluateCovariateAtPoints}
\alias{evaluateCovariateAtPixels}
\alias{evaluate2Dkernel}
\alias{ExpSmoothLog}
\alias{exactppm}
\alias{extractAtomicQtests}
\alias{fasp}
\alias{f3engine}
\alias{f3Cengine}
\alias{findbestlegendpos}
\alias{findcbind}
\alias{flatfname}
\alias{FormatFaspFormulae}
\alias{fvexprmap}
\alias{fvlabels}
\alias{fvlabels<-}
\alias{fvlabelmap}
\alias{fvlegend}
\alias{g3engine}
\alias{g3Cengine}
\alias{getSumFun}
\alias{good.correction.K}
\alias{hasenvelope}
\alias{implemented.for.K}
\alias{is.atomicQtest}
\alias{is.poisson.exactppm}
\alias{is.scov}
\alias{is.stationary.exactppm}
\alias{k3engine}
\alias{Kborder.engine}
\alias{Knone.engine}
\alias{Krect.engine}
\alias{Kount}
\alias{Kwtsum}               
\alias{localKengine}
\alias{localKmultiEngine}
\alias{localpcfengine}
\alias{localpcfmatrix}
\alias{lookup2DkernelInfo}
\alias{makefvlabel}
\alias{maskLaslett}
\alias{match2DkernelName}
\alias{mctestSigtraceEngine}
\alias{meanlistfv}
\alias{names<-.fv}
\alias{nncleanEngine}
\alias{pairs.listof}
\alias{pairs.solist}
\alias{pcf3engine}
\alias{pcfmulti.inhom}
\alias{plot.bw.frac}
\alias{plot.bw.optim}
\alias{plot.localpcfmatrix}
\alias{plot.plotpairsim}
\alias{plot.spatialcdf}
\alias{polyLaslett}
\alias{predict.exactppm}
\alias{prefixfv}
\alias{print.bw.frac}
\alias{print.bw.optim}
\alias{print.densityfun}
\alias{print.envelope}
\alias{print.exactppm}
\alias{print.fasp}       
\alias{print.fv}       
\alias{print.fvfun}       
\alias{print.hasenvelope}       
\alias{print.laslett}
\alias{print.localpcfmatrix}
\alias{print.plotpairsim}
\alias{print.quadrattest}
\alias{print.rat}
\alias{print.Smoothfun}       
\alias{print.summary.bw.optim}
\alias{print.summary.ssf}
\alias{quadrat.testEngine}
\alias{ratfv}
\alias{rebadge.as.crossfun}
\alias{rebadge.as.dotfun}
\alias{rebadge.fv}
\alias{rebadgeLabels}
\alias{reconcile.fv}
\alias{RelevantDeviation}
\alias{rename.fv}
\alias{resolve.2D.kernel}
\alias{resolve.foxall.window}
\alias{resolve.lambda}
\alias{resolve.lambda.ppp}
\alias{resolve.lambdacross}
\alias{resolve.lambdacross.ppp}
\alias{resolve.reciplambda}
\alias{resolve.reciplambda.ppp}
\alias{resolveEinfo}
\alias{rhohatEngine}
\alias{rhohatCalc}
\alias{rmax.Rigid}
\alias{rmax.rule}
\alias{rocData}
\alias{rocModel}
\alias{roseContinuous}
\alias{scanmeasure}
\alias{scanmeasure.ppp}
\alias{scanmeasure.im}
\alias{scanBinomLRTS}
\alias{scanPoisLRTS}
\alias{second.moment.calc}
\alias{second.moment.engine}
\alias{sewpcf}
\alias{sewsmod}
\alias{shift.quadrattest}              
\alias{simulrecipe}              
\alias{StieltjesCalc.fv}
\alias{Smooth.solist}
\alias{smoothcrossEngine}              
\alias{smoothpointsEngine}
\alias{spatialCDFframe}
\alias{spatialCDFtest}
\alias{spatialCDFtestCalc}
\alias{spatialCovariateEvidence}
\alias{spatialCovariateEvidence.exactppm}
\alias{spatialCovariateEvidence.ppp}
\alias{sphere.volume}
\alias{summary.envelope}
\alias{summary.bw.optim}
\alias{tweak.fv.entry}
\alias{tweak.ratfv.entry}
\alias{twostage.test}
\alias{twostage.envelope}
\alias{updateData}
\alias{updateData.default}
\alias{validate2Dkernel}
\alias{validate.angles}        
\alias{validate.weights}        
\alias{vanilla.fv}
\alias{weightedclosepairs}
\alias{X2testEngine}




%%%%%%%
\description{
  Internal spatstat.explore functions.
}
\usage{

\method{[}{localpcfmatrix}(x, i, \dots)
\method{[}{rat}(x, \dots)
adjust.ratfv(f, columns, numfactor, denfactor)
ang2rad(ang, unit, start, clockwise)
\method{as.data.frame}{bw.optim}(x, \dots)
\method{as.data.frame}{fv}(x, \dots)
assemble.plot.objects(xlim, ylim, \dots, lines, polygon)
bandwidth.is.infinite(sigma)
BartCalc(fY, fK)
bermantestCalc(fram, which, alternative, \dots)
bermantestEngine(model, covariate, which, alternative, \dots,
                 modelname, covname, dataname)
bind.ratfv(x, numerator, denominator, labl, desc, preferred,
           ratio, quotient)
bw.optim(cv, h, iopt, \dots, cvname, hname, criterion,
         optimum, warnextreme, hargnames, yexp, unitname,
         template, exponent, hword)
calc.DR(COV, z, Dim)
calc.NNIR(COV, z, pos, Dim)
calc.SAVE(COV, z, Dim)
calc.SIR(COV, z)
calc.TSE(COV, z, pos, Dim1, Dim2)
censtimeCDFest(o, cc, d, breaks, \dots,
     KM, RS, HAN, RAW, han.denom, tt, pmax, fname, fexpr)
check.testfun(f, f1, X)
circticks(R, at, unit, start, clockwise, labels)
clarkevansCalc(X, correction, clipregion, working)
\method{compatible}{rat}(A, B, \dots)
conform.ratfv(x)
CressieReadStatistic(OBS,EXP,lambda,normalise,named)
CressieReadSymbol(lambda)
CressieReadName(lambda)
cutoff2Dkernel(kernel, sigma, varcov, \dots, scalekernel, cutoff, fatal)
CVforPCF(bw, stuff)
Deviation(x, ref, leaveout, n, xi)
densitycrossEngine(Xdata, Xquery, sigma, \dots,
                    kernel, scalekernel,
                    weights, edge, varcov,
                    diggle, sorted, cutoff,
                    se, kerpow)
densitypointsEngine(x, sigma, \dots,
                    kernel, scalekernel, kerpow,
                    weights, edge, varcov,
                    leaveoneout, diggle, sorted, spill, cutoff,
                    debug)
digestCovariates(\dots, W)
digital.volume(range, nval, vside)
\method{dim}{fasp}(x)
\method{dimnames}{fasp}(x)
\method{dimnames}{fasp}(x) <- value
distributecbind(x)
ensure.listarg(x, n, singletypes, xtitle, things)
envelopeEngine(X, fun, simul,
           nsim, nrank, \dots, funargs, funYargs,
           verbose, clipdata, 
           transform, global, ginterval, use.theory,
           alternative, scale, clamp,
           savefuns, savepatterns, saveresultof,
           weights,
           nsim2, VARIANCE, nSD,
           Yname,
           maxnerr, rejectNA, silent, maxerr.action,
           internal, cl,
           envir.user, expected.arg, do.pwrong,
           foreignclass, collectrubbish)
envelopeProgressData(X, fun, \dots, exponent,
                     alternative, leaveout, scale, clamp,
                     normalize, deflate, rmin, 
                     save.envelope, savefuns, savepatterns)
envelopeTest(X, \dots, exponent, alternative,
            rinterval, leaveout, scale, clamp, tie.rule,
            interpolate, save.interpolant,
            save.envelope, savefuns, savepatterns,
            Xname, badXfatal, verbose)
\method{envelope}{hasenvelope}(Y, \dots, Yname)
\method{envelope}{matrix}(Y, \dots, argvals, rvals,
  observed, theory, funX, nsim, nsim2, jsim, jsim.mean,
  type, alternative, scale, clamp, csr, use.theory, nrank, ginterval, nSD,
  savefuns, check, Yname, argname, arg.desc,
  do.pwrong, weights, precomputed, gaveup)
evaluateCovariate(covariate, locations, \dots)
evaluateCovariateAtPixels(covariate, locations, \dots,
                          types, eps, dimyx, rule.eps)
evaluateCovariateAtPoints(covariate, locations, \dots, allow.column)
evaluate2Dkernel(kernel, x, y, sigma, varcov, \dots, scalekernel)
ExpSmoothLog(X, \dots, at, weights, se)
exactppm(X, baseline, \dots, subset, eps, dimyx, rule.eps)
extractAtomicQtests(x)
fasp(fns, which, formulae, dataname, title, rowNames, colNames, checkfv)
f3engine(x, y, z, box, vside, range, nval, correction)
f3Cengine(x, y, z, box, vside, rmax, nrval)
findbestlegendpos(\dots)
findcbind(root, depth, maxdepth)
flatfname(x)
FormatFaspFormulae(f, argname)
fvexprmap(x)
fvlabels(x, expand=FALSE)
fvlabels(x) <- value
fvlabelmap(x, dot=TRUE)
fvlegend(object, elang)
g3engine(x, y, z, box, rmax, nrval, correction)
g3Cengine(x, y, z, box, rmax, nrval)
getSumFun(abbreviation, classname, ismarked, fatal)
good.correction.K(X)
hasenvelope(X, E)
implemented.for.K(correction, windowtype, explicit)
is.atomicQtest(x)
\method{is.poisson}{exactppm}(x)
is.scov(x)
\method{is.stationary}{exactppm}(x)
k3engine(x, y, z, box, rmax, nrval, correction)
Kborder.engine(X, rmax, nr, correction, weights, ratio)
Knone.engine(X, rmax, nr, weights, ratio)
Krect.engine(X, rmax, nr, correction, weights,
             ratio, fname, use.integers)
Kount(dIJ, bI, b, breaks)
Kwtsum(dIJ, bI, wIJ, b, w, breaks, fatal)
localKengine(X, \dots, wantL, lambda, rmax, correction, verbose, rvalue)
localKmultiEngine(X, from, to,
           lambdaFrom, lambdaTo, \dots, rmax, wantL,
           correction, verbose, rvalue, sigma, varcov,
           lambdaX, update, leaveoneout,
           Iexplain, Jexplain, Ikey, Jkey)
localpcfengine(X, \dots, delta, rmax, nr, stoyan, lambda, rvalue)
localpcfmatrix(X, i, \dots, lambda, delta, rmax, nr, stoyan)
lookup2DkernelInfo(kernel)
makefvlabel(op, accent, fname, sub, argname)
maskLaslett(X, \dots, eps, dimyx, xy, rule.eps, oldX, verbose, plotit)
match2DkernelName(kernel)
mctestSigtraceEngine(R, devdata, devsim, \dots,
     interpolate, confint, alpha, exponent, unitname)
meanlistfv(z, \dots)
\method{names}{fv}(x) <- value
nncleanEngine(kthNND, k, d, \dots, tol, maxit,
              plothist, lineargs, verbose, Xname)
\method{pairs}{listof}(\dots, plot=TRUE)
\method{pairs}{solist}(\dots, plot=TRUE)
pcf3engine(x, y, z, box, rmax, nrval, correction, delta)
pcfmulti.inhom(X, I, J, lambdaI, lambdaJ, \dots,
               lambdaX, r, breaks, kernel,
               bw, adjust.bw, stoyan, correction,
               sigma, adjust.sigma, varcov, update, leaveoneout, 
               Iname, Jname)
\method{plot}{bw.frac}(x, \dots)
\method{plot}{bw.optim}(x, \dots, showopt, optargs)
\method{plot}{localpcfmatrix}(x, \dots)
\method{plot}{plotpairsim}(x, \dots)
\method{plot}{spatialcdf}(x, \dots, xlab, ylab, do.points)
polyLaslett(X, \dots, oldX, verbose, plotit)
\method{predict}{exactppm}(object, locations, \dots, eps, dimyx, rule.eps)
prefixfv(x, tagprefix, descprefix, lablprefix, whichtags)
\method{print}{bw.frac}(x, \dots)
\method{print}{bw.optim}(x, \dots)
\method{print}{densityfun}(x, \dots)
\method{print}{envelope}(x, \dots)
\method{print}{exactppm}(x, \dots)
\method{print}{fasp}(x, \dots)
\method{print}{fv}(x, \dots, tight)
\method{print}{fvfun}(x, \dots)
\method{print}{hasenvelope}(x, \dots)
\method{print}{laslett}(x, \dots)
\method{print}{localpcfmatrix}(x, \dots)
\method{print}{plotpairsim}(x, \dots)
\method{print}{quadrattest}(x, \dots)
\method{print}{rat}(x, \dots)
\method{print}{Smoothfun}(x, \dots)
\method{print}{summary.bw.optim}(x, \dots, digits)
\method{print}{summary.ssf}(x, \dots)
quadrat.testEngine(X, nx, ny, alternative, method, conditional, CR, 
     \dots, nsim, Xcount, xbreaks, ybreaks, tess, fit, df.est, Xname, fitname)
ratfv(df, numer, denom, \dots, ratio)
rebadge.as.crossfun(x, main, sub, i, j)
rebadge.as.dotfun(x, main, sub, i)
rebadge.fv(x, new.ylab, new.fname, tags, new.desc, new.labl, new.yexp,
           new.dotnames, new.preferred, new.formula, new.tags)
rebadgeLabels(x, new.fname)
reconcile.fv(\dots)
RelevantDeviation(x, alternative, clamp, scaling)
rename.fv(x, fname, ylab, yexp)
resolve.2D.kernel(\dots,
            sigma, varcov, x, mindist, adjust, bwfun, allow.zero)
resolve.foxall.window(X, Y, W, warn.trim)
resolve.lambda(X, lambda, \dots)
\method{resolve.lambda}{ppp}(X, lambda, \dots,
               sigma, varcov, leaveoneout, update, check)
resolve.lambdacross(X, I, J, lambdaI, lambdaJ, \dots)
\method{resolve.lambdacross}{ppp}(X, I, J, lambdaI, lambdaJ, \dots,
         lambdaX, sigma, varcov, leaveoneout, update, lambdaIJ,
         Iexplain, Jexplain)
resolve.reciplambda(X, lambda, reciplambda, \dots)
\method{resolve.reciplambda}{ppp}(X, lambda, reciplambda, \dots,
                    sigma, varcov, leaveoneout, update, check)
resolveEinfo(x, what, fallback, warn, atomic)
rhohatEngine(model, covariate, reference, volume, \dots,
               subset, weights, method, horvitz, smoother,
               resolution, spatCovarArgs,
               n, bw, adjust, from, to, 
               bwref, covname, covunits, confidence,
               breaks,
               modelcall, callstring)
rhohatCalc(ZX, Zvalues, lambda, denom, \dots,
           weights, lambdaX,
           method, horvitz, smoother, do.CI, 
           n, bw, adjust, from, to, 
           bwref, covname, confidence, breaks, positiveCI, markovCI,
           covunits, modelcall, callstring, savestuff)
rmax.Rigid(X, g)
rmax.rule(fun, W, lambda)
rocData(covariate, nullmodel, \dots, high, p)
rocModel(lambda, nullmodel, \dots, high, p)
roseContinuous(ang, rad, unit, \dots,
               start, clockwise, main, labels, at, do.plot)
scanmeasure(X, \dots)
\method{scanmeasure}{ppp}(X, r, \dots, method) 
\method{scanmeasure}{im}(X, r, \dots)
scanPoisLRTS(nZ, nG, muZ, muG, alternative)
scanBinomLRTS(nZ, nG, muZ, muG, alternative)
second.moment.calc(x, sigma, edge, what, \dots,
                   varcov, expand, obswin, npts, debug)
second.moment.engine(x, sigma, edge, what, \dots,
      kernel, scalekernel, kerpow,
      obswin, varcov, npts, debug, fastgauss)
sewpcf(d, w, denargs, lambda2area, divisor)
sewsmod(d, ff, wt, Ef, rvals, method="smrep", \dots, nwtsteps=500)
\method{shift}{quadrattest}(X, \dots)
simulrecipe(type, expr, envir, csr, pois, constraints)
\method{StieltjesCalc}{fv}(M, f, \dots)
\method{Smooth}{solist}(X, \dots)
smoothcrossEngine(Xdata, Xquery, values, sigma, \dots,
                    weights, varcov,
                    kernel, scalekernel, sorted, cutoff)
smoothpointsEngine(x, values, sigma, \dots,
                   kernel, scalekernel,
                   weights, varcov, leaveoneout, sorted, cutoff, debug)
spatialCDFframe(model, covariate, \dots, jitter, covariateAtPoints,
                make.quantile.function)
spatialCDFtest(model, covariate, test, \dots,
         dimyx, eps, rule.eps, interpolate, jitter,
         nsim, verbose, modelname, covname, dataname)
spatialCDFtestCalc(fra, test, \dots, details)
spatialCovariateEvidence(model, covariate, \dots)
\method{spatialCovariateEvidence}{exactppm}(model, covariate, \dots, lambdatype,
          dimyx, eps, rule.eps, interpolate, jitter, jitterfactor,
          modelname, covname, dataname, subset, clip.predict) 
\method{spatialCovariateEvidence}{ppp}(model, covariate, \dots, lambdatype,
          dimyx, eps, rule.eps, interpolate, jitter, jitterfactor,
          modelname, covname, dataname, subset, clip.predict)
sphere.volume(range, nval = 10)
\method{summary}{bw.optim}(object, \dots)
\method{summary}{envelope}(object,\dots)
tweak.fv.entry(x, current.tag, new.labl, new.desc, new.tag)
tweak.ratfv.entry(x, \dots)
twostage.test(X, \dots, exponent, nsim, nsimsub,
              alternative, reuse, leaveout, interpolate,
              savefuns, savepatterns, verbose, badXfatal, testblurb)
twostage.envelope(X, \dots, nsim, nsimsub, nrank,
                  alternative, reuse, leaveout, interpolate,
                  savefuns, savepatterns, verbose, badXfatal, testlabel)
updateData(model, X, \dots)
\method{updateData}{default}(model, X, \dots, warn)
validate2Dkernel(kernel, fatal)
validate.angles(angles, unit, guess)
validate.weights(x, recip, how, allowzero, allowinf)
vanilla.fv(x)
weightedclosepairs(X, r, correction, what)
X2testEngine(OBS, EXP, \dots, method, CR, df, nsim,
     conditional, alternative, testname, dataname)




}
\details{
  These internal \pkg{spatstat.explore} functions should not be called
  directly by the user. Their names and capabilities may change
  without warning from one version of \pkg{spatstat.explore} to the next.
}
\value{
The return values of these functions are not documented,
and may change without warning.
}
\keyword{internal}
