/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: HGGRooPdfs.cc,v 1.1 2012/02/10 15:10:48 gpetrucc Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
 
//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Power function p.d.f
// END_HTML
//
 

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>
#include "assert.h"

#include "../interface/MCGBRIntegrator.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TMath.h"
#include <Math/QuantFuncMathCore.h>
#include <Math/ProbFunc.h>
#include "../interface/MCGBRForest.h"
#include "../interface/MCGBRTreeD.h"
#include "TH1D.h"
#include "omp.h"
#include <malloc.h>
#include "../interface/GBRArrayUtils.h"
#include "../interface/GBRMath.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TH2D.h"
#include <Eigen/Dense>
#include "../interface/Faddeeva.hh"
#include "TFile.h"
#include "TTree.h"

//-----------------------------------------------------------------------------
// Testing function of Peter Lepage
//-----------------------------------------------------------------------------
double MCGBRIntegrator::Camel(int nDim, double *Xarg) const {
  double Fun1,Fun2,R1,R2;
  double pos1=1e0/3e0;
  double pos2=2e0/3e0;
  double Gam1= 0.100e0;  // as in JPC
  double Gam2= 0.100e0;  // as in JPC
  double sPi = sqrt(M_PI);
  double xn1=1e0;
  double xn2=1e0;
  int i;
  R1=0;
  R2=0;
  for(i = 0 ; i<nDim ; i++){
    R1=R1+(Xarg[i] -pos1)*(Xarg[i] -pos1);
    R2=R2+(Xarg[i] -pos2)*(Xarg[i] -pos2);
    xn1=xn1*Gam1*sPi;
    xn2=xn2*Gam2*sPi;      
  }
  R1   = sqrt(R1);
  R2   = sqrt(R2);
  Fun1 = exp(-(R1*R1)/(Gam1*Gam1))/xn1;  // Gaussian delta-like profile
  Fun2 = exp(-(R2*R2)/(Gam2*Gam2))/xn2;  // Gaussian delta-like profile
  return 0.5e0*(Fun1+ Fun2);
}

double MCGBRIntegrator::Camelrnd(int nDim, double *Xarg) const {
  
  double oval = Camel(nDim,Xarg);
  return oval;
  double logval = log(oval) + 0.5*gRandom->Gaus();
  return exp(logval);

}

// double MCGBRIntegrator::NLSkewGaus(double x, double mu) const {
//   constexpr double sigma1 = 1.;
//   constexpr double k1 = 0.5/sigma1/sigma1;
//   constexpr double k2 = 10.;
//   
// //   if ( k2*(x-mu) > 5.) {
// //     return k1*(x-mu)*(x-mu) + k2*(x-mu)*(x-mu);
// //   }
//   
//   return k1*(x-mu)*(x-mu) - log(erfc(k2*(x-mu)));
//   
// //   return -log(erfc(k2*(x-mu)));
//   
// }

double MCGBRIntegrator::NLSkewGausDMu(double x, double envelope) const {
//   const double sigma = doenv ? std::min(10.,std::max(1e-1,0.1*extx)) : 0.01;
//   const double sigma = doenv ? std::max(1e-1,0.1*x) : 0.01;
//   const double sigma = doenv ? 1. : 0.01;
//   const double alpha = doenv ? 0. : 0.;
  
  
  const double sigma = 0.01;
  const double alpha = 100.;
  
  const double k1 = 0.5/sigma/sigma;
  const double k2 = alpha/sigma/sqrt(2.);
  
//   double mu = std::abs(alpha) > 0. ? envelope - 5.*sigma/alpha : envelope;
  double mu = envelope;
  
  double g = k2*(x-mu);
  double erfcxval = 1./(sqrt(M_PI)*Faddeeva::erfcx(g));
  
  double drv = -2.*k1*(x-mu) - 2.*k2*erfcxval;
  
  return drv;
  
}

double MCGBRIntegrator::NLSkewGausD2Mu(double x, double envelope) const {
//   const double sigma = doenv ? std::min(10.,std::max(1e-1,0.1*extx)) : 0.01;
//   const double sigma = doenv ? std::max(1e-1,0.1*x) : 0.01;
//   const double sigma = doenv ? 1. : 0.01;
//   const double alpha = doenv ? 0. : 0.;
  
  const double sigma = 0.01;
  const double alpha = 100.;
  
  const double k1 = 0.5/sigma/sigma;
  const double k2 = alpha/sigma/sqrt(2.);
  
//   double mu = std::abs(alpha) > 0. ? envelope - 5.*sigma/alpha : envelope;
  double mu = envelope;
  
  double g = k2*(x-mu);
  double erfcxval = 1./(sqrt(M_PI)*Faddeeva::erfcx(g));
  
//   double drv = k1 + 4.*k2*k2*(1.-sqrt(M_PI)*g*erfcxval)/(M_PI*erfcxval*erfcxval);
//   double drv = k1 + 4.*k2*k2*(1./(M_PI*erfcxval*erfcxval)-g/(sqrt(M_PI)*erfcxval))/(M_PI*erfcxval*erfcxval);
  double drv = 2.*k1 + 4.*k2*k2*erfcxval*(erfcxval-g);
  
  if (!std::isnormal(drv)) {
    printf("x = %5e, mu = %5e, g = %5e, erfcxval = %5e\n",x,mu,g,erfcxval);
  }
  
  return drv;
  
}

double MCGBRIntegrator::NLNormDMu(double x, double mu, double sigma) const {
//   constexpr double sigma = 0.01;
  const double k = 0.5/sigma/sigma;
  
  double drv = -2.*k*(x-mu);
  return drv;
}

double MCGBRIntegrator::NLNormD2Mu(double x, double mu, double sigma) const {
//   constexpr double sigma = 0.01;
  const double k = 0.5/sigma/sigma;
  
  double drv = 2.*k;
  return drv;
}

double MCGBRIntegrator::NLLogNormDMu(double x, double mu) const {
//   const double sigmabase = 1.;
//   const double sigmasq = sigmabase*sigmabase;
//   const double sigmascale = 0.;
//   const double sigmaalt = exp(std::floor(log(x)));
//   const double sigmasq = sigmabase*sigmabase + sigmascale*sigmascale*sigmaalt*sigmaalt;
//   double xdef = x > 100. ? std::min(x,1000.) - 100. + 0.1 : 0.;
//   double xdef = x>10. ? std::min(x,100.) - 10. + 1e-2 : 0.;
//   const double xdef = exp(std::round(3.*log(x))/3.);
//   const double k = 0.5/sigmasq;
  
  constexpr double sigma = 0.1;
  constexpr double k = 0.5/sigma/sigma;
  
//   double xmod = std::max(1e-1,x);
//   double drv = -2.*k*(x-mu)/(xmod*xmod);
  
//   constexpr double xmodbase = 0.1;
//   const double xsqmod = xmodbase*xmodbase + x*x;
  double drv = -2.*k*(log(x)-log(mu))/mu;
  
//   double drv = 0 ? -2.*k*(log(x)-log(mu))/mu : -2.*k*(x-mu)/xsqmod;
  
//   double drv = x>mu ? -2.*k*(log(x)-log(mu))/mu : 0.;
//   double drv = -2.*k*(x-mu);
//   double drv = -2.*k*(x-mu);
  return drv;
}

double MCGBRIntegrator::NLLogNormD2Mu(double x, double mu) const {
//   const double sigmabase = 1.;
//   const double sigmasq = sigmabase*sigmabase;
//   const double sigmascale = 0.;
//   const double sigmaalt = exp(std::floor(log(x)));
//   const double sigmasq = sigmabase*sigmabase + sigmascale*sigmascale*sigmaalt*sigmaalt;
//   double xdef = x > 100. ? std::min(x,1000.) - 100. + 0.1 : 0.;
//   double xdef = x>10. ? std::min(x,100.) - 10. + 1e-2 : 0.;
//   const double xdef = exp(std::round(3.*log(x))/3.);
//   const double k = 0.5/sigmasq;
  
  constexpr double sigma = 0.1;
  constexpr double k = 0.5/sigma/sigma;
  
  double drv = 2.*k*(1.+log(x)-log(mu))/(mu*mu);
  
//   double drvdyn = std::max(0.5, 1.+log(x)-log(mu));
//   double drv = 2.*k*drvdyn/(mu*mu);
  
//   constexpr double xmodbase = 0.1;
//   const double xsqmod = xmodbase*xmodbase + x*x;
//   double xmod = std::max(1e-1,x);
//   double drv = 2.*k/(xmod*xmod);
  
//   double drv = 0 ? 2.*k*(1.+log(x)-log(mu))/(mu*mu) : 2.*k/xsqmod;
  
//   double drv = std::max(2.*k*(1.+log(x)-log(mu))/(mu*mu),0.);
//   double drv = x>mu ? 2.*k*(1.+log(x)-log(mu))/(mu*mu) : 0.;
//   double drv = 2.*k;
  return drv;
}


// double MCGBRIntegrator::NLSkewGausDMu(double x, double envelope, bool doenv) const {
//   const double sigma = doenv ? 1. : 0.1;
//   const double delta = 1.;
//   
//   const double k1 = 0.5/sigma/sigma;
//   const double k2 = delta/sigma;
//   
//   const double mu = envelope;
//   
//   double g = (x-mu)/sigma;
//   
//   double sign = g>0. ? 1. : -1.;
//   double drv = std::abs(g)<delta ? -2.*k1*(x-mu) : -sign*k2;
//   
//   return drv;
//   
// }
// 
// double MCGBRIntegrator::NLSkewGausD2Mu(double x, double envelope, bool doenv) const {
//   const double sigma = doenv ? 1. : 0.1;
//   const double delta = 1.;
//   
//   const double k1 = 0.5/sigma/sigma;
// //   const double k2 = delta/sigma;
//   
//   const double mu = envelope;
//   
//   double g = (x-mu)/sigma;
//   
//   double drv = std::abs(g)<delta ? 2.*k1 : 0.;
//   
//   return drv;
//   
// }

// double MCGBRIntegrator::NLSkewGausDMu(double x, double mu) const {
//   constexpr double sigma1 = 1.;
//   constexpr double k1 = 0.5/sigma1/sigma1;
//   constexpr double k2 = 10.;  
//   
//   double t = k2*(x-mu);
//   
//   if (k2*(x-mu)>10.) {
//     return -k1*(x-mu) + 2.*k2*(x-mu) + 1./(x-mu)
//   }
//   
//   constexpr double step = 1e-6;
//   constexpr double step1 = step;
//   constexpr double step2 = 0.5*step;
//   
//   double valup1 = NLSkewGaus(x,mu+step1);
//   double valdown1 = NLSkewGaus(x,mu-step1);
//   double valup2 = NLSkewGaus(x,mu+step2);
//   double valdown2 = NLSkewGaus(x,mu-step2);
//   
//   double drv1 = (valup1-valdown1)/(2.0*step1);
//   double drv2 = (valup2-valdown2)/(2.0*step2);
//   
//   double drv = (4.0*drv2 - drv1)/3.0;
//   
//   return -k1*(x-mu) + drv;
//   
// }
// 
// double MCGBRIntegrator::NLSkewGausD2Mu(double x, double mu) const {
//   constexpr double sigma1 = 1.;
//   constexpr double k1 = 0.5/sigma1/sigma1;
//   constexpr double k2 = 10.;  
//   
//   constexpr double step = 1e-6;
//   constexpr double step1 = step;
//   constexpr double step2 = 0.5*step;
//   
//   double valnom = NLSkewGaus(x,mu);
//   double valup1 = NLSkewGaus(x,mu+step1);
//   double valdown1 = NLSkewGaus(x,mu-step1);
//   double valup2 = NLSkewGaus(x,mu+step2);
//   double valdown2 = NLSkewGaus(x,mu-step2);
//   
//   double drv1 = (valup1+valdown1-2.0*valnom)/(step1*step1);
//   double drv2 = (valup2+valdown2-2.0*valnom)/(step2*step2);
//   
//   double drv = (4.0*drv2 - drv1)/3.0;
//   
//   return k1 + drv;
//   
// }



//_____________________________________________________________________________
MCGBRIntegrator::MCGBRIntegrator(const char *name, const char *title, int nevents, int neventsintial) :
  TNamed(name,title),
  fNThreads(std::max(1,omp_get_max_threads())),
//   fNThreads(1),  
  fMinEvents(-99),  
  fMinWeightTotal(-99.),
  fShrinkage(1.0),
  fNTrees(20),
//   fNQuantiles(std::min(nevents+1,std::numeric_limits<unsigned short>::max()+1)),
  fNQuantiles(std::min(nevents,std::numeric_limits<unsigned short>::max()+1)),
//   fNQuantiles(nevents+2),
//   fNQuantiles(pow(2,30)+1),
//   fNQuantiles(4*nevents),
//   fNBinsMax(128),
  fNBinsMax(128),
//   fNBinsMax(std::numeric_limits<unsigned short>::max()+1),
  fMinCutSignificance(-99.),
  fMinCutSignificanceMulti(-99.),
  fNEvents(nevents),
  fNEventsInitial(neventsintial),
  fNEventsBagged(-99),
  fMaxDepth(-1),
  fMaxNodes(-1),
  fNVars(9),
  fForest(0),
  fForestGen(0),
  fGenTree(0),
  _sepgains(0),
  _ws(0),
  fSigmaScale(1.),
  fDoEnvelope(false),
  fStagedGeneration(false),
  fSigmaRatio(1.),
  fIntegralRatio(1.)
{
  
  Eigen::initParallel();
  
  omp_set_num_threads(fNThreads);

  int nvars = fNVars;
  int ncls = 1;
  int nev = fNEventsInitial;
   
  //initialize arrays (ensure 32 byte alignment for avx vector instructions)

  _sepgains = (double*)memalign(32, nvars*sizeof(double));
  _sepgainsigs = (double*)memalign(32, nvars*sizeof(double));
  _cutvals = (float*)memalign(32, nvars*sizeof(float));
  _nlefts = (int*)memalign(32, nvars*sizeof(int));
  _nrights = (int*)memalign(32, nvars*sizeof(int));
  _sumwlefts = (double*)memalign(32, nvars*sizeof(double));
  _sumwrights = (double*)memalign(32, nvars*sizeof(double));
  _sumtgtlefts = (double*)memalign(32, nvars*sizeof(double));
  _sumtgtrights = (double*)memalign(32, nvars*sizeof(double));
  _leftvars = (float*)memalign(32, nvars*sizeof(float));
  _rightvars = (float*)memalign(32, nvars*sizeof(float));  
  _fullvars = (float*)memalign(32, nvars*sizeof(float));  
  _bestbins = (int*)memalign(32, nvars*sizeof(int));
  
  _ws = new double*[nvars];
  _ws2 = new double*[nvars];
  _wscls = new double**[nvars];
  _ns = new int*[nvars];
  _nsd = new int*[nvars];  
  _tgts = new double*[nvars];  
  _tgt2s = new double*[nvars];  
  _tgt3s = new double*[nvars];
  _tgt4s = new double*[nvars];
  _tgt5s = new double*[nvars];
  _fmins = new double*[nvars];
  _tgtmaxs = new double*[nvars];
  _sumws = new double*[nvars];
  _sumws2 = new double*[nvars];
  _sumwscls = new double**[nvars];
  _sumns = new int*[nvars];
  _sumtgts = new double*[nvars];  
  _sumtgt2s = new double*[nvars];
  _sumtgt3s = new double*[nvars];
  _sumtgt4s = new double*[nvars];
  _sumtgt5s = new double*[nvars];
  _sumfmins = new double*[nvars];
  _sumfminsr = new double*[nvars];
  _sumtgtmaxs = new double*[nvars];
  _sumtgtmaxsr = new double*[nvars];
  _varvals = new float*[nvars];    
  _varvalmaxs = new int*[nvars];  
  _varvalmins = new int*[nvars];  
  _sumvarvalmaxs = new int*[nvars];  
  _sumvarvalmins = new int*[nvars];  
  _bsepgains = new double*[nvars];
  _bsepgainsigs = new double*[nvars];
  
  _binquants  = new int*[nvars];  
  _quants  = new int*[nvars];  
  
  _clss  = (int*)memalign(32, nev*sizeof(int));
  _tgtvals  = (double*)memalign(32, nev*sizeof(double));
  _tgt2vals  = (double*)memalign(32, nev*sizeof(double));
  _tgt3vals  = (double*)memalign(32, nev*sizeof(double));
  _tgt4vals  = (double*)memalign(32, nev*sizeof(double));
  _tgt5vals  = (double*)memalign(32, nev*sizeof(double));
  _fvals  = (double*)memalign(32, nev*sizeof(double));
  _weightvals  = (double*)memalign(32, nev*sizeof(double));
  
  fQuantileMaps = new float*[nvars];  
  fQuantileMapsMin = new float*[nvars];  
  
  
  for (int ivar=0; ivar<nvars; ++ivar) {
    _ws[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _ws2[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _ns[ivar] = (int*)memalign(32, fNBinsMax*sizeof(int));
    _tgts[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _tgt2s[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));  
    _tgt3s[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _tgt4s[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));  
    _tgt5s[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));  
    _fmins[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double)); 
    _tgtmaxs[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double)); 
    _sumws[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumws2[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumns[ivar] = (int*)memalign(32, fNBinsMax*sizeof(int));
    _sumtgts[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumtgt2s[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumtgt3s[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumtgt4s[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumtgt5s[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumfmins[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumfminsr[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumtgtmaxs[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumtgtmaxsr[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _varvals[ivar] = (float*)memalign(32, fNBinsMax*sizeof(float));
    _varvalmaxs[ivar] = (int*)memalign(32, fNBinsMax*sizeof(int));
    _varvalmins[ivar] = (int*)memalign(32, fNBinsMax*sizeof(int));
    _sumvarvalmaxs[ivar] = (int*)memalign(32, fNBinsMax*sizeof(int));
    _sumvarvalmins[ivar] = (int*)memalign(32, fNBinsMax*sizeof(int));    
    _bsepgains[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _bsepgainsigs[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    
    _wscls[ivar] = new double*[ncls];
    _sumwscls[ivar] = new double*[ncls];    
    
    _binquants[ivar] = (int*)memalign(32, fNBinsMax*sizeof(int));
    _quants[ivar] = (int*)memalign(32, nev*sizeof(int));
    
    fQuantileMaps[ivar] = (float*)memalign(32, fNQuantiles*sizeof(float));
    fQuantileMapsMin[ivar] = (float*)memalign(32, fNQuantiles*sizeof(float));
    

    for (int icls=0; icls<ncls; ++icls) {
      _wscls[ivar][icls] = (double*)memalign(32, fNBinsMax*sizeof(double));
      _sumwscls[ivar][icls] = (double*)memalign(32, fNBinsMax*sizeof(double));
    }
    
  }
  
  
  fEvts.reserve(nev);
  
//   double sumw = 0.;
  double sumabsw = 0.;
  
  printf("second loop, fill events in memory\n");
  //loop over trees to fill arrays and event vector
   
  for (int ivar=0; ivar<nvars; ++ivar) {
//     fLimits.push_back(std::pair<float,float>(-std::numeric_limits<float>::max(),std::numeric_limits<float>::max()));
//     fLimits.push_back(std::pair<float,float>(-5.0,5.0));
    fLimits.push_back(std::pair<float,float>(0.,1.0));
  }

  std::vector<double> tmpvals(fNVars);
  
//   for (int iev=0; iev<nev; ++iev) {
//     fEvts.push_back(new MCGBREvent(nvars));
//     MCGBREvent &evt = *fEvts.back();
//     for (int ivar=0; ivar<nvars; ++ivar) {
//       double low = fLimits[ivar].first;
//       double high = fLimits[ivar].second;
//       double val = gRandom->Uniform(low,high);
//       evt.SetVar(ivar,val);
// //       printf("ivar = %i, val = %5f\n",ivar,val);
//     }
// //     evt.SetTarget(1.);
//     double funcval = 1.;
//     for (int ivar=0; ivar<nvars; ++ivar) {
//       double val = evt.Var(ivar);
//       funcval *= exp(-0.5*val*val);
//       tmpvals[ivar] = val;
//     }
//     funcval = Camel(fNVars,tmpvals.data());
// //     funcval = funcval + 0.05 + 0.01*gRandom->Gaus();
// //     funcval = funcval + 0.01*gRandom->Gaus();
//     evt.SetFuncVal(funcval);
//     sumabsw += std::abs(evt.Weight());
//   }
  
  int iev = 0;
  while (iev<nev) {
    
    for (int ivar=0; ivar<nvars; ++ivar) {
      double low = fLimits[ivar].first;
      double high = fLimits[ivar].second;
      double val = gRandom->Uniform(low,high);
      tmpvals[ivar] = val;
//       evt.SetVar(ivar,val);
//       printf("ivar = %i, val = %5f\n",ivar,val);
    }
//     double funcval = Camel(fNVars,tmpvals.data());
    double funcval = Camelrnd(fNVars,tmpvals.data());
    
//     double origin = 90.;
//     double target = funcval;
//     double target = 90.-funcval;
//     double origin = 510.;
//     double target = funcval;
//     double origin = 12.71;
//     double target = log(funcval) + 10.;
//     double rnd = gRandom->Uniform(0.,origin);
//     if (target < rnd) continue;
    
    fEvts.push_back(new MCGBREvent(nvars));
    MCGBREvent &evt = *fEvts.back();
//     fEvts2.push_back(&evt);
    
//     evt.SetTarget(1.);
    for (int ivar=0; ivar<nvars; ++ivar) {
      double val = tmpvals[ivar];
      evt.SetVar(ivar,val);
    }
    
//     funcval = funcval + 0.05 + 0.01*gRandom->Gaus();
//     funcval = funcval + 0.01*gRandom->Gaus();
    evt.SetFuncVal(funcval);
    evt.SetFuncValAlt(vdt::fast_log(funcval));
//     double logsigma = log(funcval) + 1.*gRandom->Gaus();
//     evt.SetFuncValAlt(exp(0.1*gRandom->Gaus()));
//     evt.SetTarget(fForestGen->InitialResponse());
//     evt.SetTargetMin(fForest->InitialResponse());
//     evt.SetWeight(1./funcval/funcval);

    sumabsw += std::abs(evt.Weight());
    ++iev;
  }  
  
  
  //second loop here
//   for (unsigned int idata=0; idata<data.size(); ++idata) {
//     std::vector<RooRealVar*> dcondvars(fCondVars.getSize());
//     std::vector<RooRealVar*> dparmvars(fParmVars.getSize());
// 
//     const RooArgSet *dset = data[idata]->get();
//     
//     for (int ivar=0; ivar<fCondVars.getSize(); ++ivar) {
//       dcondvars[ivar] = static_cast<RooRealVar*>(dset->find(fCondVars.at(ivar)->GetName()));
//     }
//       
//     for (int ivar=0; ivar<fParmVars.getSize(); ++ivar) {
//       dparmvars[ivar] = static_cast<RooRealVar*>(dset->find(fParmVars.at(ivar)->GetName()));
//     }
//       
//     for (int iev=0; iev<data[idata]->numEntries(); ++iev) {
//       data[idata]->get(iev);
//       fEvts.push_back(new MCGBREvent(nvars,fNTargets,fFullParms.getSize()));
//       MCGBREvent *evt = fEvts.back();
//       evt->SetWeight(data[idata]->weight());
//       evt->SetClass(idata);
// 
//       sumw += evt->Weight();
//       sumabsw += std::abs(evt->Weight());
//       
//       for (unsigned int ivar=0; ivar<dcondvars.size(); ++ivar) {
// 	evt->SetVar(ivar,dcondvars[ivar]->getVal());
//       }
//       for (unsigned int ivar=0; ivar<dparmvars.size(); ++ivar) {
// 	evt->SetVar(dcondvars.size() + ivar, dparmvars[ivar]->getVal());
//       }    
//     }
//   }
  
  
  printf("filled data\n");
    
  //int nevr = datau->numEntries();
  
  
  //map of input variable quantiles to values
  //fQuantileMaps.resize(nvars, std::vector<float>(fNQuantiles));
  
  //BuildQuantiles(nvars, sumw);
  BuildQuantiles(nvars, sumabsw);  
  
  
  
}

//_____________________________________________________________________________
MCGBRIntegrator::~MCGBRIntegrator() {
  
  int nvars = fNVars;
  int ncls = 1;
  
  for (int ivar=0; ivar<nvars; ++ivar) {
    
    for (int icls=0; icls<ncls; ++icls) {
      free(_wscls[ivar][icls]);
      free(_sumwscls[ivar][icls]);
    }    
    
    free(_ws[ivar]);
    free(_ws2[ivar]);
    free(_ns[ivar]);
    free(_tgts[ivar]);
    free(_tgt2s[ivar]);
    free(_sumws[ivar]);
    free(_sumws2[ivar]);
    free(_sumns[ivar]);
    free(_sumtgts[ivar]);
    free(_sumtgt2s[ivar]);
    free(_sumtgt3s[ivar]);
    free(_sumtgt4s[ivar]);
    free(_sumtgt5s[ivar]);
    free(_sumfmins[ivar]);
    free(_sumfminsr[ivar]);
    free(_sumtgtmaxs[ivar]);
    free(_sumtgtmaxsr[ivar]);
    free(_varvals[ivar]);
    free(_varvalmaxs[ivar]);
    free(_varvalmins[ivar]);
    free(_sumvarvalmaxs[ivar]);
    free(_sumvarvalmins[ivar]);    
    free(_bsepgains[ivar]);
    free(_bsepgainsigs[ivar]);
    
    delete[] _wscls[ivar];
    delete[] _sumwscls[ivar];
    
    free(_binquants[ivar]);
    free(_quants[ivar]);
    
    free(fQuantileMaps[ivar]);
    free(fQuantileMapsMin[ivar]);
  }    
  
  free(_sepgains);
  free(_sepgainsigs);
  free(_cutvals);
  free(_nlefts);
  free(_nrights);
  free(_sumwlefts);
  free(_sumwrights);
  free(_sumtgtlefts);
  free(_sumtgtrights);
  free(_leftvars);
  free(_rightvars);
  free(_fullvars);
  free(_bestbins);
  
  delete[] _ws;
  delete[] _ws2;
  delete[] _wscls;
  delete[] _ns;
  delete[] _nsd;
  delete[] _tgts;
  delete[] _tgt2s;
  delete[] _sumws;
  delete[] _sumws2;
  delete[] _sumwscls;
  delete[] _sumns;
  delete[] _sumtgts;
  delete[] _sumtgt2s;
  delete[] _sumtgt3s;
  delete[] _sumtgt4s;
  delete[] _sumtgt5s;
  delete[] _sumfmins;
  delete[] _sumfminsr;
  delete[] _sumtgtmaxs;
  delete[] _sumtgtmaxsr;
  delete[] _varvals;
  delete[] _varvalmaxs;
  delete[] _varvalmins;
  delete[] _sumvarvalmaxs;
  delete[] _sumvarvalmins;  
  delete[] _bsepgains;
  delete[] _bsepgainsigs;
  
  delete[] _binquants;
  delete[] _quants;
  
  free(_clss);
  free(_tgtvals);
  free(_tgt2vals);
  free(_tgt3vals);
  free(_tgt4vals);
  free(_tgt5vals);
  free(_fvals);
  free(_weightvals);
  
  delete[] fQuantileMaps;
  delete[] fQuantileMapsMin;


  for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
    delete fEvts[iev];
    fEvts[iev] = 0;
  }
  
  

  
}
  


//_____________________________________________________________________________
void MCGBRIntegrator::SetMinCutSignificance(double x) {
  
//   fMinCutSignificance = x;
//   fMinCutSignificance = TMath::ChisquareQuantile(TMath::Erf(x/sqrt(2)),1)/2.0;
  fMinCutSignificance = x*x/2.0;
//   fMinCutSignificanceMulti = TMath::ChisquareQuantile(TMath::Erf(x/sqrt(2)),1)/2.0;
  
  
}



void MCGBRIntegrator::BuildQuantiles(int nvars, double sumabsw) {
 
//   std::sort(fEvts.begin(),fEvts.end(),MCGBRTgtCMP());
  
  double normwq = 0;
  for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
//     normwq += std::abs(fEvts[iev]->Weight()*fEvts[iev]->FuncVal());
//     normwq += std::abs(fEvts[iev]->Weight());
    normwq += 1.;
  }
  
  //parallelize building of quantiles for each input variable
  //(sorting of event pointer vector is cpu-intensive)
  #pragma omp parallel for
  for (int ivar=0; ivar<nvars; ++ivar) {
//     printf("sorting var %i\n",ivar);
        
    
    std::map<int,float,std::greater<float> > tmpmap;
    std::map<int,float,std::greater<float> > tmpmapmin;
    std::vector<MCGBREvent*> evtsvarsort(fEvts.begin(),fEvts.end());
    
    std::sort(evtsvarsort.begin(),evtsvarsort.end(),MCGBRVarCMP(ivar));
    
    double sumwq = 0;
    int lastquant = -1;
    for (unsigned int iev=0; iev<evtsvarsort.size(); ++iev) {
//       sumwq += std::abs(evtsvarsort[iev]->Weight());
      sumwq += 1.;
//       sumwq += std::abs(evtsvarsort[iev]->Weight()*evtsvarsort[iev]->FuncVal());
      int quant = int((sumwq/normwq)*(fNQuantiles-1));
      float val = evtsvarsort[iev]->Var(ivar);
          
      //ensure that events with numerically identical values receive the same quantile
      if (iev>0 && val==evtsvarsort[iev-1]->Var(ivar)) quant = evtsvarsort[iev-1]->Quantile(ivar);
    
//       printf("iev = %i, weight = %5f, sumabsw = %5f, sumwq = %5f, ivar = %i, val = %5f, quant = %i\n",iev,evtsvarsort[iev]->Weight(),sumabsw,sumwq,ivar,val,quant);
      
      evtsvarsort[iev]->SetQuantile(ivar,quant);
    
      tmpmap[quant] = val;
      
      if (quant!=lastquant) {
        tmpmapmin[quant] = val;
      }
      
      lastquant = quant;
      
//       printf("ivar = %i, iev = %i, val = %5f, quant = %i, fNQuantiles = %i\n",ivar,iev,evtsvarsort[iev]->Var(ivar),quant,fNQuantiles);
      
// //       adjust map to put boundaries in between values for individual events
//       if (iev>0 && quant>0 && evtsvarsort[iev-1]->Quantile(ivar) < quant) {
// //         printf("adjusting quantile %i, original = %5f, new = %5f\n",int(quant-1),tmpmap[quant-1]
//         tmpmap[evtsvarsort[iev-1]->Quantile(ivar)] = 0.5*(val + evtsvarsort[iev-1]->Var(ivar));
//       }
          
    }
    

    for (int i=0; i<fNQuantiles; ++i) {
      std::map<int,float,std::greater<float> >::const_iterator mit = tmpmap.lower_bound(i);
      
      float val;
      if (mit!=tmpmap.end()) val = mit->second;
      else val = -std::numeric_limits<float>::max();
      
      fQuantileMaps[ivar][i] = val;      
      
      std::map<int,float,std::greater<float> >::const_iterator mitmin = tmpmapmin.lower_bound(i);
      float valmin;
      if (mitmin!=tmpmapmin.end()) valmin = mitmin->second;
      else valmin = -std::numeric_limits<float>::max();
      
      fQuantileMapsMin[ivar][i] = valmin;
      
    }
    
    
    
  }    
  
}


void MCGBRIntegrator::TrainForest(int ntrees, bool reuseforest) {

  if (fDoEnvelope) {
    fStagedGeneration = true;
  }
    
//   int nvars = fNVars;
//   int nvarstrain = fNVars;
  
//   double sumwd = 0.;
//   double sumwu = 0.;

  int sincereset = 0; 
//   bool havereset = false;
  
  double totvolume = 1;
  for (int ivar=0; ivar<fNVars; ++ivar) {
    double low = fLimits[ivar].first;
    double high = fLimits[ivar].second;
    totvolume *= (high-low);
  }  
  
  printf("nvars = %i, totvolume = %5f\n",fNVars,totvolume);
  
//   if (!reuseforest || !fForest) {
//     fForest = new MCGBRForest; 
//   }
  
//   double shrinkageoriginal = fShrinkage;
  
  double maxval = -std::numeric_limits<double>::max();
  
  std::vector<MCGBREvent*> evtslast;

  if (fForest) delete fForest;
  fForest = new MCGBRForest;  

  if (fForestGen) delete fForestGen;
  fForestGen = new MCGBRForest;  
  
  double rmin = std::numeric_limits<double>::max();
  double rmax = std::numeric_limits<double>::lowest();
//   for (
  
//   fForest->SetInitialResponse(1e-2);
  
  for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
    double funcval = fEvts[iev]->FuncVal();
    if (funcval>rmax) rmax = funcval;
    if (funcval<rmin) rmin = funcval;
//     fEvts[iev]->SetTarget(fForest->InitialResponse());
  }
  
  
  fForest->SetInitialResponse(-128.);
//   fForest->SetInitialResponse(0.);
//   fForest->SetInitialResponseMin(0.);
//   fForest->SetInitialResponseMin(-log(100.));
  
//   fForestGen->SetInitialResponse(exp(-128.));
  fForestGen->SetInitialResponse(0.);
//   fForestGen->SetInitialResponseMin(0.);  
  
//   const double initialraw = exp(fForest->InitialResponse());
  
  for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
    fEvts[iev]->SetTarget(fForestGen->InitialResponse());
    fEvts[iev]->SetTargetMin(fForest->InitialResponse());
  }  
  
//   double shrinkagefactor = 0.;
  shrinkagefactor = 0.;
  std::vector<double> shrinkagefactors;
  
  std::vector<double> posprobboundstrees;
  std::vector<double> posprobboundsnodes;
  std::vector<unsigned int> nodeidxs;
//   double sumposprob = std::max(0.,fForest->InitialResponse()*totvolume);
  double sumposprob = std::abs(fForestGen->InitialResponse()*totvolume);
  
  double forestintegral = fForestGen->InitialResponse()*totvolume;
  fForestIntegralNow = forestintegral;
//   double lasttreeintegral = 0.;

  printf("forestintegral = %5f\n",forestintegral);
  
  posprobboundstrees.push_back(sumposprob);
  posprobboundsnodes.push_back(sumposprob);
  nodeidxs.push_back(0);
  
  
//   int niter = ntrees;
  int niter = 1e6;
  int numtrees = 0;
  int sameiter = 0;
      
  double sumwfuncval = 0;
  double sumwfuncvalsq = 0;
  double sumwint = 0.;
  
  double sumfull = 0.;
  double sumfullsq = 0.;
  double sumfullw = 0.;
  
//   double integral2now = 0.;
  
  std::vector<float> evalv(fNVars);
  
  for (int iter=0; iter<niter; ++iter) {
    
//     if (iter==(ntrees-1)) {
//       fMinEvents = 1;
//       fMinCutSignificance = 2.;
//     }
	  
    
    double sumw = 0.;
    int nev = fEvts.size();
    
//     int nev = std::min(int(1e5),int(fEvts.size()));
    std::vector<MCGBREvent*> evtstrimmed(fEvts.end()-nev,fEvts.end());
    fEvts.swap(evtstrimmed);
    
    for (std::vector<MCGBREvent*>::iterator it=fEvts.begin(); it!=fEvts.end(); ++it) {
      sumw += (*it)->Weight();
    }  
    
//     double sumwevts2 = 0.;
//     for (std::vector<MCGBREvent*>::iterator it=fEvts2.begin(); it!=fEvts2.end(); ++it) {
//       sumwevts2 += (*it)->Weight();
//     }      

    
   std::vector<std::pair<float,float> > limits(fLimits);

    
//     if (iter==0) {
//       if (fGenTree) delete fGenTree;
//       fGenTree = new MCGBRTreeD;
//       TrainTree(fEvts,sumw,*fGenTree,0,limits,false);
//     }
    
//     fGenTree = new MCGBRTreeD;
    
  
//     printf("setting targets\n");
  
//     printf("nev = %i, sumw = %5f\n",int(nev), sumw);

//     bool regen;
    
//     double integralnow = sumwfuncval/sumwint;
//     bool useflat = forestintegral > 1.1*integralnow || forestintegral < 0.;
//     bool useflat = iter < 100;
//     useflat = false;
//     useflat = true;
    bool useflat = false;
    
    printf("useflat = %i\n", int(useflat));   
   
    bool dologtrain = (numtrees<ntrees/2);
    dologtrain = true;

    if (dologtrain) {
      fForest->Trees().emplace_back();
    }
    MCGBRTreeD &tree = fForest->Trees().back();

    fForestGen->Trees().emplace_back();
    MCGBRTreeD &gentree = fForestGen->Trees().back();

    constexpr double cutoff = exp(-128.);
    
    int nevrnd = nev;
    std::vector<MCGBREvent*> rndevts;
    
    if (fNEventsBagged>0) {
      nevrnd = std::min(int(fEvts.size()),fNEventsBagged);
      rndevts.reserve(nevrnd);
      for (int iev=0; iev<nevrnd; ++iev) {
        int evidx = gRandom->Integer(fEvts.size());
        rndevts.push_back(fEvts[evidx]);
      }
    }
    else {
      rndevts = fEvts;
    }
    
    #pragma omp parallel for simd
    for (int iev=0; iev<nevrnd; ++iev) {
      const double target = rndevts[iev]->Target();
      const double targetmin = rndevts[iev]->TargetMin();
//       const double funcval = rndevts[iev]->FuncVal();
      const double funcvalalt = rndevts[iev]->FuncValAlt();
      
      _tgtvals[iev] = vdt::fast_log(std::max(cutoff,vdt::fast_exp(targetmin)-target));
//       _tgtvals[iev] = vdt::fast_log(std::max(cutoff,funcval-target));
      _tgt2vals[iev] = funcvalalt - targetmin;
    }  
    
    for (int iev=0; iev<nevrnd; ++iev) {
      rndevts[iev]->SetArg(_tgtvals[iev]);
      rndevts[iev]->SetArgLog(_tgt2vals[iev]);
    }
    
/*    #pragma omp parallel for
    for (int iev=0; iev<nev; ++iev) {
      const double target = fEvts[iev]->Target();
      const double targetmin = fEvts[iev]->TargetMin();
//       const double funcval = fEvts[iev]->FuncVal();
      const double funcvalalt = fEvts[iev]->FuncValAlt();
      
      _tgtvals[iev] = vdt::fast_log(std::max(cutoff,vdt::fast_exp(targetmin)-target));
//       _tgtvals[iev] = vdt::fast_log(std::max(cutoff,funcval-target));
      _tgt2vals[iev] = funcvalalt - targetmin;
    }  */  
    
//     for (int iev=0; iev<nev; ++iev) {
//       fEvts[iev]->SetArg(_tgtvals[iev]);
//       fEvts[iev]->SetArgLog(_tgt2vals[iev]);
//     }
    
    printf("training tree\n");
    if (dologtrain) {
//       TrainTree(fEvts,sumw,tree,0,limits,true,false);
// #pragma omp parallel
      TrainTree(rndevts,sumw,tree,0,limits,true,false);
      
      for (MCGBREvent *evt : fEvts) {
        for (int ivar=0; ivar<fNVars; ++ivar) {
          evalv[ivar] = evt->Var(ivar);
        }
        double responsemin = tree.GetResponse(evalv.data());
        evt->SetTargetMin(evt->TargetMin() + responsemin);

      }
      
    }
//     TrainTree(fEvts,sumw,gentree,0,limits,true,true);   
// #pragma omp parallel
    TrainTree(rndevts,sumw,gentree,0,limits,true,true);   
    for (MCGBREvent *evt : fEvts) {
      for (int ivar=0; ivar<fNVars; ++ivar) {
        evalv[ivar] = evt->Var(ivar);
      }
      double response = gentree.GetResponse(evalv.data());
      evt->SetTarget(evt->Target() + response);

    }
    
    
    fGenTree = &fForest->Trees().front();

    
    double shrinkagefactordiff = fShrinkage*(1.-shrinkagefactor);
    double shrinkagesecondarydiff = fShrinkage*(shrinkagefactor - fShrinkageFactorSecondary);
    
    shrinkagefactor += shrinkagefactordiff;
    shrinkagefactors.push_back(shrinkagefactordiff);
    
    fShrinkageFactorSecondary += shrinkagesecondarydiff;
    
    
    

    printf("iter = %i, fulliter = %i, sameiter = %i, responses = %i, generesponses = %i, shrinkagefactor = %5f, fShrinkageFactorSecondary = %5f, sigmascale = %5e\n",iter,numtrees,sameiter,int(tree.Responses().size()),int(gentree.Responses().size()),shrinkagefactor,fShrinkageFactorSecondary,1./(1-fShrinkageFactorSecondary));
    
    posprobboundsnodes.reserve(posprobboundsnodes.size()+gentree.Responses().size());
    nodeidxs.reserve(nodeidxs.size()+gentree.Responses().size());
    
    std::vector<double> posprobboundsnodesnow;
    std::vector<unsigned int> nodeidxsnow;
    double sumposprobnow = 0.;
    
    posprobboundsnodesnow.reserve(gentree.Responses().size());
    nodeidxsnow.reserve(gentree.Responses().size());
    

    
    
    double treeintegral = 0.;
    for (unsigned int inode = 0; inode<gentree.Responses().size(); ++inode) {
      double responsemin = gentree.GetResponseMin(inode);
      double response = gentree.GetResponse(inode);
      double volume = 1.0;
      for (int ivar=0; ivar<fNVars; ++ivar) {
        double low = gentree.Limits()[inode][ivar].first;
        double high = gentree.Limits()[inode][ivar].second;
        volume *= (high-low);
      }
      double prob = response*volume;
      double probmin = responsemin*volume;
      double absprob = std::abs(prob);
//       double absprob = std::max(0.,prob);
      double absprobmin = std::abs(probmin);
      forestintegral += prob;
      treeintegral += prob;
      if (absprob>0.) {
        sumposprob += absprob;
        posprobboundsnodes.push_back(sumposprob);
        nodeidxs.push_back(inode);
        
     
      }
      if (absprobmin>0.) {
        sumposprobnow += absprobmin;
        posprobboundsnodesnow.push_back(sumposprobnow);
        nodeidxsnow.push_back(inode);           
      }
      
//       printf("iter = %i, inode = %i, response = %5e, volume = %5e, prob = %5f, treeintegral = %5f, forestintegral = %5f\n",iter,inode,response,volume,prob, treeintegral,forestintegral);
    }
    fForestIntegralNow = forestintegral;
//     lasttreeintegral = treeintegral;
    printf("treeintegral = %5e, forestintegral = %5e\n",treeintegral,forestintegral);
    posprobboundstrees.push_back(sumposprob);
    
    
    
    
    ++sameiter;
//     int period = 2./fShrinkage;
//     if (iter%period!=(period-1)) continue;
//     if (tree.Responses().size()>1 && sameiter<100) continue;
    
    sameiter=0;
    ++numtrees;
    
    
    ++sincereset;
    
    
//     if (numtrees>600 && numtrees<ntrees) continue;
    
//     if (iter>300 && numtrees<ntrees) continue;
//     if (numtrees>75 && numtrees<ntrees) continue;
//     if (numtrees>(ntrees-100) && numtrees<ntrees) continue;

    
//     if (iter == (niter-1)) break;
    
//     if (iter == (niter-1)) {
    if (0) {
    
      int nsample = 100;
      std::vector<float> samplevars(fNVars);
      for (unsigned int inode = 0; inode<fGenTree->Responses().size(); ++inode) {
        double maxfuncval = fGenTree->Responses()[inode];
        for (int iev=0; iev<nsample; ++iev) {
          for (int ivar=0; ivar<fNVars; ++ivar) {
            double low = fGenTree->Limits()[inode][ivar].first;
            double high = fGenTree->Limits()[inode][ivar].second;
            double val = gRandom->Uniform(low,high);
            samplevars[ivar] = val;
          }
          double funcval = 1.;
          for (int ivar=0; ivar<fNVars; ++ivar) {
            double val = samplevars[ivar];
            funcval *= exp(-0.5*val*val);
          }
          if (funcval>maxfuncval) {
            maxfuncval = funcval;
          }
          
        }
        fGenTree->Responses()[inode] = maxfuncval;
      }
    }

    
// //     std::vector<double> pdfvals(fGenTree->Responses().size());
//     std::vector<double> probbounds(fGenTree->Responses().size());
//     double sumprob = 0.;
//     double maxresponse = -std::numeric_limits<double>::max();
//     for (unsigned int inode = 0; inode<fGenTree->Responses().size(); ++inode) {
//       double response = fGenTree->GetResponse(inode);
//       maxresponse = std::max(response,maxresponse);
//       double volume = 1.0;
//       for (int ivar=0; ivar<fNVars; ++ivar) {
//         double low = fGenTree->Limits()[inode][ivar].first;
//         double high = fGenTree->Limits()[inode][ivar].second;
//         volume *= (high-low);
//       }
//       double prob = response*volume;
//       sumprob += prob;      
// //       pdfvals[inode] = prob;
//       probbounds[inode] = sumprob;      
//     }
    

    
    int nvars = fNVars;
    double sumabsw = 0.;
    
    for (int iev=0; iev<nev; ++iev) {
      sumabsw += std::abs(fEvts[iev]->Weight());
    }
  
  
  
   for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
     maxval = std::max(maxval,fEvts[iev]->FuncVal());
   }
    
        
//     int nevold = nev;

//     if (iter==0) {
//       nev += fNEvents;
//       fEvts.reserve(nev);
//     fNQuantiles = std::min(nev+1,std::numeric_limits<unsigned short>::max()+1);
//     }
//     else {
//       fEvts.clear();
//     }
    
    printf("generating events\n");
    
    int nevold = nev;

    bool dogenevents = (numtrees<ntrees/2);
//     dogenevents = (numtrees<150);
    dogenevents = true;
//     fEvts.clear();
    bool islast = false;
    if (numtrees<ntrees) {
      ;
      if (dogenevents) {
        nev += fNEvents;
      }
    }
    else {
//       nev += fNEvents;
//       nev += 100*fNEvents;
      nev += 20e3;
      islast = true;
//       nev += 2*1000*1000;
    }
    
//     if (islast) {
//       fStagedGeneration = true;
//     }
    
//     if (fDoEnvelope) {
// //       fMinEvents = std::max(10,nev/300);
//       fMinEvents = 1e3;
//     }
    
    std::vector<double> weightv;
    weightv.reserve(nev-nevold);
    
//     nev += fNEvents/100;
    fEvts.reserve(nev);    
//     fNQuantiles = nev+2;
    fNQuantiles = std::min(nev,std::numeric_limits<unsigned short>::max()+1);
    
    evtslast.clear();
    double maxratio = 0.;
    double maxratiofuncval = 0;
    double minratio = std::numeric_limits<double>::max();
    double efficiency = 0;
    double maxdiff = std::numeric_limits<double>::lowest();
//     double integral = 0.;
    
    int nattempts = 0;
    int ngen = 0;
//     double intscale = 1./double(nev);    
    
    std::vector<double> tmpvals(fNVars);
    
    double sumwnow = 0.;
    double sumw2now = 0.;
    double sumwfuncvalnow = 0;
    double sumwfuncvalsqnow = 0;    
    
    double sum1 = 0.;
    double sum1sq = 0.;
    double sumw1 = 0.;
    
    double sum1t = 0.;
    double sum1tsq = 0.;

    double sumct = 0.;
    double sumctw = 0.;
    double sumctsq = 0.;
    
    double sum2 = 0.;
    double sum2sq = 0.;
    double sumw2 = 0.;    
    
    double sumw2t = 0.;
    double sum2t = 0.;
    double sum2tsq = 0.;
    
    double suminter = 0.;
    double sumintersq = 0.;
    
    double maxsigmaratio = std::numeric_limits<double>::lowest();
    
//     double sumdiff = 0.;
//     double sumdiffsq = 0.;
//     double sumdiffw = 0.;
//     
//     for (unsigned int iev=0; iev<1e3; ++iev) {
//       for (int ivar=0; ivar<fNVars; ++ivar) {   
//         double low = fLimits[ivar].first;
//         double high = fLimits[ivar].second;             
//         double val = gRandom->Uniform(low,high);
//         evalv[ivar] = val;        
//       }     
//       
//               
//       double targetmin = fForest->GetResponse(evalv.data());
//       double target = fForestGen->GetResponse(evalv.data());
//       
//       double weight = 1.;
//       double diff = exp(targetmin) - target;
//       
//       sumdiff += weight*diff;
//       sumdiffsq += weight*diff*diff;
//       sumdiffw += weight;
//     }
//     
//     double wavgdiff = sumdiff/sumdiffw;
//     double integral2diff = wavgdiff + forestintegral;
//     
//     double sigmawdiff = sqrt(sumdiffsq/sumdiffw - wavgdiff*wavgdiff);
//     double sigmaintegral2diff = sigmawdiff/sqrt(sumdiffw);
//     
//     printf("wavgdiff = %5f, sigmawdiff = %5f, integral2diff = %5f +- %5f\n",wavgdiff,sigmawdiff,integral2diff,sigmaintegral2diff);

    
//     double scaleorigin = forestintegral > 0. && integral2now > 0. ? 1.5*integral2now/forestintegral : 1.;
//     double scaleorigin = islast ? 100. : 1.5;
//     double scaleorigin = islast ? 1.5 : 1.5;
    double scaleorigin = fStagedGeneration ? 4. : 1.;
//     if (fStagedGeneration) {
//       scaleorigin = std::min(100.,std::max(1.,std::max(4.*fIntegralRatio, 4.*fSigmaRatio*fSigmaRatio)));
//     }
    
    printf("scaleorigin = %5f, fIntegralRatio = %5f, fSigmaRatio = %5f\n",scaleorigin, fIntegralRatio,fSigmaRatio);
//     scaleorigin = 1.;
    
    while (int(fEvts.size())<nev) {
      ++nattempts;
//       fEvts.push_back(new MCGBREvent(nvars));
//       MCGBREvent &evt = *fEvts.back();
//       evtslast.push_back(&evt);
      
//       printf("finding tree\n");
//       double genpdf;
      double target;
//       double targetpos;
      double targetmin = 0.;
      double gentarget = 0.;


      double probuniform = gRandom->Uniform();
      bool douniform = probuniform > shrinkagefactor;
      douniform = false;
      
//           printf("full\n");
      double prob = gRandom->Uniform(0.,sumposprob);

      
      
      std::vector<double>::const_iterator treeit = std::upper_bound(posprobboundstrees.begin(),posprobboundstrees.end(),prob);
      unsigned int itree = treeit - posprobboundstrees.begin();  
      int treeidx = int(itree) - 1;
      
      
      
//           printf("finding node, treeidx = %i, itree = %i, prob = %5f, forest size = %i, vsize = %i\n",treeidx,itree,prob,int(fForest->Trees().size()),int(posprobboundstrees.size()));
      
      std::vector<double>::const_iterator nodeit = std::upper_bound(posprobboundsnodes.begin(),posprobboundsnodes.end(),prob);
      unsigned int nodeidx = nodeit - posprobboundsnodes.begin();
      unsigned int inode = nodeidxs[nodeidx];          
      
//           printf("inode = %i, tree size = %i\n",inode, int( fForest->Trees()[treeidx].Responses().size()));
      
      for (int ivar=0; ivar<fNVars; ++ivar) {
        double low =  treeidx >= 0 && !douniform ? fForestGen->Trees()[treeidx].Limits()[inode][ivar].first : fLimits[ivar].first;
        double high = treeidx >= 0 && !douniform ? fForestGen->Trees()[treeidx].Limits()[inode][ivar].second : fLimits[ivar].second;    
/*            double low = fLimits[ivar].first;
        double high = fLimits[ivar].second;  */           
        double val = gRandom->Uniform(low,high);
        evalv[ivar] = val;
      }       
      

      
//           double targetpos;
//           fForest->ResponseGEQ(evalv.data(),targetmin,targetpos);
//           target = fForest->GetResponse(evalv.data(),targetmin,targetpos);
      
//           double targetposnow = std::abs(tree.GetResponse(evalv.data()));
//       double targetposnow = 1.;
    
//           genpdf = fullfrac*targetpos/sumposprob + nowfrac*targetposnow/sumposprobnow;
    
//           if (false) {
//           if (numtrees==ntrees) {
    
//             printf("accept reject\n");
//           if (0) {
//             printf("targetpos = %5e, exptgtmin = %5e\n",targetpos,exp(targetmin));
//             double genorigin = initialraw;
//             gentarget = std::min(targetpos,vdt::fast_exp(targetmin));
//             gentarget = target;
//             gentarget = vdt::fast_exp(target);
//             double genorigin = 1.;

      double dummy1, dummy2;

      targetmin = fForest->GetResponse(evalv.data(),dummy1,dummy2);
      target = fForestGen->GetResponse(evalv.data(),dummy1,dummy2);
      
      
//       double genorigin = log(510.) + 10.;
//       double genorigin = islast ? 510. : 10.;
//       gentarget = exp(target);
//       gentarget = islast ? exp(target) : std::min(1.,target+10.);
//       gentarget = islast ? exp(target) : targetmin;
      
      
      
      double genorigin = scaleorigin*target;
      gentarget = genorigin;
      if (fStagedGeneration) {
        gentarget = std::min(exp(targetmin),genorigin);
      }
      double realtarget = gentarget;
      
//       double realtarget = exp(target);
//       gentarget = std::min(target,exp(targetmin));
      
      double weight = 1.;
      double ratio1 = gentarget/genorigin;
      

      sumw1 += weight;
      sum1 += weight*ratio1;
      sum1sq += weight*ratio1*ratio1;
      
      double ratio1t = exp(targetmin)/genorigin;
      sum1t += weight*ratio1t;
      sum1tsq += weight*ratio1t*ratio1t;
      
      double funcvaltest = 1.;
      for (int ivar=0; ivar<nvars; ++ivar) {
        tmpvals[ivar] = evalv[ivar];
      }
      funcvaltest = Camel(fNVars,tmpvals.data());
      
      double ratioct = forestintegral*funcvaltest/target;
      sumct += ratioct;
      sumctsq += ratioct*ratioct;
      sumctw += weight;
      
      
      double ratiointer = exp(targetmin)/genorigin;
      suminter += weight*ratiointer;
      sumintersq += weight*ratiointer*ratiointer;
      
      
// //       double sigmaratio = (exp(targetmin) - target)/exp(targetmin);
//       double sigmaratio = (exp(targetmin) - target)/exp(targetmin);
//       if (std::isnormal(sigmaratio) && sigmaratio>maxsigmaratio) {
//         maxsigmaratio = sigmaratio;
//       }
      
//       double diff1 = (gentarget - target)*forestintegral/target;
//       sumct += diff1;
//       sumctsq += diff1*diff1;
      
      double rnd = gRandom->Uniform(0.,genorigin);
      
//       printf("genorigin = %5e, gentarget = %5e, rnd = %5e\n",genorigin,gentarget,rnd);
      
//             bool pass = fForest->ResponseGEQ(evalv.data(),vdt::fast_log(rnd),target);
      
//             printf("rnd = %5e, target = %5e, pass = %i\n",rnd,target,int(pass));
      
//             if (!pass) continue;
      
//             gentarget = vdt::fast_exp(target);
//             gentarget = target;
      
//       if (gentarget < rnd) continue;
//             if (vdt::fast_exp(targetmin)<rnd) continue;        
//             if (target<rnd) continue;        
//             if (1./target < rnd) continue;
//       genpdf = target/forestintegral;
      
//             genpdf = targetmin/forestintegral;
        
        
      if (gentarget < rnd) continue;

      ++ngen;

      
//         
        
        
      
//       printf("filled vars\n");
      

      
      
      fEvts.push_back(new MCGBREvent(nvars));
      MCGBREvent &evt = *fEvts.back();
      
      for (int ivar=0; ivar<fNVars; ++ivar) {
        evt.SetVar(ivar,evalv[ivar]);
      }
      
      evt.SetTarget(target);
      evt.SetTargetMin(targetmin);
      
      


      
//       fEvts2.push_back(&evt);
      evtslast.push_back(&evt);   

      
//       evt.SetWeight(sumposprob/target);
      
//       if (numtrees==ntrees) {
//         printf("target = %5e, targetmin = %5e\n",target,targetmin);
//       }
      
//       evt.SetTargetMin(targetpos);
      
//       if (target > 1.) {
//         printf("gen target = %5f\n",target);
//       }

//       evt.SetTarget(target);

//       double weight = sumprob/(totvolume*response);
//       double weight = 1.0;
//       printf("sumprob = %5f, totvolume = %5f, probs[inode] = %5f, weight = %5f\n",sumprob,totvolume, probs[inode],weight);
//       evt.SetWeight(weight);
      double funcval = 1.;
      for (int ivar=0; ivar<nvars; ++ivar) {
        double val = evt.Var(ivar);
        funcval *= exp(-0.5*val*val);
        tmpvals[ivar] = val;
      }
      funcval = islast ? Camel(fNVars,tmpvals.data()) : Camelrnd(fNVars,tmpvals.data());
//       funcval = Camel(fNVars,tmpvals.data());
      
      
//       double sigmaratio = (exp(targetmin) - target)/exp(targetmin);
      double sigmaratio = (exp(targetmin) - target)/funcval;
      if (std::isnormal(sigmaratio) && sigmaratio>maxsigmaratio) {
        maxsigmaratio = sigmaratio;
      }      
      
      double ratio2 = funcval/gentarget;
//       double funcvalmod = funcval > 100. ? std::min(funcval,1000.) - 100. : 0.;
//       double ratio2 = funcvalmod/gentarget;
      
      sumw2 += weight;
      sum2 += weight*ratio2;
      sum2sq += weight*ratio2*ratio2;
      
      weightv.push_back(weight*ratio2);
      
      double weight2t = exp(targetmin)/gentarget;
      double ratio2t = funcval/exp(targetmin);
      sumw2t += weight2t;
      sum2t += weight2t*ratio2t;
      sum2tsq += weight2t*ratio2t*ratio2t;
      
//       double ratioct = forestintegral*funcval/target;
//       sumct += ratioct;
//       sumctsq += ratioct*ratioct;
//       sumctw += weight;
      
//       double ratiofull = forestintegral*funcval/target;
//       sumfull += weight*ratiofull;
//       sumfullsq += weight*ratiofull*ratiofull;
//       sumfullw += weight;
      
//       double wun = 1./gentarget;
//       sumctw += wun;
//       sumct += funcval*wun;
      
//       evt.SetWeight(1./funcval/funcval);

      
//       if (numtrees<ntrees) {
//       if (true) {
//         funcval = funcval + 0.05 + 0.01*gRandom->Gaus();
//       }
//       evt.SetFuncVal(funcval);

//       if (numtrees<ntrees) {
//         funcval = funcval + 0.01*gRandom->Gaus();
//       }
      evt.SetFuncVal(funcval);
      sumabsw += std::abs(evt.Weight());  
      
      evt.SetFuncValAlt(vdt::fast_log(funcval));
//       double logsigma = log(funcval) + 1.*gRandom->Gaus();
//       evt.SetFuncValAlt(exp(logsigma));
//       evt.SetFuncValAlt(exp(0.1*gRandom->Gaus()));

      
//       double fullweight = (1./totvolume)/genpdf;
//       evt.SetWeight(fullweight);
      
//       printf("fullweight = %5f, target = %5f, forestintegral = %5f\n",fullweight,target,forestintegral);
      
//       double wfunc = fullweight*totvolume*funcval;
      double intweight = 1./realtarget;
      double wfunc = funcval*intweight;
      sumwfuncval += wfunc;
      sumwfuncvalsq += wfunc*wfunc;      
      sumwfuncvalnow += wfunc;
      sumwfuncvalsqnow += wfunc*wfunc;          
      sumwint += intweight;
      

      
      double ratio = funcval/realtarget;
//       if (!std::isnormal(ratio) ) {
//         printf("funcval = %5e, gentarget = %5e, target = %5e, targetmin = %5e, exp(targetmin) = %5e\n",funcval,gentarget,target,targetmin,exp(targetmin));
//       }
//       printf("funcval = %5e, gentarget = %5e, ratio = %5e\n",funcval,gentarget,ratio);
//       double ratio = funcval/target;
      double diff = funcval - realtarget;
      if (ratio>maxratio) {
        maxratio = ratio;
        maxratiofuncval = funcval;
      }
      if (ratio<minratio) minratio = ratio;
      
      if (diff>maxdiff) {
        maxdiff = diff;
      }
      
      sumwnow += ratio;
      sumw2now += ratio*ratio;
      
      
      efficiency += ratio;

//       if (numtrees==ntrees) {
//         printf("funcval = %5f, target = %5f, ratio = %5f\n",funcval,target,ratio);
//       }

      
//       printf("nattempts = %i, fEvtsSize = %i\n",nattempts,int(fEvts.size()));
      
//       printf("iev = %i, weight = %5f, sumabsw = %5f\n",iev,weight,sumabsw);
    }    
    
//     double wavg = sumwnow/double(ngen);
//     double sigmaw = sqrt(sumw2now/double(ngen) - wavg*wavg);
    
//     double integralnow = sumwfuncvalnow/sumwint;
//     double integralerrnow = sqrt(sumwfuncvalsqnow/sumwint-integralnow*integralnow)/sqrt(sumwint);
    
    double wavg1 = sum1/sumw1;
    double sigmaw1 = sqrt(sum1sq/sumw1 - wavg1*wavg1);
    double sigmarel1 = sigmaw1/wavg1;
    
    double wavg1t = sum1t/sumw1;
    double sigmaw1t = sqrt(sum1tsq/sumw1 - wavg1t*wavg1t);
//     double sigmarel1t = sigmaw1t/wavg1t;
    
    double wavg2 = sum2/sumw2;
    double sigmaw2 = sqrt(sum2sq/sumw2 - wavg2*wavg2);
    double sigmarel2 = sigmaw2/wavg2;
    
    double integral1 = scaleorigin*forestintegral;
    double integral2 = integral1*wavg1;
//     integral2now = integral2;
    double integralfunc = integral2*wavg2;
    double sigmaintegral2 = integral1*sigmaw1/sqrt(sumw1);
    double sigmaintegralfuncpart = integral2*sigmaw2/sqrt(sumw2);
    
    double sigmaintegralfunc = sqrt(pow(sigmaintegralfuncpart/integralfunc,2)+pow(sigmaintegral2/integral2,2))*integralfunc;
    
    double wavg2t = sum2t/sumw2t;
    double sigmaw2t = sqrt(sum2tsq/sumw2t - wavg2t*wavg2t);
    double sigmarel2t = sigmaw2t/wavg2t;
    
    double integral2t = integral1*wavg1t;
    double sigmaintegral2t = integral1*sigmaw1t/sqrt(sumw1);
    
    double wavginter = suminter/sumw1;
    double sigmawinter = sqrt(sumintersq/sumw1 - wavginter*wavginter);
    double sigmarelinter = sigmawinter/wavginter;
    
    fSigmaRatio = sigmarel1/sigmarel2t;
    fIntegralRatio = integral2/forestintegral;
    
//     fSigmaScale = 1./(1-std::min(forestintegral/integral2t,shrinkagefactor));
    fSigmaScale = std::max(1.,std::min(1./maxsigmaratio,1./(1.-shrinkagefactor)));
    
//     fSigmaScale = 1./(1.-std::min(forestintegral/integralfunc,shrinkagefactor));
    
//     double ctavg = sumct/sumw1;
//     double integral2alt = ctavg + forestintegral;
//     double sigmaintegral2alt = sqrt(sumctsq/sumw1 -ctavg*ctavg)/sqrt(sumw1);
    
//     double integralfuncalt = sumct/sumctw;
    
//     double wavg3 = sum3/sumw2;
//     double sigmaw3 = sqrt(sum3sq/sumw2 - wavg3*wavg3);
//     double integral2alt = integral1*sumw2/sum3;
//     double integralfuncalt = integral1*sum2/sum3;
    
    double wavgct = sumct/sumctw;
    double sigmawct = sqrt(sumctsq/sumctw - wavgct*wavgct);
    double sigmarelct = sigmawct/wavgct;
    
//     integralfunc = wavgct;
//     sigmaintegralfunc = sigmawct/sqrt(sumctw);
    
    double weightfull = 1./sigmaintegralfunc/sigmaintegralfunc;
    sumfull += weightfull*integralfunc;
    sumfullsq += weightfull*integralfunc*integralfunc;
    sumfullw += weightfull;
    double wavgfull = sumfull/sumfullw;
//     double sigmawfull = sqrt(sumfullsq/sumfullw - wavgfull*wavgfull);
    double sigmafull = sqrt(sumfullw)/sumfullw;
    
    double sigmarelequiv = sqrt(sigmarel2*sigmarel2 + sigmarel1*sigmarel1/scaleorigin/scaleorigin);
    
//     double wavgfull = sumfull/sumfullw;
//     double sigmawfull = sqrt(sumfullsq/sumfullw - wavgfull*wavgfull);
//     double sigmafull = sigmawfull/sqrt(sumfullw);
    
    
//     fSigmaScale = 1./(1.-std::min(forestintegral/wavgfull,shrinkagefactor));
    
    const double epsilon = 5e-4;
    std::sort(weightv.begin(), weightv.end(), std::greater<double>());
    double maxwthres = epsilon*sum2;
    double runningsum2 = 0.;
    double maxweightepsilon = 0.;
    for (double weight : weightv) {
      runningsum2 += weight;
      if (runningsum2 >= maxwthres) {
        maxweightepsilon = weight;
        break;
      }
    }
    
    double effepsilon = sum2/sumw2/maxweightepsilon;
    
    
//     if (0) {
    printf("integral2t = %5e +- %5e, maxsigmaratio = %5f, fSigmaScale = %5f\n",integral2t,sigmaintegral2t,maxsigmaratio,fSigmaScale);
//     double sigmascale = 1./
//     }
//     printf("sumct = %5f, sumctw = %5f\n",sumct,sumctw);
    printf("sigmarelct = %5f, sigmarel2t = %5f, sigmarelinter = %5f, sigmarelequiv = %5f, maxweightepsilon = %5f, effepsilon = %5f\n", sigmarelct,sigmarel2t,sigmarelinter, sigmarelequiv, maxweightepsilon,effepsilon);
    printf("testgeneff = %5f, wavg1 = %5f, wavg2 = %5f, sigmaw1 = %5f, sigmarel1 = %5f, sigmaw2 = %5f, sigmarel2 = %5f, integral2 = %5f +- %5f, integralfunc = %5f +- %5f, integralfull = %5f +- %5f\n",double(ngen)/double(nattempts),wavg1, wavg2, sigmaw1,sigmarel1, sigmaw2, sigmarel2, integral2, sigmaintegral2,integralfunc, sigmaintegralfunc, wavgfull, sigmafull);
    printf("numevents = %i\n",int(fEvts.size()));
    
//     printf("testgeneff = %5f, wavg = %5f, sigmaw = %5f, sigmarel = %5f, integral = %5f += %5f\n",double(ngen)/double(nattempts),wavg,sigmaw,sigmaw/wavg,integralnow, integralerrnow);
//     printf("testgeneff = %5f\n",double(nev)/double(nattempts));
//     double testgeneff = double(ngen)/double(nattempts);
    

    
//     if (false && testgeneff<0.1 && sincereset > 5./fShrinkage) {
//     if (false && numtrees==(ntrees/2)) {      
    if (false && numtrees%300==0 && numtrees>1) {
//       havereset = true;
      sincereset = 0;
      printf("resetting target\n");
      for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
        fEvts[iev]->SetTarget(fForest->InitialResponse());
      }
      for (unsigned int itree=0; itree<fForest->Trees().size(); ++itree) {
        for (unsigned int inode=0; inode<fForest->Trees()[itree].Responses().size(); ++inode) {
          fForest->Trees()[itree].Responses()[inode] = 0.;
        }
      }
      for (unsigned int iel=1; iel<posprobboundstrees.size(); ++iel) {
        posprobboundstrees[iel] = 0.;
      }
      
//       forestintegral = 0.;
//       sumposprob = 0.;
      posprobboundsnodes.clear();
//       posprobboundstrees.clear();
      nodeidxs.clear();
      
      sumposprob = std::abs(fForest->InitialResponse()*totvolume);
      forestintegral = fForestGen->InitialResponse()*totvolume;
      fForestIntegralNow = forestintegral;

//       posprobboundstrees.push_back(sumposprob);
      posprobboundsnodes.push_back(sumposprob);
      nodeidxs.push_back(0);      
      
    }
    

    
    
    
    
//     for (int iev=fEvts.size(); iev<nev; ++iev) {
//       fEvts.push_back(new MCGBREvent(nvars));
//       MCGBREvent &evt = *fEvts.back();
//       evtslast.push_back(&evt);
//       unsigned int inode = 0;
//       while (true) {
//         double prob = gRandom->Uniform(0.,sumprob);
//         std::vector<double>::const_iterator nodeit = std::upper_bound(probbounds.begin(),probbounds.end(),prob);
//         inode = nodeit - probbounds.begin();
// //         printf("probs[inode] = %5f\n",probs[inode]);
//         if (fGenTree->GetResponse(inode)>0.) break;
//       }
//       for (int ivar=0; ivar<fNVars; ++ivar) {
//         double low = fGenTree->Limits()[inode][ivar].first;
//         double high = fGenTree->Limits()[inode][ivar].second;
//         double val = gRandom->Uniform(low,high);
//         evt.SetVar(ivar,val);
//         evalv[ivar] = val;
//       }
// //       double weight = (1./totvolume)/(probs[inode]/sumprob);
//       double response = fGenTree->GetResponse(inode);
//       double target = fForest->GetResponse(evalv.data());
//       
// //       if (target > 1.) {
// //         printf("gen target = %5f\n",target);
// //       }
// 
//       evt.SetTarget(target);
// 
// //       double weight = sumprob/(totvolume*response);
//       double weight = 1.0;
// //       printf("sumprob = %5f, totvolume = %5f, probs[inode] = %5f, weight = %5f\n",sumprob,totvolume, probs[inode],weight);
//       evt.SetWeight(weight);
//       double funcval = 1.;
//       for (int ivar=0; ivar<nvars; ++ivar) {
//         double val = evt.Var(ivar);
//         funcval *= exp(-0.5*val*val);
//       }
// //       evt.SetFuncVal(funcval);
//       evt.SetFuncVal(funcval);
//       sumabsw += std::abs(evt.Weight());  
//       
//       double ratio = funcval/response;
//       if (ratio>maxratio) {
//         maxratio = ratio;
//         maxratiofuncval = funcval;
//       }
//       if (ratio<minratio) minratio = ratio;
//       
//       efficiency += ratio;
// 
//       
// //       printf("iev = %i, weight = %5f, sumabsw = %5f\n",iev,weight,sumabsw);
//     }
    efficiency = efficiency/double(nev-nevold)/maxratio;
//        efficiency = efficiency/double(nev)/maxratio;
    
    
//     sumw = 0.;
//     double sumwfuncval = 0;
//     double sumwfuncvalsq = 0;
//     double intscale = 1./double(fEvts.size());
//     for (std::vector<MCGBREvent*>::iterator it=fEvts.begin(); it!=fEvts.end(); ++it) {
// //       sumw += (*it)->Weight();
//       double wfunc = (*it)->Weight()*totvolume*(*it)->FuncVal()*intscale;
//       sumwfuncval += wfunc;
//       sumwfuncvalsq += wfunc*wfunc;
//     }
    double integral = sumwfuncval/sumwint;
    double integralerr = sqrt(sumwfuncvalsq/sumwint-integral*integral)/sqrt(sumwint);
    printf("integral = %5f +- %5f, maxratio = %5f, maxratiofuncval = %5f, minratio = %5f, efficiency = %5f, maxdiff = %5f\n",integral, integralerr,maxratio,maxratiofuncval,minratio,efficiency,maxdiff);
    
//     if (forestintegral > 1.1*integral) {
//       double scale = 0.01;
//       sumposprob *= scale;
//       forestintegral *= scale;
//       for (unsigned int inode=0; inode<posprobboundstrees.size(); ++inode) {
//         posprobboundstrees[inode]*=scale;
//       }
//       for (unsigned int inode=0; inode<posprobboundsnodes.size(); ++inode) {
//         posprobboundsnodes[inode]*=scale;
//       }      
//       for (unsigned int itree=0; itree<fForest->Trees().size(); ++itree) {
//         MCGBRTreeD &modtree = fForest->Trees()[itree];
//         for (unsigned int inode=0; inode<modtree.Responses().size(); ++inode) {
//           modtree.Responses()[inode]*=scale;
//         }
//       }
//     }
    
//     if (maxratio<1.1) {
//       fShrinkage = 0.2;
//     }
//     else {
//       fShrinkage = shrinkageoriginal;
//     }
    
    if (numtrees==ntrees) break;
    
    
//     if (numtrees<ntrees/2) {
    if (0) {
//       const double decayrate = 0.9;
//       const double decayrate = pow(integral/forestintegral,0.5*fShrinkage);
//       const double decayrate = integral/forestintegral;
//       double decayrate = pow(1. - lasttreeintegral/forestintegral,0.95);
//       double decayrate = fShrinkage;
//       double newforestintegral = forestintegral - (1.-fShrinkage)*treeintegral;
//       newforestintegral = std::min(newforestintegral,0.9*wavgfull);
//       double decayrate  = std::max(0.5*treeintegral/forestintegral,0.5);
//       double oldforestintegral = forestintegral-treeintegral;
//       double newforestintegral = forestintegral - 0.5*treeintegral;
//       double newforestintegral = forestintegral - 2.*std::max(0.,forestintegral-wavgfull);
//       double decayrate = newforestintegral/forestintegral;
//       double decayrate = std::min((1.-fShrinkage)*integral2t/forestintegral,1.);
//       double decayrate = forestintegral > (1.+0.5*fShrinkage)*integral2t ? (1.-fShrinkage)*integral2t/forestintegral : 1.;
      double decayrate  = std::min(1.,pow(integral2t/forestintegral,2.));
//       double decayrate =  0.96;
//       decayrate = 0.5;
      printf("decaying responses, decayrate = %5f\n",decayrate);
      for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
        fEvts[iev]->SetTarget(decayrate*fEvts[iev]->Target());
//         fEvts[iev]->SetTargetMin(decayrate*fEvts[iev]->TargetMin());
      }
      sumposprob *= decayrate;
      forestintegral *= decayrate;
      shrinkagefactor *= decayrate;
      fForestGen->SetInitialResponse(decayrate*fForestGen->InitialResponse());
      for (unsigned int itree=0; itree<fForestGen->Trees().size(); ++itree) {
        for (unsigned int inode = 0; inode<fForestGen->Trees()[itree].Responses().size(); ++inode) {
          fForestGen->Trees()[itree].Responses()[inode] *= decayrate;
          fForestGen->Trees()[itree].ResponsesMin()[inode] *= decayrate;
        }
      }
      for (unsigned int iel=0; iel<posprobboundstrees.size(); ++iel) {
        posprobboundstrees[iel] *= decayrate;
      }
      for (unsigned int iel=0; iel<posprobboundsnodes.size(); ++iel) {
        posprobboundsnodes[iel] *= decayrate;
      }
      
      fForestIntegralNow = forestintegral;
    }
    
    if (0) {
//       const double decayrate = 0.9;
//       const double decayrate = pow(integral/forestintegral,0.5*fShrinkage);
//       const double decayrate = integral/forestintegral;
//       double decayrate = pow(1. - lasttreeintegral/forestintegral,0.95);
//       double decayrate = fShrinkage;
//       double newforestintegral = forestintegral - (1.-fShrinkage)*treeintegral;
//       newforestintegral = std::min(newforestintegral,0.9*wavgfull);
//       double decayrate  = std::max(0.5*treeintegral/forestintegral,0.5);
//       double oldforestintegral = forestintegral-treeintegral;
//       double newforestintegral = forestintegral - 0.5*treeintegral;
//       double newforestintegral = forestintegral - 2.*std::max(0.,forestintegral-wavgfull);
//       double decayrate = newforestintegral/forestintegral;
      double decayrate =  0.99;
      double initresponse = fForestGen->InitialResponse();
      double diffresponse = (decayrate-1.)*initresponse;
//       decayrate = 0.5;
      printf("decaying init responses\n");
      for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
        fEvts[iev]->SetTarget(fEvts[iev]->Target()+diffresponse);
//         fEvts[iev]->SetTargetMin(decayrate*fEvts[iev]->TargetMin());
      }
      sumposprob += diffresponse;
      forestintegral += diffresponse;
//       shrinkagefactor *= decayrate;
      fForestGen->SetInitialResponse(decayrate*fForestGen->InitialResponse());
      for (unsigned int iel=0; iel<posprobboundstrees.size(); ++iel) {
        posprobboundstrees[iel] += diffresponse;
      }
      for (unsigned int iel=0; iel<posprobboundsnodes.size(); ++iel) {
        posprobboundsnodes[iel] += diffresponse;
      }
      
      fForestIntegralNow = forestintegral;
    }
//     else {
//       fMinEvents=1000;
//     }
//     if (numtrees>=150) {
//       fMinEvents=1000;
//     }
    
//     bool startenvelope = (numtrees == ntrees/2);
//     if (startenvelope) {
//       fMinEvents = 10e3;
// //       fDoEnvelope = true;
//     }
    
    bool restartgen = (numtrees == ntrees/2);
//     restartgen = true;
    restartgen = false;
    if (restartgen) {
      printf("restart gentree\n");
      sumposprob = std::abs(fForestGen->InitialResponse()*totvolume);
      forestintegral = fForestGen->InitialResponse()*totvolume;
      
      posprobboundstrees.clear();
      posprobboundsnodes.clear();
      nodeidxs.clear();
      
      posprobboundstrees.push_back(sumposprob);
      posprobboundsnodes.push_back(sumposprob);
      nodeidxs.push_back(0);
      
      for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
        fEvts[iev]->SetTarget(fForestGen->InitialResponse());
      }
      
      fForestGen->Trees().clear();
      shrinkagefactor = 0.;
      fShrinkageFactorSecondary = 0.;
    
      
      fForestIntegralNow = forestintegral;
    }
    
    
//     if (numtrees>70 && numtrees<140 && numtrees%2==0) {
    if (0 && numtrees>10 && numtrees<ntrees/2 && numtrees%2==0) {
//     if (numtrees>1 && numtrees%4==0) {
      printf("removing tree\n");
      std::vector<float> evalv(nvars);
      for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
        for (int ivar=0; ivar<nvars; ++ivar) {
          evalv[ivar] = fEvts[iev]->Var(ivar);
        }
        int firstnode = fForestGen->Trees().front().TerminalIndex(evalv.data());
        double firstresponse = fForestGen->Trees().front().GetResponse(firstnode);
        double firstresponsemin = fForestGen->Trees().front().GetResponseMin(firstnode);
        fEvts[iev]->SetTarget(fEvts[iev]->Target() - firstresponse);
        fEvts[iev]->SetTargetMin(fEvts[iev]->TargetMin() - firstresponsemin);
      }
      
      
      double treeposprob = 0.;
      double treeintegral = 0.;
      for (unsigned int inode = 0; inode<fForestGen->Trees().front().Responses().size(); ++inode) {
        double response = fForestGen->Trees().front().GetResponse(inode);
        double volume = 1.0;
        for (int ivar=0; ivar<fNVars; ++ivar) {
          double low = fForestGen->Trees().front().Limits()[inode][ivar].first;
          double high = fForestGen->Trees().front().Limits()[inode][ivar].second;
          volume *= (high-low);
        }
        double prob = response*volume;
        double absprob = std::abs(prob);
//         forestintegral += prob;
        treeintegral += prob;
        if (absprob>0.) {
          treeposprob += absprob;
        }        
  //       printf("iter = %i, inode = %i, response = %5e, volume = %5e, prob = %5f, treeintegral = %5f, forestintegral = %5f\n",iter,inode,response,volume,prob, treeintegral,forestintegral);
      }
      for (unsigned int iel=1; iel<posprobboundstrees.size(); ++iel) {
        posprobboundstrees[iel] -= treeposprob;
      }
      for (unsigned int iel=1; iel<posprobboundsnodes.size(); ++iel) {
        posprobboundsnodes[iel] -= treeposprob;
      }
      posprobboundstrees.erase(posprobboundstrees.begin()+1);
      forestintegral -= treeintegral;
      sumposprob -= treeposprob;
      
      fForestGen->Trees().erase(fForestGen->Trees().begin());
      double firstshrinkagefactor = shrinkagefactors.front();
      shrinkagefactor -= firstshrinkagefactor;
      shrinkagefactors.erase(shrinkagefactors.begin());
      
      fForestIntegralNow = forestintegral;
    }
    
//     if (numtrees>=300) {
//       fMinEvents = 100;
//     }
    
//     if (iter == (niter-1)) break;
    
    if (nev>nevold) {
    
      free(_clss);
      free(_tgtvals);
      free(_tgt2vals);
      free(_tgt3vals);
      free(_tgt4vals);
      free(_tgt5vals);
      free(_fvals);
      free(_weightvals);    
      for (int ivar=0; ivar<nvars; ++ivar) {
        free(_quants[ivar]);
        free(fQuantileMaps[ivar]);
        free(fQuantileMapsMin[ivar]);
      }
      
      _clss  = (int*)memalign(32, nev*sizeof(int));
      _tgtvals  = (double*)memalign(32, nev*sizeof(double));
      _tgt2vals  = (double*)memalign(32, nev*sizeof(double));
      _tgt3vals  = (double*)memalign(32, nev*sizeof(double));
      _tgt4vals  = (double*)memalign(32, nev*sizeof(double));
      _tgt5vals  = (double*)memalign(32, nev*sizeof(double));
      _fvals  = (double*)memalign(32, nev*sizeof(double));
      _weightvals  = (double*)memalign(32, nev*sizeof(double));    
      for (int ivar=0; ivar<nvars; ++ivar) {
        _quants[ivar] = (int*)memalign(32, nev*sizeof(int));
        fQuantileMaps[ivar] = (float*)memalign(32, fNQuantiles*sizeof(float));
        fQuantileMapsMin[ivar] = (float*)memalign(32, fNQuantiles*sizeof(float));
      }
    
    }
    
    BuildQuantiles(nvars, sumabsw);
    
  
  }
  
//   evtslast = fEvts;
  
  TH1::	SetDefaultSumw2();
  
  TH1D *hevts = new TH1D("hevts","",100,-5.,5.);
  TH1D *hevtsw = new TH1D("hevtsw","",100,-5.,5.);
  TH1D *hevtstw = new TH1D("hevtstw","",100,-5.,5.);
  TH1D *htgtpull = new TH1D("htgtpull","",200,-1.,1.);
  TH1D *htgtpullrel = new TH1D("htgtpullrel","",200,-1.,1.);
  TH1D *htgtwidth = new TH1D("htgtwidth","",200,0.,0.1);
  TH2D *htgtwidth2 = new TH2D("htgtwidth2","",100,-0.05,0.1,100,-0.05,0.1);
  TH1D *htgtratiobase = new TH1D("htgtratiobase","",100,0.,2.);
  TH1D *htgtratiointer = new TH1D("htgtratiointer","",100,0.,2.);
  TH1D *htgtratiointerfull = new TH1D("htgtratiointerfull","",100,0.,2.);
  TH1D *htgtdiffinter = new TH1D("htgtdiffinter","",100,-20.,20.);
  TH1D *htgtdiffinterfull = new TH1D("htgtdiffinterfull","",100,-20.,20.);
  TH1D *htgtratio = new TH1D("htgtratio","",100,0.,2.);
  
  TH2D *htgtratio2 = new TH2D("htgtratio2","",101,-0.05,0.05,101,-0.05,0.05);
  htgtratio2->GetXaxis()->SetTitle("h");
  htgtratio2->GetYaxis()->SetTitle("r-h");
  
  double ofuncval;
  double otarget;
  double otargetmin;
  std::vector<float> ox(fNVars);
  
  TFile *fout = new TFile("ntuple.root","RECREATE");
  TTree *evttree = new TTree("evttree","");
  evttree->Branch("funcval",&ofuncval);
  evttree->Branch("target",&otarget);
  evttree->Branch("targetmin",&otargetmin);
  evttree->Branch("x",&ox);
  
  
//   double maxval = -std::numeric_limits<double>::max();
  for (unsigned int iev=0; iev<evtslast.size(); ++iev) {
    double val = evtslast[iev]->Var(0);
    double weight = evtslast[iev]->Weight();
    double funcval = evtslast[iev]->FuncVal();
//     double h = vdt::fast_exp(evtslast[iev]->Target());
//     double h = evtslast[iev]->Target();
    double h = exp(evtslast[iev]->TargetMin());
    double htgt = evtslast[iev]->TargetMin();    
    double envval = evtslast[iev]->Target();  
//     double h = std::min(evtslast[iev]->Target(),vdt::fast_exp(evtslast[iev]->TargetMin()));
//     double h = vdt::fast_exp(evtslast[iev]->TargetMin());
    double s = std::min(h,exp(evtslast[iev]->TargetMin()));
    hevts->Fill(val);
    hevtsw->Fill(val,weight);
    hevtstw->Fill(val,weight*funcval);
//     htgtpull->Fill((evtslast[iev]->FuncVal()-evtslast[iev]->Target()));
    htgtpull->Fill((funcval-h));
    htgtpullrel->Fill( log(funcval) - htgt );
    htgtwidth->Fill(s);
    htgtwidth2->Fill(h,s);
    htgtratiobase->Fill( funcval/h );
//     printf("funcval = %5e, h = %5e\n",funcval,h);
    htgtratio2->Fill(h, funcval-h );
//     htgtratio->Fill( funcval/(h+10.*std::max(s,0.)) );
    htgtratio->Fill( funcval/s );
    htgtratiointer->Fill(h/envval);
    htgtdiffinter->Fill(h-envval);
    htgtratiointerfull->Fill(funcval/envval);
    htgtdiffinterfull->Fill(funcval-envval);    
//     if (funcval>maxval) {
//       maxval = funcval;
//     }
    
    ofuncval = funcval;
    otarget = h;
    otargetmin = evtslast[iev]->TargetMin();
    for (int ivar=0; ivar<fNVars; ++ivar) {
      ox[ivar] = evtslast[iev]->Var(ivar);
    }
    evttree->Fill();
  }
  evttree->Write();
  fout->Close();
  
  hevts->SetLineColor(kBlack);
  hevtsw->SetLineColor(kRed);
  hevtstw->SetLineColor(kBlue);
  hevts->SetMarkerColor(kBlack);
  hevtsw->SetMarkerColor(kRed);
  hevtstw->SetMarkerColor(kBlue);  
  
  gStyle->SetOptStat(111111);
  
  new TCanvas;
  hevts->Draw("E");
  
//   new TCanvas;
//   hevtstw->Draw("E");
//   
//   new TCanvas;
//   hevtsw->Draw("E");
  
//   new TCanvas;
//   htgtpull->Draw("E");

  TCanvas *clogdiff = new TCanvas;
  htgtpullrel->Draw("E");  
  clogdiff->SetLogy();
//   
//   new TCanvas;
//   htgtwidth->Draw("E");  
  
//   new TCanvas;
//   htgtwidth2->Draw("COLZ");  
  
  htgtpullrel->GetXaxis()->SetTitle("ln(funcval) - target");
  htgtpullrel->GetYaxis()->SetTitle("number of events");  
  
  htgtratiobase->GetXaxis()->SetTitle("Primary weight e^{h(#bar{x})}/f(#bar{x})");
  htgtratiobase->GetYaxis()->SetTitle("number of events");
  
  TCanvas *cratiobase = new TCanvas;
  htgtratiobase->Draw("E");
  
  cratiobase->SaveAs("ratiobase.pdf");
  cratiobase->SetLogy();
  cratiobase->SaveAs("ratiobaselog.pdf");
  
  TCanvas *cratiointer = new TCanvas;
  htgtratiointer->Draw("E");    
  htgtratiointer->GetXaxis()->SetTitle("Intermediate weight e^{h(#bar{x})}/g(#bar{x})");
  htgtratiointer->GetYaxis()->SetTitle("number of events");
  cratiointer->SaveAs("ratiointer.pdf");
  cratiointer->SetLogy();
  cratiointer->SaveAs("ratiointerlog.pdf");
  
  new TCanvas;
  htgtdiffinter->Draw("E"); 

  TCanvas *cratiointerfull = new TCanvas;
  htgtratiointerfull->Draw("E");
  htgtratiointerfull->GetXaxis()->SetTitle("Secondary weight f(#bar{x})/g(#bar{x})");
  htgtratiointerfull->GetYaxis()->SetTitle("number of events");
  cratiointerfull->SaveAs("ratiointerfull.pdf");
  cratiointerfull->SetLogy();
  cratiointerfull->SaveAs("ratiointerfulllog.pdf");
    
  
  new TCanvas;
  htgtdiffinterfull->Draw("E");   
  //   
//   new TCanvas;
//   htgtratio->Draw("E");
  

  
//   printf("maxval = %5f\n",maxval);
  
  
  
  
//   fNLLVal = 0.;

//   printf("Initial fNLLVal = %5f\n",fNLLVal);
  

//   //loop over requested number of trees
//   int nunittrees = 0;
//   int nsmalltrees = 0;
//   std::vector<double> nllvals;
//   std::vector<double> dldrvals;
//   for (int itree=0; itree<ntrees; ++itree) {
//     printf("tree %i\n",itree);
//     
//     fForest->Trees().push_back(MCGBRTreeD());
//     MCGBRTreeD &tree = fForest->Trees().back();
//     TrainTree(fEvts,sumw,tree,0,limits);
//       
//     int treesize = tree.Responses().size();        
//     
//     printf("treesize = %i\n",treesize);
// 
//     if (treesize==1) {
//       ++nunittrees;
//     }
//     else {
//       nunittrees = 0.;
//     }
//     
//     if (treesize>1 && treesize<16) {
//       ++nsmalltrees;
//     }
//     else {
//       nsmalltrees = 0;
//     }
//     
//     
//   }
   
  
  
}
 
      
//_______________________________________________________________________
void MCGBRIntegrator::TrainTree(const std::vector<MCGBREvent*> &evts, double sumwtotal, MCGBRTreeD &tree, int depth, const std::vector<std::pair<float,float> > limits, bool usetarget, bool doenv) {
    
  
//   int minevents = doenv ? 10 : fMinEvents;
//   int minevents = doenv ? fMinEvents : std::max(10,int(fEvts.size()/100));
//   int minevents = fForestGen->Trees().size() > 10 ? std::max(fMinEvents,int(fEvts.size()/1e3)) : fMinEvents;
  int minevents = fMinEvents;
  double mincutsig = fMinCutSignificance;
//   int maxnodes = fForestGen->Trees().size() > 10 ? 2000 : fMaxNodes;
  int maxnodes = fMaxNodes;
//   double mincutsig = doenv ? 1e4*fMinCutSignificance : fMinCutSignificance;
//   double mincutsig = doenv ? 1e-6 : fMinCutSignificance;


  
//   printf("TrainTree\n");
  
//   const double a = 0.;
  
  int nvars = fNVars;

  int thisidx = tree.CutIndices().size();    
  
  //number of events input to node
  const int nev = evts.size();
  
  //index of best cut variable
  int bestvar = 0;
  
  double volume = 1.;
  for (std::vector<std::pair<float,float> >::const_iterator it = limits.begin(); it!=limits.end(); ++it) {
    double low = it->first;
    double high = it->second;
    volume *= (high-low);
  }
  
//   const double sigma1 = 10.0;
//   const double sigma2 = 0.01;
//   
//   const double k1 = 1./(sigma1*sqrt(2.));
//   const double k2 = 1./(sigma2*sqrt(2.));  
  
//   const double k = 1e4;
//   const double logsigma = 0.10;
//   const double k = 0.5/logsigma/logsigma;
   
//   const double kenv = 0.1;
  
//   constexpr double sigma = 0.01;
//   constexpr double k = 0.5/sigma/sigma;
//   constexpr double kenv = 100.;
  
//   double fulltgtmax = std::numeric_limits<double>::lowest();
//   for (int iev = 0; iev<nev; ++iev) {
//     double tgtval = evts[iev]->FuncVal() - evts[iev]->Target();
//     if (tgtval > fulltgtmax) {
//       fulltgtmax = tgtval;
//     }
//   }
//    const double k = 200;
//   double tgtmax = std::numeric_limits<double>::lowest();


//   double sigmaratiomax = std::numeric_limits<double>::lowest();

  if (doenv) {
// #pragma omp parallel for
    for (int iev = 0; iev<nev; ++iev) {
      const double weight = 1.;
      const double arg = evts[iev]->Arg();
      const double tgtval = -arg;
      const double tgt2val = 1.;
      
      _weightvals[iev] = weight;
      _tgtvals[iev] = tgtval;
      _tgt2vals[iev] = tgt2val;
      
      for (int ivar=0; ivar<nvars; ++ivar) {
        _quants[ivar][iev] = evts[iev]->Quantile(ivar);
      }
    }
  }
  else {
// #pragma omp parallel for
    for (int iev = 0; iev<nev; ++iev) {
      const double weight = 1.;
      const double arg = evts[iev]->ArgLog();
      const double tgtval = -arg;
      const double tgt2val = 1.;
      
      _weightvals[iev] = weight;
      _tgtvals[iev] = tgtval;
      _tgt2vals[iev] = tgt2val;
      
      for (int ivar=0; ivar<nvars; ++ivar) {
        _quants[ivar][iev] = evts[iev]->Quantile(ivar);
      }      
    }
  }

  //printf("first parallel loop\n");
  //printf("nev = %i\n",nev);
  //fill temporary array of quantiles (to allow auto-vectorization of later loops)
//   #pragma omp parallel for
//   for (int iev = 0; iev<nev; ++iev) {
// //     _tgtvals[iev] = evts[iev]->Weight()*evts[iev]->FuncVal();
//     double weight = evts[iev]->Weight();
//     double funcval = evts[iev]->FuncVal();
//     double target = evts[iev]->Target();
//     double targetmin = evts[iev]->TargetMin(); 
// //     double funcvalalt = evts[iev]->FuncValAlt();
// //     double tgtval = usetarget ? evts[iev]->FuncVal() - evts[iev]->Target() : evts[iev]->FuncVal();
// //     double tgtval = evts[iev]->FuncVal();
// //     double tgtval = evts[iev]->FuncVal() - evts[iev]->Target();
//     
//     
// //     double tgtval = doenv ? exp(targetmin) - target : log(funcval) - targetmin;
//     
// //     double tgtval = doenv ? NLSkewGausDMu(exp(targetmin),target,doenv) : log(funcval) - targetmin;
// //     double tgt2val = doenv ? NLSkewGausD2Mu(exp(targetmin),target,doenv ) : 0.;
//    
// //     double tgtval = doenv ? NLSkewGausDMu(funcval,target,doenv) : NLSkewGausDMu(log(funcval),targetmin,doenv);
// //     double tgt2val = doenv ? NLSkewGausD2Mu(funcval,target,doenv ) : NLSkewGausD2Mu(log(funcval),targetmin,doenv);
//     
//     
// //     double tgtval = doenv ? NLSkewGausDMu(exp(targetmin),target,funcval,doenv) : NLSkewGausDMu(log(funcval),targetmin,funcval,doenv);
// //     double tgt2val = doenv ? NLSkewGausD2Mu(exp(targetmin),target,funcval,doenv ) : NLSkewGausD2Mu(log(funcval),targetmin,funcval,doenv);
//     
// //     double tgtval = doenv ? NLLogNormDMu(exp(targetmin),target,exp(targetmin)-funcval) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLLogNormD2Mu(exp(targetmin),target,exp(targetmin)-funcval) : NLNormD2Mu(log(funcval),targetmin);
//     
// //     double tgtval = doenv ? NLLogNormDMu(exp(targetmin),target,funcvalalt) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLLogNormD2Mu(exp(targetmin),target,funcvalalt) : NLNormD2Mu(log(funcval),targetmin);
//     
// //     double tgtval = doenv ? targetmin : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? 0. : NLNormD2Mu(log(funcval),targetmin);
//     
// //     double tgtval = doenv ?  NLNormDMu(targetmin,0.) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLNormD2Mu(targetmin,0.) : NLNormD2Mu(log(funcval),targetmin);
//     
// //     double tgtval = doenv ?  NLLogNormDMu(exp(targetmin),target) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLLogNormD2Mu(exp(targetmin),target) : NLNormD2Mu(log(funcval),targetmin);
//     
// //     double tgtval = doenv ?  NLNormDMu(targetmin,0.) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLNormD2Mu(targetmin,0.) : NLNormD2Mu(log(funcval),targetmin);
//     
// //     double tgtval = doenv ?  NLNormDMu(log(std::max(1e-6,exp(targetmin)-target)),0.) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLNormD2Mu(log(std::max(1e-6,exp(targetmin)-target)),0.) : NLNormD2Mu(log(funcval),targetmin);
//     
// //     double sigma = doenv ? 0.01 : 0.01;
// //     double tgtval = doenv ?  NLNormDMu(log(std::max(exp(-24.),exp(targetmin)-target)),0.,sigma) : NLNormDMu(log(funcval),targetmin,sigma);
// //     double tgt2val = doenv ? NLNormD2Mu(log(std::max(exp(-24.),exp(targetmin)-target)),0.,sigma) : NLNormD2Mu(log(funcval),targetmin,sigma);
// //      double tgtval = doenv ? log(std::max(exp(-24.),exp(targetmin)-target)) : log(funcval) - targetmin;
// 
// //     double tgtval = doenv ? log(std::max(exp(-24.),exp(targetmin)-target)) : (fDoEnvelope ? NLSkewGausDMu(log(funcval),targetmin) : log(funcval) - targetmin);
// //     double tgt2val = doenv ? 0. : (fDoEnvelope ? NLSkewGausD2Mu(log(funcval),targetmin) : 0.);
//     
//     double px = doenv ? log(std::max(exp(-128.),exp(targetmin)-target)) : log(funcval);
//     double pmu = doenv ? 0. : targetmin;
// //     double sigmabase = 0.05;
// //     double psigma = doenv ? 1. : 0.3;
// //     double psigma = sigmabase;
//     double psigma = 0.05;
//     double tgtval = fDoEnvelope && !doenv ? NLSkewGausDMu(px,pmu) : NLNormDMu(px,pmu,psigma);
//     double tgt2val = fDoEnvelope && !doenv ? NLSkewGausD2Mu(px,pmu) : NLNormD2Mu(px,pmu,psigma);
//     
// //     double px = doenv ? log(std::max(exp(-128.),funcval-target)) : log(funcval);
// //     double pmu = doenv ? 0. : targetmin;
// //     double sigmabase = 0.05;
// //     double psigma = sigmabase;
// //     double tgtval = fDoEnvelope && !doenv ? NLSkewGausDMu(px,pmu) : NLNormDMu(px,pmu,psigma);
// //     double tgt2val = fDoEnvelope && !doenv ? NLSkewGausD2Mu(px,pmu) : NLNormD2Mu(px,pmu,psigma);
//     
// //     double sigmaratio = doenv ? (exp(targetmin)-target)/exp(targetmin) : 0.;
// //     
// //     if (sigmaratio>sigmaratiomax) {
// //       sigmaratiomax = sigmaratio;
// //     }
//     
// //     double tgtval = doenv ? exp(targetmin) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? exp(2.*targetmin) : NLNormD2Mu(log(funcval),targetmin);
//     
// //     double tgtval = doenv ? NLNormDMu(targetmin,target) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLNormD2Mu(targetmin,target) : NLNormD2Mu(log(funcval),targetmin);
//     
// //     double tgtval = doenv ? NLNormDMu(exp(targetmin),target) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLNormD2Mu(exp(targetmin),target) : NLNormD2Mu(log(funcval),targetmin);
//     
// //     double tgtval = doenv ? NLSkewGausDMu(funcval,target,exp(targetmin),doenv) : NLSkewGausDMu(log(funcval),targetmin,funcval,doenv);
// //     double tgt2val = doenv ? NLSkewGausD2Mu(funcval,target,exp(targetmin),doenv ) : NLSkewGausD2Mu(log(funcval),targetmin,funcval,doenv);
//     
// /*    double tgtval = doenv ? NLSkewGausDMu(funcval,target) : log(funcval) - targetmin;
//     double tgt2val = doenv ? NLSkewGausD2Mu(funcval,target) : 0.;  */  
//     
//     
// //     double tgtval = pow(funcval,0.8) - target;
// 
// //     double tgtval = pow(funcval/510.,0.9)*510. - target;
//     
// //        double tgt2val = 1./(funcval*funcval);
// //        double tgt2val = std::min(1e3,1./(funcval*funcval));
// //        double tgt2val = funcval;
// //        double tgt2val = 1.;
// //        double tgtval = (funcval-target)*tgt2val;
//     
// //     double x = log(funcval);
// //     double h = targetmin;
// //     double erfexpr = 1./Faddeeva::erfcx(k2*(x-h));
//     
// //     double tgtval = -2.*k1*k1*(x-h) - 2.*k2*erfexpr/sqrt(M_PI);
// //     double tgt2val = 2.*k1*k1 - 4.*k2*k2*k2*(x-h)*erfexpr/sqrt(M_PI) + 4.*k2*k2*erfexpr*erfexpr/M_PI;    
//     
// //     double expF = vdt::fast_exp(targetmin);
// //     double tgtval = target;
// //     double tgt2val = expF;
// //     double tgt3val = funcval/expF;
// //     double tgt4val = log(funcval) - targetmin;
//     
// //     _weightvals[iev] = weight;
// //     _tgtvals[iev] = weight*tgtval;
// //     _tgt2vals[iev] = weight*tgt2val;  
// //     _tgt3vals[iev] = weight*tgt3val;
// //     _tgt4vals[iev] = weight*tgt4val;  
//     
//     
// //     double tgtval = funcval - target;
// //     double tgt2val = funcval*vdt::fast_exp(-targetmin);
// //     double tgt3val = vdt::fast_exp(-targetmin);
// //     double tgt4val = (funcval-target)*vdt::fast_exp(-targetmin);
// //     double tgt5val = log(funcval) - targetmin;
// 
// //     double tgtval = funcval - target;
// //     double fulltgt2 = (target + fulltgtmax)*(target + fulltgtmax);
// //     double tgtval = funcval/fulltgt2;
// //     double tgt2val = funcval - target;
//     
// //     double h = target;
// //     double f = target;
// //     double x = funcval;
//     
// //     double tgtval = 1./h - (exp(k*(-1. + x/h))*k*x)/((1. + exp(k*(-1. + x/h)))*h*h);
// //     double tgt2val = 1./h/h + (2.*exp(k*(-1. + x/h))*k*x)/((1. + exp(k*(-1. + x/h)))*h*h*h) - (exp(2.*k*(-1. + x/h))*k*k*x*x)/(pow((1. + exp(k*(-1. + x/h))),2)*h*h*h*h) + (exp(k*(-1. + x/h))*k*k*x*x)/((1. + exp(k*(-1. + x/h)))*h*h*h*h);
//     
// /*    double tgtval = 1./h - (k*x)/((1. + exp(-k*(-1. + x/h)))*h*h);
//     double tgt2val = 1./h/h + (2.*k*x)/((1. + exp(-k*(-1. + x/h)))*h*h*h) - (k*k*x*x)/(pow((1. + exp(-k*(-1. + x/h))),2)*h*h*h*h) + (k*k*x*x)/((1. + exp(-k*(-1. + x/h)))*h*h*h*h); */      
//     
// //     double expexpr = 1./(1.+exp(-k*(-1. + x/h)));
// //     double tgtval = 1./h - (k*x)*expexpr/(h*h);
// //     double tgt2val = 1./h/h + (2.*k*x)*expexpr/(h*h*h) - (k*k*x*x)*expexpr*expexpr/(h*h*h*h) + (k*k*x*x)*expexpr/(h*h*h*h);     
//    
// //     double arg = k*(x/h-1.);
// //     double expnarg2 = exp(-arg*arg);
// //     double den = 1./(sqrt(M_PI)*(1.-erf(arg)));
// //     double hinv = 1./h;
// //     
// //     double tgtval = hinv - 2.*k*x*expnarg2*den*hinv*hinv;
// //     double tgt2val = 4*k*k*x*x*expnarg2*den*den*hinv*hinv*hinv*hinv - 4*k*k*k*x*x*(x*hinv-1.)*expnarg2*den*hinv*hinv*hinv*hinv + 4.*k*x*expnarg2*den*hinv*hinv*hinv - hinv*hinv;
//    
// //     printf("tgtval = %5f, tgt2val = %5f\n",tgtval,tgt2val);
//     
// /*    double expexpr = 1./(1. + exp(-k*(f*x-1)));
//     double finv = 1./f;
//     double tgtval = -finv + k*x*expexpr;
//     double tgt2val = finv*finv + k*k*x*x*expexpr - k*k*x*x*expexpr*expexpr; */   
//     
// 
// //     double tgtval = 1./target;
// //     double tgt2val = tgtval*tgtval;
// //     double tgt3val = tgt2val*tgtval;
// //     double tgt4val = tgt3val*tgtval;
// //     double tgt5val = funcval - target;
// // //     double tgt5val = -funcval + target - targetmin;
// // 
// //     _weightvals[iev] = weight;
// //     _tgtvals[iev] = weight*tgtval;
// //     _tgt2vals[iev] = weight*tgt2val;
// //     _tgt3vals[iev] = weight*tgt3val;
// //     _tgt4vals[iev] = weight*tgt4val;
// //     _tgt5vals[iev] = weight*tgt5val;
// 
// //     double rinv2 = 1./(funcval*funcval);
// //     double tgtval = funcval - target;
// //     double tgtval = funcval > 1e-12 ? log(funcval) - log(1.-exp(-1e6*funcval)) - target : -log(1e6) + 0.5*1e6*funcval - target;
// 
// 
// //     double tgtval = 1./target;
// //     double tgt2val = tgtval*tgtval;
// //     double tgt3val = tgtval*tgt2val;
// //     double tgt4val = tgtval*tgt3val;
// //     double tgtval = funcval - target;
// //     double tgtval = log(target + fulltgtmax);
// //     double tgt2val = target;
// //     double tgt5val = funcval - target;
// //     double fval = log(funcval) - targetmin;
//     
// //     double tgt5val = (1.+k)*funcval - target - k*targetmin;
// //     double fval = funcval - targetmin;   
// 
// //     double tgt5val = log(funcval) - target;
// //     double tgt5val = funcval - target;
// 
// //     double tgt5val = funcval - target + k*targetmin;
// //     double fval = funcval - target + targetmin + k*targetmin;
//     
// //     double tgtval = 1./(funcval*funcval);
// //     double tgt2val = (funcval - target)*tgtval;
// //     double tgt3val = tgt2val*tgtval;
// //     double tgt4val = tgt3val*tgtval;
// //     double tgt5val = funcval - target;
// //     double tgt5val = -funcval + target - targetmin;
// 
// //     double r = funcval;
// //     double h = target;
// //     double invh = 1./h;
// // 
// //     double tgtval = 2.*k*(log(h)-log(r))*invh;
// //     double tgt2val = 2.*k*(1. + log(r) - log(h))*invh*invh;
// 
// //     double tgtval = funcval - target;
// 
// /*    double exp2ns = vdt::fast_exp(-2.*targetmin);
//     double diff = funcval - target;
//     double tgtval = exp2ns;
//     double tgt2val = exp2ns*diff;
//     double tgt3val = tgt2val*diff;  */  
// 
// //     double tgtval = log(funcval) - targetmin;
// //     double tgt2val = funcval*exp(-targetmin);
//     
//     _weightvals[iev] = weight;
//     _tgtvals[iev] = weight*tgtval;
//     _tgt2vals[iev] = weight*tgt2val;
// //     _tgt3vals[iev] = weight*tgt3val;
// //     _tgt4vals[iev] = weight*tgt4val;
// //     _tgt5vals[iev] = weight*tgt5val;    
// //     _fvals[iev] = fval;
//     
//     
//     
// //     double tgtval = log(funcval) - target;
// //     double tgt2val = funcval*vdt::fast_exp(-target);
//     
// //     double tgtval = 1./targetmin;
// //     double tgt2val = tgtval*tgtval;
// // //     double tgt3val = tgtval*tgt2val;
// //     double tgt4val = tgtval*tgt3val;
//     
// //     double tgt3val = 1./target;
// //     double tgt4val = tgt3val*tgt3val;
//     
// //     double tgt5val = funcval - target;
// //     double fval = funcval - targetmin;
// //     double fval = funcval - target + targetmin;
//     
// //     double expF = vdt::fast_exp(targetmin);
//     
// //     double tgtval = funcval/expF;
// //     double tgt2val = expF;
// //     double tgt3val = target;
// //     double tgt4val = funcval - target;
//     //     double tgt2val = evts[iev]->TargetMin();
// //     double tgt2val = evts[iev]->FuncVal() - evts[iev]->TargetMin();
//     
//     
// //     double expS = vdt::fast_exp(targetmin);
// //     double expnF = vdt::fast_exp(-target);
// //     
// //     double tgtval = funcval*expS*expnF;
// //     double tgt2val = expS;
// //     double tgt3val = log(funcval) - target;
// //     
// //     _weightvals[iev] = weight;
// //     _tgtvals[iev] = weight*tgtval;
// //     _tgt2vals[iev] = weight*tgt2val;  
// //     _tgt3vals[iev] = weight*tgt3val;
// //     _tgt4vals[iev] = weight*tgt4val;
//     
// //     _tgt5vals[iev] = weight*tgt5val;
// //     _fvals[iev] = weight*fval;
//     
// //     _tgt2vals[iev] = weight*tgtval*tgtval;    
//     
// //     if (tgtval>tgtmax) tgtmax = tgtval;
//     
// //     printf("tgtval = %5f, funcval = %5f, tgt = %5f\n",tgtval,evts[iev]->FuncVal(),evts[iev]->Target());
//     
// //     double funcval = evts[iev]->FuncVal();
// //     double target =  evts[iev]->Target();
// //     double targetmin =  evts[iev]->TargetMin();
//     
// //     double oos = vdt::fast_exp(targetmin);
// //     double oos2 = oos*oos;
// //     double expr = target*oos2-a*oos-funcval*oos2;
// //     double expr2 = 2.*target*oos2-a*oos-2*funcval*oos2;
//     
// //     _weightvals[iev] = weight;
// //     _tgtvals[iev] = weight*expr;
// //     _tgt2vals[iev] = weight*(-1. + (target-funcval)*expr);
// //     _tgt3vals[iev] = weight*oos2;
// //     _tgt4vals[iev] = weight*(target-funcval)*expr2;
// //     _tgt5vals[iev] = weight*expr2;
//     
//     
// //     double r = evts[iev]->FuncVal();
// //     double h =  evts[iev]->Target();
// //     double s =  evts[iev]->TargetMin();    
//     
// //     double rmh = r-h;
// //     double exps = vdt::fast_exp(s);
// //     double exp2s = vdt::fast_exp(2.*s);
//     
//     
//     
// //     _weightvals[iev] = weight;
// //     _tgtvals[iev] = weight*exp2s*rmh;
// //     _tgt2vals[iev] = weight*exp2s;
// //     _tgt3vals[iev] = weight*exp2s*rmh*rmh;
// //     _tgt4vals[iev] = weight*(target-funcval)*expr2;
// //     _tgt5vals[iev] = weight*expr2;    
//     
//     
//     
// //     double r = evts[iev]->FuncVal();
// //     double h =  evts[iev]->Target();
// //     double s =  evts[iev]->TargetMin();    
// //     
// //     double rmh = r-h;
// //     double exps = vdt::fast_exp(s);
// //     double exp2s = exps*exps;
// //     
// //     
// //     
// //     _weightvals[iev] = weight;
// //     _tgtvals[iev] = weight*rmh;
// //     _tgt2vals[iev] = weight*rmh*rmh;
// //     _tgt3vals[iev] = weight*exp2s*rmh*rmh;
//     
//     
//     
//     
//     
// //     _tgt2vals[iev] = weight*tgt2val;
// //     _tgt2vals[iev] = weight*tgtmaxval;
// //     _tgt2vals[iev] = weight*tgtval*tgtval;
// //     _fvals[iev] = evts[iev]->Target();
//     for (int ivar=0; ivar<nvars; ++ivar) {
//       _quants[ivar][iev] = evts[iev]->Quantile(ivar);   
//     }
//   }  
    
//   double sigmascale = doenv ? std::max(1.,1./sigmaratiomax) : 1.;
//   double sigmascale = doenv ? (1./(1-fShrinkageFactorSecondary)) : 1.;
//   double valscale = 0 ? 1./(sigmascale*sigmascale) : 1.;
  
  constexpr double valsigma = 0.01;
  constexpr double valscale = 1./valsigma/valsigma;
//   constexpr double valscale = 1.;
//   if (doenv) printf("sigmascale = %5f\n",sigmascale);
    
  //printf("second parallel loop\n");
//   //trivial open-mp based multithreading of loop over input variables
  //The loop is thread safe since each iteration writes into its own
  //elements of the 2-d arrays
//   #pragma omp parallel for
  for (int ivar=0; ivar<nvars; ++ivar) {
             
    
    int minquant;
    int maxquant;
    
//     double lowbound = limits[ivar].first;
//     double highbound = limits[ivar].second;
//     double volumenorm = volume/(highbound-lowbound);
    
//     printf("lowbound = %5f, highbound = %5f, volumenorm = %5f, volume = %5f, nev = %i\n", lowbound,highbound,volumenorm,volume, int(evts.size()));
    
    //find max and min quantiles in the input events
    //(this loop should be vectorized by gcc with reasonable optimization options)
    GBRArrayUtils::MinMaxQuants(minquant, maxquant, _quants[ivar], nev);

    //calculate offset and scaling (powers of 2) to reduce the total number of quantiles
    //to the fNBinsMax for the search for the best split value
    int offset = minquant;
    unsigned int bincount = maxquant-minquant+1;
    unsigned int pscale = 0;
    while (bincount>fNBinsMax) {
      ++pscale;
      //bincount >>= 1;
      bincount = ((maxquant-offset)>>pscale) + 1;
    }    

    const unsigned int nbins = ((maxquant-offset)>>pscale) + 1;
    assert(nbins<=fNBinsMax);
    
//     printf("nbins = %i\n",int(nbins));
        
    //zero arrays where necessary and 
    GBRArrayUtils::InitArrays(_ns[ivar],_tgts[ivar],_tgt2s[ivar],_bsepgains[ivar],nbins);
    GBRArrayUtils::ZeroArray(_ws[ivar],nbins);
    GBRArrayUtils::ZeroArray(_tgt3s[ivar],nbins);
    GBRArrayUtils::ZeroArray(_tgt4s[ivar],nbins);
    GBRArrayUtils::ZeroArray(_tgt5s[ivar],nbins);
    
    GBRArrayUtils::MinArray(_tgtmaxs[ivar],nbins);
    GBRArrayUtils::MaxArray(_fmins[ivar],nbins);
    GBRArrayUtils::MinArray(_varvalmaxs[ivar],nbins);
    GBRArrayUtils::MaxArray(_varvalmins[ivar],nbins);
        
    //compute map between bin numbers
    //and variable cut values
    
    //printf("quant manipulation\n");
    //this loop should auto-vectorize
    GBRArrayUtils::FillBinQuants(_binquants[ivar], offset, pscale,fNQuantiles, nbins);
    
    
    //printf("done quant manipulation\n");
    
    //printf("compute bins\n");
   

      
    //printf("filling histogram-style arrays\n");
     
    //compute summed quantities differential in each bin
    //(filling 'histograms')
    //This loop is one of the most expensive in the algorithm for large training samples
    //This loop can unfortunately not be vectorized because the memory addressed 
    //are computed within the loop iteration
    //JOSH: Once the data-dependency checks are appropriately bypassed, this should actually vectorize, but only
    //for avx2 targets which support vectorized gather instructions
//    int nevd = 0;
    for (int iev=0;iev<nev;++iev) {
      int quant = _quants[ivar][iev];
      int ibin = (quant-offset)>>pscale;

//       printf("quant = %i, bin = %i, var %i = %5f\n", quant,ibin,ivar,evts[iev]->Var(ivar));
      
      ++_ns[ivar][ibin];
      
      _ws[ivar][ibin] += _weightvals[iev];                  
      _tgts[ivar][ibin] += _tgtvals[iev];
      _tgt2s[ivar][ibin] += _tgt2vals[iev];
      _tgt3s[ivar][ibin] += _tgt3vals[iev];
      _tgt4s[ivar][ibin] += _tgt4vals[iev];
      _tgt5s[ivar][ibin] += _tgt5vals[iev];
      
//       _fmins[ivar][ibin] = std::min(_fmins[ivar][ibin], _fvals[iev]);
      _tgtmaxs[ivar][ibin] = std::max(_tgtmaxs[ivar][ibin], _tgtvals[iev]);
      _fmins[ivar][ibin] = std::min(_fmins[ivar][ibin], _tgtvals[iev]);
      _varvalmaxs[ivar][ibin] = std::max(_varvalmaxs[ivar][ibin],_quants[ivar][iev]);
      _varvalmins[ivar][ibin] = std::min(_varvalmins[ivar][ibin],_quants[ivar][iev]);
    }

    
//     _varvalmaxs[ivar][0] = _binquants[ivar][0];
//     _varvalmins[ivar][0] = 0;
//     for (unsigned int ibin=1; ibin<nbins; ++ibin) {
//       _varvalmaxs[ivar][ibin] = _binquants[ivar][ibin];
//       _varvalmins[ivar][ibin] = _binquants[ivar][ibin-1]+1;
//     }


    
      
    //printf("starting split search\n");
 
    //printf("compute array integrals\n");
    //convert differential arrays to cumulative arrays by summing over
    //each element
    //loop cannot be vectorized because this is an iterative calculation
    _sumws[ivar][0] = _ws[ivar][0];
    _sumns[ivar][0] = _ns[ivar][0];
    _sumtgts[ivar][0] = _tgts[ivar][0];
    _sumtgt2s[ivar][0] = _tgt2s[ivar][0];
    _sumtgt3s[ivar][0] = _tgt3s[ivar][0];
    _sumtgt4s[ivar][0] = _tgt4s[ivar][0];
    _sumtgt5s[ivar][0] = _tgt5s[ivar][0];
    
//     _sumtgt2s[ivar][0] = _tgts[ivar][0]!=_tgt2s[ivar][0] ? 1./(_tgts[ivar][0]-_tgt2s[ivar][0]) : 1.;
//     _sumfmins[ivar][0] = _fmins[ivar][0];
//     _sumfminsr[ivar][nbins-1] = _fmins[ivar][nbins-1];
    _sumtgtmaxs[ivar][0] = _tgtmaxs[ivar][0];
    _sumtgtmaxsr[ivar][nbins-1] = _tgtmaxs[ivar][nbins-1];
    _sumfmins[ivar][0] = _fmins[ivar][0];
    _sumfminsr[ivar][nbins-1] = _fmins[ivar][nbins-1];    
    _sumvarvalmaxs[ivar][0] = _varvalmaxs[ivar][0];
    _sumvarvalmins[ivar][nbins-1] = _varvalmins[ivar][nbins-1];  
    
    for (unsigned int ibin=1; ibin<nbins; ++ibin) {      
      _sumns[ivar][ibin] = _sumns[ivar][ibin-1] + _ns[ivar][ibin];
      _sumws[ivar][ibin] = _sumws[ivar][ibin-1] + _ws[ivar][ibin];  
      _sumtgts[ivar][ibin] = _sumtgts[ivar][ibin-1] + _tgts[ivar][ibin];  
      _sumtgt2s[ivar][ibin] = _sumtgt2s[ivar][ibin-1] + _tgt2s[ivar][ibin];  
      _sumtgt3s[ivar][ibin] = _sumtgt3s[ivar][ibin-1] + _tgt3s[ivar][ibin];  
      _sumtgt4s[ivar][ibin] = _sumtgt4s[ivar][ibin-1] + _tgt4s[ivar][ibin];  
      _sumtgt5s[ivar][ibin] = _sumtgt5s[ivar][ibin-1] + _tgt5s[ivar][ibin];  
      
    }
//     for (unsigned int ibin=1; ibin<nbins; ++ibin) {      
//       if (_tgts[ivar][ibin]!=_tgt2s[ivar][ibin]) {
//         _sumtgt2s[ivar][ibin] = _sumtgt2s[ivar][ibin-1] + 1./(_tgts[ivar][ibin]-_tgt2s[ivar][ibin]);  
//       }
//       else {
//         _sumtgt2s[ivar][ibin] = _sumtgt2s[ivar][ibin-1] + 1.;  
//       }
//     }

    
    for (unsigned int ibin=1; ibin<nbins; ++ibin) {
//       _sumfmins[ivar][ibin] = std::min(_sumfmins[ivar][ibin-1], _fmins[ivar][ibin]);  
      _sumtgtmaxs[ivar][ibin] = std::max(_sumtgtmaxs[ivar][ibin-1], _tgtmaxs[ivar][ibin]);  
      _sumvarvalmaxs[ivar][ibin] = std::max(_sumvarvalmaxs[ivar][ibin-1], _varvalmaxs[ivar][ibin]);
      _sumfmins[ivar][ibin] = std::min(_sumfmins[ivar][ibin-1], _fmins[ivar][ibin]);  
    }
    if (nbins>1) {
//        _sumtgtmaxsr[ivar][nbins-2] = _tgtmaxs[ivar][nbins-1];
//        _sumvarvalmins[ivar][nbins-2] = _varvalmins[ivar][nbins-1];
      for (unsigned int ibin=nbins-1; ibin>0; --ibin) {
//         _sumtgtmaxsr[ivar][ibin-1] = std::max(_sumtgtmaxsr[ivar][ibin+1],_tgtmaxs[ivar][ibin]);
//         _sumvarvalmins[ivar][ibin-1] = std::min(_sumvarvalmins[ivar][ibin+1],_varvalmins[ivar][ibin]);
//         _sumfminsr[ivar][ibin-1] = std::min(_sumfminsr[ivar][ibin+1],_fmins[ivar][ibin]);
        _sumtgtmaxsr[ivar][ibin-1] = std::max(_sumtgtmaxsr[ivar][ibin],_tgtmaxs[ivar][ibin]);
        _sumvarvalmins[ivar][ibin-1] = std::min(_sumvarvalmins[ivar][ibin],_varvalmins[ivar][ibin]);        
        _sumfminsr[ivar][ibin-1] = std::min(_sumfminsr[ivar][ibin],_fmins[ivar][ibin]);
      }
    }

    
//     for (unsigned int ibin=0; ibin<nbins; ++ibin) {      
//       printf("ibin = %i, binquant = %i, varvalmaxs = %i, varvalmins = %i, sumvarvalmaxs = %i, sumvarvalmins = %i\n",ibin, _binquants[ivar][ibin],_varvalmaxs[ivar][ibin],_varvalmins[ivar][ibin],_sumvarvalmaxs[ivar][ibin],_sumvarvalmins[ivar][ibin]);
//     }

    
    //make sure "dead space" is attached to the larger target value for upper envelope
/*    for (unsigned int ibin=0; ibin<nbins; ++ibin) {
      if (_sumtgtmaxsr[ivar][ibin]>_sumtgtmaxs[ivar][ibin]) {
//       if (fQuantileMaps[ivar][_sumvarvalmaxs[ivar][ibin]]<0.) {
      
        int quant = _sumvarvalmaxs[ivar][ibin];
        _varvals[ivar][ibin] = fQuantileMaps[ivar][quant];
      }
      else {
        int quant = _sumvarvalmins[ivar][ibin];
        float val = fQuantileMapsMin[ivar][quant];
        //decrement by epsilon
        _varvals[ivar][ibin] = std::nextafter(val,std::numeric_limits<float>::lowest());
      }
    } */  

    for (unsigned int ibin=0; ibin<nbins; ++ibin) {      
        int quant = _sumvarvalmaxs[ivar][ibin];
        _varvals[ivar][ibin]  = fQuantileMaps[ivar][quant];
    }
    
    
/*    for (unsigned int ibin=0; ibin<nbins; ++ibin) {
      int quant = _sumvarvalmaxs[ivar][ibin];
      float val = fQuantileMaps[ivar][quant];
    
      int quantalt = _sumvarvalmins[ivar][ibin];
      float valalt = fQuantileMapsMin[ivar][quantalt];

      _varvals[ivar][ibin] = 0.5*(val+valalt);
    } */    
    
//     _varvals[ivar][nbins-1] = fQuantileMaps[ivar][_binquants[ivar][nbins-1]];    
    
//     if (doenv) {
//       for (unsigned int ibin=0; ibin<nbins; ++ibin) {      
//         _sumtgt2s[ivar][ibin] = (_varvals[ivar][ibin] - lowbound)*volumenorm;
//       }
//     }
    
    const double sumn = _sumns[ivar][nbins-1];
//     const double sumw = _sumws[ivar][nbins-1];
    const double sumtgt = _sumtgts[ivar][nbins-1];
    const double sumtgt2 = _sumtgt2s[ivar][nbins-1];
//     const double sumtgt3 = _sumtgt3s[ivar][nbins-1];
//     const double sumtgt4 = _sumtgt4s[ivar][nbins-1];
//     const double sumtgt5 = _sumtgt5s[ivar][nbins-1];
    
//     const double sumtgt2 = volume;
//     const double sumfmin = _sumfmins[ivar][nbins-1];
//     const double sumtgtmax = _sumtgtmaxs[ivar][nbins-1];
//     const double sumtgtmin = _sumfmins[ivar][nbins-1];
    
//     printf("sumtgt2 = %5f, volume = %5f\n",sumtgt2,volume);
    
    //printf("short loop\n");
    double maxsepgain = 0.;
    float cutval = 0.;
    int nleft= 0;
    int nright = 0;
    int bestbin=0;
    
    
//     Eigen::Vector2d df(sumtgt,sumtgt2);
//     Eigen::Matrix2d d2f;
//     d2f(0,0) = sumtgt3;
//     d2f(0,1) = sumtgt5;
//     d2f(1,0) = sumtgt5;
//     d2f(1,1) = sumtgt4;
//     
//     Eigen::Vector2d dx = -d2f.ldlt().solve(df);
// //     const double curval = -sumtgt*log(sumtgt/(sumtgt+sumw)) - sumw*log(sumw/(sumtgt+sumw));
// //     const double curval = -sumtgt*sumtgt/sumw;
// //     const double curval = sumtgt2*sumw/sumtgt;
// //     const double curval = -sumtgt/(sumtgt2*sumw);
// //     const double curval = sumtgt2*sumw - sumtgt;
// //     const double curval = sumw*sumtgt2;
//     double curval = sumtgt2 - sumtgt*sumtgt/sumw;
//     double curval = sumtgt2 - sumtgt*sumtgt/sumw;
//        double curval = sumw*sumtgtmax;
//     double curval = sumtgtmax;
//     double dh = sumtgtmax;
//     double ds = (sumtgt - 8.*sumtgt2)/(8.*sumw);
//     double curval = 8.*ds*sumtgt2 + 4.*sumw*ds*ds - ds*sumtgt + dh*sumtgt2 + sumw*dh*ds + 0.5*sumw*dh;
//     double curval = 0.5*sumw*dh;
//     
//        printf("curval = %5f, sumw = %5f, sumtgtmax = %5f, sumtgt = %5f, sumtgt2 = %5f\n",curval,sumw,sumtgtmax,sumtgt,sumtgt2,dh,ds);
       
//     double curval = -0.5*dx.transpose()*d2f*dx;
    
//     double dh = sumtgt/sumtgt2;
//     double ds = -0.5*log(sumtgt3/sumw);
//     double exp2ds = vdt::fast_exp(2.*ds);
//     double curval = -sumw*ds - 0.5*dh*exp2ds*sumtgt + 0.5*dh*dh*exp2ds*sumtgt2 + 0.5*exp2ds*sumtgt3;
//     double curval = sumw*(0.5-ds);
    
//     printf("curval = %5f, dh = %5f, ds = %5f, sumtgt = %5f, sumtgt2 = %5f, sumtgt3 = %5f\n",curval,dh,ds,sumtgt,sumtgt2,sumtgt3);
//     printf("curval = %5f, responseh = %5f, responsef = %5f\n",curval,dx[0],dx[1]);
    
//     double responsesimple = -sumtgt/sumtgt3;
//     printf("sumtgt = %5f, sumtgt3 = %5f, responsesimple = %5f\n",sumtgt,sumtgt3,responsesimple);
    
//     double curresponse = std::max(fShrinkage*sumtgt/sumw,-sumfmin);
//     double curval = sumtgt2 + sumw*curresponse*curresponse - 2.*curresponse*sumtgt;
//     double curval = sumw*sumtgtmax;
//     double curval = -sumtgt/sumw;
    
//     double curresponse = sumtgtmax;
//     double curval = sumtgtmax*volume;
//     double curval = (sumtgtmax-sumtgtmin)*sumtgt2;
    
//     printf("sumtgtmax = %5f, sumtgtmin = %5f, sumtgt2 = %5f, curval = %5f\n",sumtgtmax,sumtgtmin,sumtgt2,curval);
    
//     printf("tgtmax = %5f, sumtgtmax = %5f\n",tgtmax,sumtgtmax);
    
//     curval = std::max(curval,-sumtgtmin);
    
    //printf("start heavy loop\n");
    //loop over all bins and compute improvement in weighted variance of target for each split
    //This loop is relatively expensive and should auto-vectorize in the appropriate compiler/settings
//     GBRArrayUtils::FillSepGainsMC(_sumtgts[ivar],_sumtgt2s[ivar],_sumfmins[ivar], _sumfminsr[ivar], _sumws[ivar], _bsepgains[ivar], curval, sumtgt, sumtgt2, sumw, nbins, fShrinkage);
//     GBRArrayUtils::FillSepGainsMCEnvelope(_sumtgts[ivar],_sumtgt2s[ivar],_sumtgtmaxs[ivar], _sumtgtmaxsr[ivar], _sumfmins[ivar], _sumfminsr[ivar], _sumws[ivar], _bsepgains[ivar], curval, sumtgt,sumtgt2, sumtgtmax, sumtgtmin, sumw, nbins);
    
//     GBRArrayUtils::FillSepGainsMCMatrix(_sumtgts[ivar],_sumtgt2s[ivar],_sumtgt3s[ivar], _sumtgt4s[ivar],_sumtgt5s[ivar], _sumws[ivar], _bsepgains[ivar], curval, sumtgt, sumtgt2, sumtgt3,sumtgt4,sumtgt5, sumw, nbins);

//     double curval;
// //     if (fForest->Trees().size()<150) {
//     if (1) {
//       curval = sumtgt2 - sumtgt*sumtgt/sumw;
//       GBRArrayUtils::FillSepGainsMC(_sumtgts[ivar],_sumtgt2s[ivar],_sumfmins[ivar], _sumfminsr[ivar], _sumws[ivar], _bsepgains[ivar], curval, sumtgt, sumtgt2, sumw, nbins, fShrinkage);
//     }
//     else {
//       curval = 0.5*sumw*dh;
//       GBRArrayUtils::FillSepGainsMCEnvelope(_sumtgts[ivar],_sumtgt2s[ivar],_sumtgtmaxs[ivar], _sumtgtmaxsr[ivar], _sumfmins[ivar], _sumfminsr[ivar], _sumws[ivar], _bsepgains[ivar], curval, sumtgt,sumtgt2, sumtgtmax, sumtgtmin, sumw, nbins);
//     }
 
//     double curval = 2.*sumw*sumtgtmax - 0.5*vdt::fast_exp(-sumtgtmax)*sumtgt2;
//     double curval = -vdt::fast_exp(-sumtgtmax)*sumtgt2;
    
//     double dh = sumtgtmax;
//     double df = log( sumw / (dh*sumtgt2 - sumtgt) );
//     double curval = sumw*(1.-df);

//     printf("dh = %5f, df = %5f, curval = %5f\n",dh,df,curval);
  
//     double dh = sumtgtmax;
//     double expndf = sumtgt2/(sumw*dh+sumtgt3);
//     double curval = -expndf*sumtgt;
    
//     double curval = -vdt::fast_exp(-sumtgtmax)*sumtgt2;
    
//     double dh = sumtgtmax;
//     double dw = -sumfmin + dh;
//     double dhmdl = dh - dl;
//     double curval = dhmdl*sumtgt - 0.5*dhmdl*dhmdl*sumtgt2 - dh*sumtgt3 + 0.5*dh*dh*sumtgt4;
//     double curval = dw*sumtgt - 0.5*dw*dw*sumtgt2 + dw*dw*dw*sumtgt3/3. - 0.25*dw*dw*dw*dw*sumtgt4;
//     double curval = dw*sumtgt - 0.5*dw*dw*sumtgt2;
    
//     double df = sumtgtmax;
//     double ds = log( -sumw / (vdt::fast_exp(-df)*sumtgt - sumtgt2) );
//     double curval = sumw*(1.-ds);
    
//     printf("df = %5f, ds = %5f, sumtgt = %5f, sumtgt2 = %5f, sumtgt3 = %5f, sumtgt4 = %5f, curval = %5f\n",df,ds,sumtgt,sumtgt2,sumtgt3,sumtgt4,curval);
    
//     double dh = sumtgtmax;
//     double dh = sumtgt/sumw;
//     double df = sumtgtmin;
//     double curval = -2.*dh*sumtgt + dh*dh - vdt::fast_exp(-df)*sumtgt2;
//     double curval = - vdt::fast_exp(-df)*sumtgt2;
//     double curval = - vdt::fast_exp(-df)*sumtgt2 + sumtgt4 - dh*sumtgt3;
    
//     double curval = -0.5*sumtgt*sumtgt/sumtgt2;
//     double curval = sumtgtmax - 0.*sumtgt2*exp(-sumtgtmax);
    
//     const double k = 5.;
//     double dh = sumtgtmax;
// //     double dw = sumtgtmin + dh;
// //     double curval = sumtgt*dw - 0.5*sumtgt2*dw*dw + (1./3.)*sumtgt3*dw*dw*dw;
//     double curval = sumtgt*dh - 0.5*sumtgt2*dh*dh + (1./3.)*sumtgt3*dh*dh*dh -0.25*sumtgt4*dh*dh*dh*dh - k*(-sumtgt2*dh + sumtgt3*dh*dh - sumtgt4*dh*dh*dh);
//     
// //     printf("sumtgt = %5f, sumtgt2 = %5f, curval = %5f\n",sumtgt,sumtgt2,curval);
//     
//     GBRArrayUtils::FillSepGainsMCEnvelope(_sumtgts[ivar],_sumtgt2s[ivar],_sumtgt3s[ivar],_sumtgt4s[ivar],_sumtgtmaxs[ivar], _sumtgtmaxsr[ivar], _sumfmins[ivar], _sumfminsr[ivar], _sumws[ivar], _bsepgains[ivar], curval, sumtgt,sumtgt2, sumtgt3, sumtgt4, sumtgtmax, sumtgtmin, sumw, nbins);

//     const double sigma = 1e-4;
//     const double k = 0.5/sigma/sigma;
    
//     const double k = 10.;
    
//     double dh = std::max(0.,sumtgt/sumw);
//     double curval = k*(sumtgt2 - 2.*dh*sumtgt + sumw*dh*dh);
//     double curval = k*(sumtgt2 - sumtgt*sumtgt/sumw);
//     double curval = -k*sumtgt2*sumtgt2/sumtgt;
    
//     double dh = std::max(0.,sumtgtmax);
//     double curval = k*sumw*sumtgtmin + k*sumw*std::max(0.,sumtgtmax);
    
//     double dh = std::max(0.,sumtgt/sumw);
//     double dh = sumtgt/sumw;
//     double curval = k*sumw*sumtgtmin + 1e4*(sumw*dh*dh - 2.*dh*sumtgt);
    
//     double curval = k*sumw*sumtgtmin + 20.*k*sumw*sumtgtmax;    
    
//     double curval = k*sumw*sumtgtmin;
//     double curval = k*sumw*sumtgtmin + k*sumw*std::max(0.,sumtgtmax);
//     double curval = k*sumtgt;
    
    

    
//     const double factor = 1./(1.+k);
    
//     double curval = std::max(0.,sumtgtmax)*sumtgt2;
    
//     double dl = sumtgtmin;
//     double dh = std::max(0.,sumtgtmax - k*dl);    
//     double curval = sumtgt2*(factor*(dh+k*dl) - dl);
    
//     double curval = sumtgt2*std::max(0.,sumtgtmax - sumtgtmin);
    
//     double curval = k*sumw*std::max(0.,sumtgtmax);
//     double curval = k*sumw*sumtgtmax;
    
//     double curval = sumtgt2*std::max(0.,sumtgtmax);
    
//     double curval = (sumtgt < 0. && sumtgt2 > 0.) ? -0.5*sumtgt*sumtgt/sumtgt2 : 0.;
    
//     double dh = std::max(0.,sumtgt/sumw);
//     double curval = k*(dh*dh - 2.*dh*sumtgt);
    
//     double dh = doenv ? std::max(0.,sumtgtmax) : sumtgt/sumw;
//     double curval = doenv ? kenv*sumw*dh : k*(dh*dh - 2.*dh*sumtgt);
   
//     double dh = doenv ? std::max(0.,-sumtgt/sumtgt2) : sumtgt/sumw;
//     double curval = doenv ? sumtgt*dh + 0.5*sumtgt2*dh*dh : k*(dh*dh - 2.*dh*sumtgt);
   
//        double dh = doenv ? std::max(0.,-sumtgt/sumtgt2) : -sumtgt/sumtgt2;
    
//        double dh = doenv ? (sumtgt2 > 0. ? std::max(0.,-sumtgt/sumtgt2) : 0.) : -sumtgt/sumtgt2;
    
//        double dh = doenv ? std::max(0.,-sumtgt/sumtgt2) : -sumtgt/sumtgt2;
//        double curval = sumtgt*dh + 0.5*sumtgt2*dh*dh;
   
//        constexpr double sigmaenv = 0.01;
//        constexpr double kenv = 0.5/sigmaenv/sigmaenv;
//        constexpr double kenv = 10.;
       
//        double dh = 0 ? sumtgt/sumw : -sumtgt/sumtgt2;
//        double curval = 0 ? -sumw*sumtgt2/sumtgt/sumtgt : sumtgt*dh + 0.5*sumtgt2*dh*dh; 
   
//        double dh = doenv ? std::max(0.,-sumtgt/sumtgt2) : -sumtgt/sumtgt2;
//        double curval = sumtgt*dh + 0.5*sumtgt2*dh*dh; 
       
//        double dh = 0 ? std::max(0.,-sumtgt/sumtgt2) : -sumtgt/sumtgt2;
//        double curval = sumtgt*dh + 0.5*sumtgt2*dh*dh; 
    
    
//     double curval = -k*sumtgt*sumtgt/sumw;
    
//     double curval = fDoEnvelope && !doenv ? kenv*sumw*sumtgtmax : -k*sumtgt*sumtgt/sumw;
    
//     double curval = fDoEnvelope && !doenv ? -0.5*sumtgt*sumtgt/sumtgt2 : -k*sumtgt*sumtgt/sumw;
    
    double curval = -0.5*sumtgt*sumtgt/sumtgt2;
    
/*       double dh = doenv ? 0.5*(-2.*kenv*sumtgt + sqrt(4.*kenv*kenv*sumtgt*sumtgt + 8.*sumw*kenv*sumtgt2))/sumw : -sumtgt/sumtgt2;
       double curval = doenv ? -kenv*sumtgt*sumtgt/sumtgt2 + sumw*log(dh)  : sumtgt*dh + 0.5*sumtgt2*dh*dh; */   
//        double curval = doenv ? dh*sumtgt2 : sumtgt*dh + 0.5*sumtgt2*dh*dh;    
    
//     double curval = doenv ? sumtgt2*dh : k*(dh*dh - 2.*dh*sumtgt);
    
//        double dh = std::max(0., sumtgt/sumtgt2);
//        double curval = k*(dh*dh*sumtgt2 - 2.*dh*sumtgt);
    
//     double dh = std::max(0.,sumtgt2/sumtgt);
//     double dl = 0.5*vdt::fast_log((sumtgt3 - 2.*dh*sumtgt2 + dh*dh*sumtgt)/sumw);    
//     double curval = sumw*dl + 0.5*vdt::fast_exp(-2*dl)*(sumtgt3 - 2.*dh*sumtgt2 + dh*dh*sumtgt);
    
//     const double k = 1.;
    
//     double dh = std::max(0.,sumtgtmax);
//     double curval = sumtgt2*dh;
    
//     double dh = sumtgtmax;
//     double curval = k*sumw*sumtgtmax;
//     double curval = -k*exp(-sumtgtmax)*sumtgt2 + sumw*sumtgtmax;;
    
    
//     double curval = -0.5*sumtgt*sumtgt/sumtgt2;
//     printf("sumtgt = %5e, sumtgt2 = %5e, curval = %5e\n",sumtgt,sumtgt2,curval);
// #pragma omp parallel for simd
    for (unsigned int ibin=0; ibin<nbins; ++ibin) {
//       double cutvalnow = _varvals[ivar][ibin];
      
//       double leftwsum = _sumws[ivar][ibin];
//       double lefttgtmax = _sumtgtmaxs[ivar][ibin];
//       double lefttgtmin = _sumfmins[ivar][ibin];
      double lefttgtsum = _sumtgts[ivar][ibin];
      double lefttgt2sum = _sumtgt2s[ivar][ibin];
//       double lefttgt3sum = _sumtgt3s[ivar][ibin];
//       
//       double rightwsum = sumw - leftwsum;
//       double righttgtmax = _sumtgtmaxsr[ivar][ibin];
//       double righttgtmin = _sumfminsr[ivar][ibin];
      double righttgtsum = sumtgt - lefttgtsum;
      double righttgt2sum = sumtgt2 - lefttgt2sum;
//       double righttgt3sum = sumtgt3 - lefttgt3sum;
//       
      
//       double dhleft = std::max(0., lefttgtsum/lefttgt2sum);
//       double leftval = k*(dhleft*dhleft*lefttgt2sum - 2.*dhleft*lefttgtsum);
//       
//       double dhright = std::max(0., righttgtsum/righttgt2sum);
//       double rightval = k*(dhright*dhright*righttgt2sum - 2.*dhright*righttgtsum);
      
//       double leftval = -0.5*lefttgtsum*lefttgtsum/lefttgt2sum;
//       double rightval = -0.5*righttgtsum*righttgtsum/righttgt2sum;
      
//       double leftval = -k*exp(-lefttgtmax)*lefttgt2sum + leftwsum*lefttgtmax;
//       double rightval = -k*exp(-righttgtmax)*righttgt2sum + rightwsum*righttgtmax;
//       double leftval = k*leftwsum*lefttgtmax;
//       double rightval = k*rightwsum*righttgtmax;
      
//       double dhleft = std::max(0.,lefttgtmax);
//       double leftval = lefttgt2sum*dhleft;
      
//       double dhright = std::max(0.,righttgtmax);
//       double rightval = righttgt2sum*dhright;      
      
//     double dhleft = std::max(0.,lefttgt2sum/lefttgtsum);
//     double dlleft = 0.5*vdt::fast_log((lefttgt3sum - 2.*dhleft*lefttgt2sum + dhleft*dhleft*lefttgtsum)/leftwsum);    
//     double leftval = leftwsum*dlleft + 0.5*vdt::fast_exp(-2*dlleft)*(lefttgt3sum - 2.*dhleft*lefttgt2sum + dhleft*dhleft*lefttgtsum);      

//     double dhright = std::max(0.,righttgt2sum/righttgtsum);
//     double dlright = 0.5*vdt::fast_log((righttgt3sum - 2.*dhright*righttgt2sum + dhright*dhright*righttgtsum)/rightwsum);    
//     double rightval = rightwsum*dlright + 0.5*vdt::fast_exp(-2*dlright)*(righttgt3sum - 2.*dhright*righttgt2sum + dhright*dhright*righttgtsum);      
    
//       double leftval = 0.;
//       double rightval = 0.;
//       
//       for (int iev=0;iev<nev;++iev) {
//         int quant = _quants[ivar][iev];
//         int bin = (quant-offset)>>pscale;
//         unsigned int ubin = bin;
//         
//         if (ubin>ibin) {
//           rightval += k*log(_tgt2vals[iev] + righttgtmax);
//         }
//         else {
//           leftval += k*log(_tgt2vals[iev] + lefttgtmax);
//         }
//       }
      
//       double leftval = std::max(0.,lefttgtmax)*lefttgt2sum;
//       double rightval = std::max(0.,righttgtmax)*righttgt2sum;
   
/*      double dlleft = lefttgtmin;
      double dhleft = std::max(0.,lefttgtmax - k*dlleft);    
      double leftval = lefttgt2sum*(factor*(dhleft+k*dlleft) - dlleft);      

      double dlright = righttgtmin;
      double dhright = std::max(0.,righttgtmax - k*dlright);    
      double rightval = righttgt2sum*(factor*(dhright+k*dlright) - dlright); */  

//       double leftval = lefttgt2sum*std::max(0.,lefttgtmax - lefttgtmin);
//       double rightval = righttgt2sum*std::max(0.,righttgtmax - righttgtmin);

//       double leftval = k*leftwsum*std::max(0.,lefttgtmax);
//       double rightval = k*rightwsum*std::max(0.,righttgtmax);

//       double leftval = lefttgt2sum*std::max(0.,lefttgtmax);
//       double rightval = righttgt2sum*std::max(0.,righttgtmax);

//       double leftval = k*leftwsum*lefttgtmax;
//       double rightval = k*rightwsum*righttgtmax;
      
/*      double leftval = lefttgtmax*lefttgt2sum;
      double rightval = righttgtmax*righttgt2sum;   */   

//       double leftval = (lefttgtsum < 0. && lefttgt2sum > 0.) ? -0.5*lefttgtsum*lefttgtsum/lefttgt2sum : 0.;
//       double rightval = (righttgtsum < 0. && righttgt2sum > 0.) ? -0.5*righttgtsum*righttgtsum/righttgt2sum : 0.;

/*      double dhleft = std::max(0.,lefttgtsum/leftwsum);
      double leftval = k*(dhleft*dhleft - 2.*dhleft*lefttgtsum);

      double dhright = std::max(0.,righttgtsum/rightwsum);
      double rightval = k*(dhright*dhright - 2.*dhright*righttgtsum); */  

//       double dhleft = doenv ? std::max(0.,lefttgtmax) : lefttgtsum/leftwsum;
//       double leftval = doenv ? dhleft*lefttgt2sum : k*(dhleft*dhleft - 2.*dhleft*lefttgtsum);
// 
//       double dhright = doenv ? std::max(0.,righttgtmax) : righttgtsum/rightwsum;
//       double rightval = doenv ? dhright*righttgt2sum : k*(dhright*dhright - 2.*dhright*righttgtsum);

//       double dhleft = doenv ? std::max(0.,lefttgtmax) : lefttgtsum/leftwsum;
//       double leftval = doenv ? kenv*leftwsum*dhleft : k*(dhleft*dhleft - 2.*dhleft*lefttgtsum);
//       
//       double dhright = doenv ? std::max(0.,righttgtmax) : righttgtsum/rightwsum;
//       double rightval = doenv ? kenv*rightwsum*dhright : k*(dhright*dhright - 2.*dhright*righttgtsum);

//       double dhleft = doenv ? std::max(0.,-lefttgtsum/lefttgt2sum) : lefttgtsum/leftwsum;
//       double leftval = doenv ? lefttgtsum*dhleft + 0.5*lefttgt2sum*dhleft*dhleft : k*(dhleft*dhleft - 2.*dhleft*lefttgtsum);
//       
//       double dhright = doenv ? std::max(0.,-righttgtsum/righttgt2sum) : righttgtsum/rightwsum;
//       double rightval = doenv ? righttgtsum*dhright + 0.5*righttgt2sum*dhright*dhright : k*(dhright*dhright - 2.*dhright*righttgtsum);

//       double dhleft = doenv ? (lefttgt2sum > 0. ? std::max(0.,-lefttgtsum/lefttgt2sum) : 0.) : -lefttgtsum/lefttgt2sum;
//       double leftval = lefttgtsum*dhleft + 0.5*lefttgt2sum*dhleft*dhleft;
// 
//       double dhright = doenv ? (righttgt2sum > 0. ? std::max(0.,-righttgtsum/righttgt2sum) : 0.) : -righttgtsum/righttgt2sum;
//       double rightval = righttgtsum*dhright + 0.5*righttgt2sum*dhright*dhright;

//       double dhleft = doenv ? std::max(0.,-lefttgtsum/lefttgt2sum) : -lefttgtsum/lefttgt2sum;
//       double leftval = lefttgtsum*dhleft + 0.5*lefttgt2sum*dhleft*dhleft;
//       
//       double dhright = doenv ? std::max(0.,-righttgtsum/righttgt2sum) : -righttgtsum/righttgt2sum;
//       double rightval = righttgtsum*dhright + 0.5*righttgt2sum*dhright*dhright;
   
//        double dhleft = doenv ? std::max(0.,lefttgtmax) : -lefttgtsum/lefttgt2sum;
//        double leftval = doenv ? dhleft*lefttgt2sum : lefttgtsum*dhleft + 0.5*lefttgt2sum*dhleft*dhleft;  
// 
//        double dhright = doenv ? std::max(0.,righttgtmax) : -righttgtsum/righttgt2sum;
//        double rightval = doenv ? dhright*righttgt2sum : righttgtsum*dhright + 0.5*righttgt2sum*dhright*dhright;
       
//        double dhleft = doenv ? std::max(0.,lefttgtmax) : -lefttgtsum/lefttgt2sum;
//        double leftval = doenv ? dhleft*lefttgt2sum : lefttgtsum*dhleft + 0.5*lefttgt2sum*dhleft*dhleft;  
// 
//        double dhright = doenv ? std::max(0.,righttgtmax) : -righttgtsum/righttgt2sum;
//        double rightval = doenv ? dhright*righttgt2sum : righttgtsum*dhright + 0.5*righttgt2sum*dhright*dhright;  
      
/*       double dhleft = doenv ? lefttgt2sum/lefttgtsum : -sumtgt/sumtgt2;
       double leftval = doenv ? -kenv*lefttgtsum*lefttgtsum/lefttgt2sum : lefttgtsum*dhleft + 0.5*lefttgt2sum*dhleft*dhleft;

       double dhright = doenv ? righttgt2sum/righttgtsum : -sumtgt/sumtgt2;
       double rightval = doenv ? -kenv*righttgtsum*righttgtsum/righttgt2sum : righttgtsum*dhright + 0.5*righttgt2sum*dhright*dhright;   */    

//        double dhleft = 0 ? lefttgtmax : -lefttgtsum/lefttgt2sum;
//        double leftval = 0 ? kenv*leftwsum*dhleft : lefttgtsum*dhleft + 0.5*lefttgt2sum*dhleft*dhleft; 
// 
//        double dhright = 0 ? righttgtmax : -righttgtsum/righttgt2sum;
//        double rightval = 0 ? kenv*rightwsum*dhright : righttgtsum*dhright + 0.5*righttgt2sum*dhright*dhright; 
       
//        double dhleft = doenv ? 0.5*(-2.*kenv*lefttgtsum + sqrt(4.*kenv*kenv*lefttgtsum*lefttgtsum + 8.*leftwsum*kenv*lefttgt2sum))/leftwsum : -lefttgtsum/lefttgt2sum;
//        double leftval = doenv ? -kenv*lefttgtsum*lefttgtsum/lefttgt2sum + leftwsum*log(dhleft)  : lefttgtsum*dhleft + 0.5*lefttgt2sum*dhleft*dhleft;
// 
//        double dhright = doenv ? 0.5*(-2.*kenv*righttgtsum + sqrt(4.*kenv*kenv*righttgtsum*righttgtsum + 8.*rightwsum*kenv*righttgt2sum))/rightwsum : -righttgtsum/righttgt2sum;
//        double rightval = doenv ? -kenv*righttgtsum*righttgtsum/righttgt2sum + rightwsum*log(dhright)  : righttgtsum*dhright + 0.5*righttgt2sum*dhright*dhright;
       
//        double leftval = -k*lefttgtsum*lefttgtsum/leftwsum;
//        double rightval = -k*righttgtsum*righttgtsum/rightwsum;
       
//        double leftval = fDoEnvelope && !doenv ? kenv*leftwsum*lefttgtmax : -k*lefttgtsum*lefttgtsum/leftwsum;
//        double rightval = fDoEnvelope && !doenv ? kenv*rightwsum*righttgtmax : -k*righttgtsum*righttgtsum/rightwsum;
       
       
//        double leftval = fDoEnvelope && !doenv ? -0.5*lefttgtsum*lefttgtsum/lefttgt2sum : -k*lefttgtsum*lefttgtsum/leftwsum;
//        double rightval = fDoEnvelope && !doenv ? -0.5*righttgtsum*righttgtsum/righttgt2sum : -k*righttgtsum*righttgtsum/rightwsum;

       double leftval = -0.5*lefttgtsum*lefttgtsum/lefttgt2sum;
       double rightval = -0.5*righttgtsum*righttgtsum/righttgt2sum;
       
//        double dhleft = 0 ? std::max(0.,-lefttgtsum/lefttgt2sum) : -lefttgtsum/lefttgt2sum;
//        double leftval = lefttgtsum*dhleft + 0.5*lefttgt2sum*dhleft*dhleft; 
// 
//        double dhright = 0 ? std::max(0.,-righttgtsum/righttgt2sum) : -righttgtsum/righttgt2sum;
//        double rightval = righttgtsum*dhright + 0.5*righttgt2sum*dhright*dhright; 
       
//       double dhright = doenv ? std::max(0.,-righttgtsum/righttgt2sum) : -righttgtsum/righttgt2sum;
//       double rightval = righttgtsum*dhright + 0.5*righttgt2sum*dhright*dhright;
//       double dhleft = lefttgtsum/leftwsum;
//       double leftval = k*(dhleft*dhleft - 2.*dhleft*lefttgtsum);
// 
//       double dhright = righttgtsum/rightwsum;
//       double rightval = k*(dhright*dhright - 2.*dhright*righttgtsum); 
      
//       printf("ibin = %i, curval = %5e, leftval = %5e, rightval = %5e\n",ibin,curval,leftval,rightval);
      
      _bsepgains[ivar][ibin] = valscale*(curval - leftval - rightval);
    }
    
//     double curval = k*sumw*sumtgtmin + k*sumw*sumtgtmax;
//     if (0) {
//       GBRArrayUtils::FillSepGainsMCEnvelope(_sumtgts[ivar],_sumtgt2s[ivar],_sumtgt3s[ivar],_sumtgt4s[ivar],_sumtgtmaxs[ivar], _sumtgtmaxsr[ivar], _sumfmins[ivar], _sumfminsr[ivar], _sumws[ivar], _bsepgains[ivar], curval, sumtgt,sumtgt2, sumtgt3, sumtgt4, sumtgtmax, sumtgtmin, sumw, nbins);
//     }
//     for (unsigned int ibin=0; ibin<nbins; ++ibin) {
//       printf("ivar = %i, ibin = %i, sepgain = %5f\n",ivar,ibin,_bsepgains[ivar][ibin]);
//     }
    
    
    //printf("start final loop\n");
    //loop over computed variance improvements and select best split, respecting also minimum number of events per node
    //This loop cannot auto-vectorize, at least in gcc 4.6x due to the mixed type conditional, but it's relatively fast
    //in any case
    for (unsigned int ibin=0; ibin<nbins; ++ibin) {   

      if ( _bsepgains[ivar][ibin]>maxsepgain && std::isnormal(_bsepgains[ivar][ibin])) {
	
	bool passminweights = true;
	double totalweightleft = _sumws[ivar][ibin];
	double totalweightright = _sumws[ivar][nbins-1] - _sumws[ivar][ibin];
        
        int nevleft = _sumns[ivar][ibin];
        int nevright = sumn - _sumns[ivar][ibin];
	
	if (fMinWeightTotal>=0. && (totalweightleft<fMinWeightTotal || totalweightright<fMinWeightTotal) ) {
	  passminweights = false;
	}
	
	if (minevents>0 && (nevleft<minevents || nevright<minevents)) {
          passminweights = false;
        }
	
	bool passminsig = true;
	
//         double mincutsig = fMinCutSignificance;
	
	passminsig &= mincutsig<0. || (_bsepgains[ivar][ibin] > mincutsig);
	
	if (passminweights && passminsig) {
	  maxsepgain = _bsepgains[ivar][ibin];
	  bestbin = ibin;
	}
      }
      
    }
     
    cutval = _varvals[ivar][bestbin];
    nleft = _sumns[ivar][bestbin];
    nright = nev - nleft;
    
    _sepgains[ivar] = maxsepgain;
    _sepgainsigs[ivar] = _bsepgainsigs[ivar][bestbin];
    _cutvals[ivar] = cutval;
    _nlefts[ivar] = nleft;
    _nrights[ivar] = nright;
    _bestbins[ivar] = bestbin;
        
     //printf("done var %i\n",ivar);
  }
  

  
  double globalsepgain = 0.;
  for (int ivar=0; ivar<nvars; ++ivar) {
    if (_sepgains[ivar]>globalsepgain) {
      globalsepgain = _sepgains[ivar];
      bestvar = ivar;
    }
  }    
    
  //if no appropriate split found, make this node terminal
  if (globalsepgain<=0.) {
//     printf("globalsepgain<=0.\n");
    //no valid split found, making this node a leaf node
    //printf("thisidx = %i, globalsepgain = %5f, no valid split\n",thisidx, globalsepgain);
    tree.CutIndices().push_back(0);
    tree.CutVals().push_back(0);
    tree.LeftIndices().push_back(0);   
    tree.RightIndices().push_back(0);    
    
    tree.RightIndices()[thisidx] = -tree.Responses().size();
    tree.LeftIndices()[thisidx] = -tree.Responses().size();
    
    BuildLeaf(evts,tree,limits,doenv);
    return;
  }
  
  //fill vectors of event pointers for left and right nodes below this one
  std::vector<MCGBREvent*> leftevts;
  std::vector<MCGBREvent*> rightevts;
  
  leftevts.reserve(nev);
  rightevts.reserve(nev);
  
  int nleft = 0;
  int nright = 0;
  double sumwleft = 0.;
  double sumwright = 0.;
  
  for (std::vector<MCGBREvent*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
    if ((*it)->Var(bestvar)>_cutvals[bestvar]) {
//       printf("right: var = %5f, cutval = %5f\n",(*it)->Var(bestvar),_cutvals[bestvar]);
      ++nright;
      sumwright += (*it)->Weight();
      rightevts.push_back(*it);
    }
    else {
//       printf("left: var = %5f, cutval = %5f\n",(*it)->Var(bestvar),_cutvals[bestvar]);
      ++nleft;
      sumwleft += (*it)->Weight();
      leftevts.push_back(*it);
    }
  }
  
//   printf("cutval = %5f\n",_cutvals[bestvar]);
//   printf("nleft = %i, nright = %i, nlefts = %i, nrights = %i\n",nleft,nright,_nlefts[bestvar],_nrights[bestvar]);

//   printf("bestvar = %i, nleft = %i, nright = %i, nlefts[bestvar] = %i, nrights[bestvar] = %i, cutvals[bestvar] = %5f, nmismatch = %i\n",bestvar,nleft,nright,_nlefts[bestvar],_nrights[bestvar],_cutvals[bestvar],nmismatch);
  
  assert(_nlefts[bestvar]==nleft);
  assert(_nrights[bestvar]==nright);
  
  //printf("nleft = %i, nright = %i\n",nleft,nright);
  
  
  double bestcutval = _cutvals[bestvar];  
  
  //fill intermediate node
  tree.CutIndices().push_back(bestvar);
  tree.CutVals().push_back(_cutvals[bestvar]);
  tree.LeftIndices().push_back(0);   
  tree.RightIndices().push_back(0);  
  
  //check if left node is terminal
  //bool termleft = nleft<=(2*minevents) || depth==fMaxDepth;
//   bool termleft = sumwleft<=(2*minevents) || (fMaxDepth>=0 && depth==fMaxDepth) || (fMaxNodes>=0 && int(tree.Responses().size())>=fMaxNodes) ;
  bool termleft = nleft<(2*minevents) || (fMaxDepth>=0 && depth==fMaxDepth) || (maxnodes>=0 && int(tree.Responses().size())>=maxnodes) ;
  if (termleft) tree.LeftIndices()[thisidx] = -tree.Responses().size();
  else tree.LeftIndices()[thisidx] = tree.CutIndices().size();
  
  //printf("this idx = %i, termleft = %i, nleft = %i, minevents = %i\n",thisidx,  termleft,nleft,minevents);  
  
  //build left node as appropriate
  std::vector<std::pair<float,float> > limitsleft(limits);
  limitsleft[bestvar].second = bestcutval;  
  //printf("bestvar = %i, limlow = %5f, limhigh = %5f, limleftlow = %5f, limlefthigh = %5f\n",bestvar,limits[bestvar].first,limits[bestvar].second,limitsleft[bestvar].first,limitsleft[bestvar].second);
  if (termleft) {  
    BuildLeaf(leftevts,tree,limitsleft,doenv);
  }
  else {  
    TrainTree(leftevts,sumwleft,tree,depth+1,limitsleft,usetarget,doenv);  
  }
  
  //check if right node is terminal
  //bool termright = nright<=(2*minevents) || depth==fMaxDepth;
//   bool termright = sumwright<=(2*minevents) || (fMaxDepth>=0 && depth==fMaxDepth) || (fMaxNodes>=0 && int(tree.Responses().size())>=fMaxNodes);
  bool termright = nright<(2*minevents) || (fMaxDepth>=0 && depth==fMaxDepth) || (maxnodes>=0 && int(tree.Responses().size())>=maxnodes);
  if (termright) tree.RightIndices()[thisidx] = -tree.Responses().size();
  else tree.RightIndices()[thisidx] = tree.CutIndices().size();
    
  //printf("this idx = %i, termright = %i, nright = %i, minevents = %i\n",thisidx,  termright,nright,minevents);    
  
  //build right node as appropriate
  std::vector<std::pair<float,float> > limitsright(limits);
  limitsright[bestvar].first = bestcutval;  
  //printf("bestvar = %i, limlow = %5f, limhigh = %5f, limrightlow = %5f, limrighthigh = %5f\n",bestvar, limits[bestvar].first,limits[bestvar].second,limitsright[bestvar].first,limitsright[bestvar].second);  
  if (termright) {  
    BuildLeaf(rightevts,tree,limitsright,doenv);
  }
  else {      
    TrainTree(rightevts,sumwright,tree,depth+1,limitsright,usetarget,doenv);  
  }
  
}

  

  


//_______________________________________________________________________
void MCGBRIntegrator::BuildLeaf(const std::vector<MCGBREvent*> &evts, MCGBRTreeD &tree, const std::vector<std::pair<float,float> > &limits, bool doenv) {


  
//   double sumw = 0.;
  double sumtgt = 0.;
  double sumtgt2 = 0.;
//   double tgtmax = std::numeric_limits<double>::lowest();
 
  if (doenv) {
    for (std::vector<MCGBREvent*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
//       double weight = (*it)->Weight();
      const double arg = (*it)->Arg();
      const double tgtval = -arg;
      const double tgt2val = 1.;
      
      sumtgt += tgtval;
      sumtgt2 += tgt2val;
    }
  }
  else {
    for (std::vector<MCGBREvent*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
//       double weight = (*it)->Weight();
      const double arg = (*it)->ArgLog();
      const double tgtval = -arg;
      const double tgt2val = 1.;
      
      sumtgt += tgtval;
      sumtgt2 += tgt2val;
    }
  }
  
//   for (std::vector<MCGBREvent*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
//     double weight = (*it)->Weight();
//     double funcval = (*it)->FuncVal();
//     double target =  (*it)->Target();
//     double targetmin =  (*it)->TargetMin();    
// //     double funcvalalt = (*it)->FuncValAlt();
// //     double tgtval = doenv ? exp(targetmin) - target : log(funcval) - targetmin;
// 
// //     double tgtval = doenv ? NLSkewGausDMu(exp(targetmin),target) : log(funcval) - targetmin;
// //     double tgt2val = doenv ? NLSkewGausD2Mu(exp(targetmin),target) : 0.;
//     
// /*    double tgtval = doenv ? NLSkewGausDMu(funcval,target,doenv) : NLSkewGausDMu(log(funcval),targetmin,doenv);
//     double tgt2val = doenv ? NLSkewGausD2Mu(funcval,target,doenv ) : NLSkewGausD2Mu(log(funcval),targetmin,doenv);  */ 
//     
// //     double tgtval = doenv ? NLSkewGausDMu(exp(targetmin),target,doenv) : NLSkewGausDMu(log(funcval),targetmin,doenv);
// //     double tgt2val = doenv ? NLSkewGausD2Mu(exp(targetmin),target,doenv ) : NLSkewGausD2Mu(log(funcval),targetmin,doenv);
//     
// //     double tgtval = doenv ? NLSkewGausDMu(exp(targetmin),target,funcval,doenv) : NLSkewGausDMu(log(funcval),targetmin,funcval,doenv);
// //     double tgt2val = doenv ? NLSkewGausD2Mu(exp(targetmin),target,funcval,doenv ) : NLSkewGausD2Mu(log(funcval),targetmin,funcval,doenv);
// 
// //     double tgtval = doenv ? NLLogNormDMu(exp(targetmin),target) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLLogNormD2Mu(exp(targetmin),target) : NLNormD2Mu(log(funcval),targetmin);
//     
// //     double tgtval = doenv ? NLLogNormDMu(exp(targetmin),target) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLLogNormD2Mu(exp(targetmin),target) : NLNormD2Mu(log(funcval),targetmin);
//     
// //     double tgtval = doenv ? NLLogNormDMu(exp(targetmin),target,funcvalalt) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLLogNormD2Mu(exp(targetmin),target,funcvalalt) : NLNormD2Mu(log(funcval),targetmin);
// 
// //     double tgtval = doenv ? exp(targetmin) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? exp(2.*targetmin) : NLNormD2Mu(log(funcval),targetmin);
// 
// //     double tgtval = doenv ?  NLLogNormDMu(exp(targetmin),target) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLLogNormD2Mu(exp(targetmin),target) : NLNormD2Mu(log(funcval),targetmin);
// 
// //     double tgtval = doenv ?  NLNormDMu(targetmin,0.) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLNormD2Mu(targetmin,0.) : NLNormD2Mu(log(funcval),targetmin);
// 
// //     double tgtval = doenv ?  NLNormDMu(log(std::max(1e-6,exp(targetmin)-target)),0.) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLNormD2Mu(log(std::max(1e-6,exp(targetmin)-target)),0.) : NLNormD2Mu(log(funcval),targetmin);
//    
// //     double sigma = doenv ? 0.01 : 0.01;
// //     double tgtval = doenv ?  NLNormDMu(log(std::max(exp(-24.),exp(targetmin)-target)),0.,sigma) : NLNormDMu(log(funcval),targetmin,sigma);
// //     double tgt2val = doenv ? NLNormD2Mu(log(std::max(exp(-24.),exp(targetmin)-target)),0.,sigma) : NLNormD2Mu(log(funcval),targetmin,sigma);
//     
// //     double tgtval = doenv ? NLNormDMu(targetmin,target) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLNormD2Mu(targetmin,target) : NLNormD2Mu(log(funcval),targetmin);
//     
// //     printf("funcval = %5e, funcvalalt = %5e\n",funcval,funcvalalt);
//     
// //     double tgtval = doenv ? NLNormDMu(exp(targetmin),target) : NLNormDMu(log(funcval),targetmin);
// //     double tgt2val = doenv ? NLNormD2Mu(exp(targetmin),target) : NLNormD2Mu(log(funcval),targetmin);
//     
// //     double tgtval = doenv ? NLSkewGausDMu(funcval,target,funcval,doenv) : NLSkewGausDMu(log(funcval),targetmin,funcval,doenv);
// //     double tgt2val = doenv ? NLSkewGausD2Mu(funcval,target,funcval,doenv ) : NLSkewGausD2Mu(log(funcval),targetmin,funcval,doenv);
//     
// /*    double tgtval = doenv ? NLSkewGausDMu(funcval,target) : log(funcval) - targetmin;
//     double tgt2val = doenv ? NLSkewGausD2Mu(funcval,target) : 0.;  */ 
//     
// //     if (0) {
// //       if (!std::isnormal(tgtval) || !std::isnormal(tgt2val)) {
// //         printf("targetmin = %5e, exp(targetmin) = %5e, target = %5e, tgtval = %5e, tgt2val = %5e\n",targetmin,exp(targetmin),target,tgtval,tgt2val);
// //       }
// //     }
//     
// //     double tgtval = doenv ? log(std::max(exp(-24.),exp(targetmin)-target)) : log(funcval) - targetmin;
// 
// //     double tgtval = doenv ? log(std::max(exp(-24.),exp(targetmin)-target)) : (fDoEnvelope ? NLSkewGausDMu(log(funcval),targetmin) : log(funcval) - targetmin);
// //     double tgt2val = doenv ? 0. : (fDoEnvelope ? NLSkewGausD2Mu(log(funcval),targetmin) : 0.);
// 
// //     double tgtval; = doenv ? log(std::max(exp(-24.),exp(targetmin)-target)) : (fDoEnvelope ? NLSkewGausDMu(log(funcval),targetmin) : log(funcval) - targetmin);
// //     double tgt2val; = doenv ? 0. : (fDoEnvelope ? NLSkewGausD2Mu(log(funcval),targetmin) : 0.);
// 
//     double px = doenv ? log(std::max(exp(-128.),exp(targetmin)-target)) : log(funcval);
//     double pmu = doenv ? 0. : targetmin;
// //     double sigmabase = 0.05;
// //     double psigma = doenv ? 1. : 0.3;
// //     double psigma = sigmabase;
//     double psigma = 0.05;
//     double tgtval = fDoEnvelope && !doenv ? NLSkewGausDMu(px,pmu) : NLNormDMu(px,pmu,psigma);
//     double tgt2val = fDoEnvelope && !doenv ? NLSkewGausD2Mu(px,pmu) : NLNormD2Mu(px,pmu,psigma);
// 
// //     double px = doenv ? log(std::max(exp(-128.),funcval-target)) : log(funcval);
// //     double pmu = doenv ? 0. : targetmin;
// //     double sigmabase = 0.05;
// //     double psigma = sigmabase;
// //     double tgtval = fDoEnvelope && !doenv ? NLSkewGausDMu(px,pmu) : NLNormDMu(px,pmu,psigma);
// //     double tgt2val = fDoEnvelope && !doenv ? NLSkewGausD2Mu(px,pmu) : NLNormD2Mu(px,pmu,psigma);
//     
//     
// //     double tgtval = fDoEnvelope && !doenv ? 
//     
// //     if (tgtval > tgtmax) {
// //       tgtmax = tgtval;
// //     }
//     
// //     sumw += weight;
//     sumtgt += weight*tgtval;
//     sumtgt2 += weight*tgt2val;
// 
//     
//   }
  
//   if (doenv) printf("sumtgt = %5e, sumtgt2 = %5e, dh = %5e\n",sumtgt,sumtgt2,-sumtgt/sumtgt2);
  
//   double dh = doenv ? std::max(0.,tgtmax) : sumtgt/sumw;
//   double dh = doenv ? std::max(0.,-sumtgt/sumtgt2) : sumtgt/sumw;
//   double dh = doenv ? (sumtgt2 > 0. ? std::max(0.,-sumtgt/sumtgt2) : 0.) : -sumtgt/sumtgt2;
//   double dh = doenv ? std::max(0.,-sumtgt/sumtgt2) : -sumtgt/sumtgt2;
//   double dh = doenv ? std::max(0.,tgtmax) : -sumtgt/sumtgt2;
//   double response = fShrinkage*dh;
  
//   constexpr double sigmaenv = 0.01;
//   constexpr double kenv = 0.5/sigmaenv/sigmaenv;
  
//   double dh = doenv ? 0.5*(-2.*kenv*sumtgt + sqrt(4.*kenv*kenv*sumtgt*sumtgt + 8.*sumw*kenv*sumtgt2))/sumw : -sumtgt/sumtgt2;
//   double response = doenv ? dh : fShrinkage*dh;

//     double dh =  fDoEnvelope && !doenv ? -sumtgt/sumtgt2 : sumtgt/sumw;

    double dh = -sumtgt/sumtgt2;
    double response = doenv ? fShrinkage*vdt::fast_exp(dh) : fShrinkage*dh;
//     double response = doenv ? std::min(exp(fShrinkage*dh),exp(dh/fShrinkage)) : fShrinkage*dh;

//   double dh = doenv ? std::max(0.,-sumtgt/sumtgt2) : -sumtgt/sumtgt2;
//   double dh =  -sumtgt/sumtgt2;
//   double response = doenv ? fShrinkage*exp(dh) : fShrinkage*dh;
//   double response = doenv ? fShrinkage*exp(dh) : fShrinkage*dh;
  
//   if (doenv) {
//     printf("sumtgt = %5e, sumtgt2 = %5e, dh = %5e, n = %i\n",sumtgt,sumtgt2,dh,int(evts.size()));
//   }
  
//   double response = doenv ? exp(dh) : fShrinkage*dh;

  double responsemin = 0.;
  double response3 = 0.;
  
//   for (std::vector<MCGBREvent*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
//     if (doenv) {
//       (*it)->SetTarget((*it)->Target()+response);
//     }
//     else {
//       (*it)->SetTargetMin((*it)->TargetMin()+response);
//     }
// //     (*it)->SetTargetMin((*it)->TargetMin()+responsemin);
// //     (*it)->SetTarget3((*it)->Target3()+response3);
//   }
  
  tree.Responses().push_back(response);
  tree.ResponsesMin().push_back(responsemin);
  tree.Responses3().push_back(response3);
  
  
  tree.Limits().push_back(limits);
  
}

