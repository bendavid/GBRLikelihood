/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: HGGRooPdfs.h,v 1.1 2012/02/10 15:10:48 gpetrucc Exp $
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
#ifndef GBRLikelihood_MCGBRIntegrator
#define GBRLikelihood_MCGBRIntegrator

#include "TMatrixDSym.h"
#include "../interface/MCGBREvent.h"
#include "../interface/MCGBRForest.h"
#include "../interface/MCGBRTreeD.h"


class MCGBRIntegrator : public TNamed {
public:
  
  MCGBRIntegrator(const char *name, const char *title, int nevents, int neventsinitial);

  ~MCGBRIntegrator();

  
  void AddInputVar(std::string var)    { fInputVars.push_back(var); }
  void SetTargetVar(std::string var)   { fTargetVar = var;          }
  void SetMinEvents(int n)             { fMinEvents = n;            }
  void SetNEvents(int n)             { fNEvents = n;            }
  void SetShrinkage(double x)          { fShrinkage = x;            }
  void SetMinCutSignificance(double x); //  { fMinCutSignificance = x*x/2.0; }
  void SetMinWeights(const std::vector<double> &minweights) { fMinWeights = minweights; }
  void SetMinWeightTotal(double x) { fMinWeightTotal = x; }
  void SetMaxDepth(int depth) { fMaxDepth = depth; } 
  void SetMaxNodes(int max) { fMaxNodes = max; }
  void SetDoEnvelope(bool b) { fDoEnvelope = b; }
  void SetStagedGeneration(bool b) { fStagedGeneration = b; }
  void SetNEventsBagged(int n) { fNEventsBagged = n; }
//   void SetNEventsInitial(int n) { fNEventsInitial = n; }
 
  void TrainForest(int ntrees, bool reuseforest = false);  
  
  const MCGBRForest *Forest() const { return fForest; }
  const MCGBRForest *ForestGen() const { return fForestGen; }
  const MCGBRTreeD *GenTree() const { return fGenTree; }
  
  double Camel(int nDim, double *Xarg) const;
  double Camelrnd(int nDim, double *Xarg) const;

  
protected:

//   double NLSkewGaus(double x, double mu) const;
  inline double NLSkewGausDMu(double x, double envelope) const;
  inline double NLSkewGausD2Mu(double x, double envelope) const;
  
  inline double NLNormDMu(double x, double mu, double sigma) const;
  inline double NLNormD2Mu(double x, double mu, double sigma) const;

  inline double NLLogNormDMu(double x, double mu) const;
  inline double NLLogNormD2Mu(double x, double mu) const;  
  
  void BuildQuantiles(int nvars, double sumabsw);
  void FillDerivatives();
  
  void TrainTree(const std::vector<MCGBREvent*> &evts, double sumwtotal, MCGBRTreeD &tree, int depth, std::vector<std::pair<float,float> > limits, bool usetarget, bool doenv=false);      
  void BuildLeaf(const std::vector<MCGBREvent*> &evts, MCGBRTreeD &tree, const std::vector<std::pair<float,float> > &limits, bool doenv=false);
    

  int fNThreads;

  
  std::vector<double>       fTreeWeights;
  std::string               fTrainingCut;
  std::vector<std::string>  fInputVars;  
  std::string               fTargetVar;
  int                       fMinEvents;
  std::vector<double>       fMinWeights;
  double                    fMinWeightTotal;
  double                    fShrinkage;
  int                       fNTrees;
  int                       fNQuantiles;
  const unsigned int        fNBinsMax;
  double                    fMinCutSignificance;
  double                    fMinCutSignificanceMulti;
  int                       fNEvents;
  int                       fNEventsInitial;
  int                       fNEventsBagged;
  
  
  int                       fMaxDepth;
  int                       fMaxNodes;

  
  int                      fNVars;
  double                   fNLLVal;
  
  MCGBRForest             *fForest;
  MCGBRForest             *fForestGen;
  MCGBRTreeD              *fGenTree;
  std::vector<std::pair<float,float> > fLimits;
  
  
  double *_sepgains; 
  double *_sepgainsigs; 
  float *_cutvals;  
  int *_nlefts; 
  int *_nrights; 
  double *_sumwlefts; 
  double *_sumwrights;   
  double *_sumtgtlefts; 
  double *_sumtgtrights; 
  float *_leftvars; 
  float *_rightvars;       
  float *_fullvars; 
  int *_bestbins; 
  
  
  double **_ws; 
  double **_ws2; 
  double ***_wscls; 
  int **_ns; 
  int **_nsd;   
  double **_tgts; 
  double **_tgt2s;
  double **_tgt3s;
  double **_tgt4s;
  double **_tgt5s;
  
  double **_fmins;
  double **_tgtmaxs;
  double **_sumws; 
  double **_sumws2;       
  double ***_sumwscls; 
  int **_sumns; 
  int **_sumnsd;   
  double **_sumtgts; 
  double **_sumtgt2s; 
  double **_sumtgt3s; 
  double **_sumtgt4s; 
  double **_sumtgt5s; 
  
  double **_sumfmins; 
  double **_sumfminsr; 
  double **_sumtgtmaxs; 
  double **_sumtgtmaxsr;
  float **_varvals; 
  int **_varvalmaxs;
  int **_varvalmins;
  int **_sumvarvalmaxs;
  int **_sumvarvalmins;  
  double **_bsepgains; 
  double **_bsepgainsigs;       
  
  int **_quants; 
  int **_binquants; 
  
  int *_clss; 
  double *_tgtvals;
  double *_tgt2vals;
  double *_tgt3vals;
  double *_tgt4vals;
  double *_tgt5vals;
  
  double *_fvals;
  double *_weightvals;
  
  float **fQuantileMaps;
  float **fQuantileMapsMin;
  
  std::vector<int> sparserows;
  std::vector<int> sparsecols;
  std::vector<double> sparsedata;
  
  std::vector<double> fStepSizes;
  
  std::vector<MCGBREvent*> fEvts;
  std::vector<MCGBREvent*> fEvts2;
  
  double shrinkagefactor;
  double fForestIntegralNow;
  double fSigmaScale;
  double fShrinkageFactorSecondary;
  
  bool fDoEnvelope;
  bool fStagedGeneration;
  
  double fSigmaRatio;
  double fIntegralRatio;
  
  
};


#endif
