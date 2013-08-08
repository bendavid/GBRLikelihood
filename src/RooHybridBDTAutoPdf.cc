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
 
#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>

#include "../interface/RooHybridBDTAutoPdf.h"
#include "RooRealVar.h"
#include "RooAbsData.h"
#include "RooUniform.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooConstVar.h"
#include "RooDerivative.h"
#include "TFitter.h"
#include "TFitterMinuit.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TRandom.h"
#include "TMath.h"
#include <Math/QuantFuncMathCore.h>
#include <Math/ProbFunc.h>
#include "../interface/HybridGBRForest.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TDecompChol.h"
#include "TDecompBK.h"
#include "TDecompLU.h"
#include "TDecompSparse.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooMinimizer.h"
#include "RooAddition.h"
#include "TNtuple.h"
#include "TDecompQRH.h"
#include "TDecompSVD.h"
#include "RooCFunction1Binding.h"
#include "TH1D.h"
#include "omp.h"
#include "DataFormats/Math/interface/VDTMath.h"
#include <stdlib.h>
#include <malloc.h>
#include "../interface/GBRArrayUtils.h"

ClassImp(RooTreeConvert)

RooDataSet *RooTreeConvert::CreateDataSet(std::string name, TTree *tree, std::vector<std::string> vars, std::string weight) {
 
  //initialize TTreeFormulas to read variables from TTree
  std::vector<TTreeFormula*> inputforms;
  RooArgList roovars;
  
  for (std::vector<std::string>::const_iterator it = vars.begin(); 
      it != vars.end(); ++it) {
    inputforms.push_back(new TTreeFormula(it->c_str(),it->c_str(),tree));
    RooRealVar *roovar = new RooRealVar(it->c_str(),"",0.);
    roovar->setConstant(false);
    roovars.add(*roovar);
  }  
  
  RooArgSet dvars(roovars);
  
  TTreeFormula cutform(weight.c_str(),weight.c_str(),tree);  
  RooRealVar *weightvar = new RooRealVar(TString::Format("%s_weight",name.c_str()),"",1.);
  roovars.add(*weightvar);
  
  RooDataSet *dset = new RooDataSet(name.c_str(),"",roovars,RooFit::WeightVar(*weightvar));
  
  std::vector<double> minvals(vars.size(),std::numeric_limits<double>::max());
  std::vector<double> maxvals(vars.size(),-std::numeric_limits<double>::max());
  
  for (Long64_t iev=0; iev<tree->GetEntries(); ++iev) {
    if (iev%100000==0) printf("%i\n",int(iev));
    tree->LoadTree(iev);
    
    float weight = cutform.EvalInstance();
    
    if (weight==0.) continue; //skip events with 0 weight
    
    for (unsigned int ivar=0; ivar<vars.size(); ++ivar) {
      double val = inputforms[ivar]->EvalInstance();
      if (val<minvals[ivar]) minvals[ivar] = val;
      if (val>maxvals[ivar]) maxvals[ivar] = val;
      static_cast<RooRealVar*>(roovars.at(ivar))->setVal(val);
    }
    dset->add(dvars,weight);

  }

  for (std::vector<TTreeFormula*>::const_iterator it = inputforms.begin(); 
      it != inputforms.end(); ++it) {
    delete *it;
  }  
  
//   for (unsigned int ivar=0; ivar<vars.size(); ++ivar) {
//     printf("ivar = %i, min = %5f, max = %5f\n",ivar,minvals[ivar],maxvals[ivar]);
//     static_cast<RooRealVar*>(roovars.at(ivar))->setRange(minvals[ivar],maxvals[ivar]);
//   }
  
  
  return dset;
  
}

RooDataSet *RooTreeConvert::CreateDataSet(std::string name, TTree *tree, RooArgList &vars, RooRealVar &weight) {
 
  //initialize TTreeFormulas to read variables from TTree
  std::vector<TTreeFormula*> inputforms;
  
  for (int ivar=0; ivar<vars.getSize(); ++ivar) {
    RooRealVar *var = static_cast<RooRealVar*>(vars.at(ivar));
    inputforms.push_back(new TTreeFormula(var->GetTitle(),var->GetTitle(),tree));
    var->setConstant(false);
  }  
  
  RooArgList roovars(vars);
  
  TTreeFormula cutform(weight.GetTitle(),weight.GetTitle(),tree);  
  weight.setConstant(false);

  roovars.add(weight);
  
  RooDataSet *dset = new RooDataSet(name.c_str(),"",roovars,RooFit::WeightVar(weight));
  
  int currenttree = -1;
  for (Long64_t iev=0; iev<tree->GetEntries(); ++iev) {
    if (iev%100000==0) printf("%i\n",int(iev));
    tree->LoadTree(iev);
    int thistree = tree->GetTreeNumber();
    bool newtree = currenttree!=thistree;
    currenttree = thistree;    
    
    if (newtree) {
      cutform.Notify();
      for (int ivar=0; ivar<vars.getSize(); ++ivar) {
        inputforms[ivar]->Notify();      
      }
    }
    
    float weight = cutform.EvalInstance();
    
    if (weight==0.) continue; //skip events with 0 weight
    
    bool valid = true;
    for (int ivar=0; ivar<vars.getSize(); ++ivar) {
      double val = inputforms[ivar]->EvalInstance();
      RooRealVar *var = static_cast<RooRealVar*>(vars.at(ivar));
      if (val<var->getMin() || val>var->getMax()) {
        valid = false;
      }
      var->setVal(val);
    }
    if (!valid) continue;
    
    dset->add(vars,weight);

  }

  for (std::vector<TTreeFormula*>::const_iterator it = inputforms.begin(); 
      it != inputforms.end(); ++it) {
    delete *it;
  }  
  
  return dset;
  
}


ClassImp(RooRealConstraint)

RooRealConstraint::RooRealConstraint(const char *name, const char *title, RooAbsReal &real, double low, double high) :
  RooAbsReal(name,title),
  _real("real","",this,real),
  _low(low),
  _high(high),
  _scale(0.5*(_high-_low)),
  _offset(_low + 0.5*(_high-_low))
{

  RooGBRTarget *var = dynamic_cast<RooGBRTarget*>(&real);
  
  if (var) {
    double oldval = var->Var()->getVal();
    double newval = asin(2.0*(oldval-_low)/(_high-_low)-1.0);
    //double newval = atanh( (oldval - _offset)/_scale );
    //double newval = -log(_scale/(oldval-_low) - 1.0);
    //double newval = tan( (oldval - _offset)/_scale );
    var->Var()->setVal(newval);
    
    printf("oldval = %f, newval = %5f, evaluate = %5f\n",oldval,newval,evaluate());
  }
  
}
  
  
RooRealConstraint::RooRealConstraint(const RooRealConstraint& other, const char* name) :
  RooAbsReal(other,name),
  _real("real",this,other._real),
  _low(other._low),
  _high(other._high),
  _scale(other._scale),
  _offset(other._offset)
{
  
}

Double_t RooRealConstraint::evaluate() const
{
 
//   double hprd = (_high -_low);
//   double upbound = _low + hprd;
//   double lowbound = _low - hprd;
//   
//   double val = _real.arg().getVal();
//   
//   while (val>upbound) {
//     val -= 2.0*hprd;
//   }
//   while (val<lowbound) {
//     val += 2.0*hprd;
//   }
//   
//   return _low + std::abs(val-_low);
  
  //return (_low + std::abs((_real.arg().getVal() % (2.0*(_high-_low))) - _high));
  //return std::max(std::min(_real.arg().getVal(),_high),_low);
  //return std::max(std::min(_real.arg().getVal(),_high),_low);
  //return _low + 0.5*(_high-_low)*(sin(_real)+1.0);
  
  //return _offset + _scale*vdt::fast_sinf(_real);
  return _offset + _scale*vdt::fast_sin(_real);
  //return _offset + _scale*sin(_real);
  
  //return _offset + _scale*tanh(_real);
  //return _low + _scale/(1.0+exp(-_real));
  //return _offset + _scale*atan(_real);
  
}


ClassImp(RooPowerLaw)

RooPowerLaw::RooPowerLaw(const char *name, const char *title, RooAbsReal &x, RooAbsReal &p) :
  RooAbsPdf(name,title),
  _x("x","",this,x),
  _p("p","",this,p)
{

}
  
  
RooPowerLaw::RooPowerLaw(const RooPowerLaw& other, const char* name) :
  RooAbsPdf(other,name),
  _x("x",this,other._x),
  _p("p",this,other._p)
{
  
}


Double_t RooPowerLaw::evaluate() const
{
 
  //const RooArgSet* nset = _normSet ;  
  
  return pow(_x,_p);
    
}

Int_t RooPowerLaw::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
   
  if (matchArgs(allVars,analVars,_x)) return 1 ;
  //if (matchArgs(allVars,analVars,_p)) return 2 ;
  return 0;
  
}

Double_t RooPowerLaw::analyticalIntegral(Int_t code, const char* rangeName) const
{
  
  assert(code==1);
  
  double xmax = _x.max(rangeName);
  double xmin = _x.min(rangeName);
  
  double omp = 1.0 + _p;
  
  
  return  ((pow(xmax,omp)-pow(xmin,omp))/omp);
  
}







ClassImp(RooCondAddPdf)

RooCondAddPdf::RooCondAddPdf(const char *name, const char *title, const RooArgList &pdfs, const RooArgList &coeffs) :
  RooAbsPdf(name,title),
  _pdfs("pdfs","",this),
  _coeffs("coeffs","",this),
  _selfnorm(coeffs.getSize()==pdfs.getSize() ? false : true)
{
  _pdfs.add(pdfs);
  _coeffs.add(coeffs);
}
  
  
RooCondAddPdf::RooCondAddPdf(const RooCondAddPdf& other, const char* name) :
  RooAbsPdf(other,name),
  _pdfs("pdfs",this,other._pdfs),
  _coeffs("coeffs",this,other._coeffs),
  _selfnorm(other._selfnorm)
{
  
}

Double_t RooCondAddPdf::evaluate() const
{
 
  //const RooArgSet* nset = _normSet ;  
  
  
  double sumcoeff = 0.;
  double val = 0.;
  //double finalcoeff = 1.0;
  for (int i=0; i<_coeffs.getSize(); ++i) {
    double coval = static_cast<RooAbsReal*>(_coeffs.at(i))->getVal();
    val += coval*static_cast<RooAbsPdf*>(_pdfs.at(i))->getValV(_normSet);
    sumcoeff += coval;
    //finalcoeff -= coval;
  }
  
  if (_selfnorm) {
    val += (1.0-sumcoeff)*static_cast<RooAbsPdf*>(_pdfs.at(_pdfs.getSize()-1))->getValV(_normSet);
  }
  else {
    val /= sumcoeff;
  }
  
  return val;
  
}

Int_t RooCondAddPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
 
  int code = static_cast<RooAbsReal*>(_pdfs.at(0))->getAnalyticalIntegral(allVars,analVars,rangeName);
  for (int ipdf=1; ipdf<_pdfs.getSize(); ++ipdf) {
    int pdfcode = static_cast<RooAbsPdf*>(_pdfs.at(ipdf))->getAnalyticalIntegral(allVars,analVars,rangeName);
    if (pdfcode!=code) code = 0;
  }
  
  //printf("RooCondAddPdf code = %i\n",code);
  //return 1;
  return code;
  
}

Double_t RooCondAddPdf::analyticalIntegral(Int_t code, const char* rangeName) const
{
  
  return 1.0;
 
//   double finalcoeff = 1.0;
//   double integral=0.;
//   for (int i=0; i<_coeffs.getSize(); ++i) {
//     double coval = static_cast<RooAbsReal*>(_coeffs.at(i))->getVal();
//     integral += coval*static_cast<RooAbsReal*>(_pdfs.at(i))->analyticalIntegral(code,rangeName);
//     finalcoeff -= coval;
//   }
//   
//   integral += finalcoeff*static_cast<RooAbsReal*>(_pdfs.at(_pdfs.getSize()-1))->analyticalIntegral(code,rangeName);
//   
//   return integral;
  
}


ClassImp(RooPdfAddReal)

RooPdfAddReal::RooPdfAddReal(const char *name, const char *title, const RooArgList &pdfs, const RooArgList &coeffs) :
  RooAbsReal(name,title),
  _pdfs("pdfs","",this),
  _coeffs("coeffs","",this)
{
  _pdfs.add(pdfs);
  _coeffs.add(coeffs);
}
  
  
RooPdfAddReal::RooPdfAddReal(const RooPdfAddReal& other, const char* name) :
  RooAbsReal(other,name),
  _pdfs("pdfs",this,other._pdfs),
  _coeffs("coeffs",this,other._coeffs)
{
  
}

Double_t RooPdfAddReal::evaluate() const
{
 
  const RooArgSet *nset = _pdfs.nset();  
  
  //if (nset) nset->Print("V");
  
  double val = 0.;
  for (int i=0; i<_coeffs.getSize(); ++i) {
    double coval = static_cast<RooAbsReal*>(_coeffs.at(i))->getVal();
    val += coval*static_cast<RooAbsPdf*>(_pdfs.at(i))->getValV(nset);
  }
    
  //printf("val = %5f\n",val);
    
  return val;
  
}

Int_t RooPdfAddReal::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
 
  int code = static_cast<RooAbsReal*>(_pdfs.at(0))->getAnalyticalIntegral(allVars,analVars,rangeName);
  for (int ipdf=1; ipdf<_pdfs.getSize(); ++ipdf) {
    int pdfcode = static_cast<RooAbsPdf*>(_pdfs.at(ipdf))->getAnalyticalIntegral(allVars,analVars,rangeName);
    if (pdfcode!=code) code = 0;
  }
  
  //printf("RooPdfAddReal code = %i\n",code);
  //return 1;
  return code;
  
}

Double_t RooPdfAddReal::analyticalIntegral(Int_t code, const char* rangeName) const
{
  
  //return 1.0;
 
  double integral=0.;
  for (int i=0; i<_coeffs.getSize(); ++i) {
    double coval = static_cast<RooAbsReal*>(_coeffs.at(i))->getVal();
    integral += coval*static_cast<RooAbsPdf*>(_pdfs.at(i))->analyticalIntegral(code,rangeName);
  }
    
  return integral;
  
}


ClassImp(RooGBRFunction)

//_____________________________________________________________________________
RooGBRFunction::RooGBRFunction(const char *name, const char *title, const RooArgList &vars, int ntargets) :
  RooAbsReal(name,title),
  _vars("vars","",this),
  _forest(new HybridGBRForest(ntargets)),
  _eval(vars.getSize())
{
  _vars.add(vars);
}

//_____________________________________________________________________________
RooGBRFunction::RooGBRFunction(const RooGBRFunction& other, const char* name) :
  RooAbsReal(other,name), 
  _vars("vars",this,other._vars),  
  _forest(new HybridGBRForest(*other._forest)),
  _eval(other._eval)
{
  
}
  
//_____________________________________________________________________________  
RooGBRFunction::~RooGBRFunction()
{
  delete _forest;
}
  
//_____________________________________________________________________________  
float RooGBRFunction::GetResponse(int itgt) const {
  
  //printf("RooGBRFunction::GetResponse(%i)\n",itgt);
//   if (isValueDirtyAndClear()) {
//   //if (1) {
//     //printf("recomputing bdt response\n");
//     for (int ivar=0; ivar<_vars.getSize(); ++ivar) {
//       _eval[ivar] = static_cast<RooAbsReal*>(_vars.at(ivar))->getVal();
//       //printf("ivar = %i, var = %5f\n",ivar,_eval[ivar]);
//     }
//     _forest->GetResponse(&_eval[0]);
//   }
//   
//   //printf("response %i = %5f\n",itgt,_forest->GetResponse(itgt));
//   return _forest->GetResponse(itgt);
  
  
  for (int ivar=0; ivar<_vars.getSize(); ++ivar) {
    _eval[ivar] = static_cast<RooAbsReal*>(_vars.at(ivar))->getVal();
    //printf("ivar = %i, var = %5f\n",ivar,_eval[ivar]);
  }  
  return _forest->GetResponse(&_eval[0],itgt);

  
}

//_____________________________________________________________________________  
void RooGBRFunction::SetForest(HybridGBRForest *forest) {
 
  if (_forest) delete _forest;
  _forest = forest;
  setValueDirty();
  
}

ClassImp(RooGBRTarget)

//_____________________________________________________________________________
RooGBRTarget::RooGBRTarget(const char *name, const char *title, RooGBRFunction &func, int itgt, RooRealVar &var) :
  RooAbsReal(name,title),
  _func("func","",this,func,true,false),
  _itgt(itgt),
  _var("var","",this,var),
  _usefunc(false)
{
  Var()->setConstant(_usefunc);
  //unRegisterProxy(_func);
  
}

//_____________________________________________________________________________
RooGBRTarget::RooGBRTarget(const RooGBRTarget& other, const char* name) :
  RooAbsReal(other,name), 
  _func("func",this,other._func),  
  _itgt(other._itgt),
  _var("var",this,other._var),
  _usefunc(other._usefunc)  
{
  
//   if (_usefunc) {
//     unRegisterProxy(_var);
//   }
//   else {
//     unRegisterProxy(_func);
//   }
//   
//   setValueDirty();
//   setShapeDirty();
  
}

//_____________________________________________________________________________
void RooGBRTarget::SetUseFunc(bool b)
{
  
  if (_usefunc==b) return;
  
  _usefunc = b;
  Var()->setConstant(_usefunc);
  
//   if (_usefunc) {
//     registerProxy(_func);
//     unRegisterProxy(_var);
//   }
//   else {
//     unRegisterProxy(_func);
//     registerProxy(_var);
//   }  
  
  setValueDirty();
  setShapeDirty();
  
}

//ClassImp(RooHybridBDTAutoPdf) 

RooHybridBDTAutoPdf *gHybridBDTAutoPointer;


//_____________________________________________________________________________
RooHybridBDTAutoPdf::RooHybridBDTAutoPdf(const char *name, const char *title, RooGBRFunction &func, const RooArgList &tgtvars, RooAbsReal &n0, RooRealVar &r, const std::vector<RooAbsData*> &data, const std::vector<RooAbsReal*> &pdfs) :
  TNamed(name,title),
  fTgtVars(tgtvars),
  fFunc(&func),  
  //fCondVars(fFunc->Vars()),
  fResTree(0),
  fPdfs(pdfs),
  fData(data),
  fNThreads(std::max(1,omp_get_max_threads())),  
  fLambdaVal(1.0), 
  fExternal(&n0),  
  fR(&r),  
  fN0Obs(data.front()->sumEntries()),  
  fLambda(0),  
  fDrvGraph(0),
  fDrvGraphSmooth(0),
  fGraphDelta(0),
  //fLambda(new RooRealVar("lambda","",0.)),
  fMinEvents(-99),  
  fMinWeights(std::vector<double>(data.size(),1000.)),
  fShrinkage(0.5),
  fNTrees(20),
  fNQuantiles(std::numeric_limits<unsigned short>::max()+1),
  //fNQuantiles(128),
  fNBinsMax(128),
  //fNBinsMax(fNQuantiles),
  fTransitionQuantile(0.7),
  fMinCutSignificance(-99.),
  fMinCutSignificanceMulti(-99.),
  fMaxNSpurious(-99.),
  fSumWTimesNVars(0.),
  fMaxDepth(-1),
  fMaxNodes(-1),
  fNTargets(tgtvars.getSize()),  
  fPrescaleInit(-1),
  _sepgains(0),
  _ws(0)  
{
  
  //printf("fN0 = %5f\n",fN0->getVal());
  
  //create constraint term for profile likelihood gradient scan and combine with external likelihood
  fConstraintVal = new RooRealVar(TString::Format("%s_constraintval",GetName()),"",fR->getVal());
  fConstraintCoeff = new RooRealVar(TString::Format("%s_constraintcoeff",GetName()),"",0.);
  RooFormulaVar *constraint = new RooFormulaVar(TString::Format("%s_constraint",GetName()),"","@0*pow(@1-@2,2)",RooArgList(*fConstraintCoeff,*fR,*fConstraintVal));
  
  fN0 = new RooAddition(TString::Format("%s_fullexternal",GetName()),"",RooArgList(*fExternal,*constraint));
  
  fGarbageCollection.addOwned(*fConstraintVal);
  fGarbageCollection.addOwned(*fConstraintCoeff);
  fGarbageCollection.addOwned(*constraint);
  fGarbageCollection.addOwned(*fN0);
  
  
  
  
  for (int ivar=0; ivar<fFunc->Vars().getSize(); ++ivar) {
    fCondVars.add(*fFunc->Vars().at(ivar));
  }
  
  //fCondVars.Print("V");
  
  printf("filling observables\n");
  //fill observables for paramtric pdfs  
  RooArgSet *allvarss = fPdfs.front()->getObservables(*data.front());
  {
  RooArgList allvars(*allvarss);
  for (int ivar=0; ivar<allvars.getSize(); ++ivar) {
    if (!fParmVars.contains(*allvars.at(ivar)) && !fCondVars.contains(*allvars.at(ivar))) fParmVars.add(*allvars.at(ivar));
  }
  }
  delete allvarss;
  
  //fParmVars.Print("V");
  

  
  printf("creating shadow pdfs\n");
  //make static variables shadowing dynamic targets along with corresponding pdf clones
  for (int ivar=0; ivar<fTgtVars.getSize(); ++ivar) {
    //double initval = static_cast<RooAbsReal*>(fTgtVars.at(ivar))->getVal();
    //printf("ivar = %i, initval = %5f\n",ivar,initval);
    //RooRealVar *var = new RooRealVar(fTgtVars.at(ivar)->GetName(),"",initval);
    //RooRealVar *var = new RooRealVar(fTgtVars.at(ivar)->GetName(),"",1.0);
    static_cast<RooGBRTarget*>(fTgtVars.at(ivar))->SetUseFunc(false);
    RooRealVar *var = static_cast<RooGBRTarget*>(fTgtVars.at(ivar))->Var();
    //var->setConstant(false);
    //var->removeRange();
    //fStaticTgts.addOwned(*var);
    fStaticTgts.add(*var);
  }
  
  for (unsigned int ipdf=0; ipdf<fPdfs.size(); ++ipdf) {
    //fStaticPdfs.addOwned(*fPdfs[ipdf]->cloneTree());
    fStaticPdfs.add(*fPdfs[ipdf]);
  }  
  
//   for (int ipdf=0; ipdf<fStaticPdfs.getSize(); ++ipdf) {  
//     fStaticPdfs.at(ipdf)->recursiveRedirectServers(fStaticTgts);
//   }
  
  //fStaticTgts.Print("V");
  
  //targets corresponding to log likelihood ratios
  fLLRTargets.push_back(0);
  for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {  
    fLLRTargets.push_back(static_cast<RooAbsReal*>(fStaticTgts.at(ipdf-1)));
  }

  printf("first loop, count events\n");

  //fStaticPdfs.at(0)->Print("V");
  //fStaticPdfs.at(0)->getParameters(*data.front())->Print("V");
  
  
  int nev = 0;
  for (unsigned int idata=0; idata<data.size(); ++idata) {
    nev += data[idata]->numEntries();
    printf("idata = %i, sumentries = %5f, numentries = %i\n",idata,data[idata]->sumEntries(),data[idata]->numEntries());
  }
  
  int nvars = fCondVars.getSize() + fParmVars.getSize();
  
  
  //set up computation of first and second derivatives
//   fDerivatives.resize(fFullFuncs.getSize(),std::vector<RooAbsReal*>(nparms));
//   f2Derivatives.resize(fFullFuncs.getSize(),std::vector<std::vector<RooAbsReal*> >(nparms,std::vector<RooAbsReal*>(nparms)));
//   for (int ipdf=0; ipdf<fFullFuncs.getSize(); ++ipdf) {
//     for (int iparm=0; iparm<nparms; ++iparm) {
//       if (fFullFuncs.at(ipdf)->overlaps(*fFullParms.at(iparm))) {
// 	fDerivatives[ipdf][iparm] = static_cast<RooAbsReal*>(fFullFuncs.at(ipdf))->derivative(*static_cast<RooRealVar*>(fFullParms.at(iparm)),fParmVars,1);
//       }
//       else
// 	fDerivatives[ipdf][iparm] = new RooConstVar("constzero","",0.);
//     }
//     for (int iparm=0; iparm<nparms; ++iparm) {
//       for (int jparm=0; jparm<nparms; ++jparm) {
// 	if (fDerivatives[ipdf][iparm]->overlaps(*fFullParms.at(jparm))) {
// 	  f2Derivatives[ipdf][iparm][jparm] = fDerivatives[ipdf][iparm]->derivative(*static_cast<RooRealVar*>(fFullParms.at(jparm)),fParmVars,1);
// 	}
// 	else
// 	  f2Derivatives[ipdf][iparm][jparm] = new RooConstVar("constzero","",0.);
//       }
//     }    
//   }  
  


   
  
//first loop here  
  
  printf("nev = %i, nvar = %i\n",int(nev),nvars);

  int ncls = fData.size();
  
//  fEvalVector.resize(fCondVars.getSize());
  
 
   
  //initialize arrays (ensure 32 byte alignment for avx vector instructions)

  _sepgains = (float*)memalign(32, nvars*sizeof(float));
  _sepgainsigs = (float*)memalign(32, nvars*sizeof(float));
  _cutvals = (float*)memalign(32, nvars*sizeof(float));
  _nlefts = (int*)memalign(32, nvars*sizeof(int));
  _nrights = (int*)memalign(32, nvars*sizeof(int));
  _sumwlefts = (float*)memalign(32, nvars*sizeof(float));
  _sumwrights = (float*)memalign(32, nvars*sizeof(float));
  _sumtgtlefts = (float*)memalign(32, nvars*sizeof(float));
  _sumtgtrights = (float*)memalign(32, nvars*sizeof(float));
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
  _sumws = new double*[nvars];
  _sumws2 = new double*[nvars];
  _sumwscls = new double**[nvars];
  _sumns = new int*[nvars];
  _sumtgts = new double*[nvars];  
  _sumtgt2s = new double*[nvars];
  _varvals = new float*[nvars];    
  _bsepgains = new float*[nvars];
  _bsepgainsigs = new float*[nvars];
  
  _binquants  = new int*[nvars];  
  _quants  = new int*[nvars];  
  
  _clss  = (int*)memalign(32, nev*sizeof(int));
  _tgtvals  = (double*)memalign(32, nev*sizeof(double));
  _tgt2vals  = (double*)memalign(32, nev*sizeof(double));
  _weightvals  = (double*)memalign(32, nev*sizeof(double));
  
  fQuantileMaps = new float*[nvars];  
  
  
  for (int ivar=0; ivar<nvars; ++ivar) {
    _ws[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _ws2[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _ns[ivar] = (int*)memalign(32, fNBinsMax*sizeof(int));
    _tgts[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _tgt2s[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));  
    _sumws[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumws2[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumns[ivar] = (int*)memalign(32, fNBinsMax*sizeof(int));
    _sumtgts[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _sumtgt2s[ivar] = (double*)memalign(32, fNBinsMax*sizeof(double));
    _varvals[ivar] = (float*)memalign(32, fNBinsMax*sizeof(float));
    _bsepgains[ivar] = (float*)memalign(32, fNBinsMax*sizeof(float));
    _bsepgainsigs[ivar] = (float*)memalign(32, fNBinsMax*sizeof(float));
    
    _wscls[ivar] = new double*[ncls];
    _sumwscls[ivar] = new double*[ncls];    
    
    _binquants[ivar] = (int*)memalign(32, fNBinsMax*sizeof(int));
    _quants[ivar] = (int*)memalign(32, nev*sizeof(int));
    
    fQuantileMaps[ivar] = (float*)memalign(32, fNQuantiles*sizeof(float));
    

    for (int icls=0; icls<ncls; ++icls) {
      _wscls[ivar][icls] = (double*)memalign(32, fNBinsMax*sizeof(double));
      _sumwscls[ivar][icls] = (double*)memalign(32, fNBinsMax*sizeof(double));
    }
    
  }
  
/*  _sepgains = new float[nvars];
  _sepgainsigs = new float[nvars];
  _cutvals = new float[nvars];
  _nlefts = new int[nvars];
  _nrights = new int[nvars];
  _sumwlefts = new float[nvars];
  _sumwrights = new float[nvars];
  _sumtgtlefts = new float[nvars];
  _sumtgtrights = new float[nvars];
  _leftvars = new float[nvars];
  _rightvars = new float[nvars];
  _fullvars = new float[nvars];
  _bestbins = new int[nvars];

 
  
  
  
  _ws = new double*[nvars];
  _ws2 = new double*[nvars];
  _wscls = new double**[nvars];
  _ns = new int*[nvars];
  _nsd = new int*[nvars];
  _tgts = new double*[nvars];
  _tgt2s = new double*[nvars];
  _sumws = new double*[nvars];
  _sumws2 = new double*[nvars];
  _sumwscls = new double**[nvars];
  _sumns = new int*[nvars];
  _sumtgts = new double*[nvars];
  _sumtgt2s = new double*[nvars];
  _varvals = new float*[nvars];
  _bsepgains = new float*[nvars];
  _bsepgainsigs = new float*[nvars];
  
  _quants = new int*[nvars];
  _bins = new int*[nvars];
  _clss = new int[nvars];
  _tgtvals  = new double[nev];
  _tgt2vals  = new double[nev];
  _weightvals  = new double[nev];
  
  fQuantileMaps = new float*[nvars];
  
  
  for (int ivar=0; ivar<nvars; ++ivar) {
    _ws[ivar] = new double[fNBinsMax];
    _ws2[ivar] = new double[fNBinsMax];
    _ns[ivar] = new int[fNBinsMax];
    _tgts[ivar] = new double[fNBinsMax];
    _tgt2s[ivar] = new double[fNBinsMax];
    _sumws[ivar] = new double[fNBinsMax];
    _sumws2[ivar] = new double[fNBinsMax];
    _sumns[ivar] = new int[fNBinsMax];
    _sumtgts[ivar] = new double[fNBinsMax];
    _sumtgt2s[ivar] = new double[fNBinsMax];
    _varvals[ivar] = new float[fNBinsMax];
    _bsepgains[ivar] = new float[fNBinsMax];
    _bsepgainsigs[ivar] = new float[fNBinsMax];
    
    _wscls[ivar] = new double*[ncls];
    _sumwscls[ivar] = new double*[ncls];
    
    _quants[ivar] = new int[nev];
    _bins[ivar] = new int[fNBinsMax];
    
    
    fQuantileMaps[ivar] = new float[fNQuantiles];
    
    for (int icls=0; icls<ncls; ++icls) {
      _wscls[ivar][icls] = new double[fNBinsMax];
      _sumwscls[ivar][icls] = new double[fNBinsMax];
    }    
    
  }*/   
  
  
  
      


  //TrainForest(fNTrees);
//  TrainForest(fNTrees);
  
  
  fLambdaVal = 1.;
//   int nparms = fExtVars.getSize() + fNTargets*nterm; 
//   fHessian = TMatrixDSym(nparms);   
//   for (int iel=0; iel<nparms; ++iel) {
//     fHessian(iel,iel) = 1.0;
//   }
  

  
//   fExtVars.removeAll();
//   fFullParms.removeAll();
//   fFullFuncs.removeAll();
//   fOuterIndices.clear();
//   fIndices.clear();
  
  //fill list of global parameters
  std::vector<RooAbsReal*> sources;
  for (unsigned int ipdf=0; ipdf<fPdfs.size(); ++ipdf) {
    sources.push_back(fPdfs.at(ipdf));
  }
  sources.push_back(fN0);
  
  
  //printf("filling parameters\n");
  for (unsigned int isrc=0; isrc<sources.size(); ++isrc) {
    RooArgSet *allparmss = sources[isrc]->getParameters(*fData.front());
    {
    RooArgList allparms(*allparmss);
    for (int ivar=0; ivar<allparms.getSize(); ++ivar) {
      if (!fExtVars.contains(*allparms.at(ivar)) && !fStaticTgts.contains(*allparms.at(ivar)) && !allparms.at(ivar)->getAttribute("Constant")) fExtVars.add(*allparms.at(ivar));
    }
    }
    delete allparmss;
  }
  
  if (dynamic_cast<RooRealVar*>(fN0) && !fExtVars.contains(*fN0) && !fStaticTgts.contains(*fN0) && !static_cast<RooRealVar*>(fN0)->getAttribute("Constant")) fExtVars.add(*fN0);
  
  //printf("ExtVars:\n");
  //fExtVars.Print("V");  
  
  fFullParms.add(fExtVars);
  fFullParms.add(fStaticTgts);
//  int nparms = fFullParms.getSize();
  
  fFullFuncs.add(fStaticPdfs);
  fFullFuncs.add(*fN0);
  
  //fFullFuncs.Print("V");
  
  fOuterIndices.resize(fFullFuncs.getSize());
  fIndices.resize(fFullFuncs.getSize());
  
  for (int ipdf=0; ipdf<fFullFuncs.getSize(); ++ipdf) {
    for (int iparm=0; iparm<fFullParms.getSize(); ++iparm) {
      if (fFullFuncs.at(ipdf)->overlaps(*fFullParms.at(iparm))) {
	fOuterIndices[ipdf].push_back(iparm);
      }
    }
    for (unsigned int iidx=0; iidx<fOuterIndices[ipdf].size(); ++iidx) {
      for (unsigned int jidx=iidx; jidx<fOuterIndices[ipdf].size(); ++jidx) {
	fIndices[ipdf].insert(std::pair<int,int>(fOuterIndices[ipdf][iidx],fOuterIndices[ipdf][jidx]));
      }
    }
  }
  
  fCondVarsClones.resize(fNThreads);
  fParmVarsClones.resize(fNThreads);
  fStaticTgtsClones.resize(fNThreads);
  fStaticPdfsClones.resize(fNThreads);
  fFullParmsClones.resize(fNThreads);
  fExtVarsClones.resize(fNThreads);
  fParmSetClones.reserve(fNThreads);
  
  RooAddition fullpdfsum("fullpdfsum","",fFullFuncs);
  
  //fullpdfsum.getComponents()->Print("V");
  
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    printf("ithread = %i\n",ithread);
    RooAbsArg *clone = fullpdfsum.cloneTree();
    fClones.addOwned(*clone);
    
    RooArgSet *clonecomps = clone->getComponents();
    RooArgSet *clonevars = clone->getVariables();
    clonecomps->add(*clonevars);
    //clonecomps->Print("V");
    for (int ivar=0; ivar<fCondVars.getSize(); ++ivar) {
      fCondVarsClones[ithread].add(*clonecomps->find(fCondVars.at(ivar)->GetName()));
    }
    for (int ivar=0; ivar<fParmVars.getSize(); ++ivar) {
      fParmVarsClones[ithread].add(*clonecomps->find(fParmVars.at(ivar)->GetName()));
    }
    for (int ivar=0; ivar<fStaticTgts.getSize(); ++ivar) {
      fStaticTgtsClones[ithread].add(*clonecomps->find(fStaticTgts.at(ivar)->GetName()));
    }
    for (int ivar=0; ivar<fStaticPdfs.getSize(); ++ivar) {
      fStaticPdfsClones[ithread].add(*clonecomps->find(fStaticPdfs.at(ivar)->GetName()));
    }    
    for (int ivar=0; ivar<fFullParms.getSize(); ++ivar) {
      fFullParmsClones[ithread].add(*clonecomps->find(fFullParms.at(ivar)->GetName()));
    }   
    for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
      fExtVarsClones[ithread].add(*clonecomps->find(fExtVars.at(ivar)->GetName()));
    }           
    delete clonecomps;
    delete clonevars;
    
    fParmSetClones.push_back(RooArgSet(fParmVarsClones[ithread]));
  }
  
//   fCondVarsClones[0].Print("V");
//   fParmVarsClones[0].Print("V");
//   fStaticTgtsClones[0].Print("V");
//   fStaticPdfsClones[0].Print("V");
//   fFullParmsClones[0].Print("V");
  
  fEvts.reserve(nev);
  
  double sumw = 0.;
  
  printf("second loop, fill events in memory\n");
  //loop over trees to fill arrays and event vector
   
  //second loop here
  for (unsigned int idata=0; idata<data.size(); ++idata) {
    for (int iev=0; iev<data[idata]->numEntries(); ++iev) {
      const RooArgSet *ev = data[idata]->get(iev);
      fEvts.push_back(new HybridGBREvent(nvars,fNTargets,fFullParms.getSize()));
      HybridGBREvent *evt = fEvts.back();
      evt->SetWeight(data[idata]->weight());
      evt->SetClass(idata);

      sumw += evt->Weight();
      
      for (int ivar=0; ivar<fCondVars.getSize(); ++ivar) {
	evt->SetVar(ivar,static_cast<RooAbsReal*>(ev->find(fCondVars.at(ivar)->GetName()))->getVal());
      }
      for (int ivar=0; ivar<fParmVars.getSize(); ++ivar) {
	evt->SetVar(fCondVars.getSize() + ivar, static_cast<RooAbsReal*>(ev->find(fParmVars.at(ivar)->GetName()))->getVal());
      }    
    }
  }
  
  fSumWTimesNVars = sumw*fCondVars.getSize();
  
  printf("filled data\n");
    
  //int nevr = datau->numEntries();
  
  
  //map of input variable quantiles to values
  //fQuantileMaps.resize(nvars, std::vector<float>(fNQuantiles));
  
  //BuildQuantiles(nvars, sumw);
  BuildQuantiles(nvars, sumw);  
  
  
  
}

//_____________________________________________________________________________
void RooHybridBDTAutoPdf::SetMinCutSignificance(double x) {
  
  //fMinCutSignificance = TMath::ChisquareQuantile(TMath::Erf(x/sqrt(2)),1)/2.0;
  fMinCutSignificance = x*x/2.0;
  fMinCutSignificanceMulti = TMath::ChisquareQuantile(TMath::Erf(x/sqrt(2)),fNTargets)/2.0;
  
  
}

void RooHybridBDTAutoPdf::UpdateTargets(int nvars, double sumw, int itree) {
 
  //printf("UpdateTargets\n");
  
//  int tgtidx = itree%fNTargets;
  
  #pragma omp parallel for
  for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
    
    int ithread =  omp_get_thread_num();
    //int ithread =  0;
    
  //for (unsigned int iev=0; iev<1; ++iev) {
      
    //printf("iev = %i\n",iev);
//     for (int ivar=0; ivar<nvars; ++ivar) {     
//       double var = fEvts.at(iev)->Var(ivar);
//       fEvalVector[ivar] = var;      
//     }
//     fForest->GetResponse(&fEvalVector[0]);
//     
//     for (int itgt=0; itgt<fForest->NTargets(); ++itgt) {
//       fEvts.at(iev)->SetTarget(itgt,fForest->GetResponse(itgt));
//     }
    
    for (int ivar=0; ivar<fCondVars.getSize(); ++ivar) {
      static_cast<RooRealVar*>(fCondVarsClones[ithread].at(ivar))->setVal(fEvts.at(iev)->Var(ivar));
    }
    for (int ivar=0; ivar<fParmVars.getSize(); ++ivar) {
      static_cast<RooRealVar*>(fParmVarsClones[ithread].at(ivar))->setVal(fEvts.at(iev)->Var(fCondVarsClones[ithread].getSize() + ivar));
    }
    
    for (int itgt=0; itgt<fNTargets; ++itgt) {
      static_cast<RooRealVar*>(fStaticTgtsClones[ithread].at(itgt))->setVal(fEvts.at(iev)->Target(itgt));
    }
    
    
//     fStaticPdfs.at(0)->Print("V");
//     fStaticPdfs.at(0)->getComponents()->Print("V");
//     
//     fDerivatives[0][fExtVars.getSize()+0]->Print("V");
//     fDerivatives[0][fExtVars.getSize()+0]->getComponents()->Print("V");
    
    int evcls = fEvts.at(iev)->Class();
    double pdfval = fEvts.at(iev)->PdfVal();
    //double invpdf = 1.0/pdfval;
    double invpdf = vdt::fast_inv(pdfval);
    //float invpdf = vdt::fast_invf(pdfval);
    double invpdfsq = invpdf*invpdf;
    double weight = fEvts.at(iev)->Weight();
    
    //derivatives for pdfs
    //for (int itgt=0; itgt<fNTargets; ++itgt) {
      
    for (int itgt=0; itgt<fNTargets; ++itgt) {
      fEvts.at(iev)->SetTransTarget(itgt,0.);
      fEvts.at(iev)->SetTransTarget2(itgt,0.);
    }      
      
    //if (itree<fFullParms.getSize()) {
    //if (1) {
      
      
      for (unsigned int iidx=0; iidx<fOuterIndices[evcls].size(); ++iidx) {
	
	int ivar = fOuterIndices[evcls][iidx];
	int itgt = ivar - fExtVars.getSize();
	
	RooRealVar *var = static_cast<RooRealVar*>(fFullParmsClones[ithread].at(ivar));
        double startval = var->getVal();
        double step = 1e-3*var->getError();
        
        RooAbsReal *func = static_cast<RooAbsReal*>(fStaticPdfsClones[ithread].at(evcls));
        
        var->setVal(startval + step);
        double upval = func->getValV(&fParmSetClones[ithread]);
        
        var->setVal(startval - step);
        double downval = func->getValV(&fParmSetClones[ithread]);
        
        var->setVal(startval);
        
        double drvval = (upval-downval)*vdt::fast_inv(2.0*step);
        
        
	//double drvval = Derivative1Fast(static_cast<RooAbsReal*>(fStaticPdfsClones[ithread].at(evcls)),pdfval,static_cast<RooRealVar*>(fFullParmsClones[ithread].at(ivar)),&fParmSetClones[ithread],1e-3*static_cast<RooRealVar*>(fFullParmsClones[ithread].at(ivar))->getError());
	//double drvval = Derivative1(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls)),static_cast<RooRealVar*>(fFullParms.at(ivar)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(ivar))->getError());
	
	
//	double drvval = Derivative1(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls)),static_cast<RooRealVar*>(fFullParms.at(ivar)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(ivar))->getError());
	//double drvval = Derivative1Fast(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls)),pdfval,static_cast<RooRealVar*>(fFullParms.at(ivar)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(ivar))->getError());
	fEvts.at(iev)->SetDerivative(ivar,drvval);
	//fEvts.at(iev)->SetDerivative2(ivar,drv2val);
	
	if (itgt<0) continue;	
	
	//double drv2val = Derivative2Fast(static_cast<RooAbsReal*>(fStaticPdfsClones[ithread].at(evcls)),pdfval,static_cast<RooRealVar*>(fFullParmsClones[ithread].at(ivar)),&fParmSetClones[ithread],1e-3*static_cast<RooRealVar*>(fFullParmsClones[ithread].at(ivar))->getError());
	
        double drv2val = (upval + downval - 2.0*pdfval)*vdt::fast_inv(step*step);

	fEvts.at(iev)->SetDerivative2(ivar,drv2val);

//         for (unsigned int jidx=iidx+1; jidx<fOuterIndices[evcls].size(); ++jidx) {
//           int jvar = fOuterIndices[evcls][jidx];
//           double drv2ij = Derivative2Fast(static_cast<RooAbsReal*>(fStaticPdfsClones[ithread].at(evcls)),static_cast<RooRealVar*>(fFullParmsClones[ithread].at(ivar)),static_cast<RooRealVar*>(fFullParmsClones[ithread].at(jvar)),&fParmSetClones[ithread],1e-3*static_cast<RooRealVar*>(fFullParmsClones[ithread].at(ivar))->getError(),1e-3*static_cast<RooRealVar*>(fFullParmsClones[ithread].at(jvar))->getError());
//           
//           fEvts.at(iev)->ParmMatrix()(ivar,jvar) = drv2ij;
//           
//         }
        
       // if (itgt<0) continue;   


	
	//double drv2val = Derivative2(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls)),static_cast<RooRealVar*>(fFullParms.at(ivar)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(ivar))->getError());
	
	//fEvts.at(iev)->SetTransTarget(itgt,fDerivatives[evcls][fExtVars.getSize()+itgt]->getVal()/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset));
	
	//printf("evcls = %i, itgt = %i, drv = %5f, val = %5f\n",evcls,itgt, Derivative1(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls)),static_cast<RooRealVar*>(fStaticTgts.at(itgt)),&parmset),static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset));
	
  //      fEvts.at(iev)->SetTransTarget(itgt,(1.0+fLambda->getVal())*Derivative1(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls)),static_cast<RooRealVar*>(fStaticTgts.at(itgt)),&parmset,1e-3*static_cast<RooRealVar*>(fStaticTgts.at(itgt))->getError())/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset));
	
	fEvts.at(iev)->SetTransTarget(itgt,-weight*drvval*invpdf);
	fEvts.at(iev)->SetTransTarget2(itgt,-weight*drv2val*invpdf + weight*drvval*drvval*invpdfsq);

	
	//fEvts.at(iev)->SetTransTarget(itgt,1.0/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset));
      }
      
//       for (unsigned int iidx=0; iidx<fOuterIndices[evcls].size(); ++iidx) {
// 	
// 	int ivar = fOuterIndices[evcls][iidx];
// 	int itgt = ivar - fExtVars.getSize();
// 	
// 	if (itgt<0) continue;
// 	
// 	double drvi = fEvts.at(iev)->Derivative(ivar);
// 	
// 	fEvts.at(iev)->ValVector()(itgt) = -drvi*invpdf;
// 	
// 	for (unsigned int jidx=iidx; jidx<fOuterIndices[evcls].size(); ++jidx) {
// 	
// 	  int jvar = fOuterIndices[evcls][iidx];
// 	  int jtgt = jvar - fExtVars.getSize();
// 	  
// 	  if (jtgt<0) continue;
// 	  
//        	  double drvj = fEvts.at(iev)->Derivative(jvar);
// 
// 	  fEvts.at(iev)->ParmMatrix()(itgt,jtgt) = drvi*drvj*invpdfsq;
// 	}
// 	
//       }
//       
//       for (int iel=0; iel<fNTargets; ++iel) {
// 	for (int jel=0; jel<iel; ++jel) {
// 	  fEvts.at(iev)->ParmMatrix()(iel,jel) = fEvts.at(iev)->ParmMatrix()(jel,iel);
// 	}
//       }
      
      
    //}
//     else {
//       for (int iidx=0; iidx<fFullParms.getSize(); ++iidx) {
// 	valv(iidx) = pdfval;
// 	for (int jidx=0; jidx<fFullParms.getSize(); ++jidx) {
// 	  parmm(iidx,jidx) = static_cast<RooRealVar*>(fFullParms.at(jidx))->getVal();
// 	}
//       }
//       
//       parmm -= fEvts.at(iev)->ParmMatrix();
//       valv -= fEvts.at(iev)->ValVector();
//       
//       TDecompLU dbk(parmm);
//       dbk.Solve(valv);
//       
//       for (unsigned int iidx=0; iidx<fOuterIndices[evcls].size(); ++iidx) {
// 	
// 	int ivar = fOuterIndices[evcls][iidx];
// 	int itgt = ivar - fExtVars.getSize();
// 	
// 	
// 	double drvval = valv(ivar);
// 
// 	fEvts.at(iev)->SetDerivative(ivar,drvval);
// 	
// 	if (itgt<0) continue;
// 	
// 	
// 	fEvts.at(iev)->SetTransTarget(itgt,drvval/pdfval);
// 
// 	
//       }      
//       
//     }
    
/*    if (pointidx>=0) {
      fEvts.at(iev)->ValVector()(pointidx) = pdfval;
      for (int iparm=0; iparm<fFullParms.getSize(); ++iparm) {
	fEvts.at(iev)->ParmMatrix()(pointidx,iparm) = static_cast<RooRealVar*>(fFullParms.at(iparm))->getVal();
      }
    } */   


    
//     double normden = 1.0;
//     for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {
//       int itgt = ipdf - 1;
//       int ivar = fExtVars.getSize() + itgt;
//       double expF = exp(fEvts.at(iev)->Target(itgt));
//         
//       normden += expF;
//     }
//     double invnormden = 1.0/normden;
//     double invnormdensq = invnormden*invnormden;
// 
//     for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {
//       int itgt = ipdf - 1;
//       int ivar = fExtVars.getSize() + itgt;
//       double expF = exp(fEvts.at(iev)->Target(itgt));
//      
//       if (ipdf==evcls) {
//         fEvts.at(iev)->SetTransTarget(itgt,fEvts.at(iev)->TransTarget(itgt) - 1.0);
//       }
//       
//       fEvts.at(iev)->SetTransTarget(itgt,fEvts.at(iev)->TransTarget(itgt) + invnormden*expF);
//       fEvts.at(iev)->SetTransTarget2(itgt,fEvts.at(iev)->TransTarget2(itgt) + invnormden*expF - invnormdensq*expF*expF);
//     }


    //normalization terms
/*    if (evcls==0) {
      for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {
	int itgt = ipdf - 1;
	int ivar = fExtVars.getSize() + itgt;
	double nexpF = exp(-fEvts.at(iev)->Target(itgt));
// 	fEvts.at(iev)->SetTransTarget(itgt,fEvts.at(iev)->TransTarget(itgt) - nexpF/fN0Obs);
// 	fEvts.at(iev)->SetTransTarget2(itgt,fEvts.at(iev)->TransTarget2(itgt) + nexpF/fN0Obs);
	fEvts.at(iev)->SetTransTarget(itgt,fEvts.at(iev)->TransTarget(itgt) - nexpF);
	fEvts.at(iev)->SetTransTarget2(itgt,fEvts.at(iev)->TransTarget2(itgt) + nexpF);	
      }
    }
    else {
      int itgt = evcls - 1;
      int ivar = fExtVars.getSize() + itgt;
      fEvts.at(iev)->SetTransTarget(itgt,fEvts.at(iev)->TransTarget(itgt) + 1.0);
    } */   
    
    
//     for (int itgt=0; itgt<fNTargets; ++itgt) {
//       printf("itgt = %i, transtgt = %5f\n",itgt,fEvts.at(iev)->TransTarget(itgt));
//     }
    
//     //normalization terms
//     if (evcls==0) {
//       
//       
//       
//       dcoeffs[0] = fN0->getVal();
//       for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {
// 	int itgt = ipdf - 1;
//         dcoeffs[ipdf] = fR->getVal()*exp(-fEvts.at(iev)->Target(itgt));
// 	dcoeffs[0] -= dcoeffs[ipdf];
//       }
//       
//       double ntimesd = dcoeffs[0]*static_cast<RooAbsReal*>(fStaticPdfs.at(0))->getValV(&parmset));
//       for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {
// 	ntimesd + dcoeffs[ipdf]*static_cast<RooAbsReal*>(fStaticPdfs.at(ipdf))->getValV(&parmset));
//       }
//       
//       for (int itgt=fStaticPdfs.getSize()-1; itgt<fNTargets; ++itgt) {
// 	fEvts.at(iev)->SetTransTarget(itgt,0.);
// 	for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {
// 	  fEvts.at(iev)->SetTransTarget(itgt, fEvts.at(iev)->TransTarget(itgt) +  dcoeffs[ipdf]*fDerivatives[evcls][fExtVars.getSize()+itgt]->getVal()/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset));	
// 	}
//       }
// 
//       for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {
// 	int itgt = ipdf - 1;
// 	double nexpF = exp(-fEvts.at(iev)->Target(itgt));
// 	fEvts.at(iev)->SetTransTarget(itgt,fEvts.at(iev)->TransTarget(itgt) + nexpF/fN0Obs);
//       }
//     }
//     else {
//       int itgt = evcls - 1;
//       fEvts.at(iev)->SetTransTarget(itgt,fEvts.at(iev)->TransTarget(itgt) - 1.0);
//     }
    
    
//     if ( fEvts.at(iev)->Class()==0 ) {
//       //double ntimesd = ((ns/nd)*nexpF*sigpdf->getValV(&parmset)+ (1.0-(ns/nd)*nexpF)*bkgpdf->getValV(&parmset));
//        double ntimesd = ns*nexpF*sigpdf->getValV(&parmset)+ (nd-ns*nexpF)*bkgpdf->getValV(&parmset); 
// 
// 
//       //fEvts.at(iev)->SetTransTarget(0,(-ns*nexpF*fSigPdf->getValV(&parmset)+ns*nexpF*fBkgPdf->getValV(&parmset))/ntimesd + nexpF);
//       //fEvts.at(iev)->SetTransTarget(0,(-(ns)*nexpF*sigpdf->getValV(&parmset)+(ns)*nexpF*bkgpdf->getValV(&parmset))/ntimesd + (nd/fNDataObs)*nexpF);
//       fEvts.at(iev)->SetTransTarget(0,(-ns*nexpF*sigpdf->getValV(&parmset) + ns*expF*bkgpdf->getValV(&parmset))/ntimesd + (1.0/fNDataObs)*nexpF);
//       //fEvts.at(iev)->SetTransTarget(0,(-(ns)*nexpF*sigpdf->getValV(&parmset)+(ns)*nexpF*bkgpdf->getValV(&parmset))/ntimesd);
//       fEvts.at(iev)->SetTransTarget(1,ns*nexpF*dsig->getVal()/ntimesd);
//       fEvts.at(iev)->SetTransTarget(2,(nd-ns*nexpF)*dbkg->getVal()/ntimesd);
//       //fEvts.at(iev)->SetTransTarget(3,-1.0);      
//     }
//     else if ( fEvts.at(iev)->Class()==1 ) {
//       fEvts.at(iev)->SetTransTarget(0,-1.0);
//       fEvts.at(iev)->SetTransTarget(1,dsig->getVal()/sigpdf->getValV(&parmset));
//       fEvts.at(iev)->SetTransTarget(2, 0.0);
//       //fEvts.at(iev)->SetTransTarget(3,-1.0);      
//     }
  

    
  }  
  
  
}

void RooHybridBDTAutoPdf::BuildQuantiles(int nvars, double sumw) {
 
  //parallelize building of quantiles for each input variable
  //(sorting of event pointer vector is cpu-intensive)
  #pragma omp parallel for
  for (int ivar=0; ivar<nvars; ++ivar) {
    printf("sorting var %i\n",ivar);
        
    
    std::map<int,float,std::greater<float> > tmpmap;    
    std::vector<HybridGBREvent*> evtsvarsort(fEvts.begin(),fEvts.end());
    
    //std::sort(evtsvarsort.begin(),evtsvarsort.end(),GBRVarCMP(ivar));
    std::sort(evtsvarsort.begin(),evtsvarsort.end(),GBRVarCMP(ivar));
    
    double sumwq = 0;
    for (unsigned int iev=0; iev<evtsvarsort.size(); ++iev) {
      sumwq += evtsvarsort[iev]->Weight();
      int quant = int((sumwq/sumw)*(fNQuantiles-1));
      float val = evtsvarsort[iev]->Var(ivar);
    
      //ensure that events with numerically identical values receive the same quantile
      if (iev>0 && val==evtsvarsort[iev-1]->Var(ivar)) quant = evtsvarsort[iev-1]->Quantile(ivar);
    
      evtsvarsort[iev]->SetQuantile(ivar,quant);
    
      tmpmap[quant] = val;
    
    }
    

    for (int i=0; i<fNQuantiles; ++i) {
      std::map<int,float,std::greater<float> >::const_iterator mit = tmpmap.lower_bound(i);
      
      float val;
      if (mit!=tmpmap.end()) val = mit->second;
      else val = -std::numeric_limits<float>::max();
      
      fQuantileMaps[ivar][i] = val;
      
      
    }
    
    
    
  }    
  
}


const HybridGBRForest *RooHybridBDTAutoPdf::TrainForest(int ntrees, bool reuseforest) {

  for (int ithread=0; ithread<fNThreads; ++ithread) {
    static_cast<RooAbsReal*>(fClones.at(ithread))->getValV(&fParmSetClones[ithread]);
  }
  
  for (int ivar=0; ivar<fTgtVars.getSize(); ++ivar) {
    static_cast<RooGBRTarget*>(fTgtVars.at(ivar))->SetUseFunc(false);  
  }
  fFunc->setOperMode(RooAbsArg::AClean);
  
  for (unsigned int iev=0; iev<fEvts.size(); ++iev) {
    for (int itgt=0; itgt<fNTargets; ++itgt) {
      fEvts.at(iev)->SetCurrentNode(itgt,0);
    }
  }
  
  if (fResTree) delete fResTree;
  
  fResTree = new TNtuple("restree","","iter:r:nd:nllval:dldr");  
  
  
  
  int nvars = fCondVars.getSize() + fParmVars.getSize();
  int nvarstrain = fCondVars.getSize();
  
  double sumw = 0.;
  double sumwd = 0.;
  double sumwu = 0.;
  int nev = fEvts.size();
  for (std::vector<HybridGBREvent*>::iterator it=fEvts.begin(); it!=fEvts.end(); ++it) {
    sumw += (*it)->Weight();
    if ( (*it)->Class()==0 ) sumwd += (*it)->Weight();
    else if ( (*it)->Class()==1 ) sumwu += (*it)->Weight();
  }
  
  RooArgSet cvarset(fCondVars);

  
  HybridGBRForest *forest  = 0;
  if (reuseforest) {
    forest = fFunc->Forest();
  }
  else {
    forest = new HybridGBRForest(fFunc->Forest()->NTargets());
  }
  
  
  printf("setting targets\n");
 
  //printf("setting cloned vals\n");
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    //printf("ithread = %i\n",ithread);
    for (int ivar=0; ivar<fFullParms.getSize(); ++ivar) {
      //printf("ivar = %i: %s\n",ivar, fFullParms.at(ivar)->GetName());
      RooRealVar *cloneparm = static_cast<RooRealVar*>(fFullParmsClones[ithread].at(ivar));
      cloneparm->setVal(static_cast<RooRealVar*>(fFullParms.at(ivar))->getVal());
      cloneparm->setError(static_cast<RooRealVar*>(fFullParms.at(ivar))->getError());
    }
  }      
  
  if (!reuseforest) {
    {  
      gHybridBDTAutoPointer = this;
      RooAddition fullparmsum("fullparmsum","",fFullParms);
      RooCFunction1Binding<double,double> nllfunc("nllfunc","", &EvalLossNull,fullparmsum);
      
      nllfunc.Print();
      
      RooMinimizer *minim = new RooMinimizer(nllfunc);
      minim->setErrorLevel(0.5);
      minim->setStrategy(0);
      minim->minimize("Minuit2","minimize");  
      delete minim;  
      
      //if (fConstraintCoeff->getVal()>0.) fR->setError(1e3*fR->getError());
    }
    
    
    for (int itgt = 0; itgt<fStaticTgts.getSize(); ++itgt) {
      int treetgt = static_cast<RooGBRTarget*>(fTgtVars.at(itgt))->Index();
      double initF = static_cast<RooRealVar*>(fStaticTgts.at(itgt))->getVal();
      forest->SetInitialResponse(treetgt,initF);
      for (std::vector<HybridGBREvent*>::iterator it=fEvts.begin(); it!=fEvts.end(); ++it) {
	(*it)->SetTarget(itgt,initF);
      }    
    }  
    
    for (int ithread=0; ithread<fNThreads; ++ithread) {
      for (int ivar=0; ivar<fFullParms.getSize(); ++ivar) {
	RooRealVar *cloneparm = static_cast<RooRealVar*>(fFullParmsClones[ithread].at(ivar));
	cloneparm->setVal(static_cast<RooRealVar*>(fFullParms.at(ivar))->getVal());
	cloneparm->setError(static_cast<RooRealVar*>(fFullParms.at(ivar))->getError());
      }
    }
  }
  
  
  fStepSizes.resize(fFullParms.getSize());  
  for (int iparm=0; iparm<fFullParms.getSize(); ++iparm) {
    fStepSizes[iparm] = 1e-3*static_cast<RooRealVar*>(fFullParms.at(iparm))->getError();
  }
  
  
  printf("nev = %i, sumw = %5f\n",int(nev), sumw);
  
  std::vector<std::pair<float,float> > limits;
  for (int ivar=0; ivar<nvarstrain; ++ivar) {
    //limits.push_back(std::pair<float,float>(-std::numeric_limits<float>::max(),std::numeric_limits<float>::max()));
    limits.push_back(std::pair<float,float>(-10.,10.));
  }
  
  RooArgSet parmset(fParmVars);

  
  TVectorD dparnull(fFullParms.getSize());
  fNLLVal = EvalLoss(forest,0.,dparnull);
  
  //loop over requested number of trees
  int nunittrees = 0;
  std::vector<double> nllvals;
  std::vector<double> dldrvals;
  for (int itree=0; itree<ntrees; ++itree) {
    printf("tree %i\n",itree);

    UpdateTargets(nvars,sumw, itree); 

    int maxtreesize = 0;
    int maxtreeidx = 0;
    
    std::vector<int> treesizes(fNTargets);

    for (int itgt=0; itgt<fNTargets; ++itgt) {
      int treetgt = static_cast<RooGBRTarget*>(fTgtVars.at(itgt))->Index();
      forest->Trees()[treetgt].push_back(HybridGBRTree()); 
    }
    
    for (int itgt=0; itgt<fNTargets; ++itgt) {
      int treetgt = static_cast<RooGBRTarget*>(fTgtVars.at(itgt))->Index();
      HybridGBRTree &tree = forest->Trees()[treetgt].back();      
      TrainTree(fEvts,sumw,tree,nvarstrain,0.,0,limits,itgt);
      int treesize = tree.Responses().size();
      treesizes[itgt] = treesize;
      //printf("itgt = %i, treesize = %i\n",itgt,int(tree.Responses().size()));
    }
    
    for (int itgt=0; itgt<fNTargets; ++itgt) {
      int treesize = treesizes[itgt];
      if (treesize>maxtreesize) {
	maxtreesize = treesize;
        maxtreeidx = itgt;
      }      
    }
    
    
    printf("maxtreesize = %i, maxtgtidx = %i\n",maxtreesize,maxtreeidx);

    double originalshrinkage = fShrinkage;
    if (maxtreesize==1) {
      ++nunittrees;
      //fShrinkage = std::max(originalshrinkage,0.9);
    }
    else {
      nunittrees = 0.;
    }
    
    FitResponses(forest);

    fShrinkage = originalshrinkage;
    
    nllvals.push_back(fNLLVal);
    
    double dldrval = -2.0*fConstraintCoeff->getVal()*(fR->getVal()-fConstraintVal->getVal());
    dldrvals.push_back(dldrval);
    
    fResTree->Fill(itree,fR->getVal(),fN0->getVal(),fNLLVal,fdLdR);
   
    int oldnllidx = nllvals.size() - 5 - 1;
    //if (oldnllidx>=0 && (fNLLVal - nllvals[oldnllidx])>(-2e-3) && std::abs(dldrval-dldrvals[oldnllidx])<2e-1 && nunittrees>10) {
    if (oldnllidx>=0 && (fNLLVal - nllvals[oldnllidx])>(-2e-3) && nunittrees>10) {      
      break;
    }
    
    oldnllidx = nllvals.size() - 3.0/fShrinkage - 1;
    //if (oldnllidx>=0 && (fNLLVal - nllvals[oldnllidx])>(-2e-3) && std::abs(dldrval-dldrvals[oldnllidx])<2e-1) {
    if (oldnllidx>=0 && (fNLLVal - nllvals[oldnllidx])>(-2e-3)) {
      printf("breaking\n");
      break;
    }

  }
  
  fFunc->setOperMode(RooAbsArg::Auto);
  
//  printf("fNorm = %5f\n",fNorm);
  //return fully trained HybridGBRForest
  
  if (reuseforest) {
    fFunc->setValueDirty();
  }
  else {
    fFunc->SetForest(forest);
  }

  for (int ivar=0; ivar<fTgtVars.getSize(); ++ivar) {
    static_cast<RooGBRTarget*>(fTgtVars.at(ivar))->SetUseFunc(true);  
  }  
  
  return forest;  
  
}
 
    
// void RooHybridBDTAutoPdf::vtest(const float **__restrict__ ra, const float **__restrict__ rb, float **__restrict__ rc) {
//     
//   for (int iat=0; iat<128; ++iat) {
// //     const float *__restrict__ sra = (const float*)__builtin_assume_aligned(ra[iat],32);
// //     const float *__restrict__ srb = (const float*)__builtin_assume_aligned(rb[iat],32);
// //     float *__restrict__ src = (float*)__builtin_assume_aligned(rc[iat],32);    
//     //rc[iat] = ra[iat] + rb[iat];
//     for (int jat=0; jat<128; ++jat) {
//       //src[jat] = sra[jat] + srb[jat];
//       rc[iat][jat] = ra[iat][jat] + rb[iat][jat];
//     }
//   }
//        
// }  
    
  
//_______________________________________________________________________
void RooHybridBDTAutoPdf::TrainTree(const std::vector<HybridGBREvent*> &evts, double sumwtotal, HybridGBRTree &tree, const int nvars, double transition, int depth, const std::vector<std::pair<float,float> > limits, int tgtidx) {
  
  
  //alignas(32) float floats[128];
  
  //float *aarrtest = new float[128];
  
  //index of current intermediate node
  int thisidx = tree.CutIndices().size();    
  
  //number of events input to node
  const int nev = evts.size();
  const int ncls = fStaticPdfs.getSize();
  
  //index of best cut variable
  int bestvar = 0;

  
//   double sumwd = 0;
//   double sumnp = 0;
//   int numprimary = 0;
//   //int numprimary = nev;
//   for (int iev = 0; iev<nev; ++iev) {
//     //if (evts[iev]->Class()==0 && evts[iev]->IsPrimary()) {
//     if (evts[iev]->Class()==0) {  
//       sumwd += evts[iev]->Weight();
//       sumnp += evts[iev]->Weight()*exp(-evts[iev]->Target());
//     }
//     if (evts[iev]->IsPrimary()) ++numprimary;
//   }
//   
//   //double avgp = sumnp/sumwd/fNorm;
//   double avgp = 1.0;
  
  
//   const float *a = (float*)__builtin_assume_aligned(memalign(32,128*sizeof(float)),32);
//   const float *b = (float*)__builtin_assume_aligned(memalign(32,128*sizeof(float)),32);
//   float *c = (float*)__builtin_assume_aligned(memalign(32,128*sizeof(float)),32);
//    
  
//   const float *__restrict__ ar = (float*)__builtin_assume_aligned(a,32);
//   const float *__restrict__ br = (float*)__builtin_assume_aligned(b,32);
//   float *__restrict__ cr = (float*)__builtin_assume_aligned(c,32);
  

  
//   std::vector<int> threadnums(nvars);
//   std::vector<int> nbinsvar(nvars);
//   std::vector<double> timesvar(nvars);
//   
  //float *__restrict *__restrict__ sumws = _sumws;  
  
//   printf("start split search loop\n");
//   TStopwatch clockloop;
//   clockloop.Start();
  
   
  //printf("first parallel loop\n");
  //printf("nev = %i\n",nev);
  //fill temporary array of quantiles (to allow auto-vectorization of later loops)
  #pragma omp parallel for
  for (int iev = 0; iev<nev; ++iev) {
    _clss[iev] = evts[iev]->Class();
    _tgtvals[iev] = evts[iev]->TransTarget(tgtidx);
    _tgt2vals[iev] = evts[iev]->TransTarget2(tgtidx);
    _weightvals[iev] = evts[iev]->Weight();
    for (int ivar=0; ivar<nvars; ++ivar) {
      _quants[ivar][iev] = evts[iev]->Quantile(ivar);   
      //printf("quant = %i\n",_quants[ivar][iev]);
    }
  }  
  
  //printf("second parallel loop\n");
  //trivial open-mp based multithreading of loop over input variables
  //The loop is thread safe since each iteration writes into its own
  //elements of the 2-d arrays
  #pragma omp parallel for schedule(dynamic,1)
  for (int ivar=0; ivar<nvars; ++ivar) {
             
    
    int minquant;
    int maxquant;
    
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
    


    
    //zero arrays where necessary and 
    GBRArrayUtils::InitArrays(_ns[ivar],_tgts[ivar],_tgt2s[ivar],_bsepgains[ivar],nbins);
    
    //printf("touch wscls\n");
    //the inner loop here should also vectorize
    for (int icls=0; icls<ncls; ++icls) {
      GBRArrayUtils::ZeroArray(_wscls[ivar][icls],nbins);
    }    
    //printf("done wscls\n");
    
    
    //compute map between bin numbers
    //and variable cut values
    
    //printf("quant manipulation\n");
    //this loop should auto-vectorize
    GBRArrayUtils::FillBinQuants(_binquants[ivar], offset, pscale,fNQuantiles, nbins);
    
    //this loop won't auto-vectorize because it's another gather operation (maybe in avx2 and gcc>4.7)
    for (unsigned int ibin=0; ibin<nbins; ++ibin) { 
      int quant = _binquants[ivar][ibin];
      _varvals[ivar][ibin] = fQuantileMaps[ivar][quant];
    }
    //printf("done quant manipulation\n");
    
    //printf("compute bins\n");
    
    //compute reduced bin value for each event using bit-shift operations
    //This loop should auto-vectorize in appropriate compiler/settings
//     for (int iev=0;iev<nev;++iev) {
//       _bins[ivar][iev] = (_quants[ivar][iev]-offset)>>pscale;
//     }

      
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
      //if (!evts[iev]->IsPrimary()) continue;
      
      int icls = _clss[iev];
      int ibin = (_quants[ivar][iev]-offset)>>pscale;
      
      //printf("icls = %i, ibin = %i, weight = %5f, transtarget = %5f, transtarget2 = %5f\n",icls,ibin, evts[iev]->Weight(),evts[iev]->TransTarget(tgtidx),evts[iev]->TransTarget2(tgtidx));
      
      //assert(ibin<int(nbins));
      
      ++_ns[ivar][ibin];
      
      //printf("incrementing wscls: ivar = %i, icls = %i, ibin = %i, quant = %i\n",ivar,icls,ibin,_quants[ivar][iev]);
      //printf("incrementing wscls\n");
      _wscls[ivar][icls][ibin] += _weightvals[iev];            
      //printf("done incrementing wscls\n");
      
      _tgts[ivar][ibin] += _tgtvals[iev];
      _tgt2s[ivar][ibin] += _tgt2vals[iev];
      
    } 
       

      
    //printf("starting split search\n");
 
    //printf("compute array integrals\n");
    //convert differential arrays to cumulative arrays by summing over
    //each element
    //loop cannot be vectorized because this is an iterative calculation
    _sumws[ivar][0] = _ws[ivar][0];
    _sumws2[ivar][0] = _ws2[ivar][0];
    for (int icls=0; icls<ncls; ++icls) {
      _sumwscls[ivar][icls][0] = _wscls[ivar][icls][0];
    }
    _sumns[ivar][0] = _ns[ivar][0];
    _sumtgts[ivar][0] = _tgts[ivar][0];
    _sumtgt2s[ivar][0] = _tgt2s[ivar][0];    
    
    for (unsigned int ibin=1; ibin<nbins; ++ibin) {      
      _sumns[ivar][ibin] = _sumns[ivar][ibin-1] + _ns[ivar][ibin];
      _sumtgts[ivar][ibin] = _sumtgts[ivar][ibin-1] + _tgts[ivar][ibin];  
      _sumtgt2s[ivar][ibin] = _sumtgt2s[ivar][ibin-1] + _tgt2s[ivar][ibin];  
    }
    
    //printf("cls array integrals\n");
    for (int icls=0; icls<ncls; ++icls) {
      for (unsigned int ibin=1; ibin<nbins; ++ibin) {  
	_sumwscls[ivar][icls][ibin] = _sumwscls[ivar][icls][ibin-1] + _wscls[ivar][icls][ibin];
      }
    }    
   // printf("done cls array integrals\n");
    
    
    
    
    //int n = sumns[ivar][nbins-1];
//     float sumw = _sumws[ivar][nbins-1];
//     float sumw2 = _sumws2[ivar][nbins-1];
//     std::vector<float> sumwscls(fStaticTgts.getSize());
//     for (int icls=0; icls<ncls; ++icls) {
//       sumwscls[icls] = _sumwscls[ivar][nbins-1][icls];
//     }
    
    //std::vector<float> sumtgt(fNTargets);
    //float sumtgtmag2 = 0.;

    const double sumtgt = _sumtgts[ivar][nbins-1];
    //sumtgtmag2 += _sumtgts[ivar][nbins-1]*_sumtgts[ivar][nbins-1];      
    
    const double sumtgt2 = _sumtgt2s[ivar][nbins-1];      
    //if (sumtgt2<=0.) continue;
    
    
    //weighted variance of target in full dataset
    //float fullvariance = sumtgt2 - sumtgtmag2/sumw;
    //float fullvariancevar = fullvariance*fullvariance/sumw2/sumw2;
    
    //_fullvars[ivar] = fullvariance;
    
   // printf("fullrms = %5f, sumtgt2 = %5f, sumtgt = %5f, sumw = %5f\n",fullrms,sumtgt2,sumtgt,sumw);
    
    //printf("short loop\n");
    float maxsepgain = 0.;
    float cutval = 0.;
    int nleft= 0;
    int nright = 0;
//    float sumwleft=0.;
  //  float sumwright=0.;
    int bestbin=0;
    
    //const double fulldiff = std::min(0.,-0.5*sumtgt*sumtgt*vdt::fast_inv(sumtgt2));
    const double fulldiff = std::min(0.,-0.5*sumtgt*sumtgt/sumtgt2);
    

    
    //const double fulldiff = -0.5*sumtgt*sumtgt/sumtgt2;
    
    //printf("start heavy loop\n");
    //loop over all bins and compute improvement in weighted variance of target for each split
    //This loop is relatively expensive and should auto-vectorize in the appropriate compiler/settings
    GBRArrayUtils::FillSepGains(_sumtgts[ivar], _sumtgt2s[ivar], _bsepgains[ivar], fulldiff, sumtgt, sumtgt2, nbins);
     
    //printf("start final loop\n");
    //loop over computed variance improvements and select best split, respecting also minimum number of events per node
    //This loop cannot auto-vectorize, at least in gcc 4.6x due to the mixed type conditional, but it's relatively fast
    //in any case
    for (unsigned int ibin=0; ibin<nbins; ++ibin) {   
      //printf("testing split, widthleft = %5f, widthright = %5f, sepgain = %5f\n",_varvals[ivar][ibin] - limits[ivar].first,limits[ivar].second-_varvals[ivar][ibin],_bsepgains[ivar][ibin]);
//       if (_sumns[ivar][ibin]>=fMinEvents && (nev-_sumns[ivar][ibin])>=fMinEvents && _bsepgains[ivar][ibin]>maxsepgain && _bsepgainsigs[ivar][ibin]>fMinCutSignificance) {
        //if (_sumws2[ivar][ibin]>=fMinEvents && (sumw2-_sumws2[ivar][ibin])>=fMinEvents && _bsepgains[ivar][ibin]>maxsepgain && _bsepgainsigs[ivar][ibin]>fMinCutSignificance) {
      //if ( (_varvals[ivar][ibin] - limits[ivar].first)>(0.2*sigmanode) && (limits[ivar].second-_varvals[ivar][ibin])>(0.2*sigmanode) && _sumws2[ivar][ibin]>=fMinEvents && (sumw2-_sumws2[ivar][ibin])>=fMinEvents && _bsepgains[ivar][ibin]>maxsepgain && _bsepgainsigs[ivar][ibin]>fMinCutSignificance) {
      
	

	  
    // if ( (_varvals[ivar][ibin] - limits[ivar].first)>(0.2*sigmanode) && (limits[ivar].second-_varvals[ivar][ibin])>(0.2*sigmanode) && _sumns[ivar][ibin]>=fMinEvents && (numprimary-_sumns[ivar][ibin])>=fMinEvents && _bsepgains[ivar][ibin]>maxsepgain ) {	  
	
	//printf("accepted split\n");
      
        //printf("sumnleft = %i, sumnright = %i, sepgain = %5f, sepgainsig = %5f\n",_sumns[ivar][ibin], (nev-_sumns[ivar][ibin]),_bsepgains[ivar][ibin],_bsepgainsigs[ivar][ibin]);
      //if (sumns[ivar][ibin]>=fMinEvents && (nev-sumns[ivar][ibin])>=fMinEvents && bsepgains[ivar][ibin]>maxsepgain) {
	//if (_sumns[ivar][ibin]>=fMinEvents && (nev-_sumns[ivar][ibin])>=fMinEvents && _bsepgains[ivar][ibin]>maxsepgain && _bsepgainsigs[ivar][ibin]>fMinCutSignificance) {
      //if (_sumnsd[ivar][ibin]>=fMinEvents && (nevd-_sumnsd[ivar][ibin])>=fMinEvents && _bsepgains[ivar][ibin]>maxsepgain && _bsepgainsigs[ivar][ibin]>fMinCutSignificance) {
     //if ( (_varvals[ivar][ibin] - limits[ivar].first)>(1.0*sigmanode) && (limits[ivar].second-_varvals[ivar][ibin])>(1.0*sigmanode) && _sumnsd[ivar][ibin]>=fMinEvents && (nevd-_sumnsd[ivar][ibin])>=fMinEvents && (_sumns[ivar][ibin]-_sumnsd[ivar][ibin])>=fMinEvents && (nev-_sumns[ivar][ibin] - nevd + _sumnsd[ivar][ibin])>=fMinEvents && _bsepgains[ivar][ibin]>maxsepgain && _bsepgainsigs[ivar][ibin]>fMinCutSignificance) {
       
       
      //if ( (_varvals[ivar][ibin] - limits[ivar].first)>(0.2*sigmanode) && (limits[ivar].second-_varvals[ivar][ibin])>(0.2*sigmanode) && (_sumns[ivar][ibin]-_sumnsd[ivar][ibin])>=fMinEvents && (numprimary-_sumns[ivar][ibin] - nevd + _sumnsd[ivar][ibin])>=fMinEvents && _bsepgains[ivar][ibin]>maxsepgain && _bsepgainsigs[ivar][ibin]>fMinCutSignificance) {       
       
       
	//if ( (_varvals[ivar][ibin] - limits[ivar].first)>(0.2*sigmanode) && (limits[ivar].second-_varvals[ivar][ibin])>(0.2*sigmanode) && _sumws2[ivar][ibin]>=fMinEvents && (sumw2-_sumws2[ivar][ibin])>=fMinEvents && _bsepgains[ivar][ibin]>maxsepgain ) {       
       
       //if (_sumns[ivar][ibin]>=fMinEvents && (nev-_sumns[ivar][ibin])>=fMinEvents && _bsepgains[ivar][ibin]>maxsepgain && _bsepgainsigs[ivar][ibin]>fMinCutSignificance) {
       //if (_sumws[ivar][ibin]>=fMinEvents && (sumw -_sumws[ivar][ibin])>=fMinEvents && _bsepgains[ivar][ibin]>maxsepgain && _bsepgainsigs[ivar][ibin]>fMinCutSignificance) {
	 
       //printf("_bsepgains[ivar][ibin] = %5f, maxsepgain = %5f\n",_bsepgains[ivar][ibin],maxsepgain);
	 
       //if (_sumnsd[ivar][ibin]>=fMinEvents && (nevd-_sumnsd[ivar][ibin])>=fMinEvents && _bsepgains[ivar][ibin]>maxsepgain && _bsepgainsigs[ivar][ibin]>fMinCutSignificance) {
	 
          
          
// 	maxsepgain = _bsepgains[ivar][ibin];
//         bestbin = ibin;
      //}
      


      if ( _bsepgains[ivar][ibin]>maxsepgain && !std::isinf(_bsepgains[ivar][ibin]) && _sumtgt2s[ivar][ibin]>0. && (sumtgt2-_sumtgt2s[ivar][ibin])>0.) {
	
	bool passminweights = true;
	double minweights = std::numeric_limits<double>::max();
	for (int icls=0; icls<ncls; ++icls) {
	  double minweightcls = std::min(_sumwscls[ivar][icls][ibin], _sumwscls[ivar][icls][nbins-1] - _sumwscls[ivar][icls][ibin]);
	  if (minweightcls < fMinWeights[icls]) {
	    passminweights = false;
	  }
	  if (minweightcls < minweights) {
	    minweights = minweightcls;
	  }
	}
	
	
	bool passminsig = true;
	
	passminsig &= fMinCutSignificance<0. || (_bsepgains[ivar][ibin] > fMinCutSignificance);
	
        if (fMaxNSpurious >= 0.) {
          double prob = 1.0 - TMath::Erf(sqrt(_bsepgains[ivar][ibin]));
          double nspurious = fSumWTimesNVars*prob;	
          passminsig &= nspurious<fMaxNSpurious;	
        }
	
	if (passminweights && passminsig) {
	  maxsepgain = _bsepgains[ivar][ibin];
	  bestbin = ibin;
	}
      }
      
    }
     
    //printf("ivar = %i, maxsepgain = %5f, left: w0 = %5f, w1 = %5f, right: w0 = %5f, w1 = %5f\n", ivar, maxsepgain, _sumwscls[ivar][bestbin][0],_sumwscls[ivar][bestbin][1],sumwscls[0] - _sumwscls[ivar][bestbin][0],sumwscls[1] - _sumwscls[ivar][bestbin][1]);
     
     
     
    cutval = _varvals[ivar][bestbin];
    nleft = _sumns[ivar][bestbin];
    nright = nev - nleft;
   // sumwleft = _sumws[ivar][bestbin];
   // sumwright = sumw - sumwleft;        
    
    
    _sepgains[ivar] = maxsepgain;
    _sepgainsigs[ivar] = _bsepgainsigs[ivar][bestbin];
    _cutvals[ivar] = cutval;
    _nlefts[ivar] = nleft;
    _nrights[ivar] = nright;
    //_sumwlefts[ivar] = sumwleft;
    //_sumwrights[ivar] = sumwright;
    //_sumtgtlefts[ivar] = _sumtgts[ivar][bestbin];
    //_sumtgtrights[ivar] = sumtgt - _sumtgts[ivar][bestbin];
    //_leftvars[ivar] = _sumtgt2s[ivar][bestbin] - _sumtgts[ivar][bestbin]*_sumtgts[ivar][bestbin]/_sumws[ivar][bestbin];
    //_rightvars[ivar] = (sumtgt2-_sumtgt2s[ivar][bestbin]) - (sumtgt-_sumtgts[ivar][bestbin])*(sumtgt-_sumtgts[ivar][bestbin])/(sumw-_sumws[ivar][bestbin]);
    _bestbins[ivar] = bestbin;
        
     //printf("done var %i\n",ivar);
//     
//     watch.Stop();
//     timesvar[ivar] = watch.RealTime();
  }
  
//   clockloop.Stop();
//   printf("end split search loop\n");
//   
//   
//   for (int ivar=0; ivar<nvars; ++ivar) {
//     printf("ivar = %i, ithread = %i, nbins = %i, time = %5e\n",ivar,threadnums[ivar],nbinsvar[ivar],timesvar[ivar]);
//   }
//   printf("clockloop = %5f\n",clockloop.RealTime());
  
  float globalsepgain = 0.;
  for (int ivar=0; ivar<nvars; ++ivar) {
    if (_sepgains[ivar]>globalsepgain) {
      globalsepgain = _sepgains[ivar];
      bestvar = ivar;
    }
  }    
  
  
  
  //if no appropriate split found, make this node terminal
  if (globalsepgain<=0.) {
    //no valid split found, making this node a leaf node
    //printf("thisidx = %i, globalsepgain = %5f, no valid split\n",thisidx, globalsepgain);
    tree.CutIndices().push_back(0);
    tree.CutVals().push_back(0);
    tree.LeftIndices().push_back(0);   
    tree.RightIndices().push_back(0);    
    
    tree.RightIndices()[thisidx] = -tree.Responses().size();
    tree.LeftIndices()[thisidx] = -tree.Responses().size();
    
    BuildLeaf(evts,tree,tgtidx);
    return;
  }
  
  //printf("cutval = %5f\n",_cutvals[bestvar]);
  
//   TStopwatch buildclock;
//   buildclock.Start();
//   
//   printf("start filling vectors\n");
  
  //fill vectors of event pointers for left and right nodes below this one
  std::vector<HybridGBREvent*> leftevts;
  std::vector<HybridGBREvent*> rightevts;
  
  leftevts.reserve(nev);
  rightevts.reserve(nev);
  
  int nleft = 0;
  int nright = 0;
  double sumwleft = 0.;
  double sumwright = 0.;
  
  
  
  for (std::vector<HybridGBREvent*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
    if ((*it)->Var(bestvar)>_cutvals[bestvar]) {
      ++nright;
      sumwright += (*it)->Weight();
      rightevts.push_back(*it);
    }
    else {
      ++nleft;
      sumwleft += (*it)->Weight();
      leftevts.push_back(*it);
    }
    
  }
  
//   printf("done filling vectors\n");
//   
//   buildclock.Stop();
//   printf("buildclock = %5e\n",buildclock.RealTime());
  
//   for (std::vector<HybridGBREvent*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
//     if ((*it)->Var(bestvar)>_cutvals[bestvar]) {
//       ++nright;
// //       sumwright += (*it)->Weight();
// //       rightevts.push_back(*it);
//     }
//     else {
//       ++nleft;
// //       sumwleft += (*it)->Weight();
// //       leftevts.push_back(*it);
//     }
//     
//     //double sigma = 0.2;
//     //double sigma = fSigmaConsts[bestvar]*sqrt(fNorm/exp(-(*it)->Target()));
//     //double sigma = 0.2;
//     //double sigma = std::min(0.4,fSigmaConsts[bestvar]*sqrt(fNorm/exp(-(*it)->Target())));
//     //double sigmanode = fSigmaConsts[ivar]*pow(1.0/avgp,1.0/fCondVars.getSize());
//     //double sigma = std::min(2.2,fSigmaConsts[bestvar]*pow(fNorm/exp(-(*it)->Target()),1.0/fCondVars.getSize()));
//     //double sigma = 0.2;
//     //printf("sigma = %5f\n",sigma);
//     
//     double sigma = 1.0;
//     
//     double gausxcut = (_cutvals[bestvar]- (*it)->Var(bestvar))/sigma;
//     double gausxlow = (limits[bestvar].first- (*it)->Var(bestvar))/sigma;
//     double gausxhi  = (limits[bestvar].second- (*it)->Var(bestvar))/sigma;
//     
//     double leftweight;
//     double rightweight;
//     
//     if (1) {   
//       leftweight = 0.5*(1.0+TMath::Erf(gausxcut/sqrt(2))) - 0.5*(1.0+TMath::Erf(gausxlow/sqrt(2)));
//       rightweight = 0.5*(1.0+TMath::Erf(gausxhi/sqrt(2))) - 0.5*(1.0+TMath::Erf(gausxcut/sqrt(2)));
//     }
//     else {
//       if ((*it)->Var(bestvar)>_cutvals[bestvar]) {
//         rightweight = 1.0;
//         leftweight = 0.;
//       }
//       else {
//         rightweight = 0.0;
//         leftweight = 1.0;
//       }
//     }
//     
//     //printf("leftweight = %5f, rightweight = %5f, bestvar = %5f, cutval = %5f, low = %5f, hi = %5f\n",leftweight,rightweight,(*it)->Var(bestvar),_cutvals[bestvar],limits[bestvar].first,limits[bestvar].second);
//     
//     const double minweight = 1e-4;
//     
//     if (leftweight*(*it)->Weight()<minweight && sumwleft>0.0) {
//       rightweight += leftweight;
//       leftweight = 0.;
//     }
//     else if (rightweight*(*it)->Weight()<minweight && sumwright>0.0) {
//       leftweight += rightweight;
//       rightweight = 0.;
//     }
//     
//     if (leftweight>0.) {
//       HybridGBREvent *evtleft = new HybridGBREvent(nvars, fNTargets,**it);
//       evtleft->SetWeight(leftweight*evtleft->Weight());
//       if ((*it)->Var(bestvar) > _cutvals[bestvar]) evtleft->SetIsPrimary(false);
//       sumwleft += evtleft->Weight();
//       leftevts.push_back(evtleft);
//     }
//     
//     if (rightweight>0.) {
//       HybridGBREvent *evtright = new HybridGBREvent(nvars, fNTargets,**it);
//       evtright->SetWeight(rightweight*evtright->Weight());
//       if ((*it)->Var(bestvar) <= _cutvals[bestvar]) evtright->SetIsPrimary(false);
//       sumwright += evtright->Weight();
//       rightevts.push_back(evtright);
//     }
// 
//     
//   }  
  
 
//   float fullres = sqrt(_fullvars[bestvar]/sumwtotal);
//   float leftres = sqrt(_leftvars[bestvar]/sumwleft);
//   float rightres = sqrt(_rightvars[bestvar]/sumwright);
//  
//   float fullmean = (_sumtgtlefts[bestvar] + _sumtgtrights[bestvar])/sumwtotal;
//   float leftmean = _sumtgtlefts[bestvar]/sumwleft;
//   float rightmean = _sumtgtrights[bestvar]/sumwright;
  
  
  //printf("thisidx = %i, bestvar = %i, cutval = %5f, n = %i, nleft = %i, nright = %i, fullres = %5f, leftres = %5f, rightres = %5f, fullmean = %5f, leftmean = %5f, rightmrean = %5f, leftsepgain = %5f, sepgainsig = %5f\n",thisidx,bestvar,_cutvals[bestvar],nev,_nlefts[bestvar],_nrights[bestvar],fullres,leftres,rightres,fullmean, leftmean, rightmean, _sepgains[bestvar],_sepgainsigs[bestvar]);
  
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
  //bool termleft = nleft<=(2*fMinEvents) || depth==fMaxDepth;
  bool termleft = sumwleft<=(2*fMinEvents) || (fMaxDepth>=0 && depth==fMaxDepth) || (fMaxNodes>=0 && int(tree.Responses().size())>=fMaxNodes) ;
  if (termleft) tree.LeftIndices()[thisidx] = -tree.Responses().size();
  else tree.LeftIndices()[thisidx] = tree.CutIndices().size();
  
  //printf("this idx = %i, termleft = %i, nleft = %i, fMinEvents = %i\n",thisidx,  termleft,nleft,fMinEvents);  
  
  //build left node as appropriate
  std::vector<std::pair<float,float> > limitsleft(limits);
  limitsleft[bestvar].second = bestcutval;  
  //printf("bestvar = %i, limlow = %5f, limhigh = %5f, limleftlow = %5f, limlefthigh = %5f\n",bestvar,limits[bestvar].first,limits[bestvar].second,limitsleft[bestvar].first,limitsleft[bestvar].second);
  if (termleft) {  
    BuildLeaf(leftevts,tree,tgtidx);
  }
  else {  
    TrainTree(leftevts,sumwleft,tree,nvars,transition,depth+1,limitsleft,tgtidx);  
  }
  
  //check if right node is terminal
  //bool termright = nright<=(2*fMinEvents) || depth==fMaxDepth;
  bool termright = sumwright<=(2*fMinEvents) || (fMaxDepth>=0 && depth==fMaxDepth) || (fMaxNodes>=0 && int(tree.Responses().size())>=fMaxNodes);
  if (termright) tree.RightIndices()[thisidx] = -tree.Responses().size();
  else tree.RightIndices()[thisidx] = tree.CutIndices().size();
    
  //printf("this idx = %i, termright = %i, nright = %i, fMinEvents = %i\n",thisidx,  termright,nright,fMinEvents);    
  
  //build right node as appropriate
  std::vector<std::pair<float,float> > limitsright(limits);
  limitsright[bestvar].first = bestcutval;  
  //printf("bestvar = %i, limlow = %5f, limhigh = %5f, limrightlow = %5f, limrighthigh = %5f\n",bestvar, limits[bestvar].first,limits[bestvar].second,limitsright[bestvar].first,limitsright[bestvar].second);  
  if (termright) {  
    BuildLeaf(rightevts,tree,tgtidx);
  }
  else {      
    TrainTree(rightevts,sumwright,tree,nvars,transition,depth+1,limitsright, tgtidx);  
  }
  
}

  

  


//_______________________________________________________________________
void RooHybridBDTAutoPdf::BuildLeaf(const std::vector<HybridGBREvent*> &evts, HybridGBRTree &tree, int tgtidx) {

  //printf("building leaf\n");
  
  int thisidx = -tree.Responses().size();
    
  tree.Responses().push_back(0.);
  
  for (std::vector<HybridGBREvent*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {    
    (*it)->SetCurrentNode(tgtidx, -thisidx);
  }   
  
}





double RooHybridBDTAutoPdf::EvalLossNull(double dummy) {
 
  return gHybridBDTAutoPointer->EvalLossRooFit();
  
}

double RooHybridBDTAutoPdf::EvalLossRooFit() {
 
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
  RooAbsReal::clearEvalErrorLog() ;  

  for (int ithread=0; ithread<fNThreads; ++ithread) {
    for (int ivar=0; ivar<fFullParms.getSize(); ++ivar) {
      RooRealVar *cloneparm = static_cast<RooRealVar*>(fFullParmsClones[ithread].at(ivar));
      cloneparm->setVal(static_cast<RooRealVar*>(fFullParms.at(ivar))->getVal());
      cloneparm->setError(static_cast<RooRealVar*>(fFullParms.at(ivar))->getError());
    }
  } 

  double nllval = 0.;
  //RooArgSet parmset(fParmVars);
  
  std::vector<double> nllvals(fNThreads);
	
  #pragma omp parallel for
  for (unsigned int ievt=0; ievt<fEvts.size(); ++ievt) {

    //if (ievt%20!=0) continue;
    //if (ievt%100!=0) continue;
    
    if (fPrescaleInit>0 && ievt%fPrescaleInit!=0) continue;
    
    int ithread =  omp_get_thread_num();
    //int ithread =  0;

    //int termidx = fEvts[ievt]->CurrentNode();
    //if (seltermidx>=0 && termidx!=seltermidx) continue;
        
    for (int ivar=0; ivar<fCondVarsClones[ithread].getSize(); ++ivar) {
      static_cast<RooRealVar*>(fCondVarsClones[ithread].at(ivar))->setVal(fEvts[ievt]->Var(ivar));
    }
    for (int ivar=0; ivar<fParmVarsClones[ithread].getSize(); ++ivar) {
      static_cast<RooRealVar*>(fParmVarsClones[ithread].at(ivar))->setVal(fEvts[ievt]->Var(fCondVarsClones[ithread].getSize() + ivar));
    }    

    
    int evcls = fEvts[ievt]->Class(); 
    double weight = fEvts[ievt]->Weight();

    //nll value for minos constraint
    //if (evcls==0) fNLLVal += -log(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset));
    
    double pdfval = static_cast<RooAbsReal*>(fStaticPdfsClones[ithread].at(evcls))->getValV(&fParmSetClones[ithread]);
    
    //nllval += -weight*log(pdfval);
    //nllvals[ithread] += -weight*vdt::fast_logf(pdfval);
    nllvals[ithread] += -weight*log(pdfval);
    
    if (RooAbsReal::numEvalErrors()>0 || RooAbsPdf::evalError() || pdfval<0.) {
      nllvals[ithread] += std::numeric_limits<float>::max();
    }
    
/*    if (RooAbsReal::numEvalErrors()>0 || RooAbsPdf::evalError() || pdfval<0.) {
    //if (RooAbsReal::numEvalErrors()>0 || RooAbsPdf::evalError()) {  
      //printf("numEvalErrors = %i\n",RooAbsReal::numEvalErrors());
      RooAbsReal::clearEvalErrorLog();
      RooAbsPdf::clearEvalError();
      printf("eval error in EvalLoss\n");
      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
      return std::numeric_limits<double>::max();
    }   */ 
    
    
  }
  
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    nllval += nllvals[ithread];
  }
  
  int infunc = fFullFuncs.getSize()-1;
  nllval += static_cast<RooAbsReal*>(fFullFuncs.at(infunc))->getVal();

  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
  
  return nllval;  
  
}

double RooHybridBDTAutoPdf::EvalLoss(HybridGBRForest *forest, double lambda, const TVectorD &dL, int itree) {
 
 // int tgtidx = itree%fNTargets;
  
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
  RooAbsReal::clearEvalErrorLog() ;  
  
//   int nterm = tree.Responses()[0].size();
// 
  int nextvars = fExtVars.getSize();
// 
  
  int nparms = nextvars;
  std::vector<int> localidxs(fNTargets);
  for (int itgt=0; itgt<fNTargets; ++itgt) {
    localidxs[itgt] = nparms;
    int treetgt = static_cast<RooGBRTarget*>(fTgtVars.at(itgt))->Index();
    if (forest->Trees()[treetgt].size())
      nparms += forest->Trees()[treetgt].back().Responses().size();
    else 
      nparms += 1;
  }   
   
//   int nparms = fExtVars.getSize() + fNTargets*nterm; 
  
  //printf("nextvars = %i, nparms = %i\n",nextvars,nparms);
  
 // double dr = 0.;
  
  int idxglobal = 0;
  //int idxn0 = idxglobal + fExtVars.index(fN0);
  
  //printf("idxn0 = %i, fN0 = %5f\n", idxn0,fN0->getVal());

  double nllval = 0.;
  
  
  //printf("initval = %5f\n",static_cast<RooRealVar*>(fExtVars.at(1))->getVal());
  
  //global variables
//   RooArgList extbak = fExtVars;
//   
  std::vector<double> extvals(fExtVars.getSize());
  
  for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
    int iel = idxglobal + ivar;
    extvals[ivar] = static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal();
    static_cast<RooRealVar*>(fExtVars.at(ivar))->setVal(static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal() + lambda*dL(iel));   
  }
  
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
      RooRealVar *cloneparm = static_cast<RooRealVar*>(fExtVarsClones[ithread].at(ivar));
      cloneparm->setVal(static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal());
      cloneparm->setError(static_cast<RooRealVar*>(fExtVars.at(ivar))->getError());
    }
  }      
  
//   printf("intermediateval = %5f\n",static_cast<RooRealVar*>(fExtVars.at(1))->getVal());
//   
//   printf("begin loop\n");

//   int pointidx = -1;
//   if (itree>=0) {
//     pointidx = itree%fFullParms.getSize();
//   }

  std::vector<double> nllvals(fNThreads);

  #pragma omp parallel for
  for (unsigned int ievt=0; ievt<fEvts.size(); ++ievt) {

    int ithread =  omp_get_thread_num();
    //int ithread =  0;
    
    
    //int termidx = fEvts[ievt]->CurrentNode();
    //if (seltermidx>=0 && termidx!=seltermidx) continue;
    
    //int idxlocal = nextvars + fNTargets*termidx;
    
    for (int ivar=0; ivar<fCondVarsClones[ithread].getSize(); ++ivar) {
      static_cast<RooRealVar*>(fCondVarsClones[ithread].at(ivar))->setVal(fEvts[ievt]->Var(ivar));
    }
    for (int ivar=0; ivar<fParmVarsClones[ithread].getSize(); ++ivar) {
      static_cast<RooRealVar*>(fParmVarsClones[ithread].at(ivar))->setVal(fEvts[ievt]->Var(fCondVarsClones[ithread].getSize() + ivar));
    }    

    for (int itgt=0; itgt<fNTargets; ++itgt) {
      int iel = localidxs[itgt] + fEvts[ievt]->CurrentNode(itgt);
      static_cast<RooRealVar*>(fStaticTgtsClones[ithread].at(itgt))->setVal(fEvts[ievt]->Target(itgt) + float(lambda*dL[iel]));
    }
    
    int evcls = fEvts[ievt]->Class(); 
    double weight = fEvts[ievt]->Weight();

    //nll value for minos constraint
    //if (evcls==0) fNLLVal += -log(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset));
    
    double pdfval = static_cast<RooAbsReal*>(fStaticPdfsClones[ithread].at(evcls))->getValV(&fParmSetClones[ithread]);
    fEvts[ievt]->SetPdfVal(pdfval);
    
    //nllval += -weight*log(pdfval);
    //nllvals[ithread] += -weight*vdt::fast_logf(pdfval);
    nllvals[ithread] += -weight*log(pdfval); 
    
    
//     if (pointidx>=0) {
//       fEvts[ievt]->ValVector()(pointidx) = pdfval;
//       for (int iparm=0; iparm<fFullParms.getSize(); ++iparm) {
// 	fEvts[ievt]->ParmMatrix()(iparm,pointidx) = static_cast<RooRealVar*>(fFullParms.at(iparm))->getVal();
//       }
//     }
    
//     if (RooAbsReal::numEvalErrors()>0) {
//       RooAbsReal::clearEvalErrorLog();
//       return std::numeric_limits<double>::max();
//     }
    
//     double testmass = fEvts[ievt]->Var(fCondVars.getSize());
//     if (testmass>178.546 && testmass<178.548)
//       printf("evcls = %i, pdfval = %5e, logmode = %i\n", evcls,pdfval,int(RooAbsReal::evalErrorLoggingMode()));
    
    if (RooAbsReal::numEvalErrors()>0 || RooAbsPdf::evalError() || pdfval<0.) {
      nllvals[ithread] += std::numeric_limits<float>::max();
    }
    
    
/*    if (RooAbsReal::numEvalErrors()>0 || RooAbsPdf::evalError() || pdfval<0.) {
    //if (RooAbsReal::numEvalErrors()>0 || RooAbsPdf::evalError()) {  
      //printf("numEvalErrors = %i\n",RooAbsReal::numEvalErrors());
      RooAbsReal::clearEvalErrorLog();
      RooAbsPdf::clearEvalError();
      printf("eval error in EvalLoss\n");
      for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
        static_cast<RooRealVar*>(fExtVars.at(ivar))->setVal(extvals[ivar]);   
      }        
      for (int ithread=0; ithread<fNThreads; ++ithread) {
	for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
	  RooRealVar *cloneparm = static_cast<RooRealVar*>(fExtVarsClones[ithread].at(ivar));
	  cloneparm->setVal(static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal());
	  cloneparm->setError(static_cast<RooRealVar*>(fExtVars.at(ivar))->getError());
	}
      }          
      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
      return std::numeric_limits<double>::max();
    }*/    
    
    //if (std::isnan(nllval) || std::isinf(nllval) || nllval==0.) return std::numeric_limits<double>::max();
    
    //normalization terms
//     double normden = 1.0;
//     for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {
//       int itgt = ipdf - 1;
//       int iel = idxlocal + itgt;
//       double tgtval = fEvts[ievt]->Target(itgt) + lambda*dL[iel];
//       double expF = exp(tgtval);
//       normden += expF;
//       
//       if (ipdf==evcls) {
//         nllval += -weight*tgtval;
//       }
//     }
//     nllval += weight*log(normden);
    
    
/*    if (evcls==0) {
      for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {
	int itgt = ipdf - 1;
	int iel = idxlocal + itgt;
	double nexpF = exp(-fEvts[ievt]->Target(itgt) - lambda*dL[iel]);
	//nllval += weight*nexpF/fN0Obs;
	nllval += weight*nexpF;
      }
    }
    else {
      int itgt = evcls - 1;
      int iel = idxlocal + itgt;
      nllval += weight*(fEvts[ievt]->Target(itgt) + lambda*dL[iel]);      
    }*/    
    
  }
  
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    nllval += nllvals[ithread];
  }
  
  //global terms
  //if (seltermidx<0) {
  int infunc = fFullFuncs.getSize()-1;
  nllval += static_cast<RooAbsReal*>(fFullFuncs.at(infunc))->getVal();
  //nllval += static_cast<RooAbsReal*>(fFullFuncs.at(infunc))->getVal() - fN0Obs*log(static_cast<RooAbsReal*>(fFullFuncs.at(infunc))->getVal());  
  //}
  
  //fExtVars = extbak;
  
  for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
    static_cast<RooRealVar*>(fExtVars.at(ivar))->setVal(extvals[ivar]);   
  }  
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
      RooRealVar *cloneparm = static_cast<RooRealVar*>(fExtVarsClones[ithread].at(ivar));
      cloneparm->setVal(static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal());
      cloneparm->setError(static_cast<RooRealVar*>(fExtVars.at(ivar))->getError());
    }
  }      
  
  //printf("finalval = %5f\n",static_cast<RooRealVar*>(fExtVars.at(1))->getVal());

  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
  
  return nllval;
      
      
}


double RooHybridBDTAutoPdf::EvalLossAvg() {
 
  //double nllval = 0.;
  
  RooArgSet parmset(fParmVars);
	
  
  
  TH1D *hmass = new TH1D("hmass","",400,0.,100.);

    
  
  //std::vector<double> clsvals(fData.size());  
  
  for (std::vector<HybridGBREvent*>::const_iterator it = fEvts.begin(); it!=fEvts.end(); ++it) {
    
    int evcls = (*it)->Class(); 
    if (evcls!=0) continue;
    
    double weight = (*it)->Weight();    
    double pdfval = (*it)->PdfVal();
    
    double expF = exp((*it)->Target(0));
    double mpdfval = (1.0-expF)*pdfval;
    
    double mval = (*it)->Var(fCondVars.getSize() + 0);
    
    hmass->Fill(mval,weight*mpdfval);
    
    //clsvals[evcls] += clsweights[evcls]*weight*static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset);
        
    //nllval += -weight*log(pdfval); 
         
  }  
  
  
  hmass->Scale(1.0/hmass->GetSumOfWeights());
  
  double nllval = 0.;
  for (std::vector<HybridGBREvent*>::const_iterator it = fEvts.begin(); it!=fEvts.end(); ++it) {
    int evcls = (*it)->Class(); 
    if (evcls!=0) continue;
    
    double weight = (*it)->Weight();     
    
    double mval = (*it)->Var(fCondVars.getSize() + 0);

    
    double pdfval = hmass->GetBinContent(hmass->GetXaxis()->FindFixBin(mval));
    nllval += -weight*log(pdfval);
  }
  
  //delete hmass;
  
  new TCanvas;
  hmass->Draw();
  
  return nllval;

  
  
//   for (std::vector<HybridGBREvent*>::const_iterator it = fEvts.begin(); it!=fEvts.end(); ++it) {
//     
// 
//     for (int ivar=0; ivar<fParmVars.getSize(); ++ivar) {
//       static_cast<RooRealVar*>(fParmVars.at(ivar))->setVal((*it)->Var(fCondVars.getSize() + ivar));
//     }
//         
//     int evcls = (*it)->Class(); 
//     double weight = (*it)->Weight();    
//     
//     double pdfval = 0.;
//     
//     for (std::vector<HybridGBREvent*>::const_iterator jt = fEvts.begin(); jt!=fEvts.end(); ++jt) { 
//       int jcls = (*jt)->Class();
//       
//       if (jcls!=evcls) continue;
//       
//       int jweight = (*jt)->Weight();
//       
//       for (int ivar=0; ivar<fCondVars.getSize(); ++ivar) {
// 	static_cast<RooRealVar*>(fCondVars.at(ivar))->setVal((*jt)->Var(ivar));
//       }    
//       
//       for (int itgt=0; itgt<fNTargets; ++itgt) {
// 	static_cast<RooRealVar*>(fStaticTgts.at(itgt))->setVal((*jt)->Target(itgt));
//       }
//       
//       pdfval += clsweights[evcls]*jweight*static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset);
//       
//     }
//    
//     nllval += -weight*log(pdfval); 
//          
//   }
  
//   int infunc = fFullFuncs.getSize()-1;
//   nllval += static_cast<RooAbsReal*>(fFullFuncs.at(infunc))->getVal();
//   
//   double sumweight = 0.;
//   double sumnll = 0.;
//   double sumnllsq = 0.;
//   for (std::vector<HybridGBREvent*>::const_iterator it = fEvts.begin(); it!=fEvts.end(); ++it) {
//     double weight = (*it)->Weight();
//     double pdfval = (*it)->PdfVal();
//     
//     double modnll = nllval + weight*log(pdfval);
//     
//     sumweight += weight;
//     sumnll += modnll;
//     sumnllsq += modnll*modnll;
//     
// 
//   }
//   
//   
//   return sqrt(sumnllsq/sumweight - sumnll*sumnll/sumweight/sumweight);
      
      
}









TMatrixD RooHybridBDTAutoPdf::vmultT(const TVectorD &v, const TVectorD &vT) const
{
  int nrows = v.GetNrows();
  int ncols = vT.GetNrows();
  TMatrixD m(nrows,ncols);
  
  for (int irow=0; irow<nrows; ++irow) {
    for (int icol=0; icol<ncols; ++icol) {
      m(irow,icol) = v(irow)*vT(icol);
    }
  }
  
  return m;
  
}
  
double RooHybridBDTAutoPdf::vmult(const TVectorD &vT, const TVectorD &v) const
{
 
  assert(vT.GetNrows()==v.GetNrows());
  
  double val = 0.;
  for (int iel=0; iel<vT.GetNrows(); ++iel) {
    val += vT(iel)*v(iel);
  }
  
  return val;
  
}

void RooHybridBDTAutoPdf::FillDerivatives() {
 
  

//  int nextvars = fExtVars.getSize();
  
  int nparms = fFullParms.getSize();
  
  //printf("nextvars = %i, nparms = %i\n",nextvars,nparms);
  
//  double dr = 0.;
  
//  int idxglobal = 0;
//  int idxlocal = nextvars;
  //int idxn0 = idxglobal + fExtVars.index(fN0);
  
  //printf("idxn0 = %i, fN0 = %5f\n", idxn0,fN0->getVal());

  double nmc = 0.;
  double nmcalt = 0.;
  
    
  int msize = nparms;
  
  
  //double lambda = 1.0;
//    double nll = 0.;
  fNLLVal = 0.; 



  RooArgSet parmset(fParmVars);
	
  //printf("begin loop\n");
  
  for (std::vector<HybridGBREvent*>::const_iterator it = fEvts.begin(); it!=fEvts.end(); ++it) {

    int ievt(it-fEvts.begin());
    
    TVectorD &dL = fGradients[ievt];
    TMatrixD &d2L = fHessians[ievt];
    
    //int termidx = (*it)->CurrentNode();
    //int idxlocal = nextvars + fNTargets*termidx;
    
    
    for (int ivar=0; ivar<fCondVars.getSize(); ++ivar) {
      static_cast<RooRealVar*>(fCondVars.find(*fCondVars.at(ivar)))->setVal((*it)->Var(ivar));
    }
    for (int ivar=0; ivar<fParmVars.getSize(); ++ivar) {
      static_cast<RooRealVar*>(fParmVars.find(*fParmVars.at(ivar)))->setVal((*it)->Var(fCondVars.getSize() + ivar));
    }    

    for (int itgt=0; itgt<fNTargets; ++itgt) {
      static_cast<RooRealVar*>(fStaticTgts.at(itgt))->setVal((*it)->Target(itgt));
    }
    
    int evcls = (*it)->Class(); 

    //nll value for minos constraint
    fNLLVal += -log(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset));
    
    
    //derivatives for pdfs
    //for (int ivar=0; ivar<fFullParms.getSize(); ++ivar) {
    for (unsigned int iidx=0; iidx<fOuterIndices[evcls].size(); ++iidx) {
      
      int ivar = fOuterIndices[evcls][iidx];
      int idrv = ivar;
//      int itgt = ivar - fExtVars.getSize();
//      int iel;
//      if (itgt>=0) iel = idxlocal + itgt;
//      else iel = idxglobal + ivar;
      
      double drvi = Derivative1(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls)),static_cast<RooRealVar*>(fFullParms.at(idrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(idrv))->getError());
      //double drvialt = fDerivatives[evcls][idrv]->getVal();
      
      //printf("drv%i = %5f, alt = %5f\n",idrv,drvi,drvialt);
      
      
      //dL[iel] += -fDerivatives[evcls][idrv]->getVal()/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset);
      dL[ivar] = -drvi/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset);
      
      
// 	if (ivar==0 && std::abs(fDerivatives[evcls][idrv]->getVal())>0.01) printf("derivative for class %i, drv %i = %5f\n",evcls,idrv,fDerivatives[evcls][idrv]->getVal());
// 	if (static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset)<1e-4) printf("pdf %i val = %5f\n",evcls,static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset));
      //if (ivar==0) printf("dL[0] = %5f, fN0 = %5f\n", dL[0],fN0->getVal());
      //for (int jvar=0; jvar<fFullParms.getSize(); ++jvar) {
	
      for (unsigned int jidx=iidx; jidx<fOuterIndices[evcls].size(); ++jidx) {
	
	int jvar = fOuterIndices[evcls][jidx];
	int jdrv = jvar;
//	int jtgt = jvar - fExtVars.getSize();
	//int jel;
	//if (jtgt>=0) jel = idxlocal + jtgt;
	//else jel = idxglobal + jvar;
	
	double drvj; 
	double drv2ij;
	if (idrv==jdrv) {
	  drvj = drvi;
	  drv2ij = Derivative2(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls)),static_cast<RooRealVar*>(fFullParms.at(idrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(idrv))->getError());
	}
	else {
	  drvj = Derivative1(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls)),static_cast<RooRealVar*>(fFullParms.at(jdrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(jdrv))->getError());
	  drv2ij = Derivative2(static_cast<RooAbsReal*>(fStaticPdfs.at(evcls)),static_cast<RooRealVar*>(fFullParms.at(idrv)),static_cast<RooRealVar*>(fFullParms.at(jdrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(idrv))->getError(),1e-3*static_cast<RooRealVar*>(fFullParms.at(jdrv))->getError());
	}
	
	d2L[ivar][jvar] = -drv2ij/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset) + drvi*drvj/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset)/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset);	  
	
	//d2L[iel][jel] += -f2Derivatives[evcls][idrv][jdrv]->getVal()/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset) + fDerivatives[evcls][idrv]->getVal()*fDerivatives[evcls][jdrv]->getVal()/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset)/static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset);

	  
	  
      }
    }
    
    for (int i=0; i<msize; ++i) {
      for (int j=0; j<i; ++j) {
	//assert(d2L(i,j)==0.);
	d2L(i,j) = d2L(j,i);
      }
    }  
    
    //normalization terms
    if (evcls==0) {
      for (int ipdf=1; ipdf<fStaticPdfs.getSize(); ++ipdf) {
	int itgt = ipdf - 1;
//	int iel = idxlocal + itgt;
	double nexpF = exp(-(*it)->Target(itgt));
	fNLLVal += nexpF/fN0Obs;
	nmc += nexpF/fN0Obs;
	
// 	dL[iel] += -nexpF/fN0Obs;
// 	d2L[iel][iel] += nexpF/fN0Obs;
      }
    }
    else {
      int itgt = evcls - 1;
//      int iel = idxlocal + itgt;
      fNLLVal += (*it)->Target(itgt);
      double expF = exp((*it)->Target(itgt));
      nmcalt += expF;
      
      //dL[iel] += 1.0;
    }

  }
  
  int infunc = fFullFuncs.getSize()-1;
  fNLLVal += static_cast<RooAbsReal*>(fFullFuncs.at(infunc))->getVal() - fN0Obs*log(static_cast<RooAbsReal*>(fFullFuncs.at(infunc))->getVal());  
  
  
}


void RooHybridBDTAutoPdf::FitResponses(HybridGBRForest *forest) {

  printf("fit responses\n");
  //int tgtidx = itree%(fNTargets);

  
  //RooRealVar &minosarg = *fR;
  //bool doconstraint = findError && std::abs(fNLLVal-threshold)<10.;
  
  //int nterm = tree.Responses()[0].size();

  int nextvars = fExtVars.getSize();
  
  //printf("build indexes\n");
  int nparms = nextvars;
  std::vector<int> localidxs(fNTargets);
  for (int itgt=0; itgt<fNTargets; ++itgt) {
    localidxs[itgt] = nparms;
    int treetgt = static_cast<RooGBRTarget*>(fTgtVars.at(itgt))->Index();
    nparms += forest->Trees()[treetgt].back().Responses().size();
  }
  //printf("done build indexes\n");
  
//  double dr = 0.;
  
  int idxglobal = 0;
  //int idxn0 = idxglobal + fExtVars.index(fN0);
  
  //printf("idxn0 = %i, fN0 = %5f\n", idxn0,fN0->getVal());

  double nmc = 0.;
  double nmcalt = 0.;
  //double dldr = 0.;
  
  
      
  int msize = nparms;
  
  //bool usematrix = true;
  bool usematrix = msize<8000;
  //bool usematrix = false;
  
  //double lambda = 1.0;
//    double nll = 0.;
  //fNLLVal = 0.; 

  printf("allocating Hessian with size %i\n",msize);
  
  TMatrixDSym d2L;
  if (usematrix) {
    d2L.ResizeTo(msize,msize);
  }
  //std::map<std::pair<int,int>,double> d2Lmap;
  TVectorD dL(msize);

  RooArgSet parmset(fParmVars);
	
  //printf("begin loop\n");
  
  //std::vector<TMatrixDSym> d2Ls(fNThreads,TMatrixDSym(msize));
  //std::vector<std::map<std::pair<int,int>,double> > d2Lmaps(fNThreads);
  std::vector<TVectorD> dLs(fNThreads,TVectorD(msize));
  
  
  std::vector<double> nmcs(fNThreads);
  
  //printf("FitResponses start loop\n");
  
  #pragma omp parallel for
  for (unsigned int iev=0; iev<fEvts.size(); ++iev) {

    int ithread = omp_get_thread_num();
      
//    int termidx = fEvts[iev]->CurrentNode();
//    int idxlocal = nextvars + fNTargets*termidx;
    
    int evcls = fEvts[iev]->Class(); 
    double weight = fEvts[iev]->Weight();
    
    //double pdfval = static_cast<RooAbsReal*>(fStaticPdfs.at(evcls))->getValV(&parmset);
    double pdfval = fEvts[iev]->PdfVal();
    
    //double invpdf = 1.0/pdfval;
    double invpdf = vdt::fast_inv(pdfval);
    //double invpdf = vdt::fast_invf(pdfval);
    double invpdfsq = invpdf*invpdf;
    
    //double fval = fEvts[iev]->Target(0);
    //double flim = 0. + 0.5*(1.0-0.)*(sin(fval)+1.0);
    //double flim = 0. + 0.5*(1.0-0.)*(atan(fval)+1.0)*2.0/TMath::Pi();
    
    //double flim = TMath::Pi()/2.0 + atan(fval)/TMath::Pi();
    //double flim = 0.5 + atan(fval)/TMath::Pi();
    
    //double flim = 0. + 0.5*(1.0-0.)*(tanh(fval)+1.0);
    //double flim = 0.5 + atan(fval)/TMath::Pi();    
    //nmcs[ithread] += 1.0 - fval;
    //nmcs[ithread] += weight*1.0/(1.0+exp(fval));
    //nmcs[ithread] += weight*(1.0-flim);
    //nmcs[ithread] += weight*(1.0-fval);
      
    for (unsigned int iidx=0; iidx<fOuterIndices[evcls].size(); ++iidx) {
	
      int ivar = fOuterIndices[evcls][iidx];
      int idrv = ivar;
      int itgt = ivar - fExtVars.getSize(); 
      int iel;
      if (itgt>=0) iel = localidxs[itgt] + fEvts[iev]->CurrentNode(itgt);
      else iel = idxglobal + ivar;

      
      double drvi = fEvts[iev]->Derivative(idrv);
      
      dLs[ithread][iel] += -weight*drvi*invpdf;
      
      if (!usematrix) continue;
      
      //double drv2i = fEvts[iev]->Derivative2(idrv);
      //d2Lmaps[ithread][std::pair<int,int>(iel,iel)] += -weight*drv2i*invpdf;
      
      //d2Ls[ithread][iel][iel] += -weight*drv2i*invpdf + weight*drvi*drvi*invpdfsq;
      //d2Ls[ithread][iel][iel] += -weight*drv2i*invpdf;
      

      for (unsigned int jidx=iidx; jidx<fOuterIndices[evcls].size(); ++jidx) {
      
	
	int jvar = fOuterIndices[evcls][jidx];
	int jdrv = jvar;
	int jtgt = jvar - fExtVars.getSize();
	int jel;      
	if (jtgt>=0) jel = localidxs[jtgt] + fEvts[iev]->CurrentNode(jtgt);
	else jel = idxglobal + jvar;    
	
	
	
	double drvj = fEvts[iev]->Derivative(jdrv);
        //double drv2ij = fEvts[iev]->ParmMatrix()(idrv,jdrv);
        
        //double updval = d2L(iel,jel) + weight*drvi*drvj*invpdfsq;
        
        #pragma omp atomic
        d2L(iel,jel) += weight*drvi*drvj*invpdfsq;
        
        //d2Lmaps[ithread][std::pair<int,int>(iel,jel)] += weight*drvi*drvj*invpdfsq;
	//d2Ls[ithread][iel][jel] += weight*drvi*drvj*invpdfsq;
        //d2Ls[ithread][iel][jel] +=  -weight*drv2ij*invpdf + weight*drvi*drvj*invpdfsq;

      
      }
	    
    }

  }
  
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    dL += dLs[ithread];
    //d2L += d2Ls[ithread];
    nmc += nmcs[ithread];
    
    if (!usematrix) continue;

    
//     for (std::map<std::pair<int, int>, double>::const_iterator it=d2Lmaps[ithread].begin(); it!=d2Lmaps[ithread].end(); ++it) {
//       d2Lmap[it->first] += it->second;
//     }
  }

  
  int infunc = fFullFuncs.getSize()-1;
   
  for (unsigned int iidx=0; iidx<fOuterIndices[infunc].size(); ++iidx) {
	int ivar = fOuterIndices[infunc][iidx];
	int idrv = ivar;      
	int iel = idxglobal + ivar;

	double drvi = Derivative1(static_cast<RooAbsReal*>(fFullFuncs.at(infunc)),static_cast<RooRealVar*>(fFullParms.at(idrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(idrv))->getError());
	
	dL[iel] += drvi;
      
        if (!usematrix) continue;

	
    for (unsigned int jidx=iidx; jidx<fOuterIndices[infunc].size(); ++jidx) {
	
	int jvar = fOuterIndices[infunc][jidx];
	int jdrv = jvar; 
	int jel = idxglobal + jvar;
	
	//double drvj;
	double drv2ij;
	if (idrv==jdrv) {
	  //drvj = drvi;
	  drv2ij = Derivative2(static_cast<RooAbsReal*>(fFullFuncs.at(infunc)),static_cast<RooRealVar*>(fFullParms.at(idrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(idrv))->getError());
	}
	else {
	  //drvj = Derivative1(static_cast<RooAbsReal*>(fFullFuncs.at(infunc)),static_cast<RooRealVar*>(fFullParms.at(jdrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(jdrv))->getError());
	  drv2ij = Derivative2(static_cast<RooAbsReal*>(fFullFuncs.at(infunc)),static_cast<RooRealVar*>(fFullParms.at(idrv)),static_cast<RooRealVar*>(fFullParms.at(jdrv)),&parmset,1e-3*static_cast<RooRealVar*>(fFullParms.at(idrv))->getError(),1e-3*static_cast<RooRealVar*>(fFullParms.at(jdrv))->getError());
	}
	
	d2L[iel][jel] += drv2ij;
	//d2Lmap[std::pair<int,int>(iel,jel)] += drv2ij;
	
    }
	
  }
  
 // printf("FitRespones end loop\n");
  
//  double sumweight = 0.;
  
  //symmetrize matrix
  if (usematrix) {
    #pragma omp parallel for
    for (int i=0; i<msize; ++i) {
      for (int j=0; j<i; ++j) {
	assert(d2L(i,j)==0.);
	d2L(i,j) = d2L(j,i);
      }
    }
  }

  bool solved = false;
  TVectorD dpar(msize);    
  
//solve
  if (usematrix) {
    TVectorD dLc(dL);
    printf("start matrix decomposition\n");  
    //TDecompLU dbk(d2L);
    TDecompBK dbk(d2L);
    //TDecompChol dbk(d2L);
    solved = dbk.Solve(dLc);
    printf("done matrix decomposition\n");  
    dpar = -1.0*dLc;
  }

  

//   bool solved = false;
//   TVectorD dpar(msize);  
  
//   double maxscale = 0.;
//   for (int iel=0; iel<msize; ++iel) {
//     int iparm;
//     if (iel<fExtVars.getSize()) {
//       iparm = iel;
//     }
//     else {
//       iparm = (iel-fExtVars.getSize())%fNTargets;
//     }
//     double scale = dL[iel]/fStepSizes[iparm];
//     if (scale>maxscale) {
//       maxscale = scale;
//     }
//   }
//   
//   double nlldrv=0;
//   for (int iel=0; iel<msize; ++iel) {
//     nlldrv += dL[iel];
//   }
//   
//   
//   double drvstep = 1.0/maxscale;
//   
//   double upnllval = EvalLoss(forest,drvstep,dL);
//   double downnllval = EvalLoss(forest,-drvstep,dL);
//   
//   double nlldrv2 = (upnllval + downnllval - 2.0*fNLLVal)/drvstep/drvstep;
//   
//   dpar = (-nlldrv*nlldrv/nlldrv2)*dL;
//   solved = true;
  
  
  
  
  
  
//   if (usematrix) {
//     //symmetrize matrix
//     for (std::map<std::pair<int, int>, double>::const_iterator it=d2Lmap.begin(); it!=d2Lmap.end(); ++it) {
//       d2Lmap[std::pair<int,int>(it->first.second,it->first.first)] = it->second;
//     }
//     
//     printf("create sparse matrix, msize = %i, nnzr = %i\n",msize,int(d2Lmap.size()));;
//     
//     
//     sparserows.resize(d2Lmap.size());
//     sparsecols.resize(d2Lmap.size());
//     sparsedata.resize(d2Lmap.size());
//     
//     {
//       int isparse = 0;
//       for (std::map<std::pair<int, int>, double>::const_iterator it=d2Lmap.begin(); it!=d2Lmap.end(); ++it, ++isparse) {
// 	sparserows[isparse] = it->first.first;
// 	sparsecols[isparse] = it->first.second;
// 	sparsedata[isparse] = it->second;
//       }
//     }
// 
//     TMatrixDSparse d2L(msize,msize);
//     d2L.SetMatrixArray(sparsedata.size(), &sparserows[0], &sparsecols[0], &sparsedata[0]);
//   
//     printf("start matrix decomposition\n");
//     TVectorD dLc(dL);
//     TDecompSparse dsp(d2L,0);
//     solved = dsp.Solve(dLc);
//     dpar = -1.0*dLc;
//     printf("done matrix decomposition\n");
//     
//     double drv=0.;
//     for (int ipar=0; ipar<msize; ++ipar) {
//       drv += dpar[ipar]*dL[ipar];
//       if (std::isnan(dL[ipar]) || std::isinf(dL[ipar])) {
// 	solved = false;
//       }
//     }
//     
//     if (drv>=0) {
//       solved = false;
//     }
//     
//   }
  
 // printf("FitRespones done invert\n");

  
  double step = fShrinkage;

  double drvstep = 1e-3;
  if (!solved) {
 // if (1) {  
    printf("fallback to gradient descent\n");
    dpar = -1.0*dL;
    drvstep = 1e-9;
  }
  
//   double maxscale = 0.;
//   for (int iel=0; iel<msize; ++iel) {
//     int iparm;
//     if (iel<fExtVars.getSize()) {
//       iparm = iel;
//     }
//     else {
//       //iel = localidxs[itgt] + termidx
//       //iparm = (iel-fExtVars.getSize())%fNTargets;
//       int tgt = 0;
//       for (int itgt=0; itgt<fNTargets; ++itgt) {
//         if (iel>=localidxs[itgt]) {
//           tgt = itgt;
//           break;
//         }        
//       }
//       iparm = fExtVars.getSize() + tgt;
//     }
//     double scale = dpar[iel]/fStepSizes[iparm];
//     if (scale>maxscale) {
//       maxscale = scale;
//     }
//   }
  
//   double nlldrv=0;
//   for (int iel=0; iel<msize; ++iel) {
//     nlldrv += dpar[iel]*dL[iel];
//   }
  
  
  //double drvstep = 1.0/maxscale;
  if (fShrinkage<0.85)
  {
  
    double upnllval = EvalLoss(forest,drvstep,dpar);
    double downnllval = EvalLoss(forest,-drvstep,dpar);
    
    double nlldrv = (upnllval - downnllval)/(2.0*drvstep);
    double nlldrv2 = (upnllval + downnllval - 2.0*fNLLVal)/drvstep/drvstep;
    
    double maxscale = 0.;
    step = -fShrinkage*nlldrv/nlldrv2;
    
    //step = -0./0.;
    
    printf("upnllval = %5f, downnllval = %5f, fNLLVal = %5f, maxscale = %5f, drvstep = %5f, nlldrv = %5f, nlldrv2 = %5f, step = %5f\n",upnllval,downnllval,fNLLVal,maxscale,drvstep,nlldrv,nlldrv2,step);
    
    //if (nlldrv>=0. || nlldrv2<=0. || !std::isfinite(step) ) step = 0.1*fShrinkage;
    if (nlldrv>=0. || nlldrv2<=0. || !std::isnormal(step) || step<(0.1*fShrinkage) ) step = 0.1*fShrinkage;
    //if (nlldrv>=0. || nlldrv2<=0. || std::isinf(step) || std::isnan(step) ) step = 0.1*fShrinkage;
        
    
    
  }
  else {
    step = fShrinkage;
  }
   
//   if (!solved) {
//     printf("fallback to gradient descent\n");
//     
//     dpar = -1.0*dL;
//     drvstep = 1e-9;
//     
//     double upnllval = EvalLoss(forest,drvstep,dpar);
//     double downnllval = EvalLoss(forest,-drvstep,dpar);
//     
//     double nlldrv = (upnllval - downnllval)/(2.0*drvstep);
//     double nlldrv2 = (upnllval + downnllval - 2.0*fNLLVal)/drvstep/drvstep;
//     
//     double maxscale = 0.;
//     step = -fShrinkage*nlldrv/nlldrv2;
//     
//     printf("upnllval = %5f, downnllval = %5f, fNLLVal = %5f, maxscale = %5f, drvstep = %5f, nlldrv = %5f, nlldrv2 = %5f, step = %5f\n",upnllval,downnllval,fNLLVal,maxscale,drvstep,nlldrv,nlldrv2,step);    
//     
//   }
  //step = fShrinkage;
  
//   for (int iel=0; iel<msize; ++iel) {
//     dpar[iel] = -dL[iel]/d2L(iel,iel);
//   }

  double nllval = fNLLVal;
  int stepiter = 0;
  do {

//     if (stepiter==1) {
//       printf("fallback to gradient descent\n");
//       dpar = -1.0*dL;
//       step = 1e-8;      
//     }
    
    nllval = EvalLoss(forest,step,dpar);
    printf("step = %5f, nllval = %5f, fNLLVal = %5f\n",step,nllval,fNLLVal);
    if ( (nllval-fNLLVal)<1e-3 ) {
      break; 
    }
    else {
      step /= 2.0;
    }
    ++stepiter;
  }
  while (stepiter<50);

  if (!((nllval-fNLLVal)<1e-3)) {
    step = 0.;
    nllval = EvalLoss(forest, step,dpar);
  }
  
  fNLLVal = nllval;
  
  
//   for (int itgt=0; itgt<fNTargets; ++itgt) {
//     double currentval = (*it)->Target(itgt);
//     double newval = currentval + step*dpar(ivar);
//     double minval = static_cast<RooRealVar*>(fStaticTgts.at(itgt))->getMin();
//     double maxval = static_cast<RooRealVar*>(fStaticTgts.at(itgt))->getMax();
//     if (newval<minval) {
//       dpar(ivar) = (minval-currentval)/fShrinkage;
//     }
//     if (newval>maxval) {
//       dpar(ivar) = (maxval-currentval)/fShrinkage;
//     }
//   }
  
  
  for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
    static_cast<RooRealVar*>(fExtVars.at(ivar))->setVal(static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal() + step*dpar(ivar));   
  }
  
  for (int ithread=0; ithread<fNThreads; ++ithread) {
    for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
      RooRealVar *cloneparm = static_cast<RooRealVar*>(fExtVarsClones[ithread].at(ivar));
      cloneparm->setVal(static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal());
      cloneparm->setError(static_cast<RooRealVar*>(fExtVars.at(ivar))->getError());
    }
  }    
  
  
  for (int itgt=0; itgt<fNTargets; ++itgt) {
    int treetgt = static_cast<RooGBRTarget*>(fTgtVars.at(itgt))->Index();
    int nterm = forest->Trees()[treetgt].back().Responses().size();
    for (int iterm=0; iterm<nterm; ++iterm) {
      forest->Trees()[treetgt].back().Responses()[iterm] = forest->Trees()[treetgt].back().Responses()[iterm] + step*dpar(localidxs[itgt] + iterm);
    }
  }
  
  
  for (std::vector<HybridGBREvent*>::const_iterator it = fEvts.begin(); it!=fEvts.end(); ++it) {
    for (int itgt=0; itgt<fNTargets; ++itgt) {
      int termidx = (*it)->CurrentNode(itgt);
      (*it)->SetTarget(itgt,(*it)->Target(itgt)+step*dpar(localidxs[itgt] + termidx));
    }
  }  

  double dldr = -Derivative1(fN0,fR,&parmset,1e-3*fR->getError());
  fdLdR = dldr;
  
  double etermval = fN0->getVal();
  

  
  printf("r = %5f +- %5f, dldr = %5f, dL[0] = %5e,, n0 = %5f, ndobs = %5f, fNLLVal = %5f, bnll = %5f, nmc = %5f, nmcalt = %5f\n",fR->getVal(),fR->getError(),dldr,dL[0],fN0->getVal(),fN0Obs,fNLLVal,fNLLVal-etermval,nmc,nmcalt);
    
  
  
}

void RooHybridBDTAutoPdf::RecomputeTargets() {
 
  std::vector<std::vector<float> > evalvectors(fNThreads, std::vector<float>(fCondVars.getSize()));
  
  #pragma omp parallel for
  for (unsigned int ievt=0; ievt<fEvts.size(); ++ievt) {
    int ithread =  omp_get_thread_num();
    
    for (int ivar=0; ivar<fCondVars.getSize(); ++ivar) {
      evalvectors[ithread][ivar] = fEvts.at(ievt)->Var(ivar);
    }
    
    for (int itgt=0; itgt<fStaticTgts.getSize(); ++itgt) {
      int treetgt = static_cast<RooGBRTarget*>(fTgtVars.at(itgt))->Index();
      fEvts[ievt]->SetTarget(itgt, fFunc->Forest()->GetResponse(&evalvectors[ithread][0],treetgt));
    }
    
  }
   
  
}

void RooHybridBDTAutoPdf::GradientMinos() {
 
//     SetMinCutSignificance(1.0);
//     TrainForest(1e6,true);
//     return; 
  
    //save initial state
    HybridGBRForest *origforest = new HybridGBRForest(*fFunc->Forest());
    std::vector<double> extvals(fExtVars.getSize());
    for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
      extvals[ivar] = static_cast<RooRealVar*>(fExtVars.at(ivar))->getVal();
    }

  
    double rval = fR->getVal(); 
    double rerr = fR->getError();
        
    //double errscale = 1e-2;
    double errscale = 1e-3;
    
    int nstep = 3;
    double stepsize = 1.2*rerr*(1.0+errscale*errscale)/(double)nstep;
    
    std::vector<std::pair<double,double> > drvvals;
    drvvals.push_back(std::pair<double,double>(rval,0.));
    
    double coeffval = 0.5/pow(errscale*rerr,2);
    double errval = errscale*rerr;
    
    fConstraintCoeff->setVal(coeffval);  
    fR->setError(errval);     
    
    //upper uncertainty
    double rstep = rval;
    //fConstraintVal->setVal(rstep);
    //fConstraintCoeff->setVal(0.);  
    //TrainForest(0,false);
    //fConstraintCoeff->setVal(coeffval);  
    //fR->setError(errval);
    //TrainForest(1e6,true);
/*    TrainForest(1e6,false);
    double drvval = -2.0*fConstraintCoeff->getVal()*(fR->getVal()-fConstraintVal->getVal());
    drvvals.push_back(std::pair<double,double>(fR->getVal(),drvval));*/    
    
    for (int istep=0; istep<nstep; ++istep) {
      rstep += stepsize;
      fConstraintVal->setVal(rstep);
      //fR->setVal(rstep);
//       fConstraintCoeff->setVal(0.);  
//       TrainForest(0,false);
//       fConstraintCoeff->setVal(coeffval);  
//       fR->setError(errval); 
//       TrainForest(1e6,true);
      TrainForest(1e6,false);
      double drvval = -2.0*fConstraintCoeff->getVal()*(fR->getVal()-fConstraintVal->getVal());
      drvvals.push_back(std::pair<double,double>(fR->getVal(),drvval));
    }   
    
//     //restore initial state
//     fFunc->SetForest(new HybridGBRForest(*origforest));
//     for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
//       static_cast<RooRealVar*>(fExtVars.at(ivar))->setVal(extvals[ivar]);
//     }
//     RecomputeTargets();
    
 
    
    //lower uncertainty
    rstep = rval;
    for (int istep=0; istep<nstep; ++istep) {
      rstep -= stepsize;
      fConstraintVal->setVal(rstep);
      //fR->setVal(rstep);
//       fConstraintCoeff->setVal(0.);  
//       TrainForest(0,false);
//       fConstraintCoeff->setVal(coeffval);  
//       fR->setError(errval);
//       TrainForest(1e6,true);
      TrainForest(1e6,false);
      double drvval = -2.0*fConstraintCoeff->getVal()*(fR->getVal()-fConstraintVal->getVal());
      drvvals.push_back(std::pair<double,double>(fR->getVal(),drvval));
    }    
    
    std::sort(drvvals.begin(),drvvals.end());
    
    if (fDrvGraph) delete fDrvGraph;
    fDrvGraph = new TGraph(drvvals.size());		
    for (unsigned int istep=0; istep<drvvals.size(); ++istep) {
      fDrvGraph->SetPoint(istep,drvvals[istep].first,2.0*drvvals[istep].second);
    }
    
    int nstepint = 1001;
    double stepsizeint = (drvvals.back().first-drvvals.front().first)/(double)(nstepint-1);
    rstep = drvvals.front().first;
    if (fDrvGraphSmooth) delete fDrvGraphSmooth;
    fDrvGraphSmooth = new TGraph(nstepint);
   
    TGraph *intgraph = new TGraph(nstepint);
    double intval = 0.; 
    
    for (int istep = 0; istep<nstepint; ++istep) {

      double drvval = fDrvGraph->Eval(rstep,0,"S");
      intval += drvval*stepsizeint;
      
      double rstepint = rstep + 0.5*stepsizeint;
      
      fDrvGraphSmooth->SetPoint(istep,rstep,drvval);
      intgraph->SetPoint(istep,rstepint,intval);
      
      rstep += stepsizeint;
    }
    
    double minval = std::numeric_limits<double>::max();
    double minr = 0.;
    int minstep = 0;
    for (int istep=0; istep<nstepint; ++istep) {
      double val;
      intgraph->GetPoint(istep,rstep,val);
      if (val<minval) {
	minval = val;
	minr = rstep;
	minstep = istep;
      }
    }
    
    double rlow = 0.;
    double rhigh = 0.;
    
    for (int istep=minstep; istep>=0; --istep) {
      double val;
      intgraph->GetPoint(istep,rstep,val);
      if ( (val-minval)>1.0 ) {
	rlow = rstep;
	break;
      }
    }

    for (int istep=minstep; istep<nstepint; ++istep) {
      double val;
      intgraph->GetPoint(istep,rstep,val);
      if ( (val-minval)>1.0 ) {
	rhigh = rstep;
	break;
      }
    }
    
    if (fGraphDelta) delete fGraphDelta;
    fGraphDelta = new TGraph(nstepint);			
    for (int istep=0; istep<nstepint; ++istep) {
      double val;
      intgraph->GetPoint(istep,rstep,val);
      fGraphDelta->SetPoint(istep,rstep,(val-minval));
    }
    delete intgraph;
    
    fRMin = minr;
    fRHigh = rhigh;
    fRLow = rlow;
    
    printf("rfit = %5f, rerr = %5f, minr = %5f, rlow = %5f, rhigh = %5f, sigmu = %5f\n",rval,rerr,minr,rlow,rhigh,0.5*(rhigh-rlow));
    printf("nomfit: r = %5f +%5f -%5f\n",rval,rhigh-rval,rval-rlow);
    printf("minr  : r = %5f +%5f -%5f\n",minr,rhigh-minr,minr-rlow);      
    
   FILE * fp;
   fp = fopen ("confint.txt", "w");
   fprintf(fp,"rfit = %5f, rerr = %5f, minr = %5f, rlow = %5f, rhigh = %5f, sigmu = %5f,\n",rval,rerr,minr,rlow,rhigh,0.5*(rhigh-rlow));
   fprintf(fp,"nomfit: r = %5f +%5f -%5f\n",rval,rhigh-rval,rval-rlow);
   fprintf(fp,"minr  : r = %5f +%5f -%5f\n",minr,rhigh-minr,minr-rlow);    
   fclose(fp);     
    
    //restore initial state
    fFunc->SetForest(new HybridGBRForest(*origforest));
    delete origforest;
    for (int ivar=0; ivar<fExtVars.getSize(); ++ivar) {
      static_cast<RooRealVar*>(fExtVars.at(ivar))->setVal(extvals[ivar]);
    }
    fR->setError(rerr);
    
    return;
  
}

void RooHybridBDTAutoPdf::fitWithMinosFast() {
  TrainForest(fNTrees);
  //TrainForest(fNTrees);
  
  //return;

  
  
  //fLambda->setVal(1.0);
  fLambda->setVal(1.);
//  double thres = fNLLVal + 0.5;
  //fNLLVal = 0.;
  //double thres = -245e3;
//  TrainForest(10*fNTrees,true,thres,fR->getVal(),fR->getError());
  //TrainForest(1,true,fNLLVal+10e3,fR->getVal(),fR->getError());
  fLambda->setVal(0.);
}


void RooHybridBDTAutoPdf::fitWithMinos() {
  
    RooRealVar &r = *fR;
    bool isconstant = r.getAttribute("Constant");
    
    
    TrainForest(fNTrees);
    double rerr = r.getError();
  
    //RooAbsReal *nll = this->createNLL(data,RooFit::ConditionalObservables(fCondVars),RooFit::Optimize(false),RooFit::Extended());  

    double r0 = r.getVal();
    double rMax = r.getMax();
    double rMin = r.getMin();
    
    r.setConstant();
    
//    int ndim = 1;
    //double delta68 = 0.5*ROOT::Math::chisquared_quantile_c(1-0.68,ndim);
    double delta68 = 0.5;
    double nll0 = fNLLVal;
    double threshold68 = nll0 + delta68; 
    
    double hi68 = findCrossing(r, threshold68, r0,   rMax, rerr);
    
    fNLLVal = nll0;
    //double hi95 = do95_ ? findCrossing(minim2, *nll, r, threshold95, std::isnan(hi68) ? r0 : hi68, std::max(rMax, std::isnan(hi68*2-r0) ? r0 : hi68*2-r0)) : r0;
    // low error 
    double lo68 = findCrossing(r, threshold68, r0,   rMin, rerr); 
    //double lo95 = do95_ ? findCrossing(minim2, *nll, r, threshold95, std::isnan(lo68) ? r0 : lo68, rMin) : r0;

    fNLLVal = nll0;
    
    r.setVal(r0);
    r.setError(rerr);
    r.setAsymError(!std::isnan(lo68) ? lo68 - r0 : 0, !std::isnan(hi68) ? hi68 - r0 : 0);    
    
    r.setConstant(isconstant);
    
    
    
    printf("r = %5f +%5f -%5f\n",r.getVal(),std::abs(r.getAsymErrorHi()),std::abs(r.getAsymErrorLo()));
    
    //delete nll;
}
 
double RooHybridBDTAutoPdf::findCrossing(RooRealVar &r, double level, double rStart, double rBound, double rerr) {

    r.setVal(rStart);
    //double rInc = 0.1*(rBound - rStart);
    double rInc = rerr*(rBound-rStart)/std::abs(rBound-rStart);
    
    int verbose = 9;
    
    double here = fNLLVal;
    do {
        rStart += rInc;
        if (rInc*(rStart - rBound) > 0) { // went beyond bounds
            rStart -= rInc;
            rInc    = 0.5*(rBound-rStart);
        }
        r.setVal(rStart);
        //nll.clearEvalErrorLog(); nll.getVal();
        TrainForest(fNTrees);
        double there = here;
        here = fNLLVal;
        if (verbose > 0) { printf("%f    %+.5f  %+.5f    %f\n", rStart, level-here, level-there, rInc); fflush(stdout); }
        if ( fabs(here - level) < 4*1e-4 ) {
            // set to the right point with interpolation
            r.setVal(rStart + (level-here)*(level-there)/(here-there));
            return r.getVal();
        } else if (here > level) {
            // I'm above the level that I wanted, this means I stepped too long
            // First I take back all the step
            rStart -= rInc; 
            // Then I try to figure out a better step
            if (1) {
                if (fabs(there - level) > 0.05) { // If started from far away, I still have to step carefully
                    double rIncFactor = std::max(0.2, std::min(0.7, 0.75*(level-there)/(here-there)));
                    //printf("\t\t\t\t\tCase A1: level-there = %f,  here-there = %f,   rInc(Old) = %f,  rInFactor = %f,  rInc(New) = %f\n", level-there, here-there, rInc, rIncFactor, rInc*rIncFactor);
                    rInc *= rIncFactor;
                } else { // close enough to jump straight to target
                    double rIncFactor = std::max(0.05, std::min(0.95, 0.95*(level-there)/(here-there)));
                    //printf("\t\t\t\t\tCase A2: level-there = %f,  here-there = %f,   rInc(Old) = %f,  rInFactor = %f,  rInc(New) = %f\n", level-there, here-there, rInc, rIncFactor, rInc*rIncFactor);
                    rInc *= rIncFactor;
                }
            } else {
                rInc *= 0.3;
            }
            //if (allpars.get() == 0) allpars.reset(nll.getParameters((const RooArgSet *)0));
            //RooArgSet oldparams(checkpoint->floatParsFinal());
            //*allpars = oldparams;
        } else if ((here-there)*(level-there) < 0 && // went wrong
                   fabs(here-there) > 0.1) {         // by more than roundoff
//            if (allpars.get() == 0) allpars.reset(nll.getParameters((const RooArgSet *)0));
            //RooArgSet oldparams(checkpoint->floatParsFinal());
            //*allpars = oldparams;
            rStart -= rInc; rInc *= 0.5;
        } else {
            // I did this step, and I'm not there yet
            if (1) {
                if (fabs(here - level) > 0.05) { // we still have to step carefully
                    if ((here-there)*(level-there) > 0) { // if we went in the right direction
                        // then optimize step size
                        double rIncFactor = std::max(0.2, std::min(2.0, 0.75*(level-there)/(here-there)));
                        //printf("\t\t\t\t\tCase B1: level-there = %f,  here-there = %f,   rInc(Old) = %f,  rInFactor = %f,  rInc(New) = %f\n", level-there, here-there, rInc, rIncFactor, rInc*rIncFactor);
                        rInc *= rIncFactor;
                    } //else printf("\t\t\t\t\tCase B3: level-there = %f,  here-there = %f,   rInc(Old) = %f\n", level-there, here-there, rInc);
                } else { // close enough to jump straight to target
                    double rIncFactor = std::max(0.05, std::min(4.0, 0.95*(level-there)/(here-there)));
                    //printf("\t\t\t\t\tCase B2: level-there = %f,  here-there = %f,   rInc(Old) = %f,  rInFactor = %f,  rInc(New) = %f\n", level-there, here-there, rInc, rIncFactor, rInc*rIncFactor);
                    rInc *= rIncFactor;
                }
            } else {
                //nothing?
            }
            //checkpoint.reset(minim.save());
        }
    } while (fabs(rInc) > 1e-4*0.1*std::max(1.0,rBound-rStart));
    return r.getVal();
//     if (fabs(here - level) > 0.01) {
//         std::cout << "Error: closed range without finding crossing." << std::endl;
//         return NAN;
//     } else {
//         return r.getVal();
//     } 
  
}



double RooHybridBDTAutoPdf::Derivative1Fast(RooAbsReal *function, double currentval, RooRealVar *var, RooArgSet *nset, double step) {
 
  //printf("var = %s, step = %5e\n",var->GetName(),step);
  
  double startval = var->getVal();
  
  var->setVal(startval+step);
  double valup = function->getValV(nset);

  var->setVal(startval-step);
  double valdown = function->getValV(nset);  
  
  double drv = (valup-valdown)*vdt::fast_inv(2.0*step);
  
  var->setVal(startval);
  
  return drv;
  
}

double RooHybridBDTAutoPdf::Derivative2Fast(RooAbsReal *function, double currentval, RooRealVar *var, RooArgSet *nset, double step) {
 
  double startval = var->getVal();
  
  double valnom = currentval;
  
  var->setVal(startval+step);
  double valup1 = function->getValV(nset);
  
  var->setVal(startval-step);
  double valdown1 = function->getValV(nset);  

  double drv1 = (valup1+valdown1-2.0*valnom)*vdt::fast_inv(step*step);
  
  var->setVal(startval);
  
  return drv1;
  
}

double RooHybridBDTAutoPdf::Derivative2Fast(RooAbsReal *function, RooRealVar *var1, RooRealVar *var2, RooArgSet *nset, double step1, double step2) {

  double startval1 = var1->getVal();
  double startval2 = var2->getVal();
    
  //double valnom = function->getValV(nset);
  
  double stepa1 = step1;
  double stepb1 = step2;
  
  var1->setVal(startval1+stepa1);
  var2->setVal(startval2+stepb1);
  double valupup1 = function->getValV(nset);

  var1->setVal(startval1+stepa1);
  var2->setVal(startval2-stepb1);
  double valupdown1 = function->getValV(nset);  
  
  var1->setVal(startval1-stepa1);
  var2->setVal(startval2+stepb1);
  double valdownup1 = function->getValV(nset);    
  
  var1->setVal(startval1-stepa1);
  var2->setVal(startval2-stepb1);
  double valdowndown1 = function->getValV(nset);    
  
  double drv1 = (valupup1+valdowndown1-valupdown1-valdownup1)/(4.0*stepa1*stepb1);
  
  var1->setVal(startval1);
  var2->setVal(startval2);
  
  return drv1;
  
}



double RooHybridBDTAutoPdf::Derivative1(RooAbsReal *function, RooRealVar *var, RooArgSet *nset, double step) {
 
  //printf("var = %s, step = %5e\n",var->GetName(),step);
  
  double startval = var->getVal();
  
  var->setVal(startval+step);
  double valup1 = function->getValV(nset);
  
  var->setVal(startval-step);
  double valdown1 = function->getValV(nset);  

  double drv1 = (valup1-valdown1)/(2.0*step);
  
  double step2 = 0.5*step;
  
  var->setVal(startval+step2);
  double valup2 = function->getValV(nset);
  
  var->setVal(startval-step2);
  double valdown2 = function->getValV(nset);    
  
  double drv2 = (valup2-valdown2)/(2.0*step2);
  
  double drv = (4.0*drv2 - drv1)/3.0;
  
  var->setVal(startval);
  
  return drv;
  
}

double RooHybridBDTAutoPdf::Derivative2(RooAbsReal *function, RooRealVar *var, RooArgSet *nset, double step) {
 
  double startval = var->getVal();
  
  double valnom = function->getValV(nset);
  
  var->setVal(startval+step);
  double valup1 = function->getValV(nset);
  
  var->setVal(startval-step);
  double valdown1 = function->getValV(nset);  

  double drv1 = (valup1+valdown1-2.0*valnom)/(step*step);
  
  double step2 = 0.5*step;  
  
  var->setVal(startval+step2);
  double valup2 = function->getValV(nset);
  
  var->setVal(startval-step2);
  double valdown2 = function->getValV(nset);    
  
  double drv2 = (valup2+valdown2-2.0*valnom)/(step2*step2);
  
  double drv = (4.0*drv2 - drv1)/3.0;
  
  var->setVal(startval);
  
  return drv;
  
}

double RooHybridBDTAutoPdf::Derivative2(RooAbsReal *function, RooRealVar *var1, RooRealVar *var2, RooArgSet *nset, double stepa, double stepb) {
 
  //double step = std::min(stepa,stepb);
  
  double startval1 = var1->getVal();
  double startval2 = var2->getVal();
    
  //double valnom = function->getValV(nset);
  
  double stepa1 = stepa;
  double stepb1 = stepb;
  
  var1->setVal(startval1+stepa1);
  var2->setVal(startval2+stepb1);
  double valupup1 = function->getValV(nset);

  var1->setVal(startval1+stepa1);
  var2->setVal(startval2-stepb1);
  double valupdown1 = function->getValV(nset);  
  
  var1->setVal(startval1-stepa1);
  var2->setVal(startval2+stepb1);
  double valdownup1 = function->getValV(nset);    
  
  var1->setVal(startval1-stepa1);
  var2->setVal(startval2-stepb1);
  double valdowndown1 = function->getValV(nset);    
  
  double drv1 = (valupup1+valdowndown1-valupdown1-valdownup1)/(4.0*stepa1*stepb1);
  
  
      
  
  double stepa2 = 0.5*stepa;
  double stepb2 = 0.5*stepb;  
  
  var1->setVal(startval1+stepa2);
  var2->setVal(startval2+stepb2);
  double valupup2 = function->getValV(nset);

  var1->setVal(startval1+stepa2);
  var2->setVal(startval2-stepb2);
  double valupdown2 = function->getValV(nset);  
  
  var1->setVal(startval1-stepa2);
  var2->setVal(startval2+stepb2);
  double valdownup2 = function->getValV(nset);    
  
  var1->setVal(startval1-stepa2);
  var2->setVal(startval2-stepb2);
  double valdowndown2 = function->getValV(nset);    
  
  double drv2 = (valupup2+valdowndown2-valupdown2-valdownup2)/(4.0*stepa2*stepb2);
  
  double drv = (4.0*drv2 - drv1)/3.0;
  
  var1->setVal(startval1);
  var2->setVal(startval2);
  
  return drv;
  
}
