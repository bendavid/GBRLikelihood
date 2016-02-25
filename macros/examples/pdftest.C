//Example macro for training a simple BDT classifier using the GBRLikelihood machinery
//J.Bendavid

#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooHybridBDTAutoPdf.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooUniform.h"
#include "TRandom.h"
#include "TGraph.h"
#include "RooAddPdf.h"
#include "RooNDKeysPdf.h"
#include "RooExtendPdf.h"
#include "RooMinimizer.h"
#include "TFile.h"
#include "TNtuple.h"
#include "HybridGBRForest.h"
#include "HybridGBRForestFlex.h"
#include "RooProduct.h"
#include "RooGenericPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "TH1D.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TLegend.h"
#include "RooRandom.h"
#include "RooAddition.h"
#include "TSystem.h"
 
using namespace RooFit;
  
void pdftest() {
     
  int nvars = 8;

  std::vector<RooRealVar*> xvars;
  std::vector<RooAbsPdf*> gpdfs;
  std::vector<RooAbsReal*> absgauss;
  
  RooArgList condvars;
  RooArgList gpdflist;
  RooArgList gauslist;
  
  for (int ivar=0; ivar<nvars; ++ivar) {
    RooRealVar *x = new RooRealVar(TString::Format("x_%i",ivar),"",0.,-5.,5.);
    RooGaussian *gpdf = new RooGaussian(TString::Format("gpdf_%i",ivar),"",*x,RooConst(0.0),RooConst(1.0));
    RooFormulaVar *absgaus = new RooFormulaVar(TString::Format("absgaus_%i",ivar),"","exp(-0.5*@0*@0)",*x);
    
    xvars.push_back(x);
    gpdfs.push_back(gpdf);
    absgauss.push_back(absgaus);
    
    condvars.add(*x);  
    gpdflist.add(*gpdf);
    gauslist.add(*absgaus);
    
  }
  
  RooProdPdf *prodpdf = new RooProdPdf("prodpdf","",gpdflist);
  RooProduct *prodgaus = new RooProduct("prodgaus","",gauslist);
  
  
  RooRealVar *x = xvars[0];
  RooAbsReal *absgaus = absgauss[0];
  
  RooUniform *updf = new RooUniform("updf","",condvars);
  
  RooDataSet *udata = updf->generate(condvars,1e5);
  RooRealVar *absgausvar = (RooRealVar*)udata->addColumn(*prodgaus);
  
  RooDataSet *gdata = new RooDataSet("gdata","",*udata->get(),Import(*udata),WeightVar(*absgausvar));
  

  
//   return;
  
  
  //RooRealVars corresponding to regressed parameters (in the simple classifier case this will just be the ln s/b)
  RooRealVar logsbvar("logsbvar","",log(gdata->sumEntries()/udata->sumEntries()));
  logsbvar.setConstant(false);
  
  //define non-parametric functions for each regressed parameter
  RooGBRFunctionFlex *logsbfunc = new RooGBRFunctionFlex("logsbfunc","");
  
  //define mapping of input variables to non-parametric functions (in this case trivial since all 4 functions depend on the same inputs, but this is not a requirement)
  RooGBRTargetFlex *logsbtgt = new RooGBRTargetFlex("logsbtgt","",*logsbfunc,logsbvar,condvars);  

  //define list of mapped functions to regress
  RooArgList tgts;
  tgts.add(*logsbtgt); 
  


  //define log-likelihoods for signal and background, which for the simple classifier case lead to the "logistic regression" likelihood used for example in TMVA gradient boosted BDT classifiers
  //the s/b ratio is appropriately rescaled to compensate for possible different training sample sizes between signal and background
  //For the likelihood defined below, the minimized loss function is
  //L = -sum_signal w*ln(S*exp(F)/(B + S*exp(F))) - sum_background w*ln(B/(B + S*exp(F)))
  // where:
  //S is the sum of weights for the signal training sample
  //B is the sum of weights for the background training sample
  //w is the per-event weight if applicable (otherwise 1)
  //F is the regression (classifier) output such that for input variables xbar and probability density functons for signal and background p_s(xbar), p_b(xbar), then exp(F) = p_s(xbar)/p_b(xbar)
  
  RooConstVar *sumwsig = new RooConstVar("sumwsig","",gdata->sumEntries());
  RooConstVar *sumwbkg = new RooConstVar("sumwbkg","",udata->sumEntries());
  
  RooFormulaVar *sb = new RooFormulaVar("sb","","exp(@0)",RooArgList(*logsbtgt));
  RooFormulaVar *bkgll = new RooFormulaVar("bkgll","","1./(1.+@0)",RooArgList(*sb));
  
  RooProduct *sigll = new RooProduct("sigll","",RooArgList(*sb,*bkgll));
  
  RooFormulaVar *regnorm = new RooFormulaVar("regnorm","","@0/@1",RooArgList(RooConst(10.),*sumwbkg));
  
  //dummy variable
  RooConstVar etermconst("etermconst","",0.);  
   
  //dummy variable
  RooRealVar r("r","",1.);
  r.setConstant();

  //define list of pdfs
  std::vector<RooAbsReal*> vpdf;
  vpdf.push_back(sigll);
  vpdf.push_back(bkgll);

  //define list of training datasets
  std::vector<RooAbsData*> vdata;
  vdata.push_back(gdata);
  vdata.push_back(udata);
  
  
  //temp output file
  TFile *fres = new TFile("fres.root","RECREATE");

  //run training
  //n.b if using MinCutSignificance then the training will converge and terminate on its own, so the number of trees should be set very large
  //if MinCutSignificance is not set, then the number of trees needs to be explicitly set to a finite number
  {
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff","",tgts,etermconst,r,vdata,vpdf);
    bdtpdfdiff.SetMinCutSignificance(0.01);  //minimum statistical significance of estimated likelihood gain for a tree split
//     bdtpdfdiff.SetMaxDepth(5); //max decision tree depth
    bdtpdfdiff.SetShrinkage(0.1); //shrinkage parameter by which the current tree response is de-weighted at each iteration ("learning-rate")
//     bdtpdfdiff.SetMinWeightTotal(30.); //minimum sum of weights on a terminal node
    bdtpdfdiff.SetMaxNodes(750); //maximum terminal nodes per tree, needed to keep matrix inversion CPU time under control at each iteration
    bdtpdfdiff.TrainForest(1e6);   //this is the number of trees
  }
  
  RooFormulaVar *sbreal = new RooFormulaVar("sbreal","","exp(@0)",RooArgList(*logsbtgt));
  
//   RooGenericPdf *sbpdf = new RooGenericPdf("sbpdf","","exp(@0)",*logsbtgt);
//   RooDataSet *sbdata = sbpdf->generate(*x,100000);
  
  RooPlot *plot = x->frame(-5.,5.,100);
//   sbdata->plotOn(plot,MarkerColor(kGreen));
//   udata->plotOn(plot,MarkerColor(kRed));
//   gdata->plotOn(plot,MarkerColor(kBlue));
  sbreal->plotOn(plot);
//   absgaus->plotOn(plot,LineColor(kRed),LineStyle(kDashed));
  plot->Draw();  
  
  
  const HybridGBRForestFlex *forest = logsbtgt->Forest();
  
  int nbins = 100;
  double low = -5.;
  double high = 5.;
  double step = (high-low)/double(nbins);
  TGraph *hfunc = new TGraph(nbins);
  std::vector<float> vals(nvars,0.);
  for (int ibin=0; ibin<nbins; ++ibin) {
    float x = low + ibin*step;
    vals[0] = x;
    double yval = forest->GetResponse(vals.data());
    yval = exp(yval);
    printf("ibin = %i, x = %5f, yval = %5f\n",ibin,x,yval);
    hfunc->SetPoint(ibin,x,yval);
  }  
  
//   TF1 *gaus = new TF1("gaus","exp(-0.5*x*x)",low,high);
  new TCanvas;
//   gaus->Draw();
//   hfunc->Draw("LPSAME");  
  hfunc->Draw("ALP");
  
  
  return;
     
  //create workspace and output to file
  RooWorkspace *weclass = new RooWorkspace("wclass");
  weclass->import(*sigll,RooFit::RecycleConflictNodes());
  weclass->import(*bkgll,RooFit::RecycleConflictNodes());
  
  weclass->writeToFile("wepdf.root");    
  
  return;
  
  
}