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
#include "RooGausDoubleExp.h"
#include "RooDoubleCBFast.h"
#include "RooGaussianFast.h"
 
using namespace RooFit;
  
void pdfregtest() {
     

  RooRealVar *x = new RooRealVar("x","",0.,-5.,5.);
  
  RooGaussian *gpdf = new RooGaussian("gpdf","",*x,RooConst(0.0),RooConst(1.0));
  
  RooFormulaVar *absgaus = new RooFormulaVar("absgaus","","exp(-0.5*@0*@0)",*x);
  
  RooArgList condvars;
  condvars.add(*x);
  
  RooUniform *updf = new RooUniform("updf","",*x);
  
  RooDataSet *udata = updf->generate(*x,1000);
  RooRealVar *absgausvar = (RooRealVar*)udata->addColumn(*absgaus);
  
  RooDataSet *gdata = new RooDataSet("gdata","",*udata->get(),Import(*udata),WeightVar(*absgausvar));
  

  
//   return;
  
  //RooRealVars corresponding to regressed parameters (sigma, mean, left tail parameter, right tail parameter)
  RooRealVar sigwidthtvar("sigwidthtvar","",0.01);
  sigwidthtvar.setConstant(false);
  
  RooRealVar sigmeantvar("sigmeantvar","",0.5);
  sigmeantvar.setConstant(false); 
  
  RooRealVar signvar("signvar","",3.);
  signvar.setConstant(false);       
  
  RooRealVar sign2var("sign2var","",3.);
  sign2var.setConstant(false);     

  //define non-parametric functions for each regressed parameter
  RooGBRFunctionFlex *sigwidthtfunc = new RooGBRFunctionFlex("sigwidthtfunc","");
  RooGBRFunctionFlex *sigmeantfunc = new RooGBRFunctionFlex("sigmeantfunc","");
  RooGBRFunctionFlex *signfunc = new RooGBRFunctionFlex("signfunc","");
  RooGBRFunctionFlex *sign2func = new RooGBRFunctionFlex("sign2func","");

  //define mapping of input variables to non-parametric functions (in this case trivial since all 4 functions depend on the same inputs, but this is not a requirement)
  RooGBRTargetFlex *sigwidtht = new RooGBRTargetFlex("sigwidtht","",*sigwidthtfunc,sigwidthtvar,condvars);  
  RooGBRTargetFlex *sigmeant = new RooGBRTargetFlex("sigmeant","",*sigmeantfunc,sigmeantvar,condvars);  
  RooGBRTargetFlex *signt = new RooGBRTargetFlex("signt","",*signfunc,signvar,condvars);  
  RooGBRTargetFlex *sign2t = new RooGBRTargetFlex("sign2t","",*sign2func,sign2var,condvars);  

  //define list of mapped functions to regress
  RooArgList tgts;
//   tgts.add(*sigwidtht);
  tgts.add(*sigmeant);
//   tgts.add(*signt);
//   tgts.add(*sign2t);  
  
  //define transformations corresponding to parameter bounds for non-parametric outputs  
  RooRealConstraint sigwidthlim("sigwidthlim","",*sigwidtht,0.0002,0.5);
  RooRealConstraint sigmeanlim("sigmeanlim","",*sigmeant,0.,1.0);
  RooRealConstraint signlim("signlim","",*signt,1.01,5000.); 
  RooRealConstraint sign2lim("sign2lim","",*sign2t,1.01,5000.); 

  //define pdf, which depends on transformed outputs (and is intended to be treated as a conditional pdf over the
  //regression inputs in this case)
  //The actual pdf below is a double crystal ball, with crossover points alpha_1 and alpha_2 set constant, but all other
  //parameters regressed
//   RooDoubleCBFast sigpdf("sigpdf","",*absgausvar,sigmeanlim,sigwidthlim,RooConst(1.),RooConst(100.),RooConst(1.),RooConst(100.));  
//   RooGausDoubleExp sigpdf("sigpdf","",*absgausvar,sigmeanlim,sigwidthlim,RooConst(1.),RooConst(1.));

  RooGaussianFast sigpdf("sigpdf","",*absgausvar,sigmeanlim,RooConst(1.));  
  
 
  
  
  //dummy variable
  RooConstVar etermconst("etermconst","",0.);  
   
  //dummy variable
  RooRealVar r("r","",1.);
  r.setConstant();

  //define list of pdfs
  std::vector<RooAbsReal*> vpdf;
//   vpdf.push_back(sigll);
//   vpdf.push_back(bkgll);
  vpdf.push_back(&sigpdf);

  //define list of training datasets
  std::vector<RooAbsData*> vdata;
//   vdata.push_back(gdata);
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
  
//   RooFormulaVar *sbreal = new RooFormulaVar("sbreal","","exp(@0)",RooArgList(*logsbtgt));
  
//   RooGenericPdf *sbpdf = new RooGenericPdf("sbpdf","","exp(@0)",*logsbtgt);
//   RooDataSet *sbdata = sbpdf->generate(*x,100000);
  
  RooPlot *plot = absgausvar->frame(0.,1.,200);
//   sbdata->plotOn(plot,MarkerColor(kGreen));
  udata->plotOn(plot);
//   gdata->plotOn(plot,MarkerColor(kBlue));
//   sbreal->plotOn(plot);
  sigpdf.plotOn(plot,ProjWData(*udata));
  plot->Draw();  
  
  
  new TCanvas;
  RooRealVar *meanvar = (RooRealVar*)udata->addColumn(sigmeanlim);
  RooPlot *meanplot = meanvar->frame(0.,1.,200);
  udata->plotOn(meanplot);
  meanplot->Draw();
  
  new TCanvas;
  RooRealVar *sigmavar = (RooRealVar*)udata->addColumn(sigwidthlim);
  RooPlot *sigmaplot = sigmavar->frame(0.,1.,200);
  udata->plotOn(sigmaplot);
  sigmaplot->Draw();  
  
  
  return;
     
  //create workspace and output to file
//   RooWorkspace *weclass = new RooWorkspace("wclass");
//   weclass->import(*sigll,RooFit::RecycleConflictNodes());
//   weclass->import(*bkgll,RooFit::RecycleConflictNodes());
//   
//   weclass->writeToFile("wepdf.root");    
  
  return;
  
  
}