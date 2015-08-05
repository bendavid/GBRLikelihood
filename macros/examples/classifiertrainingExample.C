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
 
using namespace RooFit;
  
void classifiertrainingExample() {
     
  //build vectors with list of input variables
  std::vector<std::string> *varslist = new std::vector<std::string>;
  varslist->push_back("ph.scrawe");
  varslist->push_back("ph.sceta");
  varslist->push_back("ph.scphi");
  varslist->push_back("ph.r9");  
  varslist->push_back("ph.scetawidth");
  varslist->push_back("ph.scphiwidth");  
  varslist->push_back("ph.scnclusters");
  varslist->push_back("ph.hoveretower");
  varslist->push_back("rho");
  varslist->push_back("nVtx");  
   
  //create RooRealVars for each input variable
  RooArgList vars;
  for (unsigned int ivar=0; ivar<varslist->size(); ++ivar) {
    RooRealVar *var = new RooRealVar(TString::Format("var_%i",ivar),varslist->at(ivar).c_str(),0.);
    vars.addOwned(*var);
  }
  
  //make list of input variable RooRealVars
  RooArgList condvars(vars);
      
  //RooRealVar for event weight 
  RooRealVar weightvar("weightvar","",1.);

  TString treeloc = "RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonIDModPreselInvert/PhotonTreeWriterSingleInvert/hPhotonTreeSingle";
  TChain *tree = new TChain(treeloc);
  tree->Add("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-zllm50-v7n_noskim.root");
  
  //training selection cuts
  TCut basesel = "ph.genpt>16. && ph.isbarrel && ph.ispromptgen";
  TCut sigcut = "ph.scrawe/ph.gene > 0.95";
  TCut bkgcut = !sigcut;
  
  TCut prescale10 = "(evt%10==0)";
  TCut prescale20 = "(evt%20==0)";
  TCut prescale25 = "(evt%25==0)";
  TCut prescale50 = "(evt%50==0)";
  TCut prescale100 = "(evt%100==0)";  
  TCut prescale1000 = "(evt%1000==0)";  
  TCut evenevents = "(evt%2==0)";
  TCut oddevents = "(evt%2==1)";  


  //create RooDataSet from TChain
  weightvar.SetTitle(prescale100*evenevents*basesel*sigcut);
  RooDataSet *hdatasig = RooTreeConvert::CreateDataSet("hdatasig",tree,vars,weightvar);   

  weightvar.SetTitle(prescale100*evenevents*basesel*bkgcut);
  RooDataSet *hdatabkg = RooTreeConvert::CreateDataSet("hdatabkg",tree,vars,weightvar);     
  
  //RooRealVars corresponding to regressed parameters (in the simple classifier case this will just be the ln s/b)
  RooRealVar logsbvar("logsbvar","",0.);
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
  //L = -sum_signal w*ln(S*exp(F)/(B + w*S*exp(F))) - sum_background w*ln(B/(B + w*S*exp(F)))
  // where:
  //S is the sum of weights for the signal training sample
  //B is the sum of weights for the background training sample
  //w is the per-event weight if applicable (otherwise 1)
  //F is the regression (classifier) output such that for input variables xbar and probability density functons for signal and background p_s(xbar), p_b(xbar), then exp(F) = p_s(xbar)/p_b(xbar)
  
  RooConstVar *sumwsig = new RooConstVar("sumwsig","",hdatasig->sumEntries());
  RooConstVar *sumwbkg = new RooConstVar("sumwbkg","",hdatabkg->sumEntries());
  
  RooFormulaVar *sumwsigsb = new RooFormulaVar("sumwsigsb","","@0*exp(@1)",RooArgList(*sumwsig,*logsbtgt));
  RooFormulaVar *denom = new RooFormulaVar("denom","","1./(@0+@1)",RooArgList(*sumwbkg,*sumwsigsb));
  
  RooProduct *sigll = new RooProduct("sigll","",RooArgList(*sumwsigsb,*denom));
  RooProduct *bkgll = new RooProduct("bkgll","",RooArgList(*sumwbkg,*denom));
  
  
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
  vdata.push_back(hdatasig);
  vdata.push_back(hdatabkg);
  
  
  //temp output file
  TFile *fres = new TFile("fres.root","RECREATE");

  //run training
  //n.b if using MinCutSignificance then the training will converge and terminate on its own, so the number of trees should be set very large
  //if MinCutSignificance is not set, then the number of trees needs to be explicitly set to a finite number
  {
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff","",tgts,etermconst,r,vdata,vpdf);
    bdtpdfdiff.SetMinCutSignificance(3.);  //minimum statistical significance of estimated likelihood gain for a tree split
    bdtpdfdiff.SetMaxDepth(5); //max decision tree depth
    bdtpdfdiff.SetShrinkage(0.1); //shrinkage parameter by which the current tree response is de-weighted at each iteration ("learning-rate")
    bdtpdfdiff.SetMinWeightTotal(1000.); //minimum sum of weights on a terminal node
    bdtpdfdiff.SetMaxNodes(750); //maximum terminal nodes per tree, needed to keep matrix inversion CPU time under control at each iteration
    bdtpdfdiff.TrainForest(1e6);   //this is the number of trees
  }
     
  //create workspace and output to file
  RooWorkspace *weclass = new RooWorkspace("wclass");
  weclass->import(*sigll,RooFit::RecycleConflictNodes());
  weclass->import(*bkgll,RooFit::RecycleConflictNodes());
  
  weclass->writeToFile("weclass.root");    
  
  return;
  
  
}