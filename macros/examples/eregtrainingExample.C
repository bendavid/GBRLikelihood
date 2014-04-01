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
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
//#include "HZZ2L2QRooPdfs.h"
#include "RooDoubleCBFast.h"
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
#include "RooLinearVar.h"
#include "RooCBExp.h"
#include "RooCBFast.h"
#include "RooGaussianFast.h"

 
using namespace RooFit;
  
double getweight(TFile *file, double xsec) {
 
  TDirectory *dir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  TH1D *hallevts = (TH1D*)dir->Get("hDAllEvents");
  
  return xsec/hallevts->GetSumOfWeights();
  
}

float xsecweights[50];
float xsecweight(int procidx=0) {
  return xsecweights[procidx];
}

void initweights(TChain *chain, float *xsecs, float lumi) {
 
  TObjArray *files = chain->GetListOfFiles();
  for (int i=0; i<files->GetEntries(); ++i) {    
    TFile *file = TFile::Open(files->At(i)->GetTitle(),"READ");
    
    xsecweights[i] = getweight(file,lumi*xsecs[i]);
    
    file->Close();    
  } 
  
  chain->SetAlias("procidx","This->GetTreeNumber()");
  
}

void eregtrainingExample(bool dobarrel=true, bool doele=true) {
   
  //output dir
  TString dirname = "/data/bendavid/eregexampletest/"; 
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);  
  
  
  //build vectors with list of input variables
  std::vector<std::string> *varsf = new std::vector<std::string>;
  varsf->push_back("ph.scrawe");
  varsf->push_back("ph.sceta");
  varsf->push_back("ph.scphi");
  varsf->push_back("ph.r9");  
  varsf->push_back("ph.scetawidth");
  varsf->push_back("ph.scphiwidth");  
  varsf->push_back("ph.scnclusters");
  varsf->push_back("ph.hoveretower");
  varsf->push_back("rho");
  varsf->push_back("nVtx");  
 
  varsf->push_back("ph.etaseed-ph.sceta");
  varsf->push_back("atan2(sin(ph.phiseed-ph.scphi),cos(ph.phiseed-ph.scphi))");
  varsf->push_back("ph.eseed/ph.scrawe");
  
  varsf->push_back("ph.e3x3seed/ph.e5x5seed");
  varsf->push_back("ph.sigietaietaseed");   
  varsf->push_back("ph.sigiphiphiseed");   
  varsf->push_back("ph.covietaiphiseed");
  varsf->push_back("ph.emaxseed/ph.e5x5seed");
  varsf->push_back("ph.e2ndseed/ph.e5x5seed");
  varsf->push_back("ph.etopseed/ph.e5x5seed");
  varsf->push_back("ph.ebottomseed/ph.e5x5seed");
  varsf->push_back("ph.eleftseed/ph.e5x5seed");
  varsf->push_back("ph.erightseed/ph.e5x5seed");
  varsf->push_back("ph.e2x5maxseed/ph.e5x5seed");
  varsf->push_back("ph.e2x5topseed/ph.e5x5seed");
  varsf->push_back("ph.e2x5bottomseed/ph.e5x5seed");
  varsf->push_back("ph.e2x5leftseed/ph.e5x5seed");
  varsf->push_back("ph.e2x5rightseed/ph.e5x5seed");
  
  std::vector<std::string> *varseb = new std::vector<std::string>(*varsf);
  std::vector<std::string> *varsee = new std::vector<std::string>(*varsf);
  
  varseb->push_back("ph.e5x5seed/ph.eseed");
  
  varseb->push_back("ph.ietaseed");
  varseb->push_back("ph.iphiseed");
  varseb->push_back("(ph.ietaseed-1*abs(ph.ietaseed)/ph.ietaseed)%5");
  varseb->push_back("(ph.iphiseed-1)%2");       
  varseb->push_back("(abs(ph.ietaseed)<=25)*((ph.ietaseed-1*abs(ph.ietaseed)/ph.ietaseed)%25) + (abs(ph.ietaseed)>25)*((ph.ietaseed-26*abs(ph.ietaseed)/ph.ietaseed)%20)");
  varseb->push_back("(ph.iphiseed-1)%20"); 
  varseb->push_back("ph.etacryseed");
  varseb->push_back("ph.phicryseed");

  varsee->push_back("ph.scpse/ph.scrawe");
    
  //select appropriate input list for barrel or endcap
  std::vector<std::string> *varslist;
  if (dobarrel) varslist = varseb;
  else varslist = varsee;
  
  //create RooRealVars for each input variable
  RooArgList vars;
  for (unsigned int ivar=0; ivar<varslist->size(); ++ivar) {
    RooRealVar *var = new RooRealVar(TString::Format("var_%i",ivar),varslist->at(ivar).c_str(),0.);
    vars.addOwned(*var);
  }
  
  //make list of input variable RooRealVars
  RooArgList condvars(vars);
  
  //create RooRealVar for target
  RooRealVar *tgtvar = new RooRealVar("tgtvar","ph.gene/ph.scrawe",1.);
  if (!dobarrel) tgtvar->SetTitle("ph.gene/(ph.scrawe + ph.scpse)");  
  
  //add target to full list
  vars.addOwned(*tgtvar);
    
  //RooRealVar for event weight 
  RooRealVar weightvar("weightvar","",1.);

  //Initialize TChains with event weights if needed
  TString treeloc;
  if (doele) {
    treeloc = "RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonIDModPreselInvert/PhotonTreeWriterSingleInvert/hPhotonTreeSingle";
  }
  else {
    treeloc = "RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonIDModPresel/PhotonTreeWriterSingle/hPhotonTreeSingle";
  }

  TChain *tree;
  float xsecs[50];

      
  if (doele) {
    tree = new TChain("RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonIDModPreselInvert/PhotonTreeWriterSingleInvert/hPhotonTreeSingle");
    //tree->Add("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-zllm50-v7n_noskim.root");
    tree->Add("/data/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-zllm50-v7n_noskim.root");
    
    xsecs[0] = 1.;
    initweights(tree,xsecs,1.);      
    
    xsecweights[0] = 1.0;
    
  }
  else {
    tree = new TChain("RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonIDModPresel/PhotonTreeWriterSingle/hPhotonTreeSingle");
//     tree->Add("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-pj20_40-2em-v7n_noskim.root");
//     tree->Add("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-pj40-2em-v7n_noskim.root");
    tree->Add("/data/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-pj20_40-2em-v7n_noskim.root");
    tree->Add("/data/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-pj40-2em-v7n_noskim.root");    
    
    xsecs[0] = 0.001835*81930.0;
    xsecs[1] = 0.05387*8884.0;    
    initweights(tree,xsecs,1.);  
    
    double weightscale = xsecweights[1];
    xsecweights[0] /= weightscale;
    xsecweights[1] /= weightscale;
  }
  
  //training selection cut
  TCut selcut;
  if (dobarrel) {
    selcut = "ph.genpt>16. && ph.isbarrel && ph.ispromptgen"; 
  }
  else {
    selcut = "ph.genpt>16. && !ph.isbarrel && ph.ispromptgen";     
  }
  

  
  TCut selweight = "xsecweight(procidx)";
  TCut prescale10 = "(evt%10==0)";
  TCut prescale20 = "(evt%20==0)";
  TCut prescale25 = "(evt%25==0)";
  TCut prescale50 = "(evt%50==0)";
  TCut prescale100 = "(evt%100==0)";  
  TCut prescale1000 = "(evt%1000==0)";  
  TCut evenevents = "(evt%2==0)";
  TCut oddevents = "(evt%2==1)";  


  //weightvar title used for per-event weights and selection cuts
  if (doele) {
    weightvar.SetTitle(prescale100*evenevents*selcut);
  }
  else {
    weightvar.SetTitle(prescale100*selweight*selcut);
  }
  //create RooDataSet from TChain
  RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",tree,vars,weightvar);   
  
  
  //RooRealVars corresponding to regressed parameters (sigma, mean, left tail parameter, right tail parameter)
  RooRealVar sigwidthtvar("sigwidthtvar","",0.01);
  sigwidthtvar.setConstant(false);
  
  RooRealVar sigmeantvar("sigmeantvar","",1.);
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
  tgts.add(*sigwidtht);
  tgts.add(*sigmeant);
  tgts.add(*signt);
  tgts.add(*sign2t);  
  
  //define transformations corresponding to parameter bounds for non-parametric outputs  
  RooRealConstraint sigwidthlim("sigwidthlim","",*sigwidtht,0.0002,0.5);
  RooRealConstraint sigmeanlim("sigmeanlim","",*sigmeant,0.2,2.0);
  RooRealConstraint signlim("signlim","",*signt,1.01,5000.); 
  RooRealConstraint sign2lim("sign2lim","",*sign2t,1.01,5000.); 

  //define pdf, which depends on transformed outputs (and is intended to be treated as a conditional pdf over the
  //regression inputs in this case)
  //The actual pdf below is a double crystal ball, with crossover points alpha_1 and alpha_2 set constant, but all other
  //parameters regressed
  RooDoubleCBFast sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,RooConst(2.),signlim,RooConst(1.),sign2lim);
  
  //dummy variable
  RooConstVar etermconst("etermconst","",0.);  
   
  //dummy variable
  RooRealVar r("r","",1.);
  r.setConstant();

  //define list of pdfs
  std::vector<RooAbsReal*> vpdf;
  vpdf.push_back(&sigpdf);  

  //define list of training datasets
  std::vector<RooAbsData*> vdata;
  vdata.push_back(hdata);     
  
  //define minimum event weight per tree node
  double minweight = 200;
  std::vector<double> minweights;
  minweights.push_back(minweight);
  
  //temp output file
  TFile *fres = new TFile("fres.root","RECREATE");

  //run training
  if (1) {
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff","",tgts,etermconst,r,vdata,vpdf);
    bdtpdfdiff.SetMinCutSignificance(5.);
    //bdtpdfdiff.SetPrescaleInit(100);
    bdtpdfdiff.SetShrinkage(0.1);
    bdtpdfdiff.SetMinWeights(minweights);
    bdtpdfdiff.SetMaxNodes(750);
    bdtpdfdiff.TrainForest(1e6);   
  }
     
  //create workspace and output to file
  RooWorkspace *wereg = new RooWorkspace("wereg");
  wereg->import(sigpdf);
  
  if (doele && dobarrel)
    wereg->writeToFile("wereg_ele_eb.root");    
  else if (doele && !dobarrel) 
    wereg->writeToFile("wereg_ele_ee.root");    
  else if (!doele && dobarrel)
    wereg->writeToFile("wereg_ph_eb.root");    
  else if (!doele && !dobarrel)
    wereg->writeToFile("wereg_ph_ee.root");    
  
  
  return;
  
  
}