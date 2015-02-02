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

#include "Cintex/Cintex.h"

#include "HybridGBRForestFlex.h"

#ifndef __CINT__
  #include "CondFormats/EgammaObjects/interface/GBRForestD.h"
#endif

using namespace RooFit;

void regcmsformatExample(bool dobarrel=true, bool doele=true) {
  
  //output dir
//   TString dirname = "/data/bendavid/eregexampletest/eregexampletest_test/"; 
//   gSystem->mkdir(dirname,true);
//   gSystem->cd(dirname);    
  
  //read workspace from training
  TString fname;
  if (doele && dobarrel) 
    fname = "wereg_ele_eb.root";
  else if (doele && !dobarrel) 
    fname = "wereg_ele_ee.root";
  else if (!doele && dobarrel) 
    fname = "wereg_ph_eb.root";
  else if (!doele && !dobarrel) 
    fname = "wereg_ph_ee.root";
  
  TString infile = TString::Format("%s",fname.Data());
  
  TFile *fws = TFile::Open(infile); 
  RooWorkspace *ws = (RooWorkspace*)fws->Get("wereg");
  
  //read variables from workspace
  RooGBRTargetFlex *meantgt = static_cast<RooGBRTargetFlex*>(ws->arg("sigmeant"));  
  RooRealVar *tgtvar = ws->var("tgtvar");
  
  
  RooArgList vars;
  vars.add(meantgt->FuncVars());
  vars.add(*tgtvar);
  
  ROOT::Cintex::Cintex::Enable();
  
  TFile *fout = new TFile("wereg_cms.root","RECREATE");
  GBRForestD *forestout = new GBRForestD(*meantgt->Forest());
  fout->WriteObject(forestout,"regmean");
  fout->Close();
  
  
}