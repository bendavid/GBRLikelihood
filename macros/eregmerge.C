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


using namespace RooFit;
 


RooAbsArg *cloneRecursiveRename(RooAbsArg *arg, const char *postfix) {
    
    RooAbsArg *clone = arg->cloneTree();
    
    RooArgSet *clonecomps = clone->getComponents();
    RooArgSet *clonevars = clone->getVariables();
    
    RooArgList cloneargs;
    cloneargs.add(*clonecomps);
    cloneargs.add(*clonevars);
    delete clonecomps;
    delete clonevars;
    
    for (int iarg=0; iarg<cloneargs.getSize(); ++iarg) {
      cloneargs.at(iarg)->SetName(TString::Format("%s_%s",cloneargs.at(iarg)->GetName(),postfix));
    }    
    
    return clone;
    
}

void eregmerge(bool doele) {
  
  TString dirname = "/afs/cern.ch/user/b/bendavid/CMSSWhgg/CMSSW_5_3_11_patch5/src/HiggsAnalysis/GBRLikelihoodEGTools/data/"; 
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);    
  
  TString fnameeb;
  TString fnameee;
  if (doele) { 
    fnameeb = "wereg_ele_eb.root";
    fnameee = "wereg_ele_ee.root";
  }
  else if (!doele) {
    fnameeb = "wereg_ph_eb.root";
    fnameee = "wereg_ph_ee.root";
  }    
    
   
  TString infileeb = TString::Format("/afs/cern.ch/work/b/bendavid/bare/eregAug10RCalphafix/%s",fnameeb.Data());
  TString infileee = TString::Format("/afs/cern.ch/work/b/bendavid/bare/eregAug10RCalphafix/%s",fnameee.Data());
  
  TFile *fwseb = TFile::Open(infileeb); 
  TFile *fwsee = TFile::Open(infileee); 
  
  RooWorkspace *wseb = (RooWorkspace*)fwseb->Get("wereg");
  RooWorkspace *wsee = (RooWorkspace*)fwsee->Get("wereg");
  
  RooAbsPdf *sigpdfeborig = wseb->pdf("sigpdf");
  RooAbsPdf *sigpdfeeorig = wsee->pdf("sigpdf");
  
  RooAbsPdf *sigpdfeb = static_cast<RooAbsPdf*>(cloneRecursiveRename(sigpdfeborig,"EB"));
  RooAbsPdf *sigpdfee = static_cast<RooAbsPdf*>(cloneRecursiveRename(sigpdfeeorig,"EE"));
    
  RooWorkspace *wsout = new RooWorkspace("EGRegressionWorkspace");
  wsout->import(*sigpdfeb);
  wsout->import(*sigpdfee);
  
  TString outname;
  if (doele) outname = "regweights_v4_ele.root";
  else outname = "regweights_v4_ph.root";
  
  wsout->writeToFile(outname);
  
  RooArgList pdfeblist;
  RooArgSet *pdfebcomps = sigpdfeb->getComponents();
  RooArgSet *pdfebvars = sigpdfeb->getVariables();
  pdfeblist.add(*pdfebcomps);
  pdfeblist.add(*pdfebvars);
  delete pdfebcomps;
  delete pdfebvars;
  
  
  RooArgList pdfeelist;
  RooArgSet *pdfeecomps = sigpdfee->getComponents();
  RooArgSet *pdfeevars = sigpdfee->getVariables();
  pdfeelist.add(*pdfeecomps);
  pdfeelist.add(*pdfeevars);
  delete pdfeecomps;
  delete pdfeevars;  
  
  
//   RooArgList components(ws->components());
//   for (int iarg=0; iarg<components.getSize(); ++iarg) {
//     components.at(iarg)->SetName(TString::Format("%s_1",components.at(iarg)->GetName()));
//   }
  
  RooGBRFunction *funceb = static_cast<RooGBRFunction*>(pdfeblist.find("func_EB"));
  RooGBRFunction *funcee = static_cast<RooGBRFunction*>(pdfeelist.find("func_EE"));
  
//   funceb->Vars().Print("V");
//   funcee->Vars().Print("V");

  for (int ivar=0; ivar<funceb->Vars().getSize(); ++ivar) {
    printf("%i: %s, %s\n",ivar,funceb->Vars().at(ivar)->GetName(),funceb->Vars().at(ivar)->GetTitle());
  }
  
  for (int ivar=0; ivar<funcee->Vars().getSize(); ++ivar) {
    printf("%i: %s, %s\n",ivar,funcee->Vars().at(ivar)->GetName(),funcee->Vars().at(ivar)->GetTitle());
  }
  
  TString outnameforest;
  if (doele) outnameforest = "regweights_v4_forest_ele.root";
  else outnameforest = "regweights_v4_forest_ph.root";  
  
  TFile *fforest = new TFile(outnameforest,"RECREATE");
  fforest->WriteObject(funceb->Forest(),"EGRegressionForest_EB");
  fforest->WriteObject(funcee->Forest(),"EGRegressionForest_EE");
  fforest->Close();
  
}