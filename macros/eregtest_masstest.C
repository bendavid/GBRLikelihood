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

void eregtest_masstest() {   
 
  TString dirname = "/scratch/bendavid/root/bare/regflexmasstesting_test/"; 
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);    
  
  TFile *fin = TFile::Open("/scratch/bendavid/root/bare/regflexmasstesting_mod100/weregmass.root");
  RooWorkspace *ws = (RooWorkspace*)fin->Get("weregmass");
  
  const RooArgSet *varsset = ws->set("vars");
  RooArgList vars(*varsset);
  
  RooRealVar weightvar("weightvar","",1.);
  
  TString treeloc = "RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/JetPub/JetCorrectionMod/SeparatePileUpMod/ElectronIDMod/MuonIDMod/PhotonPairSelectorPreselInvertEleVetoNoSmear/PhotonTreeWriterPreselInvertEleVetoNoSmear/hPhotonTree";
     
  TChain *tree = new TChain(treeloc);
  tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV/merged/hgg-2013Final8TeV_s12-zllm50-v7n_noskim.root");

  TChain *treedata = new TChain(treeloc);
  treedata->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV/merged/hgg-2013Final8TeV_r12a-pho-j22-v1_noskim.root");  
  treedata->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV/merged/hgg-2013Final8TeV_r12b-dph-j22-v1_noskim.root");  
  treedata->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV/merged/hgg-2013Final8TeV_r12c-dph-j22-v1_noskim.root");  
  treedata->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV/merged/hgg-2013Final8TeV_r12d-dph-j22-v1_noskim.root");   


  
  TCut selcut = "mass>70. && mass<110. && (ph1.pt/mass) > (20./70.) && (ph2.pt/mass) > (20./70.) && evt%4==1 && ph1.isbarrel && ph2.isbarrel";

  
  TCut selweight = "xsecweight(procidx)";
  TCut prescale10 = "(evt%10==0)";
  TCut prescale20 = "(evt%20==0)";
  TCut prescale25 = "(evt%25==0)";
  TCut prescale50 = "(evt%50==0)";
  TCut prescale100 = "(evt%100==0)";  
  TCut prescale1000 = "(evt%1000==0)";  
  TCut evenevents = "(evt%2==0)";
  TCut oddevents = "(evt%2==1)";
  

  
  
  weightvar.SetTitle(selcut);
  RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",tree,vars,weightvar);   
  RooDataSet *hdataD = RooTreeConvert::CreateDataSet("hdataD",treedata,vars,weightvar);
  
  
  
  
  RooAbsReal *scalelim_ph1 = ws->function("scalelim_ph1");
  RooAbsReal *scalelim_ph2 = ws->function("scalelim_ph2");
  RooAbsReal *smearlim_ph1 = ws->function("smearlim_ph1");
  RooAbsReal *smearlim_ph2 = ws->function("smearlim_ph2");  
  
  RooRealVar *mass = ws->var("mass");
  
  RooFormulaVar *scaledmass = new RooFormulaVar("scaledmass","","sqrt(@0*@1)*@2",RooArgList(*scalelim_ph1,*scalelim_ph2,*mass));
  RooRealVar *scaledmassvar = (RooRealVar*)hdataD->addColumn(*scaledmass);  
  
  RooFormulaVar *smearedmass = new RooFormulaVar("smearedmass","","@0*(1.0 + 0.5*sqrt(@1+@2)*RooRandom::gaussian())",RooArgList(*mass,*smearlim_ph1,*smearlim_ph2));
  RooRealVar *smearedmassvar = (RooRealVar*)hdata->addColumn(*smearedmass);  
  
  TH1 *hmass = hdataD->createHistogram("hmass",*mass,Binning(80,70.,110.));
  TH1 *hscaledmass = hdataD->createHistogram("hscaledmass",*scaledmassvar,Binning(80,70.,110.));

  TH1 *hmassMC = hdata->createHistogram("hmassMC",*mass,Binning(80,70.,110.));
  TH1 *hsmearedmass = hdata->createHistogram("hsmearedmass",*smearedmassvar,Binning(80,70.,110.));  
  
  hmassMC->Scale(hmass->GetSumOfWeights()/hmassMC->GetSumOfWeights());
  hsmearedmass->Scale(hscaledmass->GetSumOfWeights()/hsmearedmass->GetSumOfWeights());
  
  hmass->SetLineColor(kRed);
  hscaledmass->SetLineColor(kBlue);
  hmass->SetMarkerColor(kRed);
  hscaledmass->SetMarkerColor(kBlue);
  
  hmassMC->SetLineColor(kRed);
  hsmearedmass->SetLineColor(kBlue);
  
  
  RooRealVar *smearvar1 = (RooRealVar*)hdata->addColumn(*smearlim_ph1);
  RooRealVar *smearvar2 = (RooRealVar*)hdata->addColumn(*smearlim_ph2);
  
  new TCanvas;
  RooPlot *plotsmear1 = smearvar1->frame(-0.1*0.1,0.2*0.2,100);
  hdata->plotOn(plotsmear1);
  plotsmear1->Draw();

  new TCanvas;
  RooPlot *plotsmear2 = smearvar2->frame(-0.1*0.1,0.2*0.2,100);
  hdata->plotOn(plotsmear2);
  plotsmear2->Draw();  
  
  new TCanvas;
  hmassMC->Draw("HIST");
  hsmearedmass->Draw("HISTSAME");
  hscaledmass->Draw("ESAME");
  hmass->Draw("ESAME");
  
}
