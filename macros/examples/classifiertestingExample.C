//Example macro for testing a simple BDT classifier using the GBRLikelihood machinery
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
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGraphErrors.h"

using namespace RooFit;
 
void classifiertestingExample() {
    
  TString infile = "weclass.root";
  
  TFile *fws = TFile::Open(infile); 
  RooWorkspace *ws = (RooWorkspace*)fws->Get("wclass");
  
  //read variables from workspace
  RooGBRTargetFlex *logsbtgt = static_cast<RooGBRTargetFlex*>(ws->arg("logsbtgt"));  
    
  RooArgList vars;
  vars.add(logsbtgt->FuncVars());
   
  //read testing dataset from TTree
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
  TCut prescale100test = "(evt%100==1)";  


  //create RooDataSet from TChain
  weightvar.SetTitle(prescale100test*basesel*sigcut);
  RooDataSet *hdatasig = RooTreeConvert::CreateDataSet("hdatasig",tree,vars,weightvar);   

  weightvar.SetTitle(prescale100test*basesel*bkgcut);
  RooDataSet *hdatabkg = RooTreeConvert::CreateDataSet("hdatabkg",tree,vars,weightvar);     
 
  //S/(S+B) defined on interval 0,1
  RooFormulaVar *ssb = new RooFormulaVar("ssb","","1./(1.+exp(-@0))",RooArgList(*logsbtgt));
  RooRealVar *ssbvar = (RooRealVar*)hdatasig->addColumn(*ssb);
  hdatabkg->addColumn(*ssb);
  ssbvar->setRange(0.,1.);
  ssbvar->setBins(100);  
  
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();


  //plot signal and background distributions for predicted S/(S+B) (ie the classifier output)
  TH1 *hsig = hdatasig->createHistogram("hsig",*ssbvar);
  TH1 *hbkg = hdatabkg->createHistogram("hbkg",*ssbvar);
    
  hsig->SetLineColor(kBlue);
  hbkg->SetLineColor(kRed);
  
  hsig->Scale(1./hsig->GetSumOfWeights());
  hbkg->Scale(1./hbkg->GetSumOfWeights());
  
  hsig->GetXaxis()->SetTitle("Predicted S/(S+B)");
  
  TCanvas *cssb = new TCanvas;
  hsig->Draw("HIST");
  hbkg->Draw("HISTSAME");
  cssb->SaveAs("ssb.pdf");


  //plot measured vs predicted S/(S+B)
  TF1 *fssb = new TF1("fssb","x",0.,1.); 
  fssb->SetLineColor(kRed);
  fssb->SetLineStyle(9);
  
  TH1D *hsigbkg = new TH1D(*(TH1D*)hsig);
  hsigbkg->Add(hbkg);    
  
  TH1D *hssb = new TH1D(*(TH1D*)hsig);
  int nbins = hssb->GetNbinsX();
  for (int ibin=1; ibin<=nbins; ++ibin) {
    double sval = hsig->GetBinContent(ibin);
    double serr = hsig->GetBinError(ibin);
    double bval = hbkg->GetBinContent(ibin);
    double berr = hbkg->GetBinError(ibin);
    
    double sbval = sval+bval;
    
    double ssb = sval/(sval+bval);
    double ssberr = pow(sval+bval,-2)*sqrt(bval*bval*serr*serr + sval*sval*berr*berr);
    
    if (sbval<=0.) {
      ssb = 0.;
      ssberr = 0.;
    }
    
    hssb->SetBinContent(ibin,ssb);
    hssb->SetBinError(ibin,ssberr);
  }
  hssb->GetXaxis()->SetTitle("Predicted S/(S+B)");
  hssb->GetYaxis()->SetTitle("Measured S/(S+B)");
  
  
  TCanvas *cssbcompare = new TCanvas;
  hssb->Draw("E");
  fssb->Draw("SAME");
  cssbcompare->SaveAs("ssbcompare.pdf");



  int nbinsroc = 1000;
  
  TH1 *hsigroc = hdatasig->createHistogram("hsigroc",*ssbvar,RooFit::Binning(nbinsroc,0.,1.));
  TH1 *hbkgroc = hdatabkg->createHistogram("hbkgroc",*ssbvar,RooFit::Binning(nbinsroc,0.,1.));  
  //plot ROC curve with (some approximation of) error bands
  TGraphErrors *hroc = new TGraphErrors(nbinsroc);
  for (int ibin=1; ibin<=nbinsroc; ++ibin) {
    double stot = hsigroc->Integral(0,nbinsroc+1);
    double sacc = hsigroc->Integral(ibin,nbinsroc+1);
    double srej = stot-sacc;
    double seff = sacc/stot;
    
    double btot = hbkgroc->Integral(0,nbinsroc+1);
    double bacc = hbkgroc->Integral(ibin,nbinsroc+1);
    double brej = btot-bacc;
    double beff = bacc/btot;
    
    double srejerrsq = 0.;
    for (int jbin=0; jbin<ibin; ++jbin) {
      double binerr = hsigroc->GetBinError(jbin);
      srejerrsq += binerr*binerr;
    }
    double saccerrsq = 0.;
    for (int jbin=ibin; jbin<=(nbinsroc+1); ++jbin) {
      double binerr = hsigroc->GetBinError(jbin);
      saccerrsq += binerr*binerr;
    }
    double brejerrsq = 0.;
    for (int jbin=0; jbin<ibin; ++jbin) {
      double binerr = hbkgroc->GetBinError(jbin);
      brejerrsq += binerr*binerr;
    }
    double baccerrsq = 0.;
    for (int jbin=ibin; jbin<=(nbinsroc+1); ++jbin) {
      double binerr = hbkgroc->GetBinError(jbin);
      baccerrsq += binerr*binerr;
    }

    double sefferr = pow(stot,-2)*sqrt(srej*srej*saccerrsq + sacc*sacc*srejerrsq);
    double befferr = pow(btot,-2)*sqrt(brej*brej*baccerrsq + bacc*bacc*brejerrsq);
    
    if (seff<=0.) {
      seff = 0.;
      sefferr = 0.;
    }
    if (beff<=0.) {
      beff = 0.;
      befferr = 0.;
    }

    hroc->SetPoint(ibin-1,seff,beff);
    hroc->SetPointError(ibin-1,sefferr,befferr);
    
  }
  hroc->SetFillStyle(3244);
  hroc->SetFillColor(kBlue);
  
  TCanvas *croc = new TCanvas;
  hroc->Draw("AL2");
  hroc->GetXaxis()->SetTitle("Signal Efficiency");
  hroc->GetYaxis()->SetTitle("Background Efficiency");  
  croc->SaveAs("roc.pdf");
  
  
}