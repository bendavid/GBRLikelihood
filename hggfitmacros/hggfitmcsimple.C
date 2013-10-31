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
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooCBShape.h"
//#include "RooCBShapeModified.h"
//#include "RooBernsteinFast.h"
#include "RooWorkspace.h"
#include "TH1D.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TLegend.h"
#include "RooRandom.h"
#include "RooAddition.h"
#include "TSystem.h"
#include "RooCBFast.h"
#include "RooDoubleCBFast.h"
#include "RooGaussianFast.h"
#include "RooBernsteinFast.h"
#include "TProfile.h"
#include "TStyle.h"
#include "RooRealBinding.h"

 using namespace RooFit;
 
 
RooAbsArg *cloneRecursiveRename(RooAbsArg *arg, const char *postfix, const char *titlein, const char *titleout) {
    
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
      
      TString oldtitle = cloneargs.at(iarg)->GetTitle();
      TString newtitle = oldtitle.ReplaceAll(titlein,titleout);
      cloneargs.at(iarg)->SetTitle(newtitle);
    }    
    
    return clone;
    
}
 
 
TH1D *getpuweights(TFile *file, TH1D *target) {
  
  TDirectory *dirmcpv = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  TH1D *hnpu = (TH1D*)dirmcpv->Get("hNPU");
  TH1D *hpumc = (TH1D*)hnpu->Clone();
  
  hpumc->Sumw2();
  hpumc->Scale(1.0/hpumc->Integral(0,hpumc->GetNbinsX()+1));
  
  
  TH1D *htargettmp = new TH1D("htargettmp","", hpumc->GetNbinsX(), hpumc->GetXaxis()->GetXmin(), hpumc->GetXaxis()->GetXmax());
  htargettmp->Sumw2();
  for (int ibin = 0; ibin<=(htargettmp->GetNbinsX()+1); ++ibin) {
    htargettmp->Fill(htargettmp->GetBinCenter(ibin),target->GetBinContent(target->FindFixBin(htargettmp->GetBinCenter(ibin))));
  }
  htargettmp->Scale(1.0/htargettmp->Integral(0,htargettmp->GetNbinsX()+1));
  
  TH1D *puweights = new TH1D((*htargettmp)/(*hpumc));
    
  delete htargettmp;
  
  return puweights;
    
}

double getweight(TFile *file, double xsec) {
 
  TDirectory *dir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  TH1D *hallevts = (TH1D*)dir->Get("hDTotalMCWeight");
  
  return xsec/hallevts->GetSumOfWeights();
  
}


TH1D *puweights[50];
float puweight(float npu, int wset=0) {
  if (!puweights[wset]) return 1.0;
  if (npu<0) return 1.0;
  return puweights[wset]->GetBinContent(puweights[wset]->FindFixBin(npu));
}

double xsecweights[50];
double xsecweight(int procidx=0) {
  return xsecweights[procidx];
}

float kfact(int procidx, bool ispromptgen1, bool ispromptgen2) {
 
  if (procidx==0 || procidx==9) return 0.975; //gf
  
  if (procidx>=4 && procidx<=8) {
    int npromptgen = int(ispromptgen1) + int(ispromptgen2);
    
    if (npromptgen==0) return 1.0; //fake fake
    else if (npromptgen==1) return 1.3; //prompt-fake
    //else if (npromptgen==1) return 1.0; //prompt-fake with sherpa
    else if (npromptgen==2) return 1.0; //prompt-prompt (for sherpa diphoj)
  }
  
  return 1.0;
  
}



void initweights(TChain *chain, TH1D *target, double *xsecs, double lumi) {
 
  TObjArray *files = chain->GetListOfFiles();
  for (int i=0; i<files->GetEntries(); ++i) {    
    TFile *file = TFile::Open(files->At(i)->GetTitle(),"READ");
    
    puweights[i] = getpuweights(file,target);
    xsecweights[i] = getweight(file,lumi*xsecs[i]);
    
    file->Close();    
  } 
  
  chain->SetAlias("procidx","This->GetTreeNumber()");
  
}
 
 
void hggfitmcsimple(double nommass=123., double tgtr=1., int ijob=0) {
    
  //gSystem->cd("/scratch/bendavid/root/bare/fitplotsJun10test/");
  
  int seed = 65539+ijob+1; 
  
  //TString dirname = "/scratch/bendavid/root/bare/hggfitmctestSep17_sherpa_badratio_PPpow/";
  TString dirname = "/scratch/bendavid/root/bare/gbrMultiClassOct31_dijet/";
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);
  
  
  //nommass=150.;
 // gSystem->cd("/scratch/bendavid/root/bare/fitplotsJun8_150_2x/");
  
  gRandom->SetSeed(seed);
  RooRandom::randomGenerator()->SetSeed(seed);    
  
  
//   TFile *fin = TFile::Open("/home/mingyang/cms/hist_approval/hgg-2013Moriond/merged/hgg-2013Moriond_s12-h150gg-gf-v7a_noskim.root");
//   TDirectory *hdir = (TDirectory*)fin->FindObjectAny("PhotonTreeWriterPresel");
//   TTree *htree = (TTree*)hdir->Get("hPhotonTree");

//   TFile *fdin = TFile::Open("/home/mingyang/cms/hist/hgg-2013Moriond/merged/hgg-2013Moriond_r12_ABCD.root");
//   TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("PhotonTreeWriterPresel");
//   TTree *dtree = (TTree*)ddir->Get("hPhotonTree");  
  
  //TCut selcut = "(ph1.pt > (mass/3.0) && ph2.pt > (mass/4.0) && mass>100. && mass<180. && ph1.idmva>-0.2 && ph2.idmva>-0.2)";
  //TCut selcutnoid = "(ph1.pt > (mass/3.0) && ph2.pt > (mass/4.0) && mass>100. && mass<180.)";
  
  TCut selcutbase = "(ph1.pt > (mass/3.0) && ph2.pt > (mass/4.0) && mass>100. && mass<180. && ph1.idmva>-0.2 && ph2.idmva>-0.2)";
  TCut jetcut = "(jet1pt>30. && jet2pt>20.)";
  TCut selcut = selcutbase*jetcut;
  
  
  //TCut selweight = "xsecweight(procidx)*puweight(numPU,procidx)";
  TCut selweight = "xsecweight(procidx)*mcweight*kfact(procidx,ph1.ispromptgen,ph2.ispromptgen)";
  
  TCut sigGFcut = "(procidx==0)";
  TCut sigVBFcut = "(procidx==1)";
  TCut sigVHcut = "(procidx==2)";
  TCut sigTTHcut = "(procidx==3)";
  
  TCut bkgPPcut = "(procidx==4)";
  TCut bkgPFcut = "(procidx==5 || procidx==6)";
  TCut bkgFFcut = "(procidx==7 || procidx==8)";
  
  
  
  TCut bkgcut = "(procidx>3)";
  TCut bkgcutnoq = "(procidx>3 && procidx<7)";
  
  TCut prescalenone = "(1==1)";
  TCut evenevents = "(evt%2==0)";
  TCut oddevents =  "(evt%2==1)";
  TCut prescale10 = "(evt%10==0)";
  TCut prescale25 = "(evt%25==0)";
  TCut prescale50 = "(evt%50==0)";
  TCut prescale100 = "(evt%100==0)";
  
  TCut trainingsel = "(procidx>=0 && procidx<=3) || (evt%2==0 && (procidx>=4 && procidx<=8))";
  //TCut trainingsel = "(procidx>=0 && procidx<=3) || (evt%2==1 && (procidx>=4 && procidx<=8))";
  
  
  
  TCut fcut = prescale50;  
  
  
  double xsecs[50];

  
  TChain *tree = new TChain("RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/JetPub/JetCorrectionMod/SeparatePileUpMod/ElectronIDMod/MuonIDMod/PhotonPairSelectorPresel/PhotonTreeWriterPresel/hPhotonTree");

  
  if (1) {
    tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_0/merged/hgg-2013Final8TeV_s12-h123gg-gf-v7n_noskim.root");
    tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_0/merged/hgg-2013Final8TeV_s12-h123gg-vbf-v7n_noskim.root");
    tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_0/merged/hgg-2013Final8TeV_s12-h123gg-vh-v7n_noskim.root");
    tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_0/merged/hgg-2013Final8TeV_s12-h123gg-tt-v7n_noskim.root");
    tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_0/merged/hgg-2013Final8TeV_s12-diphoj-m60-v7n_noskim.root");
    tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_0/merged/hgg-2013Final8TeV_s12-pj20_40-2em-v7n_noskim.root");
    tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_0/merged/hgg-2013Final8TeV_s12-pj40-2em-v7n_noskim.root");    
    tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_0/merged/hgg-2013Final8TeV_s12-qcd30_40-2em-v7n_noskim.root");
    tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_0/merged/hgg-2013Final8TeV_s12-qcd40-2em-v7n_noskim.root"); 
    
    tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_0/merged/hgg-2013Final8TeV_s12-h124gg-gf-v7n_noskim.root");
    tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_0/merged/hgg-2013Final8TeV_s12-h124gg-vbf-v7n_noskim.root");
    tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_0/merged/hgg-2013Final8TeV_s12-h124gg-vh-v7n_noskim.root");
    tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_0/merged/hgg-2013Final8TeV_s12-h124gg-tt-v7n_noskim.root");    
    
    xsecs[0] = 19.88*2.27e-03;
    xsecs[1] = 1.608*2.27e-03;
    xsecs[2] = (0.7434+0.4366)*2.27e-03;
    xsecs[3] = 0.1356*2.27e-03;      
    
    xsecs[4] = 120.354;  //sherpa diphoj
    xsecs[5] = 0.001835*81930.0; //gamma+jets
    xsecs[6] = 0.05387*8884.0;   //gamma+jets
    xsecs[7] = 0.000235*5.195e+07; //dijet
    xsecs[8] = 0.002175*2.365e+07; //dijet
    
    xsecs[9] = 19.57*2.27e-03;
    xsecs[10] = 1.595*2.27e-03;
    xsecs[11] = (0.7239+0.4260)*2.27e-03;
    xsecs[12] = 0.1324*2.27e-03;         
  }
  
  tree->SetCacheSize(64*1024*1024);
  
  
  const double lumi = 19620.;  
  TFile *filepuest = new TFile("/scratch/bendavid/root/Moriond_pileup_190456-208686_minBiasXsec69400_corr_observed.root","READ");
  TH1D *hpuest = (TH1D*) filepuest->Get("pileup");
  
  initweights(tree,hpuest,xsecs,lumi);  
  
//   new TCanvas;
//   tree->Draw("mcweight",selcutbase*TCut("procidx==4"));
//   
//   
//   new TCanvas;
//   tree->Draw("mcweight",selweight*selcutbase*TCut("procidx==4"));    
//   
//   
//   new TCanvas;
//   tree->Draw("mcweight",selcut*TCut("procidx==4"));
//   
//   
//   new TCanvas;
//   tree->Draw("mcweight",selweight*selcut*TCut("procidx==4"));  
//   
//   return;
  
//   new TCanvas;
//   tree->Draw("vtxprob",selweight*selcut*TCut("procidx==1"));  
//   
//   new TCanvas;
//   tree->Draw("vtxprob",selweight*selcut*TCut("procidx==1 && jet1pt>30. && jet2pt>25."));
// 
//   new TCanvas;
//   tree->Draw("vtxprob",selweight*selcut*TCut("procidx==1 && jet1pt>30. && jet2pt>25. && dijetmass>250."));  
//   
//   
//   new TCanvas;
//   tree->Draw("abs(vtxZ-genHiggsZ)<1.0",selweight*selcut*TCut("procidx==1"));  
//   
//   new TCanvas;
//   tree->Draw("abs(vtxZ-genHiggsZ)<1.0",selweight*selcut*TCut("procidx==1 && jet1pt>30. && jet2pt>25."));
// 
//   new TCanvas;
//   tree->Draw("abs(vtxZ-genHiggsZ)<1.0",selweight*selcut*TCut("procidx==1 && jet1pt>30. && jet2pt>25. && dijetmass>250."));    
//   
//   return;

  RooRealVar masserr("masserr","masserrsmeared/mass",0.);
  RooRealVar masserrwrong("masserrwrong","masserrsmearedwrongvtx/mass",0.);
  RooRealVar vtxprob("vtxprob","vtxprob",0.);
  RooRealVar pt1("pt1","ph1.pt/mass",0.);
  RooRealVar pt2("pt2","ph2.pt/mass",0.);
  RooRealVar eta1("eta1","ph1.eta",0.);
  RooRealVar eta2("eta2","ph2.eta",0.);  
  RooRealVar dphi("dphi","TMath::Cos(ph1.phi-ph2.phi)",0.);
  RooRealVar idmva1("idmva1","ph1.idmva",0.);
  RooRealVar idmva2("idmva2","ph2.idmva",0.);
  RooRealVar mass("mass","mass",100.,100.,180.);

  RooRealVar jeteta1("jeteta1","(jet1pt>30. && jet2pt>20.)*jet1eta",0.);
  RooRealVar jeteta2("jeteta2","(jet1pt>30. && jet2pt>20.)*jet2eta",0.);  
  RooRealVar jetpt1("jetpt1","(jet1pt>30. && jet2pt>20.)*jet1pt",0.);
  RooRealVar jetpt2("jetpt2","(jet1pt>30. && jet2pt>20.)*jet2pt",0.);
  RooRealVar zeppenfeld("zeppenfeld","(jet1pt>30. && jet2pt>20.)*zeppenfeld",0.);
  RooRealVar dphidijetgg("dphidijetgg","(jet1pt>30. && jet2pt>20.)*TMath::Min(dphidijetgg,2.916)",0.);  
  RooRealVar dijetmass("dijetmass","(jet1pt>30. && jet2pt>20.)*dijetmass",0.);
  RooRealVar ptgg("ptgg","ptgg/mass",0.);

  
  RooRealVar procidx("procidx","procidx",0.);
  RooRealVar evt("evt","evt",0.);
  
  RooArgList vars;
  vars.add(masserr);
  vars.add(masserrwrong);
  vars.add(vtxprob);
  vars.add(pt1);
  vars.add(pt2);
  vars.add(eta1);
  vars.add(eta2);
  vars.add(dphi);
  vars.add(idmva1);
  vars.add(idmva2);

  RooArgList condvars(vars);
  
  vars.add(mass);
  
  vars.add(jeteta1);
  vars.add(jeteta2);  
  vars.add(jetpt1);
  vars.add(jetpt2);
  vars.add(zeppenfeld);  
  vars.add(dphidijetgg);
  vars.add(dijetmass);
  vars.add(ptgg);  
  
  vars.add(procidx);
  vars.add(evt);
  
  RooArgList condvarsdijet;
  condvarsdijet.add(jeteta1);
  condvarsdijet.add(jeteta2);  
  condvarsdijet.add(jetpt1);
  condvarsdijet.add(jetpt2);
  condvarsdijet.add(zeppenfeld);  
  condvarsdijet.add(dphidijetgg);
  condvarsdijet.add(dijetmass);
  condvarsdijet.add(ptgg);  
  
  RooArgList condvarsfull(condvars);
  condvarsfull.add(jeteta1);
  condvarsfull.add(jeteta2);  
  condvarsfull.add(jetpt1);
  condvarsfull.add(jetpt2);
  condvarsfull.add(zeppenfeld);  
  condvarsfull.add(dphidijetgg);
  condvarsfull.add(dijetmass);
  
  //RooArgList trainingvars = condvars;
  //RooArgList trainingvars = condvarsfull;
  RooArgList trainingvars = condvarsdijet;
  
  vars.Print("V");
  
  condvars.Print("V");
  
  RooRealVar weightvar("weightvar","",1.);

  weightvar.SetTitle(selweight*selcut);
  RooDataSet *hdataAllW = RooTreeConvert::CreateDataSet("hdataAllW",tree,vars,weightvar);  
  
//   RooDataSet *hdataBkgPPtmp = (RooDataSet*)hdataAllW->reduce(bkgPPcut);
//   RooDataSet *hdataBkgPFtmp = (RooDataSet*)hdataAllW->reduce(bkgPFcut);
//   RooDataSet *hdataBkgFFtmp = (RooDataSet*)hdataAllW->reduce(bkgFFcut);    
  
//   tree->SetAlias("sbmass",TString::Format("1.0*%5f",sbmass.getVal()));
//   tree->SetAlias("bkgPPpeakpdf",TString::Format("1.0*%5f",bkgPPpdft.getVal()));
//   tree->SetAlias("bkgPFpeakpdf",TString::Format("1.0*%5f",bkgPFpdft.getVal()));
//   tree->SetAlias("bkgFFpeakpdf",TString::Format("1.0*%5f",bkgFFpdft.getVal()));
//   
//   TCut sbpeakweight = "((procidx>=0 && procidx<=3) || (procidx>=9 && procidx<=12))*(vtxprob/(sbmass*masserrsmeared*sqrt(2.0*TMath::Pi())/mass)+(1.0-vtxprob)/(sbmass*masserrsmearedwrongvtx*sqrt(2.0*TMath::Pi())/mass)) + (procidx==4)*bkgPPpeakpdf + (procidx==5 || procidx==6)*bkgPFpeakpdf + (procidx==7 || procidx==8)*bkgFFpeakpdf";
  
//   tree->SetAlias("sbmass",TString::Format("1.0*%5f",sbmass.getVal()));
//   tree->SetAlias("bkgPPpeakpdf",TString::Format("1.0*%5f",bkgPPpdft.getVal()));
//   tree->SetAlias("bkgPFpeakpdf",TString::Format("1.0*%5f",bkgPFpdft.getVal()));
//   tree->SetAlias("bkgFFpeakpdf",TString::Format("1.0*%5f",bkgFFpdft.getVal()));
  
  TCut sbpeakweight = "((procidx>=0 && procidx<=3) || (procidx>=9 && procidx<=12))*(vtxprob*mass/masserrsmeared+(1.0-vtxprob)*mass/masserrsmearedwrongvtx) + (procidx>=4 && procidx<=8)*1.0";
  
  weightvar.SetTitle(selweight*sbpeakweight*selcut);
  RooDataSet *hdataAllWPeak = RooTreeConvert::CreateDataSet("hdataAllWpeak",tree,vars,weightvar);  
  
  RooConstVar sigGFscale("sigGFscale","",xsecweights[0]);
  xsecweights[0]/=sigGFscale.getVal();
  
  RooConstVar sigVBFscale("sigVBFscale","",xsecweights[1]);
  xsecweights[1]/=sigVBFscale.getVal();

  RooConstVar sigVHscale("sigVHscale","",xsecweights[2]);
  xsecweights[2]/=sigVHscale.getVal();  

  RooConstVar sigTTHscale("sigTTHscale","",xsecweights[3]);
  xsecweights[3]/=sigTTHscale.getVal();    
  
  RooConstVar bkgPPscale("bkgPPscale","",xsecweights[4]);  
  xsecweights[4]/=bkgPPscale.getVal();

  RooConstVar bkgPFscale("bkgPFscale","",xsecweights[6]);  
  xsecweights[5]/=bkgPFscale.getVal();  
  xsecweights[6]/=bkgPFscale.getVal();  

  RooConstVar bkgFFscale("bkgFFscale","",xsecweights[8]);  
  xsecweights[7]/=bkgFFscale.getVal();  
  xsecweights[8]/=bkgFFscale.getVal();      
  

  weightvar.SetTitle(trainingsel*selweight*selcut);  
  RooDataSet *hdataAll = RooTreeConvert::CreateDataSet("hdataAll",tree,vars,weightvar);  
  
  RooDataSet *hdataSigGF = (RooDataSet*)hdataAll->reduce(sigGFcut);
  RooDataSet *hdataSigVBF = (RooDataSet*)hdataAll->reduce(sigVBFcut);
  RooDataSet *hdataSigVH = (RooDataSet*)hdataAll->reduce(sigVHcut);
  RooDataSet *hdataSigTTH = (RooDataSet*)hdataAll->reduce(sigTTHcut);
  RooDataSet *hdataBkgPP = (RooDataSet*)hdataAll->reduce(bkgPPcut);
  RooDataSet *hdataBkgPF = (RooDataSet*)hdataAll->reduce(bkgPFcut);
  RooDataSet *hdataBkgFF = (RooDataSet*)hdataAll->reduce(bkgFFcut);
  
  RooRealVar bkgPPp0var("bkgPPp0var","",0.1);
  bkgPPp0var.setConstant(false);   
  
  RooRealVar bkgPFp0var("bkgPFp0var","",0.);
  bkgPFp0var.setConstant(false);   
  
  RooRealVar bkgFFp0var("bkgFFp0var","",0.);
  bkgFFp0var.setConstant(false);       
  
  RooPowerLaw bkgPPpdft("bkgPPpdft","",mass,bkgPPp0var);
  RooPowerLaw bkgPFpdft("bkgPFpdft","",mass,bkgPFp0var);   
  RooPowerLaw bkgFFpdft("bkgFFpdft","",mass,bkgFFp0var);     
  
  bkgPPpdft.fitTo(*hdataBkgPP,RooFit::Strategy(0),RooFit::ConditionalObservables(condvars),RooFit::NumCPU(16),RooFit::Minimizer("Minuit2","minimize"));
  bkgPFpdft.fitTo(*hdataBkgPF,RooFit::Strategy(0),RooFit::ConditionalObservables(condvars),RooFit::NumCPU(16),RooFit::Minimizer("Minuit2","minimize"));
  bkgFFpdft.fitTo(*hdataBkgFF,RooFit::Strategy(0),RooFit::ConditionalObservables(condvars),RooFit::NumCPU(16),RooFit::Minimizer("Minuit2","minimize"));

  RooConstVar sbmass("sbmass","",125.);
  mass.setVal(sbmass.getVal());
  
  RooConstVar bkgPPpeakpdf("bkgPPpeakpdf","",bkgPPpdft.getVal(mass));
  RooConstVar bkgPFpeakpdf("bkgPFpeakpdf","",bkgPFpdft.getVal(mass));
  RooConstVar bkgFFpeakpdf("bkgFFpeakpdf","",bkgFFpdft.getVal(mass)); 
  
  
  RooRealVar sigGFllrvar("sigGFllrvar","",0.);
  sigGFllrvar.setConstant(false);  

  RooRealVar sigVBFllrvar("sigVBFllrvar","",0.);
  sigVBFllrvar.setConstant(false);  

  RooRealVar sigVHllrvar("sigVHllrvar","",0.);
  sigVHllrvar.setConstant(false);    
  
  RooRealVar sigTTHllrvar("sigTTHllrvar","",0.);
  sigTTHllrvar.setConstant(false);    
  
  RooRealVar bkgPPllrvar("bkgPPllrvar","",0.);
  bkgPPllrvar.setConstant(false);   

  RooRealVar bkgPFllrvar("bkgPFllrvar","",0.);
  bkgPFllrvar.setConstant(false);     
  
  RooRealVar bkgFFllrvar("bkgFFllrvar","",0.);
  bkgFFllrvar.setConstant(false);     
     
  RooGBRFunctionFlex sigGFllrfunc("sigGFllrfunc","");
  RooGBRFunctionFlex sigVBFllrfunc("sigVBFllrfunc","");
  RooGBRFunctionFlex sigVHllrfunc("sigVHllrfunc","");
  RooGBRFunctionFlex sigTTHllrfunc("sigTTHllrfunc","");
  RooGBRFunctionFlex bkgPFllrfunc("bkgPFllrfunc","");
  RooGBRFunctionFlex bkgFFllrfunc("bkgFFllrfunc","");
  RooGBRFunctionFlex bkgPPllrfunc("bkgPPllrfunc","");
  
  RooGBRTargetFlex sigGFllr("sigGFllr","",sigGFllrfunc,sigGFllrvar,trainingvars);
  RooGBRTargetFlex sigVBFllr("sigVBFllr","",sigVBFllrfunc,sigVBFllrvar,trainingvars);
  RooGBRTargetFlex sigVHllr("sigVHllr","",sigVHllrfunc,sigVHllrvar,trainingvars);
  RooGBRTargetFlex sigTTHllr("sigTTHllr","",sigTTHllrfunc,sigTTHllrvar,trainingvars);
  RooGBRTargetFlex bkgPFllr("bkgPFllr","",bkgPFllrfunc,bkgPFllrvar,trainingvars);
  RooGBRTargetFlex bkgFFllr("bkgFFllr","",bkgFFllrfunc,bkgFFllrvar,trainingvars);
  RooGBRTargetFlex bkgPPllr("bkgPPllr","",bkgPPllrfunc,bkgPPllrvar,trainingvars);
   
  RooArgList tgts;  
  tgts.add(sigGFllr);
  tgts.add(sigVBFllr);
  tgts.add(sigVHllr);
  tgts.add(sigTTHllr);
  tgts.add(bkgPFllr);
  tgts.add(bkgFFllr);
  
  double totalEntries = hdataAll->sumEntries();
  
  sigGFllrvar.setVal( sqrt(hdataSigGF->sumEntries()/hdataBkgPP->sumEntries()) );
  RooFormulaVar sigGFllrlim("sigGFllrlim","","@0*@0",RooArgList(sigGFllr));
  
  sigVBFllrvar.setVal( sqrt(hdataSigVBF->sumEntries()/hdataBkgPP->sumEntries()) );
  RooFormulaVar sigVBFllrlim("sigVBFllrlim","","@0*@0",RooArgList(sigVBFllr));  

  sigVHllrvar.setVal( sqrt(hdataSigVH->sumEntries()/hdataBkgPP->sumEntries()) );
  RooFormulaVar sigVHllrlim("sigVHllrlim","","@0*@0",RooArgList(sigVHllr));    
  
  sigTTHllrvar.setVal( sqrt(hdataSigTTH->sumEntries()/hdataBkgPP->sumEntries()) );
  RooFormulaVar sigTTHllrlim("sigTTHllrlim","","@0*@0",RooArgList(sigTTHllr));    
  
  bkgPFllrvar.setVal( sqrt(hdataBkgPF->sumEntries()/hdataBkgPP->sumEntries()) );
  RooFormulaVar bkgPFllrlim("bkgPFllrlim","","@0*@0",RooArgList(bkgPFllr));  

  bkgFFllrvar.setVal( sqrt(hdataBkgFF->sumEntries()/hdataBkgPP->sumEntries()) );
  RooFormulaVar bkgFFllrlim("bkgFFllrlim","","@0*@0",RooArgList(bkgFFllr));    

  RooConstVar bkgPPllrlim("bkgPPllrlim","",1.);

  RooFormulaVar llrden("llrden","","1.0/(1.0 + @0 + @1 + @2 + @3 + @4 + @5)",RooArgList(sigGFllrlim,sigVBFllrlim,sigVHllrlim,sigTTHllrlim,bkgPFllrlim,bkgFFllrlim));
  
  RooFormulaVar bkgPPreal("bkgPPreal","","@0*@1",RooArgList(llrden,bkgPPllrlim));
  RooFormulaVar bkgPFreal("bkgPFreal","","@0*@1",RooArgList(llrden,bkgPFllrlim));
  RooFormulaVar bkgFFreal("bkgFFreal","","@0*@1",RooArgList(llrden,bkgFFllrlim));
  RooFormulaVar sigGFreal("sigGFreal","","@0*@1",RooArgList(llrden,sigGFllrlim));
  RooFormulaVar sigVBFreal("sigVBFreal","","@0*@1",RooArgList(llrden,sigVBFllrlim));  
  RooFormulaVar sigVHreal("sigVHreal","","@0*@1",RooArgList(llrden,sigVHllrlim)); 
  RooFormulaVar sigTTHreal("sigTTHreal","","@0*@1",RooArgList(llrden,sigTTHllrlim));    
  
  RooFormulaVar sigpeakpdf("sigpeakpdf","","@0/(@1*@2*sqrt(2.0*TMath::Pi())) + (1.0-@0)/(@1*@3*sqrt(2.0*TMath::Pi()))",RooArgList(vtxprob,sbmass,masserr,masserrwrong));
  
  RooFormulaVar sigGFscaled("sigGFscaled","","@0*@1*@2",RooArgList(sigGFscale,sigGFllrlim,sigpeakpdf));
  RooFormulaVar sigVBFscaled("sigVBFscaled","","@0*@1*@2",RooArgList(sigVBFscale,sigVBFllrlim,sigpeakpdf));
  RooFormulaVar sigVHscaled("sigVHscaled","","@0*@1*@2",RooArgList(sigVHscale,sigVHllrlim,sigpeakpdf));
  RooFormulaVar sigTTHscaled("sigTTHscaled","","@0*@1*@2",RooArgList(sigTTHscale,sigTTHllrlim,sigpeakpdf));
  RooFormulaVar bkgPPscaled("bkgPPscaled","","@0*@1*@2",RooArgList(bkgPPscale,bkgPPllrlim,bkgPPpeakpdf));
  RooFormulaVar bkgPFscaled("bkgPFscaled","","@0*@1*@2",RooArgList(bkgPFscale,bkgPFllrlim,bkgPFpeakpdf));
  RooFormulaVar bkgFFscaled("bkgFFscaled","","@0*@1*@2",RooArgList(bkgFFscale,bkgFFllrlim,bkgFFpeakpdf));
  
  RooAddition denscaled("denscaled","",RooArgList(sigGFscaled,sigVBFscaled,sigVHscaled,sigTTHscaled,bkgPPscaled,bkgPFscaled,bkgFFscaled));
  RooFormulaVar sigxxb("sigxxb","","log(@0+@1+@2+@3)-log(@4)",RooArgList(sigGFscaled,sigVBFscaled,sigVHscaled,sigTTHscaled,denscaled));
  RooFormulaVar sigFxxb("sigFxxb","","log(@0+@1)-log(@2)",RooArgList(sigGFscaled,sigTTHscaled,denscaled));
  RooFormulaVar sigVxxb("sigVxxb","","log(@0+@1)-log(@2)",RooArgList(sigVBFscaled,sigVHscaled,denscaled));
  RooFormulaVar sigGFxxb("sigGFxxb","","log(@0)-log(@1)",RooArgList(sigGFscaled,denscaled));
  RooFormulaVar sigVBFxxb("sigVBFxxb","","log(@0)-log(@1)",RooArgList(sigVBFscaled,denscaled));
  RooFormulaVar sigVHxxb("sigVHxxb","","log(@0)-log(@1)",RooArgList(sigVHscaled,denscaled));
  RooFormulaVar sigTTHxxb("sigTTHxxb","","log(@0)-log(@1)",RooArgList(sigTTHscaled,denscaled));
  
  
  RooConstVar eterm("eterm","",0.);  
  RooRealVar dummy("dummy","",1.0);  
  
  std::vector<RooAbsData*> vdata;
  vdata.push_back(hdataBkgPP);
  vdata.push_back(hdataBkgPF);
  vdata.push_back(hdataBkgFF);
  vdata.push_back(hdataSigGF);
  vdata.push_back(hdataSigVBF);
  vdata.push_back(hdataSigVH);
  vdata.push_back(hdataSigTTH);
  
  std::vector<RooAbsReal*> vpdf;
  vpdf.push_back(&bkgPPreal);
  vpdf.push_back(&bkgPFreal);
  vpdf.push_back(&bkgFFreal);
  vpdf.push_back(&sigGFreal);
  vpdf.push_back(&sigVBFreal);
  vpdf.push_back(&sigVHreal);
  vpdf.push_back(&sigTTHreal);

  
  {
    RooHybridBDTAutoPdf bdtpdf("bdtpdf","",tgts,eterm,dummy,vdata,vpdf);
    bdtpdf.SetPrescaleInit(10);
    bdtpdf.SetMinCutSignificance(5.0);
    bdtpdf.SetShrinkage(0.1);
    //bdtpdf.SetMinWeightTotal(100.);
    bdtpdf.SetMaxNodes(500.);
    bdtpdf.TrainForest(1e6);
  }
    //return;
    
  RooWorkspace *wsout = new RooWorkspace("wsfitmcfull");
  wsout->import(*hdataAllW);
  wsout->import(*hdataAllWPeak);
  wsout->import(*hdataAll);
  
  wsout->import(sigxxb,RecycleConflictNodes());
  wsout->import(sigFxxb,RecycleConflictNodes());
  wsout->import(sigVxxb,RecycleConflictNodes());
  wsout->import(sigGFxxb,RecycleConflictNodes());
  wsout->import(sigVBFxxb,RecycleConflictNodes());
  wsout->import(sigVHxxb,RecycleConflictNodes());
  wsout->import(sigTTHxxb,RecycleConflictNodes());
  wsout->defineSet("vars",vars,true);
  wsout->defineSet("condvars",condvars,true);
  wsout->defineSet("condvarsdijet",condvarsdijet,true);
  wsout->defineSet("condvarsfull",condvarsfull,true);
  wsout->defineSet("trainingvars",trainingvars,true);
  
    
  wsout->writeToFile("hggfitmcfull.root");

  //wsout->components().Print("V");
  
  RooArgList allcomps(wsout->components());
  for (int iarg=0; iarg<allcomps.getSize(); ++iarg) {
    allcomps.at(iarg)->setOperMode(RooAbsArg::ADirty,false);
  }
  
  RooWorkspace *wsoutsmall = new RooWorkspace("wsfitmc");  
//   for (int iarg=0; iarg<allcomps.getSize(); ++iarg) {
//     wsoutsmall->import(*allcomps.at(iarg),RecycleConflictNodes());
//   }
//   wsoutsmall->defineSet("vars",*wsout->set("vars"),false);
//   wsoutsmall->defineSet("condvars",*wsout->set("condvars"),false);
//   wsoutsmall->defineSet("condvarsdijet",*wsout->set("condvarsdijet"),false);
//   wsoutsmall->defineSet("condvarsfull",*wsout->set("condvarsfull"),false);
//   wsoutsmall->defineSet("trainingvars",*wsout->set("trainingvars"),false);
  
  wsoutsmall->import(sigxxb,RecycleConflictNodes());
  wsoutsmall->import(sigFxxb,RecycleConflictNodes());
  wsoutsmall->import(sigVxxb,RecycleConflictNodes());
  wsoutsmall->import(sigGFxxb,RecycleConflictNodes());
  wsoutsmall->import(sigVBFxxb,RecycleConflictNodes());
  wsoutsmall->import(sigVHxxb,RecycleConflictNodes());
  wsoutsmall->import(sigTTHxxb,RecycleConflictNodes());
  wsoutsmall->defineSet("vars",vars,true);
  wsoutsmall->defineSet("condvars",condvars,true);
  wsoutsmall->defineSet("condvarsdijet",condvarsdijet,true);
  wsoutsmall->defineSet("condvarsfull",condvarsfull,true);
  wsoutsmall->defineSet("trainingvars",trainingvars,true);  
  wsoutsmall->writeToFile("hggfitmc.root");
  
  
}