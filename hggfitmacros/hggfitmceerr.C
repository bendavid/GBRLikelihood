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
#include "RooRevCBFast.h"
#include "RooBifurGauss.h"

 using namespace RooFit;
 
 
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

float xsecweights[50];
float xsecweight(int procidx=0) {
  return xsecweights[procidx];
}

float kfact(int procidx, bool ispromptgen1, bool ispromptgen2) {
 
  if (procidx==0) return 0.975; //gf
  
  if (procidx>=4 && procidx<=8) {
    int npromptgen = int(ispromptgen1) + int(ispromptgen2);
    
    if (npromptgen==0) return 1.0; //fake fake
    //else if (npromptgen==1) return 1.3; //prompt-fake
    else if (npromptgen==1) return 1.0; //prompt-fake with sherpa
    else if (npromptgen==2) return 1.0; //prompt-prompt (for sherpa diphoj)
  }
  
  return 1.0;
  
}



void initweights(TChain *chain, float *xsecs, float lumi) {
 
  TObjArray *files = chain->GetListOfFiles();
  for (int i=0; i<files->GetEntries(); ++i) {    
    TFile *file = TFile::Open(files->At(i)->GetTitle(),"READ");
    
    //puweights[i] = getpuweights(file,target);
    xsecweights[i] = getweight(file,lumi*xsecs[i]);
    
    file->Close();    
  } 
  
  chain->SetAlias("procidx","This->GetTreeNumber()");
  
}
 
 
void hggfitmceerr(double nommass=123., double tgtr=1., int ijob=0) {
    
  //gSystem->cd("/scratch/bendavid/root/bare/fitplotsJun10test/");
  
  int seed = 65539+ijob+1; 
  
  TString dirname = "/scratch/bendavid/root/bare/hggfiteerrtestall_large2/";
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
  TCut selcut = "(ph1.pt > (mass/3.0) && ph2.pt > (mass/4.0) && mass>100. && mass<180. && ph1.idmva>-0.2 && ph2.idmva>-0.2)";
  //TCut selweight = "xsecweight(procidx)*puweight(numPU,procidx)";
  TCut selweight = "xsecweight(procidx)*mcweight*kfact(procidx,ph1.ispromptgen,ph2.ispromptgen)";
  
  TCut sigFcut = "(procidx==0 || procidx==3)";
  TCut sigVcut = "(procidx==1 || procidx==2)";
  
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
  
  
  
  TCut fcut = prescale50;  
  
  
  float xsecs[50];

  
  //TCut selcutsingle = "ph.pt>25. && ph.isbarrel && ph.ispromptgen";
  //TCut selcutsingle = "ph.pt>25.&& ph.ispromptgen";
  TCut selcutsingle = "ph.genpt>16.&& ph.ispromptgen";
  TCut selweightsingle = "xsecweight(procidx)";
  
  
//   TChain *tree = new TChain("RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/JetPub/JetCorrectionMod/SeparatePileUpMod/ElectronIDMod/MuonIDMod/PhotonPairSelectorPresel/PhotonTreeWriterPresel/hPhotonTreeSingle");
//   tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV/merged/hgg-2013Final8TeV_s12-diphoj-v7n_noskim.root");

  TChain *tree = new TChain("RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/PhotonIDModPresel/PhotonTreeWriterSingle/hPhotonTreeSingle");
  tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_reg_trans/merged/hgg-2013Final8TeV_reg_trans_s12-pj20_40-2em-v7n_noskim.root");
  tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_reg_trans/merged/hgg-2013Final8TeV_reg_trans_s12-pj40-2em-v7n_noskim.root");
  
  xsecs[0] = 0.001835*81930.0;
  xsecs[1] = 0.05387*8884.0;
  initweights(tree,xsecs,1.);
  
  double weightscale = xsecweights[1];
  xsecweights[0] /= weightscale;
  xsecweights[1] /= weightscale;  
 
  
  
  tree->SetCacheSize(64*1024*1024);
  
  

  RooRealVar energy("energy","ph.e",0);
  RooRealVar sceta("sceta","ph.sceta",0.);
  RooRealVar idmva("idmva","ph.idmva",0.,-1.,1.);  
  RooRealVar eerr("eerr","(ph.isbarrel + 0.5*!ph.isbarrel)*ph.eerr/ph.e",0.);
  RooRealVar evt("evt","evt",0.);
  
  
  RooArgList vars;
  vars.add(energy);
  vars.add(sceta);
  //vars.add(idmva);
  
  
  RooArgList condvars(vars);
  
  vars.add(eerr);  
  
  RooArgList condvarsid(vars);
  
  vars.add(idmva);
  vars.add(evt);
  

//   RooPowerLaw("testpow","",pt1,pt2);
//   return;
  

//    new TCanvas;
//    tree->Draw("mass>>htmpall(80,100.,180.)",bkgcut*selcut*selweight,"HIST");
// 
//    new TCanvas;
//    tree->Draw("mass>>htmpallid(80,100.,180.)",idcut*bkgcut*selcut*selweight,"HIST");   
//    
//    new TCanvas;
//    tree->Draw("mass>>htmp(80,100.,180.)",bkgcutnoq*selcut*selweight,"HIST");   
// //   
//    return;
    
  
  
  //RooRealVar weightvar("weightvar","(ph1.pt > (mass/3.0) && ph2.pt > (mass/4.0) && mass>100. && mass<180. && ph1.idmva>-0.2 && ph2.idmva>-0.2 && evt%100!=0)",1.);
  RooRealVar weightvar("weightvar","",1.);
  //RooRealVar weightvar("weightvar","(ph1.pt > (mass/3.0) && ph2.pt > (mass/4.0) && mass>100. && mass<180. && ph1.idmva>-0.2 && ph2.idmva>-0.2)",1.);
  
  weightvar.SetTitle(selcutsingle*selweightsingle);
  RooDataSet *hdataSingle = RooTreeConvert::CreateDataSet("hdataSingle",tree,vars,weightvar);    
  
  int ngauseerr = 4;
  //int nparmseerr = 3*ngauseerr + 2;
  //int nparmseerr = 3*ngauseerr + 2;
  int nparmseerr = 5*ngauseerr;
  
  RooArgList tgtseerr;
  RooGBRFunction funceerr("funceerr","",condvars,nparmseerr);
 
  int iparmeerr = 0;  
  
  
  RooArgList eerrgauspdfs;
  RooArgList eerrgauscoeffs;
  
  double stepeerr = 0.07/double(std::max(1,ngauseerr-1));
  
//    RooRealVar *gmeanvar = new RooRealVar(TString::Format("gmeanvar_eerr_%i",0),"",0.007+stepeerr*0);
//    RooRealVar *gsigmavar = new RooRealVar(TString::Format("gsigmavar_eerr_%i",0),"",0.01);
// //   //RooRealVar *gsigmaRvar = new RooRealVar(TString::Format("gsigmaRvar_eerr_%i",0),"",0.02);
// //   
// //   //if (0==0) gmeanvar->setVal(0.007);
// //   
//    gmeanvar->setConstant(false);
//    gsigmavar->setConstant(false);
// //   //gsigmaRvar->setConstant(false);
// //   
// //   
//    RooGBRTarget *gmean = new RooGBRTarget(TString::Format("gmean_eerr_%i",0),"",funceerr,iparmeerr++,*gmeanvar);
//    RooGBRTarget *gsigma = new RooGBRTarget(TString::Format("gsigma_eerr_%i",0),"",funceerr,iparmeerr++,*gsigmavar);
// //   //RooGBRTarget *gsigmaR = new RooGBRTarget(TString::Format("gsigmaR_eerr_%i",0),"",funceerr,iparmeerr++,*gsigmaRvar);
// // 
//    RooRealConstraint *gmeanlim = new RooRealConstraint(TString::Format("gmeanlim_eerr_%i",0),"",*gmean,0.,0.5);   
//    RooRealConstraint *gsigmalim = new RooRealConstraint(TString::Format("gsigmalim_eerr_%i",0),"",*gsigma,1e-7,0.5);
//   //RooRealConstraint *gsigmaRlim = new RooRealConstraint(TString::Format("gsigmaRlim_eerr_%i",0),"",*gsigmaR,1e-7,0.2);  
//   
//   tgtseerr.add(*gmean);
//   tgtseerr.add(*gsigma);
  
  for (int igaus=0; igaus<ngauseerr; ++igaus) {
    RooRealVar *gmeanvar = new RooRealVar(TString::Format("gmeanvar_eerr_%i",igaus),"",0.007+stepeerr*igaus);
    RooRealVar *gsigmavar = new RooRealVar(TString::Format("gsigmavar_eerr_%i",igaus),"",0.01);
    RooRealVar *galphavar = new RooRealVar(TString::Format("galphavar_eerr_%i",igaus),"",1.0);
    RooRealVar *gnvar = new RooRealVar(TString::Format("gnvar_eerr_%i",igaus),"",2.);    
    RooRealVar *gfracvar = new RooRealVar(TString::Format("gfracvar_eerr_%i",igaus),"",1.0);  
    
    //if (igaus==0) gmeanvar->setVal(0.007);
    
    gmeanvar->setConstant(false);
    gsigmavar->setConstant(false);
    galphavar->setConstant(false);
    gnvar->setConstant(false);    
    gfracvar->setConstant(false);
    
    
    RooGBRTarget *gmean = new RooGBRTarget(TString::Format("gmean_eerr_%i",igaus),"",funceerr,iparmeerr++,*gmeanvar);
    RooGBRTarget *gsigma = new RooGBRTarget(TString::Format("gsigma_eerr_%i",igaus),"",funceerr,iparmeerr++,*gsigmavar);
    RooGBRTarget *galpha = new RooGBRTarget(TString::Format("galpha_eerr_%i",igaus),"",funceerr,iparmeerr++,*galphavar);
    RooGBRTarget *gn = new RooGBRTarget(TString::Format("gn_eerr_%i",igaus),"",funceerr,iparmeerr++,*gnvar);    
    RooGBRTarget *gfrac = new RooGBRTarget(TString::Format("gfrac_eerr_%i",igaus),"",funceerr,iparmeerr++,*gfracvar);

    RooRealConstraint *gmeanlim = new RooRealConstraint(TString::Format("gmeanlim_eerr_%i",igaus),"",*gmean,0.,0.5);   
    RooRealConstraint *gsigmalim = new RooRealConstraint(TString::Format("gsigmalim_eerr_%i",igaus),"",*gsigma,1e-5,0.1);
    RooRealConstraint *galphalim = new RooRealConstraint(TString::Format("galphalim_eerr_%i",igaus),"",*galpha,0.05,8.);
    RooRealConstraint *gnlim = new RooRealConstraint(TString::Format("gnlim_eerr_%i",igaus),"",*gn,1.01,5000.);
    //RooRealConstraint *gfraclim = new RooRealConstraint(TString::Format("gfraclim_eerr_%i",igaus),"",*gfrac,0.,1.);
    RooAbsReal *gfraclim = new RooProduct(TString::Format("gfraclim_eerr_%i",igaus),"",RooArgList(*gfrac,*gfrac));
 
    
    if (igaus==0) {
      gfraclim = new RooConstVar(TString::Format("gfraclimconst_eerr_%i",igaus),"",1.);
    }
    else {
      tgtseerr.add(*gfrac);   
    }
    
    //RooGaussianFast *gpdf = new RooGaussianFast(TString::Format("gdf_eerr_%i",igaus),"",eerr,*gmeanlim,*gsigmalim);
    //RooBifurGauss *gpdf = new RooBifurGauss(TString::Format("gdf_eerr_%i",igaus),"",eerr,*gmeanlim,*gsigmalim,*galphalim);
    
    
    if (igaus==0) {
      RooRevCBFast *gpdf = new RooRevCBFast(TString::Format("gdf_eerr_%i",igaus),"",eerr,*gmeanlim,*gsigmalim,*galphalim, *gnlim);
    
      tgtseerr.add(*gmean);
      tgtseerr.add(*gsigma);
      tgtseerr.add(*galpha);
      tgtseerr.add(*gn);
      
      eerrgauspdfs.add(*gpdf);      
    
    }
    else {
      RooGaussianFast *gpdf = new RooGaussianFast(TString::Format("gdf_eerr_%i",igaus),"",eerr,*gmeanlim,*gsigmalim);
    
      tgtseerr.add(*gmean);
      tgtseerr.add(*gsigma);
      
      eerrgauspdfs.add(*gpdf);

    }      
      
    
    eerrgauscoeffs.add(*gfraclim);    
    
  }
  RooCondAddPdf eerrpdf("eerrpdf","",eerrgauspdfs,eerrgauscoeffs);  
  
  
  RooAbsPdf *pdf0 = static_cast<RooAbsPdf*>(eerrgauspdfs.at(0));
  
  
  int ngaus = 6;
  int nparms = 4*ngaus;
  
  RooArgList tgtsid;
  RooGBRFunction funcid("funcid","",condvarsid,nparms);  
  
  RooArgList gauspdfs;
  RooArgList gauscoeffs;
  
  double step = 0.5/double(std::max(1,ngaus-1));
  
  int iparm = 0;
  for (int igaus=0; igaus<ngaus; ++igaus) {
    RooRealVar *gmeanvar = new RooRealVar(TString::Format("gmeanvar_%i",igaus),"",-0.2+step*igaus);
    RooRealVar *gsigmavar = new RooRealVar(TString::Format("gsigmavar_%i",igaus),"",0.1);
    RooRealVar *gsigmaRvar = new RooRealVar(TString::Format("gsigmaRvar_%i",igaus),"",0.1);
    RooRealVar *gfracvar = new RooRealVar(TString::Format("gfracvar_%i",igaus),"",1.0);  
    
    gmeanvar->setConstant(false);
    gsigmavar->setConstant(false);
    gsigmaRvar->setConstant(false);
    gfracvar->setConstant(false);
    
    RooGBRTarget *gmean = new RooGBRTarget(TString::Format("gmean_%i",igaus),"",funcid,iparm++,*gmeanvar);
    RooGBRTarget *gsigma = new RooGBRTarget(TString::Format("gsigma_%i",igaus),"",funcid,iparm++,*gsigmavar);
    RooGBRTarget *gsigmaR = new RooGBRTarget(TString::Format("gsigmaR_%i",igaus),"",funcid,iparm++,*gsigmaRvar);
    RooGBRTarget *gfrac = new RooGBRTarget(TString::Format("gfrac_%i",igaus),"",funcid,iparm++,*gfracvar);

    RooRealConstraint *gmeanlim = new RooRealConstraint(TString::Format("gmeanlim_%i",igaus),"",*gmean,-1.,1.);   
    RooRealConstraint *gsigmalim = new RooRealConstraint(TString::Format("gsigmalim_%i",igaus),"",*gsigma,1e-4,2.);
    RooRealConstraint *gsigmaRlim = new RooRealConstraint(TString::Format("gsigmaRlim_%i",igaus),"",*gsigmaR,1e-4,2.);
    //RooRealConstraint *gfraclim = new RooRealConstraint(TString::Format("gfraclim_%i",igaus),"",*gfrac,0.,1.);
    
    RooAbsReal *gfraclim = new RooProduct(TString::Format("gfraclim_%i",igaus),"",RooArgList(*gfrac,*gfrac));
 
    
    if (igaus==0) {
      gfraclim = new RooConstVar(TString::Format("gfraclimconst_%i",igaus),"",1.);
    }
    else {
      tgtsid.add(*gfrac);   
    }    
    
    
    RooGaussianFast *gpdf = new RooGaussianFast(TString::Format("gdf_%i",igaus),"",idmva,*gmeanlim,*gsigmalim);
    //RooBifurGauss *gpdf = new RooBifurGauss(TString::Format("gdf_%i",igaus),"",idmva,*gmeanlim,*gsigmalim,*gsigmaRlim);
    
    gauspdfs.add(*gpdf);
    gauscoeffs.add(*gfraclim);
    
    tgtsid.add(*gmean);
    tgtsid.add(*gsigma);
    //tgtsid.add(*gsigmaR);
    //tgtsid.add(*gfrac);    
  }
  RooCondAddPdf idpdf("idpdf","",gauspdfs,gauscoeffs);
  
  RooConstVar etermconst("etermconst","",0.);  
  RooAbsReal &eterm = etermconst;
  RooRealVar dummy("dummy","",1.0);
   

  std::vector<RooAbsData*> vdata;
  vdata.push_back(hdataSingle);

  
  std::vector<RooAbsReal*> vpdf;
  vpdf.push_back(&eerrpdf);
  //vpdf.push_back(pdf0);

  std::vector<RooAbsReal*> vpdfid;
  vpdfid.push_back(&idpdf);  
 
  
  RooHybridBDTAutoPdf bdtpdf("bdtpdf","",funceerr,tgtseerr,eterm,dummy,vdata,vpdf);
  bdtpdf.SetPrescaleInit(100);
  bdtpdf.SetMinCutSignificance(5.0);
  bdtpdf.SetShrinkage(0.1);
  bdtpdf.SetMinWeightTotal(200.);
  bdtpdf.SetMaxNodes(200);
  bdtpdf.TrainForest(1e6);   
  
  RooHybridBDTAutoPdf bdtpdfid("bdtpdfid","",funcid,tgtsid,eterm,dummy,vdata,vpdfid);
  bdtpdfid.SetPrescaleInit(100);
  bdtpdfid.SetMinCutSignificance(5.0);
  bdtpdfid.SetShrinkage(0.1);
  bdtpdfid.SetMinWeightTotal(200.);
  bdtpdfid.SetMaxNodes(200);
  bdtpdfid.TrainForest(1e6);    
  
  
  RooAbsReal *finalcdferr = eerrpdf.createCDF(eerr);
  
  RooFormulaVar transerr("transerr","","sqrt(2.)*TMath::ErfInverse(2.*@0-1.)",*finalcdferr);
  


  RooAbsReal *finalcdfid = idpdf.createCDF(idmva);
  
  RooFormulaVar transid("transid","","sqrt(2.)*TMath::ErfInverse(2.*@0-1.)",*finalcdfid);
  
  
  RooWorkspace *wsout = new RooWorkspace("wsfiteerr");
  wsout->import(*hdataSingle);
  
  wsout->import(eerrpdf,RecycleConflictNodes());
  wsout->import(idpdf,RecycleConflictNodes());
//   wsout->import(transerr,RecycleConflictNodes());
//   wsout->import(transid,RecycleConflictNodes());
  
  wsout->defineSet("datavars",vars,true);
    
  wsout->writeToFile("hggfiteerr.root");  
  
  
    
  
  RooRealVar *cdfidvar = (RooRealVar*)hdataSingle->addColumn(*finalcdfid);    
  RooRealVar *transidvar = (RooRealVar*)hdataSingle->addColumn(transid);
    
  RooGaussianFast unormpdfid("unormpdfid","",*transidvar,RooConst(0.),RooConst(1.));    

  RooRealVar *cdferrvar = (RooRealVar*)hdataSingle->addColumn(*finalcdferr);    
  RooRealVar *transerrvar = (RooRealVar*)hdataSingle->addColumn(transerr);
    
  RooGaussianFast unormpdferr("unormpdferr","",*transerrvar,RooConst(0.),RooConst(1.));    
  
  
  //RooDataSet *testdata = (RooDataSet*)hdataSingle->reduce("abs(sceta)>1.3 && abs(sceta)<1.4");
  RooDataSet *testdata = hdataSingle;
    
  
  new TCanvas;
  RooPlot *eerrplot = eerr.frame(0.,0.1,200);
  testdata->plotOn(eerrplot);
  eerrpdf.plotOn(eerrplot,ProjWData(*testdata));
  eerrplot->Draw();    
  

  
  new TCanvas;
  RooPlot *transplot = transerrvar->frame(-5.,5.,100);
  hdataSingle->plotOn(transplot);
  unormpdferr.plotOn(transplot);
  transplot->Draw();
  //return;
  
  new TCanvas;
  RooPlot *cdfploterr = cdferrvar->frame(0.,1.,100);
  hdataSingle->plotOn(cdfploterr);
  //unormpdf.plotOn(transplot);
  cdfploterr->Draw();
  //return;    
  
  
  
  new TCanvas;
  RooPlot *idplot = idmva.frame(-1.,1.,200);
  testdata->plotOn(idplot);
  idpdf.plotOn(idplot,ProjWData(*testdata));
  idplot->Draw();  
  

  
  new TCanvas;
  RooPlot *transplotid = transidvar->frame(-5.,5.,100);
  testdata->plotOn(transplotid);
  unormpdfid.plotOn(transplotid);
  transplotid->Draw();
  //return;
  
  new TCanvas;
  RooPlot *cdfplotid = cdfidvar->frame(0.,1.,100);
  testdata->plotOn(cdfplotid);
  //unormpdf.plotOn(transplot);
  cdfplotid->Draw();
  //return;        
  

  TH1 *herrid = testdata->createHistogram("herrid",eerr,Binning(30,0.,0.1), YVar(idmva,Binning(30,-0.5,0.6)));
  TH1 *herre = testdata->createHistogram("herre",energy,Binning(30,0.,200.), YVar(eerr,Binning(30,0.,0.1)));
  TH1 *hideta = testdata->createHistogram("hideta",sceta,Binning(40,-2.5,2.5), YVar(idmva,Binning(30,-0.5,0.6)));

  TH1 *herridtrans = testdata->createHistogram("herridtrans",*transerrvar,Binning(30,-5.,5.), YVar(*transidvar,Binning(30,-5.,5.)));
  TH1 *herrtranse = testdata->createHistogram("herrtranse",energy,Binning(30,0.,200.), YVar(*transerrvar,Binning(30,-5.,5.)));
  TH1 *hidtranseta = testdata->createHistogram("hidtranseta",sceta,Binning(40,-2.5,2.5), YVar(*transidvar,Binning(30,-5.,5.)));  

  new TCanvas;
  herrid->Draw("COLZ");

  new TCanvas;
  herre->Draw("COLZ");
  
  new TCanvas;
  hideta->Draw("COLZ");    
  
  
  new TCanvas;
  herridtrans->Draw("COLZ");

  new TCanvas;
  herrtranse->Draw("COLZ");
  
  new TCanvas;
  hidtranseta->Draw("COLZ");  
  

  
//   new TCanvas;
//   RooRealVar *meanvar = (RooRealVar*)hdataSingle->addColumn(eerrmeanlim);
//   RooPlot *meanplot = meanvar->frame(0.,0.1,200);
//   hdataSingle->plotOn(meanplot);
//   meanplot->Draw();
  

  return;
  
}