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
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "TStyle.h"
#include "TH2.h"
#include "TProfile.h"
#include "TMath.h"
 
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
 
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetCanvasDefW(800);
  gStyle->SetCanvasDefH(800);
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.4*gStyle->GetTitleYOffset());
  
  bool useraw = false;
  
  TString dirname = "/data/bendavid/regmassplotsOct14_EE/"; 
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);    
  
  TFile *fin = TFile::Open("/data/bendavid/regflexmasstesting_mod100_small_even/weregmass.root");
  RooWorkspace *ws = (RooWorkspace*)fin->Get("weregmass");
  
  const RooArgSet *varsset = ws->set("vars");
  RooArgList vars(*varsset);
  
  RooRealVar *eratio1 = new RooRealVar("eratio1","(ph1.scrawe+ph1.scpse)/ph1.gene",1.);
  RooRealVar *eratio2 = new RooRealVar("eratio2","(ph2.scrawe+ph2.scpse)/ph2.gene",1.);
  
/*  RooRealVar *eratio1 = new RooRealVar("eratio1","ph1.e/ph1.gene",1.);
  RooRealVar *eratio2 = new RooRealVar("eratio2","ph2.e/ph2.gene",1.); */ 
  
  vars.add(*eratio1);
  vars.add(*eratio2);
  
  RooRealVar weightvar("weightvar","",1.);
  
  TString treeloc = "RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/JetPub/JetCorrectionMod/SeparatePileUpMod/ElectronIDMod/MuonIDMod/PhotonPairSelectorPreselInvertEleVetoNoSmear/PhotonTreeWriterPreselInvertEleVetoNoSmear/hPhotonTree";
     
  TChain *tree = new TChain(treeloc);
  tree->Add("/data/bendavid/diphoTrees8TeVOct6/hgg-2013Final8TeV_s12-zllm50-v7n_noskim.root");

  TChain *treedata = new TChain(treeloc);
  treedata->Add("/data/bendavid/diphoTrees8TeVOct6/hgg-2013Final8TeV_r12a-pho-j22-v1_noskim.root");  
  treedata->Add("/data/bendavid/diphoTrees8TeVOct6/hgg-2013Final8TeV_r12b-dph-j22-v1_noskim.root");  
  treedata->Add("/data/bendavid/diphoTrees8TeVOct6/hgg-2013Final8TeV_r12c-dph-j22-v1_noskim.root");  
  treedata->Add("/data/bendavid/diphoTrees8TeVOct6/hgg-2013Final8TeV_r12d-dph-j22-v1_noskim.root");   

/*  TString treelocs = "RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/JetPub/JetCorrectionMod/SeparatePileUpMod/ElectronIDMod/MuonIDMod/PhotonPairSelectorPreselInvertEleVeto/PhotonTreeWriterPreselInvertEleVeto/hPhotonTree";
     

  TChain *treedatas = new TChain(treelocs);
  treedatas->Add("/data/bendavid/diphoTrees8TeVOct6/hgg-2013Final8TeV_r12a-pho-j22-v1_noskim.root");  
  treedatas->Add("/data/bendavid/diphoTrees8TeVOct6/hgg-2013Final8TeV_r12b-dph-j22-v1_noskim.root");  
  treedatas->Add("/data/bendavid/diphoTrees8TeVOct6/hgg-2013Final8TeV_r12c-dph-j22-v1_noskim.root");  
  treedatas->Add("/data/bendavid/diphoTrees8TeVOct6/hgg-2013Final8TeV_r12d-dph-j22-v1_noskim.root");  */   
  

  
//   new TCanvas;
//   tree->Draw("ph1.ispromptgen+ph2.ispromptgen");
//   
//   new TCanvas;
//   tree->Draw("ph1.gene");
// 
//   new TCanvas;
//   tree->Draw("genmass");  
//   return;
  
  
  TCut selcut = "mass>70. && mass<110. && ph1.pt > 28. && ph2.pt > 20. && evt%2==1 && (!ph1.isbarrel && !ph2.isbarrel)";
  TCut gencut = "ph1.ispromptgen && ph2.ispromptgen";
  
  //selcut = selcut*gencut;
  
  TCut selweight = "xsecweight(procidx)";
  TCut prescale10 = "(evt%10==0)";
  TCut prescale20 = "(evt%20==0)";
  TCut prescale25 = "(evt%25==0)";
  TCut prescale50 = "(evt%50==0)";
  TCut prescale100 = "(evt%100==0)";  
  TCut prescale1000 = "(evt%1000==0)";  
  TCut evenevents = "(evt%2==0)";
  TCut oddevents = "(evt%2==1)";
  

  
  
//   weightvar.SetTitle(selcut);
//   RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",tree,vars,weightvar);   
//   weightvar.SetTitle(selcut);
//   RooDataSet *hdataD = RooTreeConvert::CreateDataSet("hdataD",treedata,vars,weightvar);
  RooRealVar *logmass = ws->var("logmass");
  RooRealVar *mass = ws->var("mass");
  
  
  //RooRealVar *mass = ws->var("mass");
  
  weightvar.SetTitle(selcut);
  logmass->SetTitle("log(mass)-log(91.1876)");
  mass->SetTitle("mass");
  eratio1->SetTitle("ph1.e/ph1.gene");
  eratio2->SetTitle("ph2.e/ph2.gene");
  RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",tree,vars,weightvar);   
  //RooDataSet *hdataD = RooTreeConvert::CreateDataSet("hdataD",treedata,vars,weightvar);
  
  RooDataSet *hdataD = 0;
  
  if (useraw) {
    logmass->SetTitle("0.5*(log(2.0) + log(ph1.scrawe+ph1.scpse) + log(ph2.scrawe+ph2.scpse) + log(1.0-costheta)) - log(91.1876)");
    mass->SetTitle("sqrt(2.0*(ph1.scrawe+ph1.scpse)*(ph2.scrawe+ph2.scpse)*(1.0-costheta))");
    weightvar.SetTitle(selcut);
    eratio1->SetTitle("(ph1.scrawe+ph1.scpse)/ph1.gene");
    eratio2->SetTitle("(ph2.scrawe+ph2.scpse)/ph2.gene");
    hdataD = RooTreeConvert::CreateDataSet("hdataD",tree,vars,weightvar);
  }
  else {
    //RooDataSet *hdataD = RooTreeConvert::CreateDataSet("hdataD",tree,vars,weightvar);
    hdataD = RooTreeConvert::CreateDataSet("hdataD",treedata,vars,weightvar);  
  }
  
  
  RooAbsReal *scalelim_ph1 = ws->function("scalelim_ph1");
  RooAbsReal *scalelim_ph2 = ws->function("scalelim_ph2");
  RooAbsReal *smearlim_ph1 = ws->function("smearlim_ph1");
  RooAbsReal *smearlim_ph2 = ws->function("smearlim_ph2");  
  
  
  RooFormulaVar *scaledmass = new RooFormulaVar("scaledmass","","sqrt(@0*@1)*@2",RooArgList(*scalelim_ph1,*scalelim_ph2,*mass));
  RooRealVar *scaledmassvar = (RooRealVar*)hdataD->addColumn(*scaledmass);  
  
  RooRealVar *e1 = ws->var("var_ph1_0");
  RooRealVar *e2 = ws->var("var_ph2_0");

  RooRealVar *relerr1 = ws->var("var_ph1_1");  
  RooRealVar *relerr2 = ws->var("var_ph2_1");  
  
  RooRealVar *gene_ph1 = ws->var("gene_ph1");
  RooRealVar *gene_ph2 = ws->var("gene_ph2");
  RooRealVar *costheta = ws->var("costheta");
  
  RooRealVar *sceta1 = ws->var("var_ph1_3");
  RooRealVar *sceta2 = ws->var("var_ph2_3");  

  RooRealVar *scphi1 = ws->var("var_ph1_4");
  RooRealVar *scphi2 = ws->var("var_ph2_4");   
  
  RooRealVar *localeta1 = ws->var("var_ph1_38");
  RooRealVar *localeta2 = ws->var("var_ph2_38");    

  RooRealVar *localphi1 = ws->var("var_ph1_39");
  RooRealVar *localphi2 = ws->var("var_ph2_39");      

  RooRealVar *r91 = ws->var("var_ph1_5");
  RooRealVar *r92 = ws->var("var_ph2_5");  
  
//   RooFormulaVar *smearede1 = new RooFormulaVar("smearede1","","(@1>=0.)*(@0*(1.0+sqrt(@1)*RooRandom::gaussian())) + (@1<0.)*(@2 + (@0-@2)*sqrt(@3*@3+@1)/(@3*@3))",RooArgList(*e1,*smearlim_ph1,*gene_ph1,*relerr1);
//   
//   RooFormulaVar *smearede2 = new RooFormulaVar("smearede2","","(@1>=0.)*(@0*(1.0+sqrt(@1)*RooRandom::gaussian())) + (@1<0.)*(@2 + (@0-@2)*sqrt(@3*@3+@1)/(@3*@3))",RooArgList(*e2,*smearlim_ph2,*gene_ph2,*relerr2);
  
  RooFormulaVar *smearedmass = new RooFormulaVar("smearedmass","","@0*(1.0 + 0.5*sqrt(TMath::Max(0.,@1+@2))*RooRandom::gaussian())",RooArgList(*mass,*smearlim_ph1,*smearlim_ph2));
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
  
  RooFormulaVar *ssmear1 = new RooFormulaVar("ssmear1","","sqrt(abs(@0))*@0/abs(@0)",*smearlim_ph1);
  RooFormulaVar *ssmear2 = new RooFormulaVar("ssmear2","","sqrt(abs(@0))*@0/abs(@0)",*smearlim_ph2);
  
  RooRealVar *smearvar1 = (RooRealVar*)hdata->addColumn(*ssmear1);
  RooRealVar *smearvar2 = (RooRealVar*)hdata->addColumn(*ssmear2);

  RooRealVar *scalevar1 = (RooRealVar*)hdataD->addColumn(*scalelim_ph1);
  RooRealVar *scalevar2 = (RooRealVar*)hdataD->addColumn(*scalelim_ph2);  
  
  RooFormulaVar *eratiocor1 = new RooFormulaVar("eratiocor1","","@0*@1",RooArgList(*scalelim_ph1,*eratio1));
  RooFormulaVar *eratiocor2 = new RooFormulaVar("eratiocor2","","@0*@1",RooArgList(*scalelim_ph2,*eratio2));

  RooRealVar *eratiocorvar1 = (RooRealVar*)hdataD->addColumn(*eratiocor1);
  RooRealVar *eratiocorvar2 = (RooRealVar*)hdataD->addColumn(*eratiocor2);  
  
  RooFormulaVar *smearederatio1 = new RooFormulaVar("smearederatio1","","@0*(1.0 + sqrt(TMath::Max(0.,@1))*RooRandom::gaussian())",RooArgList(*eratio1,*smearlim_ph1));  
  RooFormulaVar *smearederatio2 = new RooFormulaVar("smearederatio2","","@0*(1.0 + sqrt(TMath::Max(0.,@1))*RooRandom::gaussian())",RooArgList(*eratio2,*smearlim_ph2));  
  
  RooRealVar *smearederatiovar1 = (RooRealVar*)hdata->addColumn(*smearederatio1);
  RooRealVar *smearederatiovar2 = (RooRealVar*)hdata->addColumn(*smearederatio2);  
  
//   TH1D *hmassStd = new TH1D("hmassStd","",80,70.,110.);
//   treedatas->Draw("mass>>hmassStd",selcut,"goff");
//   hmassStd->SetMarkerColor(kMagenta);
//   hmassStd->SetLineColor(kMagenta);
  
  
  
  hmassMC->GetYaxis()->SetTitle("# of events / 0.5 GeV");
  hmassMC->GetXaxis()->SetTitle("m_{#gamma #gamma} (GeV)");
  
  hscaledmass->GetYaxis()->SetTitle("# of events / 0.5 GeV");
  hscaledmass->GetXaxis()->SetTitle("m_{#gamma #gamma} (GeV)");    
  
  new TCanvas;
  RooPlot *plotsmear1 = smearvar1->frame(-0.1,0.2,100);
  hdata->plotOn(plotsmear1);
  plotsmear1->Draw();

  new TCanvas;
  RooPlot *plotsmear2 = smearvar2->frame(-0.1,0.2,100);
  hdata->plotOn(plotsmear2);
  plotsmear2->Draw();  
  
  TCanvas *cmass = new TCanvas;
  hmassMC->Draw("HIST");
  hsmearedmass->Draw("HISTSAME");
  hscaledmass->Draw("ESAME");
  hmass->Draw("ESAME");
  //hmassStd->Draw("ESAME");
  
  TLegend *legmass = new TLegend(0.55,0.68,0.88,0.88);
  legmass->SetFillStyle(0);
  legmass->AddEntry(hmass, "Data (Regression)","LP");
  legmass->AddEntry(hscaledmass, "Data (Regression + Residual)","LP");
  legmass->AddEntry(hmassMC, "MC (Regression)","L");
  legmass->AddEntry(hscaledmass, "MC (Regression + Smearing)","L");
  legmass->Draw();  
  
  cmass->SaveAs("massall.eps");
  
  
  TCanvas *cmassdata = new TCanvas;
  hscaledmass->Draw("E");
  hmass->Draw("ESAME");
  TLegend *legmassdata = new TLegend(0.55,0.68,0.88,0.88);
  legmassdata->SetFillStyle(0);
  legmassdata->AddEntry(hmass, "Data (Regression)","LP");
  legmassdata->AddEntry(hscaledmass, "Data (Regression + Residual)","LP");
  legmassdata->Draw();    
  
  cmassdata->SaveAs("massdata.eps");  
  
  TH1D *he1 = (TH1D*)hdataD->createHistogram("he1",*eratio1,Binning(80,0.6,1.2));
  TH1D *he2 = (TH1D*)hdataD->createHistogram("he2",*eratio2,Binning(80,0.6,1.2));

  TH1D *hecor1 = (TH1D*)hdataD->createHistogram("hecor1",*eratiocorvar1,Binning(80,0.6,1.2));
  TH1D *hecor2 = (TH1D*)hdataD->createHistogram("hecor2",*eratiocorvar2,Binning(80,0.6,1.2));   
  
  TH1D *he1mc = (TH1D*)hdata->createHistogram("he1mc",*eratio1,Binning(80,0.6,1.2));
  TH1D *he2mc = (TH1D*)hdata->createHistogram("he2mc",*eratio2,Binning(80,0.6,1.2));

  TH1D *hesmear1 = (TH1D*)hdata->createHistogram("hesmear1",*smearederatiovar1,Binning(80,0.6,1.2));
  TH1D *hesmear2 = (TH1D*)hdata->createHistogram("hesmear2",*smearederatiovar2,Binning(80,0.6,1.2));   
  
  TH1D *he = new TH1D( (*he1) + (*he2) );
  TH1D *hecor = new TH1D( (*hecor1) + (*hecor2) );
  TH1D *hemc = new TH1D( (*he1mc) + (*he2mc) );
  TH1D *hesmear = new TH1D( (*hesmear1) + (*hesmear2) );
  
  
  hemc->SetLineColor(kRed);
  hesmear->SetLineColor(kBlue);
//   TH1D *he = he1;
//   TH1D *hecor = hecor1;

  TCanvas *ceratio = new TCanvas;
  he->SetLineColor(kRed);
  hecor->SetLineColor(kBlue);
  he->SetMarkerColor(kRed);
  hecor->SetMarkerColor(kBlue);
//   hemc->GetYaxis()->SetTitle("# of events / 0.5 GeV");
//   hemc->GetXaxis()->SetTitle("m_{#gamma #gamma} (GeV)");
//   
//   hecor->GetYaxis()->SetTitle("# of events / 0.5 GeV");
//   hecor->GetXaxis()->SetTitle("m_{#gamma #gamma} (GeV)");  
  
  hemc->Draw("HIST");
  hesmear->Draw("HISTSAME");
  hecor->Draw("ESAME");
  he->Draw("ESAME");
  
  legmass->Draw();
  ceratio->SaveAs("eratio.eps");
  
  

  
  {
    TCanvas *ceeta = new TCanvas;
    TProfile *prof = static_cast<TH2*>(hdataD->createHistogram("profeta",*sceta1,Binning(100,-2.5,2.5), YVar(*scalevar1,Binning(5000, 0.5,2.0))))->ProfileX("profeta_profx",1,-1,"s");
    prof->SetMinimum(0.96);
    prof->SetMaximum(1.02);
    prof->GetYaxis()->SetTitle("Residual Energy Correction");
    prof->Draw();
    ceeta->SaveAs("eeta.eps");
  }

  {
    TCanvas *ceeta = new TCanvas;
    TProfile *prof = static_cast<TH2*>(hdataD->createHistogram("profphi",*scphi1,Binning(360,-TMath::Pi(),TMath::Pi()), YVar(*scalevar1,Binning(5000, 0.5,2.0))))->ProfileX("profphi_profx",1,-1,"s");
    prof->SetMinimum(0.96);
    prof->SetMaximum(1.02);
    prof->GetYaxis()->SetTitle("Residual Energy Correction");
    prof->Draw();
    ceeta->SaveAs("ephi.eps");
  }  
  
  {
    TCanvas *ceeta = new TCanvas;
    TProfile *prof = static_cast<TH2*>(hdataD->createHistogram("proflocaleta",*localeta1,Binning(100,-0.6,0.6), YVar(*scalevar1,Binning(5000, 0.5,2.0))))->ProfileX("proflocaleta_profx",1,-1,"s");
    prof->SetMinimum(0.96);
    prof->SetMaximum(1.02);
    prof->GetYaxis()->SetTitle("Residual Energy Correction");
    prof->Draw();
    ceeta->SaveAs("elocaleta.eps");
  }
  
  {
    TCanvas *ceeta = new TCanvas;
    TProfile *prof = static_cast<TH2*>(hdataD->createHistogram("proflocalphi",*localphi1,Binning(100,-0.6,0.6), YVar(*scalevar1,Binning(5000, 0.5,2.0))))->ProfileX("proflocalphi_profx",1,-1,"s");
    prof->SetMinimum(0.96);
    prof->SetMaximum(1.02);
    prof->GetYaxis()->SetTitle("Residual Energy Correction");
    prof->Draw();
    ceeta->SaveAs("elocalphi.eps");
  }      
  
  {
    TCanvas *ceeta = new TCanvas;
    TProfile *prof = static_cast<TH2*>(hdataD->createHistogram("profr9",*r91,Binning(100,0.,1.), YVar(*scalevar1,Binning(5000, 0.5,2.0))))->ProfileX("profr9_profx",1,-1,"s");
    prof->SetMinimum(0.96);
    prof->SetMaximum(1.02);
    prof->GetYaxis()->SetTitle("Residual Energy Correction");
    prof->Draw();
    ceeta->SaveAs("er9.eps");
  }   
  
  RooRealVar mfit("mfit","",1.);
  mfit.setRange(70.,110.);
  mfit.setBins(320);
  mfit.setBins(10e3,"cache");
  
  RooRealVar mucb("mucb","",0.);
  RooRealVar sigmacb("sigmacb","",2.);
  RooRealVar ncb("ncb","",2.);
  RooRealVar alphacb("alphacb","",1.);
  
  mucb.setConstant(false);
  sigmacb.setConstant(false);
  ncb.setConstant(false);
  alphacb.setConstant(false);
  
  const double mpdg = 91.1876;
  const double wpdg = 2.4952;
  
  RooBreitWigner bwpdf("bwpdf","",mfit,RooConst(mpdg),RooConst(wpdg));
  RooCBShape cbpdf("cbpdf","",mfit,mucb, sigmacb, alphacb,ncb);
  RooFFTConvPdf zpdf("zpdf","",mfit, bwpdf,cbpdf);
  
  
  TH1 *hmassfine = hdataD->createHistogram("hmassfine",*mass,Binning(320,70.,110.));
  TH1 *hscaledmassfine = hdataD->createHistogram("hscaledmassfine",*scaledmassvar,Binning(320,70.,110.));

  TH1 *hmassMCfine = hdata->createHistogram("hmassMCfine",*mass,Binning(320,70.,110.));
  TH1 *hsmearedmassfine = hdata->createHistogram("hsmearedmassfine",*smearedmassvar,Binning(320,70.,110.));    
  
  RooDataHist *dhistmc = new RooDataHist("dhistmc","",mfit,hmassMCfine);
  RooDataHist *dhistd = new RooDataHist("dhistd","",mfit,hmassfine);
  RooDataHist *dhistdcor = new RooDataHist("dhistdcor","",mfit,hscaledmassfine);
  
  zpdf.fitTo(*dhistmc);
  double relerrmc = sigmacb.getVal()/(mpdg+mucb.getVal());  
  
  alphacb.setConstant();
  ncb.setConstant();
  
  {
    new TCanvas;
    RooPlot *plotmfit = mfit.frame(80);
    dhistmc->plotOn(plotmfit);
    zpdf.plotOn(plotmfit);
    plotmfit->Draw();
  }
  
  zpdf.fitTo(*dhistdcor);  
  double relerrdcor = sigmacb.getVal()/(mpdg+mucb.getVal());

  {
    new TCanvas;
    RooPlot *plotmfit = mfit.frame(80);
    dhistdcor->plotOn(plotmfit);
    zpdf.plotOn(plotmfit);
    plotmfit->Draw();
  }  
  
  zpdf.fitTo(*dhistd);  
  double relerrd = sigmacb.getVal()/(mpdg+mucb.getVal());  
  
  {
    new TCanvas;
    RooPlot *plotmfit = mfit.frame(80);
    dhistd->plotOn(plotmfit);
    zpdf.plotOn(plotmfit);
    plotmfit->Draw();
  }  
  
  printf("relerrmc = %5f, relerrdcor = %5f, relerrd = %5f\n",relerrmc,relerrdcor,relerrd);
  
  
  
}
