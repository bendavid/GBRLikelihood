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

void eregtraining_masstest(bool dobarrel, bool doele) {
   
//   gSystem->Setenv("OMP_WAIT_POLICY","PASSIVE");
  
  //candidate to set fixed alpha values (0.9,3.8)
  //TString dirname = TString::Format("/afs/cern.ch/work/b/bendavid/bare/eregtesteleJul30_sig5_01_alphafloat5_%i/",int(minevents)); 
  
  //TString dirname = "/afs/cern.ch/work/b/bendavid/bare/eregAug23test/";
  TString dirname = "/scratch/bendavid/root/bare/regflexmasstesting_mod100/"; 
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);  
  
  std::vector<std::string> *varsf = new std::vector<std::string>;
  varsf->push_back("ph.e");
  varsf->push_back("ph.eerr/ph.e");
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
  
  varseb->push_back("ph.e/ph.scrawe");  
  varseb->push_back("ph.e5x5seed/ph.eseed");
  
  varseb->push_back("ph.ietaseed");
  varseb->push_back("ph.iphiseed");
  varseb->push_back("(ph.ietaseed-1*abs(ph.ietaseed)/ph.ietaseed)%5");
  varseb->push_back("(ph.iphiseed-1)%2");       
  varseb->push_back("(abs(ph.ietaseed)<=25)*((ph.ietaseed-1*abs(ph.ietaseed)/ph.ietaseed)%25) + (abs(ph.ietaseed)>25)*((ph.ietaseed-26*abs(ph.ietaseed)/ph.ietaseed)%20)");
  varseb->push_back("(ph.iphiseed-1)%20"); 
  varseb->push_back("ph.etacryseed");
  varseb->push_back("ph.phicryseed");

  varsee->push_back("ph.e/(ph.scrawe + ph.scpse)");  
  varsee->push_back("ph.scpse/ph.scrawe");
    
  std::vector<std::string> *varslist;
  if (dobarrel) varslist = varseb;
  else varslist = varsee;
  
  RooArgList vars_ph1;
  RooArgList vars_ph2;  
  RooArgList vars_all;  
  //RooArgList vars_common;
  for (unsigned int ivar=0; ivar<varslist->size(); ++ivar) {
    TString sname = varslist->at(ivar).c_str();
    TString title_ph1(sname);
    title_ph1.ReplaceAll("ph.","ph1.");
    TString title_ph2(sname);
    title_ph2.ReplaceAll("ph.","ph2.");
    
    printf("ivar = %i, sname = %s, t1 = %s, t2 = %s\n",ivar,sname.Data(),title_ph1.Data(),title_ph2.Data());
    if (title_ph1==title_ph2) {
      RooRealVar *var = new RooRealVar(TString::Format("var_common_%i",ivar), sname ,0.);
      vars_ph1.add(*var);
      vars_ph2.add(*var);
      vars_all.add(*var);
    }
    else {
      RooRealVar *var_ph1 = new RooRealVar(TString::Format("var_ph1_%i",ivar), title_ph1 ,0.);
      RooRealVar *var_ph2 = new RooRealVar(TString::Format("var_ph2_%i",ivar), title_ph2 ,0.);
      vars_ph1.add(*var_ph1);
      vars_ph2.add(*var_ph2);
      vars_all.add(*var_ph1);
      vars_all.add(*var_ph2);      
    }
  }
  
  vars_ph1.Print("V");
  vars_ph2.Print("V");
  vars_all.Print("V");
  
//  RooArgList condvars(vars);
  

  
  
//   RooRealVar *tgtvar = new RooRealVar("tgtvar","ph.scrawe/ph.gene",1.);
//   if (!dobarrel) tgtvar->SetTitle("(ph.scrawe + ph.scpse)/ph.gene");
  
/*  RooRealVar *tgtvar = new RooRealVar("tgtvar","sqrt(2.0*ph1.scrawe*ph2.scrawe*(1.0-costheta)",1.);
  if (!dobarrel) tgtvar->SetTitle("sqrt(2.0*(ph1.scrawe + ph1.scpse)*(ph2.scrawe + ph2.scpse)*(1.0-costheta))"); */ 

  RooRealVar *tgtvar = new RooRealVar("tgtvar","mass",1.);

  const double pdgzmass = 91.1876;
  const double logpdgzmass = log(pdgzmass);
  
  RooRealVar *mass = new RooRealVar("mass","mass",90.,70.,110.);
  
  RooRealVar *logmass = new RooRealVar("logmass",TString::Format("log(mass)-%5f",logpdgzmass),0.,log(70.)-logpdgzmass,log(110.)-logpdgzmass);
  RooRealVar *relpt_ph1 = new RooRealVar("relpt_ph1","ph1.pt/mass",1.);
  RooRealVar *relpt_ph2 = new RooRealVar("relpt_ph2","ph2.pt/mass",1.);
  RooRealVar *eta_ph1 = new RooRealVar("eta_ph1","ph1.eta",1.);
  RooRealVar *eta_ph2 = new RooRealVar("eta_ph2","ph2.eta",1.);
  RooRealVar *cosdphi = new RooRealVar("cosdphi","TMath::Cos(ph1.phi-ph2.phi)",1.);
  
/*  RooRealVar *tgtvar = new RooRealVar("tgtvar","log(ph.gene/ph.scrawe)",1.);
  if (!dobarrel) tgtvar->SetTitle("log(ph.gene/(ph.scrawe + ph.scpse))");  */  
  
  //tgtvar->setRange(0.8,2.);
  
  RooArgList condvarsmass;
  condvarsmass.add(*relpt_ph1);
  condvarsmass.add(*relpt_ph2);
  condvarsmass.add(*eta_ph1);
  condvarsmass.add(*eta_ph2);
  condvarsmass.add(*cosdphi);

  RooArgList condvars;
  condvars.add(vars_all);
  condvars.add(condvarsmass);  
  
  RooArgList condvarssingle;
  condvarssingle.add(vars_ph1);
  condvarssingle.add(condvarsmass);   
  
  RooArgList vars;
  vars.add(vars_all);
  vars.add(*logmass);
  vars.add(*mass);
  vars.add(condvarsmass);
  
  condvarsmass.add(*vars_ph1.at(1));
  condvarsmass.add(*vars_ph2.at(1));


  
 
  
  //varstest.add(*tgtvar);
    
  RooRealVar weightvar("weightvar","",1.);

  //TFile *fdin = TFile::Open("/home/mingyang/cms/hist/hgg-2013Moriond/merged/hgg-2013Moriond_s12-diphoj-3-v7a_noskim.root");
//   TFile *fdin = TFile::Open("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/trainingtreesJul1/hgg-2013Final8TeV_s12-zllm50-v7n_noskim.root");
//   TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("PhotonTreeWriterSingleInvert");
//   TTree *dtree = (TTree*)ddir->Get("hPhotonTreeSingle");    
  
/*  TFile *fdinsig = TFile::Open("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/trainingtreesJul1/hgg-2013Moriond_s12-h125gg-gf-v7a_noskim.root");
  TDirectory *ddirsig = (TDirectory*)fdinsig->FindObjectAny("PhotonTreeWriterPreselNoSmear");
  TTree *dtreesig = (TTree*)ddirsig->Get("hPhotonTreeSingle");     */ 
  
  
  
  TString treeloc = "RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/JetPub/JetCorrectionMod/SeparatePileUpMod/ElectronIDMod/MuonIDMod/PhotonPairSelectorPreselInvertEleVetoNoSmear/PhotonTreeWriterPreselInvertEleVetoNoSmear/hPhotonTree";
     
  TChain *tree = new TChain(treeloc);
  tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV/merged/hgg-2013Final8TeV_s12-zllm50-v7n_noskim.root");

  TChain *treedata = new TChain(treeloc);
  treedata->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV/merged/hgg-2013Final8TeV_r12a-pho-j22-v1_noskim.root");  
  treedata->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV/merged/hgg-2013Final8TeV_r12b-dph-j22-v1_noskim.root");  
  treedata->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV/merged/hgg-2013Final8TeV_r12c-dph-j22-v1_noskim.root");  
  treedata->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV/merged/hgg-2013Final8TeV_r12d-dph-j22-v1_noskim.root");   

    
  
//   new TCanvas;
//   tree->Draw("ph1.pt>>htmp1(100,0.,100.)");
// 
//   new TCanvas;
//   tree->Draw("ph2.pt>>htmp2(100,0.,100.)");  
// 
//   new TCanvas;
//   tree->Draw("ph1.pt/mass>>htmp3(100,0.,2.)");    
// 
//   new TCanvas;
//   tree->Draw("ph2.pt/mass>>htmp4(100,0.,2.)");      
//   
//   return;
  

  
  TCut selcut = "mass>70. && mass<110. && (ph1.pt/mass) > (20./70.) && (ph2.pt/mass) > (20./70.) && ph1.isbarrel && ph2.isbarrel && evt%100==0";

  
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
  
  
  RooRealVar *scalevar_ph1 = new RooRealVar("scalevar_ph1","",1.);
  RooRealVar *scalevar_ph2 = new RooRealVar("scalevar_ph2","",1.);
  RooRealVar *smearvar_ph1 = new RooRealVar("smearvar_ph1","",pow(0.005,2));
  RooRealVar *smearvar_ph2 = new RooRealVar("smearvar_ph2","",pow(0.005,2));
  
  scalevar_ph1->setConstant(false);
  scalevar_ph2->setConstant(false);
  smearvar_ph1->setConstant(false);
  smearvar_ph2->setConstant(false);
  
  RooGBRFunctionFlex *scalefunc = new RooGBRFunctionFlex("scalefunc","");
  RooGBRFunctionFlex *smearfunc = new RooGBRFunctionFlex("smearfunc","");
  
  RooGBRTargetFlex *scale_ph1 = new RooGBRTargetFlex("scale_ph1","",*scalefunc,*scalevar_ph1,vars_ph1);
  RooGBRTargetFlex *scale_ph2 = new RooGBRTargetFlex("scale_ph2","",*scalefunc,*scalevar_ph2,vars_ph2);
  RooGBRTargetFlex *smear_ph1 = new RooGBRTargetFlex("smear_ph1","",*smearfunc,*smearvar_ph1,vars_ph1);
  RooGBRTargetFlex *smear_ph2 = new RooGBRTargetFlex("smear_ph2","",*smearfunc,*smearvar_ph2,vars_ph2);  
  
  RooRealConstraint *scalelim_ph1 = new RooRealConstraint("scalelim_ph1","",*scale_ph1,0.7,1.3);
  RooRealConstraint *scalelim_ph2 = new RooRealConstraint("scalelim_ph2","",*scale_ph2,0.7,1.3);  
  RooRealConstraint *smearlim_ph1 = new RooRealConstraint("smearlim_ph1","",*smear_ph1,0.,pow(0.2,2));
  RooRealConstraint *smearlim_ph2 = new RooRealConstraint("smearlim_ph2","",*smear_ph2,0.,pow(0.2,2));    
  
  //RooFormulaVar *mscale = new RooFormulaVar("mscale","","sqrt(@0*@1)",RooArgList(*scalelim_ph1,*scalelim_ph2));
  //RooLinearVar *scaledmass = new RooLinearVar("scaledmass","",*mass,*mscale,RooConst(0.));

  //RooLinearVar *scaledmasssingle = new RooLinearVar("scaledmasssingle","",*mass,*scalelim_ph1,RooConst(0.));  
  
  int ngaus = 4;
  
  RooArgList tgtsmass;
  
  RooArgList gauspdfs;
  RooArgList gauspdfsdata;
  RooArgList gauspdfsdatasingle;
  RooArgList gauscoeffs;
  
  
  double step = 0.2/double(std::max(1,ngaus-1));
  
  for (int igaus=0; igaus<ngaus; ++igaus) {
    RooRealVar *gmeanvar = new RooRealVar(TString::Format("gmeanvar_%i",igaus),"",log(80.)-logpdgzmass + step*igaus);
    RooRealVar *gsigmavar = new RooRealVar(TString::Format("gsigmavar_%i",igaus),"",0.02);
    RooRealVar *gsigmaRvar = new RooRealVar(TString::Format("gsigmaRvar_%i",igaus),"",0.02);
    RooRealVar *gfracvar = new RooRealVar(TString::Format("gfracvar_%i",igaus),"",1.0);
    
    gmeanvar->setConstant(false);
    gsigmavar->setConstant(false);
    gsigmaRvar->setConstant(false);
    gfracvar->setConstant(false);
    
    RooGBRFunctionFlex *gmeanfunc = new RooGBRFunctionFlex(TString::Format("gmeanfunc_%i",igaus),"");
    RooGBRFunctionFlex *gsigmafunc = new RooGBRFunctionFlex(TString::Format("gsigmafunc_%i",igaus),"");
    RooGBRFunctionFlex *gsigmaRfunc = new RooGBRFunctionFlex(TString::Format("gsigmaRfunc_%i",igaus),"");
    RooGBRFunctionFlex *gfracfunc = new RooGBRFunctionFlex(TString::Format("gfracfunc_%i",igaus),"");
    
    RooGBRTargetFlex *gmean = new RooGBRTargetFlex(TString::Format("gmean_%i",igaus),"",*gmeanfunc,*gmeanvar,condvarsmass);
    RooGBRTargetFlex *gsigma = new RooGBRTargetFlex(TString::Format("gsigma_%i",igaus),"",*gsigmafunc,*gsigmavar,condvarsmass);
    RooGBRTargetFlex *gsigmaR = new RooGBRTargetFlex(TString::Format("gsigmaR_%i",igaus),"",*gsigmaRfunc,*gsigmaRvar,condvarsmass);
    RooGBRTargetFlex *gfrac = new RooGBRTargetFlex(TString::Format("gfrac_%i",igaus),"",*gfracfunc,*gfracvar,condvarsmass);

    RooRealConstraint *gmeanlim = new RooRealConstraint(TString::Format("gmeanlim_%i",igaus),"",*gmean,log(75.)-logpdgzmass,log(105.)-logpdgzmass);
    RooRealConstraint *gsigmalim = new RooRealConstraint(TString::Format("gsigmalim_%i",igaus),"",*gsigma,1e-5,0.3);
    RooRealConstraint *gsigmaRlim = new RooRealConstraint(TString::Format("gsigmaRlim_%i",igaus),"",*gsigmaR,1e-5,0.3);
    //RooRealConstraint *gfraclim = new RooRealConstraint(TString::Format("gfraclim_%i",igaus),"",*gfrac,0.,1.);
    
    RooAbsReal *gfraclim = new RooProduct(TString::Format("gfraclim_%i",igaus),"",RooArgList(*gfrac,*gfrac));
     
    RooFormulaVar *gmeanscale = new RooFormulaVar(TString::Format("gmeanscale_%i",igaus),"","@0 - 0.5*log(@1*@2)",RooArgList(*gmeanlim,*scalelim_ph1,*scalelim_ph2));
    RooFormulaVar *gsigmasmear = new RooFormulaVar(TString::Format("gsigmasmear_%i",igaus),"","sqrt(@0*@0 + 0.25*(@1+@2))",RooArgList(*gsigmalim,*smearlim_ph1,*smearlim_ph2));

    RooFormulaVar *gmeanscalesingle = new RooFormulaVar(TString::Format("gmeanscalesingle_%i",igaus),"","@0 - log(@1)",RooArgList(*gmeanlim,*scalelim_ph1));    
    RooFormulaVar *gsigmasmearsingle = new RooFormulaVar(TString::Format("gsigmasmearsingle_%i",igaus),"","sqrt(@0*@0 + 0.5*@1)",RooArgList(*gsigmalim,*smearlim_ph1));    
    
    if (igaus==0) {
      gfraclim = new RooConstVar(TString::Format("gfraclimconst_%i",igaus),"",1.);
    }
    else {
      tgtsmass.add(*gfrac);
    }
    
    RooGaussianFast *gpdf = new RooGaussianFast(TString::Format("gdf_%i",igaus),"",*logmass,*gmeanlim,*gsigmalim);
    RooGaussianFast *gpdfdata = new RooGaussianFast(TString::Format("gdfdata_%i",igaus),"",*logmass,*gmeanscale,*gsigmasmear);
    RooGaussianFast *gpdfdatasingle = new RooGaussianFast(TString::Format("gpdfdatasingle_%i",igaus),"",*logmass,*gmeanscalesingle,*gsigmasmearsingle);
    
    //RooBifurGauss *gpdf = new RooBifurGauss(TString::Format("gdf_%i",igaus),"",idmva,*gmeanlim,*gsigmalim,*gsigmaRlim);
    
    gauspdfs.add(*gpdf);
    gauscoeffs.add(*gfraclim);
    
    gauspdfsdata.add(*gpdfdata);
    gauspdfsdatasingle.add(*gpdfdatasingle);
    
    tgtsmass.add(*gmean);
    tgtsmass.add(*gsigma);
    //tgtsid.add(*gsigmaR);
    //tgtsid.add(*gfrac);
    
    
  }  
  RooCondAddPdf masspdf("masspdf","",gauspdfs,gauscoeffs);  
  RooCondAddPdf masspdfdata("masspdfdata","",gauspdfsdata,gauscoeffs);  
  RooCondAddPdf masspdfdatasingle("masspdfdatasingle","",gauspdfsdatasingle,gauscoeffs);  
  
  
  
  

  
   
  RooArgList tgts;
  tgts.add(tgtsmass);
  tgts.add(*scale_ph1);
  tgts.add(*scale_ph2);
  tgts.add(*smear_ph1);
  tgts.add(*smear_ph2);

  RooArgList tgtssingle;
  tgtssingle.add(tgtsmass);
  tgtssingle.add(*scale_ph1);
  tgtssingle.add(*smear_ph1);
  
  RooConstVar etermconst("etermconst","",0.);  
   
  RooRealVar r("r","",1.);
  r.setConstant();



  double minweight = 200;

  std::vector<double> minweights;
  minweights.push_back(200.);
  minweights.push_back(0.);
  
  //ntot.setConstant();

  TFile *fres = new TFile("fres.root","RECREATE");

  
//   masspdf.fitTo(*hdata,ConditionalObservables(condvarsmass),NumCPU(16),Minimizer("Minuit2","minimize"));
//   
//   RooArgSet *gausvars = masspdf.getParameters(hdata);
//   gausvars->setAttribAll("Constant",true);  
//   masspdfdatasingle.fitTo(*hdataD,ConditionalObservables(condvarssingle),NumCPU(16),Minimizer("Minuit2","minimize"));
//   gausvars->setAttribAll("Constant",false);  
  
  
  if (1) {
    std::vector<RooAbsData*> vdatainit;
    vdatainit.push_back(hdata);
    
    std::vector<RooAbsReal*> vpdfinit;
    vpdfinit.push_back(&masspdf);
    
    RooHybridBDTAutoPdf bdtpdfinit("bdtpdfinit","",tgtsmass,etermconst,r,vdatainit,vpdfinit);
    bdtpdfinit.SetPrescaleInit(1);    
    bdtpdfinit.TrainForest(0);
    
  }
  

  
  if(0) {
    RooArgSet *gausvars = masspdf.getParameters(hdata);
    printf("gausvars:\n");
    gausvars->Print("V");
    //return;
    
    gausvars->setAttribAll("Constant",true);

    std::vector<RooAbsData*> vdatainitd;
    vdatainitd.push_back(hdataD);        
    
    std::vector<RooAbsReal*> vpdfinitd;
    vpdfinitd.push_back(&masspdfdatasingle);     
    
    RooHybridBDTAutoPdf bdtpdfinitd("bdtpdfinitd","",tgtssingle,etermconst,r,vdatainitd,vpdfinitd);
    bdtpdfinitd.SetPrescaleInit(1);    
    bdtpdfinitd.TrainForest(0);
    
    scalevar_ph2->setVal(scalevar_ph1->getVal());
    scalevar_ph2->setError(scalevar_ph1->getError());

    smearvar_ph2->setVal(smearvar_ph1->getVal());
    smearvar_ph2->setError(smearvar_ph1->getError());     
    
    gausvars->setAttribAll("Constant",false);
  }
  
  if(1) {    

    std::vector<RooAbsData*> vdatasingle;
    vdatasingle.push_back(hdata);    
    vdatasingle.push_back(hdataD);      
    
    std::vector<RooAbsReal*> vpdfsingle;
    vpdfsingle.push_back(&masspdf);  
    vpdfsingle.push_back(&masspdfdatasingle);     
    
    RooHybridBDTAutoPdf bdtpdfsingle("bdtpdfsingle","",tgtssingle,etermconst,r,vdatasingle,vpdfsingle);
    bdtpdfsingle.SetPrescaleInit(1);    
    bdtpdfsingle.TrainForest(0);
    
    scalevar_ph2->setVal(scalevar_ph1->getVal());
    scalevar_ph2->setError(scalevar_ph1->getError());

    smearvar_ph2->setVal(smearvar_ph1->getVal());
    smearvar_ph2->setError(smearvar_ph1->getError());     
  }  
  
  if (1) {  

    std::vector<RooAbsData*> vdata;
    vdata.push_back(hdata);    
    vdata.push_back(hdataD);        
    
    std::vector<RooAbsReal*> vpdf;
    vpdf.push_back(&masspdf);  
    vpdf.push_back(&masspdfdata);       
    
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff","",tgts,etermconst,r,vdata,vpdf);
    bdtpdfdiff.SetMinCutSignificance(5.);
    //bdtpdfdiff.SetPrescaleInit(10);
    //bdtpdfdiff.SetPrescaleInit(100);
   // bdtpdfdiff.SetPrescaleInit(10);
    //bdtpdfdiff.SetMaxNSpurious(300.);
    //bdtpdfdiff.SetMaxNSpurious(2400.);
    //bdtpdfdiff.SetShrinkage(0.1);
    //bdtpdfdiff.SetShrinkage(0.1);
    bdtpdfdiff.SetDoInitialFit(false);
    bdtpdfdiff.SetShrinkage(0.1);
    bdtpdfdiff.SetMinWeights(minweights);
    //bdtpdfdiff.SetMinWeightTotal(100.);
    //bdtpdfdiff.SetMaxNodes(270);
    bdtpdfdiff.SetMaxNodes(750);
    //bdtpdfdiff.SetMaxNodes(600);
    //bdtpdfdiff.SetMaxNodes(500);
    //bdtpdfdiff.SetMaxDepth(8);
    //bdtpdfdiff.TrainForest(1e6);
    bdtpdfdiff.TrainForest(500);  
    
  }   
     
  RooWorkspace *wereg = new RooWorkspace("weregmass");
  wereg->import(masspdf,RecycleConflictNodes());
  wereg->import(masspdfdata,RecycleConflictNodes());
  wereg->defineSet("vars",vars,true);
  
  wereg->writeToFile("weregmass.root");    

     
  new TCanvas;
  RooPlot *plot = logmass->frame(100);
  hdata->plotOn(plot);
  masspdf.plotOn(plot,ProjWData(*hdata));
  plot->Draw();

  new TCanvas;
  RooPlot *plotD = logmass->frame(100);
  hdataD->plotOn(plotD);
  masspdfdata.plotOn(plotD,ProjWData(*hdataD));
  plotD->Draw();
  
  
//   RooFormulaVar *scaledmass = new RooFormulaVar("scaledmass","","sqrt(@0*@1)*@2",RooArgList(*scalelim_ph1,*scalelim_ph2,*mass));
//   RooRealVar *scaledmassvar = (RooRealVar*)hdataD->addColumn(*scaledmass);
  
//   new TCanvas;
//   RooPlot *plotm = mass->frame
  
  
  
  RooAbsReal *condnll = masspdf.createNLL(*hdata,ConditionalObservables(condvarsmass),NumCPU(16));
  RooAbsReal *condnllD = masspdfdata.createNLL(*hdataD,ConditionalObservables(condvars),NumCPU(16));
  
  printf("condnll = %5f, condnllD = %5f, totalnll = %5f\n",condnll->getVal(),condnllD->getVal(), condnll->getVal()+condnllD->getVal());
  
  
  
  //return;
  
   
  
  
  return;
  
  
}