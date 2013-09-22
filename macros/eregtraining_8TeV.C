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

void eregtraining_8TeV(bool dobarrel, bool doele) {
   
//   gSystem->Setenv("OMP_WAIT_POLICY","PASSIVE");
  
  //candidate to set fixed alpha values (0.9,3.8)
  //TString dirname = TString::Format("/afs/cern.ch/work/b/bendavid/bare/eregtesteleJul30_sig5_01_alphafloat5_%i/",int(minevents)); 
  
  //TString dirname = "/afs/cern.ch/work/b/bendavid/bare/eregAug23test/";
  TString dirname = "/data/bendavid/ereg8TeV_Sept16_full/"; 
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);  
  
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
    
  std::vector<std::string> *varslist;
  if (dobarrel) varslist = varseb;
  else varslist = varsee;
  
  RooArgList vars;
  for (unsigned int ivar=0; ivar<varslist->size(); ++ivar) {
    RooRealVar *var = new RooRealVar(TString::Format("var_%i",ivar),varslist->at(ivar).c_str(),0.);
    vars.addOwned(*var);
  }
  
  RooArgList condvars(vars);
  
//   RooRealVar *tgtvar = new RooRealVar("tgtvar","ph.scrawe/ph.gene",1.);
//   if (!dobarrel) tgtvar->SetTitle("(ph.scrawe + ph.scpse)/ph.gene");
  
  RooRealVar *tgtvar = new RooRealVar("tgtvar","ph.gene/ph.scrawe",1.);
  if (!dobarrel) tgtvar->SetTitle("ph.gene/(ph.scrawe + ph.scpse)");  
  
/*  RooRealVar *tgtvar = new RooRealVar("tgtvar","log(ph.gene/ph.scrawe)",1.);
  if (!dobarrel) tgtvar->SetTitle("log(ph.gene/(ph.scrawe + ph.scpse))");  */  
  
  //tgtvar->setRange(0.8,2.);
  
  vars.addOwned(*tgtvar);


  
  //varstest.add(*tgtvar);
    
  RooRealVar weightvar("weightvar","",1.);

  //TFile *fdin = TFile::Open("/home/mingyang/cms/hist/hgg-2013Moriond/merged/hgg-2013Moriond_s12-diphoj-3-v7a_noskim.root");
//   TFile *fdin = TFile::Open("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/trainingtreesJul1/hgg-2013Final8TeV_s12-zllm50-v7n_noskim.root");
//   TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("PhotonTreeWriterSingleInvert");
//   TTree *dtree = (TTree*)ddir->Get("hPhotonTreeSingle");    
  
/*  TFile *fdinsig = TFile::Open("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/trainingtreesJul1/hgg-2013Moriond_s12-h125gg-gf-v7a_noskim.root");
  TDirectory *ddirsig = (TDirectory*)fdinsig->FindObjectAny("PhotonTreeWriterPreselNoSmear");
  TTree *dtreesig = (TTree*)ddirsig->Get("hPhotonTreeSingle");     */ 
  
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
  //TCut oddevents = prescale100; 
  
  //weightvar.SetTitle(prescale10*selcut);
  
/*  new TCanvas;
  tree->Draw("ph.genpt>>hpt(200,0.,100.)",selweight*selcut);

  return;*/  
  
  if (doele) {
    weightvar.SetTitle(evenevents*selcut);
  }
  else {
    weightvar.SetTitle(selweight*selcut);
  }
  RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",tree,vars,weightvar);   
  
//   weightvar.SetTitle(prescale1000*selcut);
//   RooDataSet *hdatasig = RooTreeConvert::CreateDataSet("hdatasig",dtree,vars,weightvar);   
//   RooDataSet *hdatasigtest = RooTreeConvert::CreateDataSet("hdatasigtest",dtree,varstest,weightvar); 
  
  RooDataSet *hdatasig = 0;
  RooDataSet *hdatasigtest = 0;
  
//   weightvar.SetTitle(prescale10*selcut);
//   RooDataSet *hdatasigsmall = RooTreeConvert::CreateDataSet("hdatasigsmall",dtreesig,vars,weightvar);   
  
  RooRealVar sigwidthtvar("sigwidthtvar","",0.01);
  sigwidthtvar.setConstant(false);
  
  RooRealVar sigmeantvar("sigmeantvar","",1.);
  sigmeantvar.setConstant(false); 

  RooRealVar sigalphavar("sigalphavar","",1.);
  sigalphavar.setConstant(false);   
  
  RooRealVar signvar("signvar","",3.);
  signvar.setConstant(false);     

  RooRealVar sigalpha2var("sigalpha2var","",1.);
  sigalpha2var.setConstant(false);   
  
  RooRealVar sign2var("sign2var","",3.);
  sign2var.setConstant(false);     
  
  RooRealVar sigwidth2var("sigwidth2var","",0.5);
  sigwidth2var.setConstant(false);
  
  RooRealVar sigmean2var("sigmean2var","",1.);
  sigmean2var.setConstant(false);   
  
  RooRealVar sigfracvar("sigfracvar","",0.9);
  sigfracvar.setConstant(false);       
  
  
   
  RooArgList tgts;
  RooGBRFunction func("func","",condvars,4);
  RooGBRTarget sigwidtht("sigwidtht","",func,0,sigwidthtvar);
  RooGBRTarget sigmeant("sigmeant","",func,1,sigmeantvar);
  RooGBRTarget signt("signt","",func,2,signvar);
  RooGBRTarget sign2t("sign2t","",func,3,sign2var);
  
    
  tgts.add(sigwidtht);
  tgts.add(sigmeant);
   //tgts.add(sigalpha);
  tgts.add(signt);
   //tgts.add(sigalpha2);
  tgts.add(sign2t);
//   tgts.add(sigmean2);
//   tgts.add(sigwidth2);
//   tgts.add(sigfrac);

  double lowboundalpha = 1e-2;
//   if (dobarrel) lowboundalpha = 0.05;
//   else lowboundalpha = 0.05;
  
  RooRealConstraint sigwidthlim("sigwidthlim","",sigwidtht,0.0002,0.5);
  RooRealConstraint sigmeanlim("sigmeanlim","",sigmeant,0.2,2.0);
  //RooRealConstraint sigmeanlim("sigmeanlim","",sigmeant,-2.0,2.0);
  //RooRealConstraint sigmeanlim("sigmeanlim","",sigmeant,-2.0,-0.2); 
  
  RooRealConstraint signlim("signlim","",signt,1.01,5000.); 
  //RooRealConstraint sigalphalim("sigalphalim","",sigalpha,lowboundalpha,8.0);

  RooRealConstraint sign2lim("sign2lim","",sign2t,1.01,5000.); 
  //RooRealConstraint sigalpha2lim("sigalpha2lim","",sigalpha2,lowboundalpha,8.0);  
  
  
/*  RooRealConstraint signlim("signlim","",signvar,1.01,500.); 
  RooRealConstraint sigalphalim("sigalphalim","",sigalphavar,lowboundalpha,8.0);

  RooRealConstraint sign2lim("sign2lim","",sign2var,1.01,500.); 
  RooRealConstraint sigalpha2lim("sigalpha2lim","",sigalpha2var,lowboundalpha,8.0);  */  
  
  
//   RooRealConstraint sigwidth2lim("sigwidth2lim","",sigwidth2,0.0002,10.);
//   RooRealConstraint sigmean2lim("sigmean2lim","",sigmean2,0.,10.0);  
//   
//   RooRealConstraint sigfraclim("sigfraclim","",sigfrac,0.,1.0);  
  
  
  //RooLinearVar tgtscaled("tgtscaled","",*tgtvar,sigmeanlim,RooConst(0.));
  
  //printf("tgtvar: min = %5e, max = %5e\n",tgtvar->getMin(),tgtvar->getMax());
  //printf("tgtscaled: min = %5e, max = %5e\n",tgtscaled.getMin(),tgtscaled.getMax());
  
 // RooDoubleCBFast sigpdf("sigpdf","",tgtscaled,RooConst(1.),sigwidthlim,sigalphalim,signlim,sigalpha2lim,sign2lim);
  
  //RooDoubleCBFast sigpdf("sigpdf","",tgtscaled,RooConst(1.),sigwidthlim,RooConst(2.),signlim,sigalpha2lim,sign2lim);
 
  
  RooDoubleCBFast sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,RooConst(2.),signlim,RooConst(1.),sign2lim);
  //RooDoubleCBFast sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,sigalphalim,signlim,sigalpha2lim,sign2lim);
  
  

  //RooDoubleCBFast sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,sigalphalim,signlim,sigalpha2lim,sign2lim);
  //RooCBFast sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,sigalpha2lim,sign2lim);
  
 // RooDoubleCBFast sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,RooConst(5.),RooConst(2.),sigalpha2lim,sign2lim);
  
 // RooDoubleCBFast sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,RooConst(5.),RooConst(2.),RooConst(5.),RooConst(2.));
  
   //RooDoubleCBFast sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,sigalphalim,RooConst(2.),sigalpha2lim,RooConst(2.));

//   RooDoubleCBFast sigg1("sigg1","",*tgtvar,sigmeanlim,sigwidthlim,RooConst(5.),RooConst(2.),RooConst(5.),RooConst(2.));
//   //RooGaussian sigg1("sigg1","",*tgtvar,sigmeanlim,sigwidthlim);
//   RooGaussianFast sigg2("sigg2","",*tgtvar,sigmean2lim,sigwidth2lim);
//   RooCondAddPdf sigpdf("sigpdf","",RooArgList(sigg1,sigg2),sigfraclim);
  
  
//   RooProduct sigwidthscaled("sigwidthscaled","",RooArgList(sigmeanlim,sigwidthlim));
//   RooDoubleCBFast sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthscaled,sigalphalim,signlim,sigalpha2lim,sign2lim);
  
  
  
  //RooDoubleCBFast sigpdf("sigpdf","",tgtscaled,RooConst(1.),sigwidthlim,RooConst(2.),signlim,RooConst(1.),sign2lim);
  
  //RooDoubleCBFast sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,sigalphalim,signlim,sigalpha2lim,sign2lim);
  
  //RooDoubleCBFast sigpdf("sigpdf","",tgtscaled,RooConst(1.),sigwidthlim,RooConst(2.0),signlim,RooConst(1.0),sign2lim);
  
  //RooCBExp sigpdf("sigpdf","",tgtscaled,RooConst(-1.),sigwidthlim,sigalpha2lim,sign2lim,sigalphalim);
  
  //RooDoubleCBFast sigpdf("sigpdf","",tgtscaled,RooConst(1.),sigwidthlim,RooConst(100.),RooConst(100.),sigalpha2lim,sign2lim);
  //RooDoubleCBFast sigpdf("sigpdf","",tgtscaled,RooConst(1.),sigwidthlim,sigalphalim,signlim,RooConst(3.),sign2lim);
  //RooCBShape sigpdf("sigpdf","",tgtscaled,RooConst(1.),sigwidthlim,sigalphalim,signlim);
  
  RooConstVar etermconst("etermconst","",0.);  
  //RooFormulaVar etermconst("etermconst","","1000.*(@0-1.)*(@0-1.)",RooArgList(tgtscaled));
   
  RooRealVar r("r","",1.);
  r.setConstant();

  std::vector<RooAbsReal*> vpdf;
  vpdf.push_back(&sigpdf);  

  double minweight = 200;

  std::vector<double> minweights;
  minweights.push_back(minweight);
  
  //ntot.setConstant();

  TFile *fres = new TFile("fres.root","RECREATE");

  if (1) {  
    std::vector<RooAbsData*> vdata;
    vdata.push_back(hdata);    
    
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff","",func,tgts,etermconst,r,vdata,vpdf);
    bdtpdfdiff.SetMinCutSignificance(5.);
    //bdtpdfdiff.SetPrescaleInit(5);
    bdtpdfdiff.SetPrescaleInit(100);
   // bdtpdfdiff.SetPrescaleInit(10);
    //bdtpdfdiff.SetMaxNSpurious(300.);
    //bdtpdfdiff.SetMaxNSpurious(2400.);
    //bdtpdfdiff.SetShrinkage(0.1);
    //bdtpdfdiff.SetShrinkage(0.1);
    bdtpdfdiff.SetShrinkage(0.1);
    bdtpdfdiff.SetMinWeights(minweights);
    //bdtpdfdiff.SetMaxNodes(270);
    bdtpdfdiff.SetMaxNodes(750);
    //bdtpdfdiff.SetMaxNodes(600);
    //bdtpdfdiff.SetMaxNodes(500);
    //bdtpdfdiff.SetMaxDepth(8);
    //bdtpdfdiff.TrainForest(1e6);
    bdtpdfdiff.TrainForest(1e6);  
    
  }   
     
  
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