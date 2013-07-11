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
 
void eregtraining(double minevents=200) {
  
//   gSystem->Setenv("OMP_WAIT_POLICY","PASSIVE");
  
  TString dirname = TString::Format("/afs/cern.ch/work/b/bendavid/bare/eregtesteleJul7_%i/",int(minevents)); 
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);  
  
  std::vector<std::string> *varsf = new std::vector<std::string>;
  varsf->push_back("ph.scrawe");
  varsf->push_back("ph.sceta");
  varsf->push_back("ph.scphi");
  varsf->push_back("ph.r9");  
  varsf->push_back("ph.e5x5/ph.scrawe");  
  varsf->push_back("ph.scetawidth");
  varsf->push_back("ph.scphiwidth");  
  varsf->push_back("ph.scnclusters");
  varsf->push_back("ph.hoveretower");
  varsf->push_back("rho");
  varsf->push_back("nVtx");  
 
  varsf->push_back("ph.etaseed-ph.sceta");
  varsf->push_back("atan2(sin(ph.phiseed-ph.scphi),cos(ph.phiseed-ph.scphi))");
  varsf->push_back("ph.eseed/ph.scrawe");
  varsf->push_back("ph.e3x3seed/ph.eseed");
  varsf->push_back("ph.e5x5seed/ph.eseed");
  varsf->push_back("ph.sigietaietaseed");   
  varsf->push_back("ph.sigiphiphiseed");   
  varsf->push_back("ph.covietaiphiseed");
  varsf->push_back("ph.emaxseed/ph.eseed");
  varsf->push_back("ph.e2ndseed/ph.eseed");
  varsf->push_back("ph.etopseed/ph.eseed");
  varsf->push_back("ph.ebottomseed/ph.eseed");
  varsf->push_back("ph.eleftseed/ph.eseed");
  varsf->push_back("ph.erightseed/ph.eseed");
  varsf->push_back("ph.e2x5maxseed/ph.eseed");
  varsf->push_back("ph.e2x5topseed/ph.eseed");
  varsf->push_back("ph.e2x5bottomseed/ph.eseed");
  varsf->push_back("ph.e2x5leftseed/ph.eseed");
  varsf->push_back("ph.e2x5rightseed/ph.eseed");
  
  std::vector<std::string> *varseb = new std::vector<std::string>(*varsf);
  std::vector<std::string> *varsee = new std::vector<std::string>(*varsf);
  
  varseb->push_back("ph.ietaseed");
  varseb->push_back("ph.iphiseed");
  varseb->push_back("ph.ietaseed%5");
  varseb->push_back("ph.iphiseed%2");       
  varseb->push_back("(abs(ph.ietaseed)<=25)*(ph.ietaseed%25) + (abs(ph.ietaseed)>25)*((ph.ietaseed-25*abs(ph.ietaseed)/ph.ietaseed)%20)");
  varseb->push_back("ph.iphiseed%20"); 
  varseb->push_back("ph.etacryseed");
  varseb->push_back("ph.phicryseed");

  varsee->push_back("ph.scpse/ph.scrawe");
    
  RooArgList vars;
  for (unsigned int ivar=0; ivar<varseb->size(); ++ivar) {
    RooRealVar *var = new RooRealVar(TString::Format("var_%i",ivar),varseb->at(ivar).c_str(),0.);
    vars.addOwned(*var);
  }
  
  RooArgList condvars(vars);
  
  RooRealVar *tgtvar = new RooRealVar("tgtvar","ph.scrawe/ph.gene",1.);
  //RooRealVar *tgtvar = new RooRealVar("tgtvar","ph.gene/ph.scrawe",1.,0.,2.);
  //tgtvar->setBins(800);
  vars.addOwned(*tgtvar);

  RooRealVar testvar("testvar","ph.e/ph.gene",1.,0.,2.);
  testvar.setBins(800);
  RooArgList varstest;
  varstest.add(testvar);
  
  //varstest.add(*tgtvar);
    
  RooRealVar weightvar("weightvar","",1.);

  //TFile *fdin = TFile::Open("/home/mingyang/cms/hist/hgg-2013Moriond/merged/hgg-2013Moriond_s12-diphoj-3-v7a_noskim.root");
  TFile *fdin = TFile::Open("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/trainingtreesJul1/hgg-2013Final8TeV_s12-zllm50-v7n_noskim.root");
  TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("PhotonTreeWriterSingleInvert");
  TTree *dtree = (TTree*)ddir->Get("hPhotonTreeSingle");    
  
/*  TFile *fdinsig = TFile::Open("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/trainingtreesJul1/hgg-2013Moriond_s12-h125gg-gf-v7a_noskim.root");
  TDirectory *ddirsig = (TDirectory*)fdinsig->FindObjectAny("PhotonTreeWriterPreselNoSmear");
  TTree *dtreesig = (TTree*)ddirsig->Get("hPhotonTreeSingle");     */ 
  
  TCut selcut = "ph.pt>25. && ph.isbarrel && ph.ispromptgen"; 
  //TCut selcut = "ph.pt>25. && ph.isbarrel && (ph.scrawe/ph.gene)>0. && (ph.scrawe/ph.gene)<2. && ph.ispromptgen";
  //TCut selcut = "ph.pt>25. && ph.isbarrel && (ph.gene/ph.scrawe)>0. && (ph.gene/ph.scrawe)<2.";
  TCut selweight = "xsecweight(procidx)*puweight(numPU,procidx)";
  TCut prescale10 = "(evt%10==0)";
  TCut prescale25 = "(evt%25==0)";
  TCut prescale100 = "(evt%100==0)";  
  TCut prescale1000 = "(evt%1000==0)";  
  TCut evenevents = "(evt%2==0)";
  TCut oddevents = "(evt%2==1)";
  //TCut oddevents = prescale100; 
  
  weightvar.SetTitle(prescale100*selcut);
  RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",dtree,vars,weightvar);   
  
//   weightvar.SetTitle(prescale1000*selcut);
//   RooDataSet *hdatasig = RooTreeConvert::CreateDataSet("hdatasig",dtree,vars,weightvar);   
//   RooDataSet *hdatasigtest = RooTreeConvert::CreateDataSet("hdatasigtest",dtree,varstest,weightvar); 
  
  RooDataSet *hdatasig = 0;
  RooDataSet *hdatasigtest = 0;
  
//   weightvar.SetTitle(prescale10*selcut);
//   RooDataSet *hdatasigsmall = RooTreeConvert::CreateDataSet("hdatasigsmall",dtreesig,vars,weightvar);   
  
  RooRealVar sigwidthtvar("sigwidthtvar","",0.008);
  sigwidthtvar.setConstant(false);
  
  RooRealVar sigmeantvar("sigmeantvar","",1.);
  sigmeantvar.setConstant(false); 

  RooRealVar sigalphavar("sigalphavar","",1.0);
  sigalphavar.setConstant(false);   
  
  RooRealVar signvar("signvar","",1.);
  signvar.setConstant(false);     

  RooRealVar sigalpha2var("sigalpha2var","",1.0);
  sigalpha2var.setConstant(false);   
  
  RooRealVar sign2var("sign2var","",1.);
  sign2var.setConstant(false);     
  
  
  RooRealVar sigmean2var("sigmean2var","",1.);
  sigmean2var.setConstant(false);       

  RooRealVar sigwidth2var("sigwidth2var","",0.04);
  sigwidth2var.setConstant(false);
  
  RooRealVar sigmean3var("sigmean3var","",1.);
  sigmean3var.setConstant(false);       

  RooRealVar sigwidth3var("sigwidth3var","",8.0);
  sigwidth3var.setConstant(false);
   
  RooRealVar sigfracvar("sigfracvar","",0.8);
  sigfracvar.setConstant(false);  

  RooRealVar sigfrac2var("sigfrac2var","",0.8);
  sigfrac2var.setConstant(false);    
  
  RooArgList tgts;
  RooGBRFunction func("func","",condvars,9);
  RooGBRTarget sigwidtht("sigwidtht","",func,0,sigwidthtvar);
  RooGBRTarget sigmeant("sigmeant","",func,1,sigmeantvar);
  RooGBRTarget sigalpha("sigalpha","",func,2,sigalphavar);
  RooGBRTarget signt("signt","",func,3,signvar);
  RooGBRTarget sigmean2("sigmean2","",func,4,sigmean2var);
  RooGBRTarget sigwidth2("sigwidth2","",func,5,sigwidth2var);  
  RooGBRTarget sigfrac("sigfrac","",func,6,sigfracvar);
  RooGBRTarget sigalpha2("sigalpha2","",func,7,sigalpha2var);
  RooGBRTarget sign2t("sign2t","",func,8,sign2var);
   
  
 tgts.add(sigwidtht);
  tgts.add(sigmeant);
  tgts.add(sigalpha);
  tgts.add(signt);
  tgts.add(sigalpha2);
  tgts.add(sign2t);
//    tgts.add(sigmean2);
//   tgts.add(sigwidth2);
//   tgts.add(sigfrac);


   
  
  RooRealConstraint sigwidthlim("sigwidthlim","",sigwidtht,0.0002,0.5);
  RooRealConstraint sigmeanlim("sigmeanlim","",sigmeant,0.2,2.0); 
  
  RooRealConstraint signlim("signlim","",signt,0.,20.); 
  RooRealConstraint sigalphalim("sigalphalim","",sigalpha,0.,5.);

  RooRealConstraint sign2lim("sign2lim","",sign2t,0.,20.); 
  RooRealConstraint sigalpha2lim("sigalph2alim","",sigalpha2,0.,5.);  
  
  //RooRealConstraint sigalphalim("sigalphalim","",sigalpha,-20.,0.); 
  
//   RooAbsReal &sigmeanlim = sigmeant;
//   RooAbsReal &sigwidthlim = sigwidtht;
  
  RooRealConstraint sigwidth2lim("sigwidth2lim","",sigwidth2,0.002,0.5);
  RooRealConstraint sigmean2lim("sigmean2lim","",sigmean2,0.05,2.);      
  
  RooRealConstraint sigfraclim("sigfraclim","",sigfrac,0.5,1.0);     
  
  RooLinearVar tgtscaled("tgtscaled","",*tgtvar,sigmeanlim,RooConst(0.));
  //RooCBShape sigpdf("sigpdf","",tgtscaled,RooConst(1.),sigwidthlim,sigalphalim,signlim);
 // RooCBShape sigpdf("sigpdf","",tgtscaled,RooConst(1.),sigwidthlim,sigalphalim,signlim);
  //RooCBShape sigpdf("sigpdf","",tgtscaled,RooConst(1.),RooConst(0.008),sigalphalim,signlim);
  //RooCBShape sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,sigalphalim,signlim);
 // RooGaussian siggaus("siggaus","",tgtscaled,sigmean2lim,sigwidth2lim);
  
  RooDoubleCBFast sigpdf("sigpdf","",tgtscaled,RooConst(1.),sigwidthlim,sigalphalim,signlim,sigalpha2lim,sign2lim);
  //RooDoubleCBSlow sigpdf("sigpdf","",tgtscaled,RooConst(1.),sigwidthlim,sigalphalim,signlim,sigalpha2lim,sign2lim);
  
  
  //RooDoubleCB sigpdf("sigpdf","",tgtscaled,RooConst(1.),sigwidthlim,RooConst(2.),signlim,RooConst(2.),sign2lim);
   //RooDoubleCB sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,sigalphalim,signlim,sigalpha2lim,sign2lim);
  
 // RooGaussian sigpdf("sigpdf","",tgtscaled,RooConst(1.),RooConst(0.005));
  
 // RooFormulaVar sigreal("sigreal","","@1*exp(-100.*(@0-1.)*(@0-1.))",RooArgList(tgtscaled,sigpdf));
  
 //  sigpdf.fitTo(*hdata,ConditionalObservables(condvars),NumCPU(16));
//   return;
  

   //RooCondAddPdf sigpdf("sigpdf","",RooArgList(sigcb,siggaus),RooArgList(sigfraclim));
  //RooCBShape &sigpdf = sigcb;
  //RooDoubleCB sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,sigalphalim,signlim,sigalpha2lim,sign2lim);
  //RooGaussian sigpdf("sigpdf","",tgtscaled,RooConst(1.),sigwidthlim);
  
  //RooFormulaVar sreal("sreal","","@0*exp(-@1*@1)",RooArgList(sigpdf,sigwidthlim));
  
  RooConstVar etermconst("etermconst","",0.);  
  //RooFormulaVar etermconst("etermconst","","1000.*(@0-1.)*(@0-1.)",RooArgList(tgtscaled));
  
  RooRealVar r("r","",1.);
  r.setConstant();

  std::vector<RooAbsReal*> vpdf;
  vpdf.push_back(&sigpdf);  

  std::vector<double> minweights;
  minweights.push_back(minevents);
  
  //ntot.setConstant();

  TFile *fres = new TFile("fres.root","RECREATE");

  if (1) {  
    std::vector<RooAbsData*> vdata;
    vdata.push_back(hdata);    
    
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff","",func,tgts,etermconst,r,vdata,vpdf);
    bdtpdfdiff.SetMinCutSignificance(7.);
    bdtpdfdiff.SetShrinkage(0.3);
    bdtpdfdiff.SetMinWeights(minweights);
    //bdtpdfdiff.SetMaxDepth(8);
    bdtpdfdiff.TrainForest(1e6);  
    
  }   
  
  
  RooWorkspace *wereg = new RooWorkspace("wereg");
  wereg->import(sigpdf);
//   wereg->import(*hdata);
  //wereg->import(*hdatasigsmall);
//   wereg->import(*hdatasigtest);
  wereg->writeToFile("wereg.root");    
  
  
  return;
  
  RooFormulaVar ecor("ecor","","@0*@1",RooArgList(*tgtvar,sigmeanlim));
  RooRealVar *ecorvar = (RooRealVar*)hdatasig->addColumn(ecor);
  ecorvar->setRange(0.,2.);
  ecorvar->setBins(800);

/*  RooFormulaVar eraw("eraw","","@0",RooArgList(*tgtvar));
  RooRealVar *erawvar = (RooRealVar*)hdatasig->addColumn(eraw);
  erawvar->setRange(0.,2.);
  erawvar->setBins(400); */ 
  
  RooDataSet *hdataclone = new RooDataSet(*hdata,"hdataclone");
  RooRealVar *meanvar = (RooRealVar*)hdataclone->addColumn(sigmeanlim);
  RooRealVar *widthvar = (RooRealVar*)hdataclone->addColumn(sigwidthlim);
  
  new TCanvas;
  RooPlot *plot = tgtvar->frame(0.6,1.2,100);
  hdata->plotOn(plot);
  sigpdf.plotOn(plot,ProjWData(*hdata));
  plot->Draw();
  
/*  new TCanvas;
  RooPlot *plotsig = tgtvar->frame(0.6,1.2,100);
  hdatasig->plotOn(plotsig);
  sigpdf.plotOn(plotsig,ProjWData(*hdatasig));
  plotsig->Draw(); */ 
  
  new TCanvas;
  RooPlot *plotmean = meanvar->frame(0.6,1.2,100);
  hdataclone->plotOn(plotmean);
  plotmean->Draw();  
  
  new TCanvas;
  RooPlot *plotwidth = widthvar->frame(0.,0.05,100);
  hdataclone->plotOn(plotwidth);
  plotwidth->Draw();    
  
  TH1 *heold = hdatasigtest->createHistogram("heold",testvar);
  TH1 *heraw = hdatasig->createHistogram("heraw",*tgtvar,Binning(800,0.,2.));
  TH1 *hecor = hdatasig->createHistogram("hecor",*ecorvar);
  
  heold->SetLineColor(kRed);
  hecor->SetLineColor(kBlue);
  heraw->SetLineColor(kMagenta);
  
  hecor->GetXaxis()->SetRangeUser(0.6,1.2);
  heold->GetXaxis()->SetRangeUser(0.6,1.2);
  
  new TCanvas;
  
  hecor->Draw("HIST");
  heold->Draw("HISTSAME");
  heraw->Draw("HISTSAME");
  
/*  new TCanvas;
  RooPlot *ploteold = testvar.frame(0.6,1.2,100);
  hdatasigtest->plotOn(ploteold);
  ploteold->Draw();    
  
  new TCanvas;
  RooPlot *plotecor = ecorvar->frame(0.6,1.2,100);
  hdatasig->plotOn(plotecor);
  plotecor->Draw(); */   
  
  
}