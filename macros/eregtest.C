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
 
//effsigma function from Chris
Double_t effSigma(TH1 * hist)
{

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
//   if(total < 100.) {
//     cout << "effsigma: Too few entries " << total << endl;
//     return 0.;
//   }
  Int_t ierr=0;
  Int_t ismin=999;
  
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }   
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
  
  return widmin;
  
}

void eregtest(bool dobarrel, bool doele) {
  
  TString dirname = "/afs/cern.ch/work/b/bendavid/bare/eregtestoutAug2b/"; 
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);    
  
  
  
  TString fname;
  if (doele && dobarrel) 
    fname = "wereg_ele_eb.root";
  else if (doele && !dobarrel) 
    fname = "wereg_ele_ee.root";
  else if (!doele && dobarrel) 
    fname = "wereg_ph_eb.root";
  else if (!doele && !dobarrel) 
    fname = "wereg_ph_ee.root";
  
  TString infile = TString::Format("/afs/cern.ch/work/b/bendavid/bare/eregRC4Aug6alt2/%s",fname.Data());
  
  TFile *fws = TFile::Open(infile); 
  RooWorkspace *ws = (RooWorkspace*)fws->Get("wereg");
  
  RooGBRFunction *func = static_cast<RooGBRFunction*>(ws->arg("func"));
  RooRealVar *tgtvar = ws->var("tgtvar");
  
  
  RooArgList vars;
  vars.add(func->Vars());
  vars.add(*tgtvar);
   
  
  RooRealVar weightvar("weightvar","",1.);

  TTree *dtree;
  
  if (doele) {
    TFile *fdin = TFile::Open("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-zllm50-v7n_noskim.root");
    TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("PhotonTreeWriterSingleInvert");
    dtree = (TTree*)ddir->Get("hPhotonTreeSingle");       
  }
  else {
    TFile *fdin = TFile::Open("root://eoscms.cern.ch///eos/cms/store/cmst3/user/bendavid/idTreesAug1/hgg-2013Final8TeV_ID_s12-h124gg-gf-v7n_noskim.root");
    TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("PhotonTreeWriterPreselNoSmear");
    dtree = (TTree*)ddir->Get("hPhotonTreeSingle");       
  }
  
  
//   //TFile *fdin = TFile::Open("/home/mingyang/cms/hist/hgg-2013Moriond/merged/hgg-2013Moriond_s12-diphoj-3-v7a_noskim.root");
//   //TFile *fdin = TFile::Open("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/trainingtreesJul1/hgg-2013Final8TeV_s12-zllm50-v7n_noskim.root");
//   TFile *fdin = TFile::Open("root://eoscms.cern.ch///eos/cms/store/cmst3/user/bendavid/idTreesAug1/hgg-2013Final8TeV_ID_s12-h124gg-gf-v7n_noskim.root");
//   //TFile *fdin = TFile::Open("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-zllm50-v7n_noskim.root");
//   //TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("PhotonTreeWriterSingleInvert");
//   TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("PhotonTreeWriterPreselNoSmear");
//   TTree *dtree = (TTree*)ddir->Get("hPhotonTreeSingle");    
  
/*  TFile *fdinsig = TFile::Open("/home/mingyang/cms/hist/hgg-2013Moriond/merged/hgg-2013Moriond_s12-h125gg-gf-v7a_noskim.root");
  TDirectory *ddirsig = (TDirectory*)fdinsig->FindObjectAny("PhotonTreeWriterPreselNoSmear");
  TTree *dtreesig = (TTree*)ddirsig->Get("hPhotonTreeSingle"); */     
  
  TCut selcut;
  if (dobarrel) 
    selcut = "ph.genpt>25. && ph.isbarrel && ph.ispromptgen"; 
  else
    selcut = "ph.genpt>25. && !ph.isbarrel && ph.ispromptgen"; 
  
//  TCut selcut = "ph.pt>25. && ph.isbarrel && ph.ispromptgen && abs(ph.sceta)<1.0"; 
  //TCut selcut = "ph.pt>25. && ph.isbarrel && (ph.scrawe/ph.gene)>0. && (ph.scrawe/ph.gene)<2. && ph.ispromptgen";
  //TCut selcut = "ph.pt>25. && ph.isbarrel && (ph.gene/ph.scrawe)>0. && (ph.gene/ph.scrawe)<2.";
  TCut selweight = "xsecweight(procidx)*puweight(numPU,procidx)";
  TCut prescale10 = "(evt%10==0)";
  TCut prescale10alt = "(evt%10==1)";
  TCut prescale25 = "(evt%25==0)";
  TCut prescale100 = "(evt%100==0)";  
  TCut prescale1000 = "(evt%1000==0)";  
  TCut evenevents = "(evt%2==0)";
  TCut oddevents = "(evt%2==1)";
  TCut prescale100alt = "(evt%100==1)";
  TCut prescale1000alt = "(evt%1000==1)";
  TCut prescale50alt = "(evt%50==1)";
  //TCut oddevents = prescale100;
  
  if (doele) 
    weightvar.SetTitle(prescale100alt*selcut);
  else
    weightvar.SetTitle(selcut);
  
  RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",dtree,vars,weightvar);   

  if (doele) 
    weightvar.SetTitle(prescale1000alt*selcut);
  else
    weightvar.SetTitle(prescale10alt*selcut);
  
  RooDataSet *hdatasmall = RooTreeConvert::CreateDataSet("hdatasmall",dtree,vars,weightvar);     
  
    
  const HybridGBRForest *forest = func->Forest();
  for (unsigned int itgt=0; itgt<forest->Trees().size(); ++itgt) {
    int ntrees = 0;
    for (unsigned int itree = 0; itree<forest->Trees().at(itgt).size(); ++itree) {
      if (forest->Trees()[itgt][itree].Responses().size()>1) ++ntrees;
    }
    printf("itgt = %i, ntrees = %i\n", int(itgt),ntrees);
  }
  
  
  RooAbsPdf *sigpdf = ws->pdf("sigpdf");
  
  RooRealVar *scetavar = ws->var("var_1");
  
  RooAbsReal *sigmeanlim = ws->function("sigmeanlim");
  RooAbsReal *sigwidthlim = ws->function("sigwidthlim");
  RooAbsReal *signlim = ws->function("signlim");
  RooAbsReal *sign2lim = ws->function("sign2lim");
  RooAbsReal *alphalim = ws->function("sigalphalim");
  RooAbsReal *alpha2lim = ws->function("sigalpha2lim");  

  RooFormulaVar ecor("ecor","","1./(@0*@1)",RooArgList(*tgtvar,*sigmeanlim));
  //RooFormulaVar ecor("ecor","","@1/@0",RooArgList(*tgtvar,*sigmeanlim));
  //RooFormulaVar ecor("ecor","","@0/@1",RooArgList(*tgtvar,*sigmeanlim));
  //RooFormulaVar ecor("ecor","","@0",RooArgList(*tgtvar));
  RooRealVar *ecorvar = (RooRealVar*)hdata->addColumn(ecor);
  ecorvar->setRange(0.,2.);
  ecorvar->setBins(800);
  
   RooFormulaVar raw("raw","","1./@0",RooArgList(*tgtvar));
   RooRealVar *rawvar = (RooRealVar*)hdata->addColumn(raw);
   rawvar->setRange(0.,2.);
   rawvar->setBins(800);

/*  RooFormulaVar eraw("eraw","","@0",RooArgList(*tgtvar));
  RooRealVar *erawvar = (RooRealVar*)hdatasig->addColumn(eraw);
  erawvar->setRange(0.,2.);
  erawvar->setBins(400); */ 
  
  RooDataSet *hdataclone = new RooDataSet(*hdata,"hdataclone");
  RooRealVar *meanvar = (RooRealVar*)hdataclone->addColumn(*sigmeanlim);
  RooRealVar *widthvar = (RooRealVar*)hdataclone->addColumn(*sigwidthlim);
  RooRealVar *nvar = (RooRealVar*)hdataclone->addColumn(*signlim);
  RooRealVar *n2var = (RooRealVar*)hdataclone->addColumn(*sign2lim);
  
  RooRealVar *alphavar = (RooRealVar*)hdataclone->addColumn(*alphalim);
  RooRealVar *alpha2var = (RooRealVar*)hdataclone->addColumn(*alpha2lim);
  
 // hdataclone = (RooDataSet*)hdataclone->reduce("sigwidthlim>0.017");
  
  
  TCanvas *craw = new TCanvas;
  //RooPlot *plot = tgtvar->frame(0.6,1.2,100);
  RooPlot *plot = tgtvar->frame(0.6,2.0,100);
  hdata->plotOn(plot);
  sigpdf->plotOn(plot,ProjWData(*hdatasmall));
  plot->Draw();
  craw->SaveAs("RawE.eps");
  craw->SetLogy();
  plot->SetMinimum(0.1);
  craw->SaveAs("RawElog.eps");
  
/*  new TCanvas;
  RooPlot *plotsig = tgtvar->frame(0.6,1.2,100);
  hdatasig->plotOn(plotsig);
  sigpdf.plotOn(plotsig,ProjWData(*hdatasig));
  plotsig->Draw(); */ 
  
  TCanvas *cmean = new TCanvas;
  RooPlot *plotmean = meanvar->frame(0.8,2.0,100);
  hdataclone->plotOn(plotmean);
  plotmean->Draw();
  cmean->SaveAs("mean.eps");
  
  
  TCanvas *cwidth = new TCanvas;
  RooPlot *plotwidth = widthvar->frame(0.,0.05,100);
  hdataclone->plotOn(plotwidth);
  plotwidth->Draw();
  cwidth->SaveAs("width.eps");
  
  TCanvas *cn = new TCanvas;
  RooPlot *plotn = nvar->frame(0.,111.,200);
  hdataclone->plotOn(plotn);
  plotn->Draw();
  cn->SaveAs("n.eps");

  TCanvas *cn2 = new TCanvas;
  RooPlot *plotn2 = n2var->frame(0.,111.,100);
  hdataclone->plotOn(plotn2);
  plotn2->Draw();
  cn2->SaveAs("n2.eps");
  
  
  TCanvas *calpha = new TCanvas;
  RooPlot *plotalpha = alphavar->frame(0.,6.,200);
  hdataclone->plotOn(plotalpha);
  plotalpha->Draw();    
  calpha->SaveAs("alpha.eps");
  calpha->SetLogy();
  plotalpha->SetMinimum(0.1);
  
  TCanvas *calpha2 = new TCanvas;
  RooPlot *plotalpha2 = alpha2var->frame(0.,6.,200);
  hdataclone->plotOn(plotalpha2);
  plotalpha2->Draw();      
  calpha2->SaveAs("alpha2.eps");
  
  TCanvas *ceta = new TCanvas;
  RooPlot *ploteta = scetavar->frame(-2.6,2.6,200);
  hdataclone->plotOn(ploteta);
  ploteta->Draw();      
  ceta->SaveAs("eta.eps");  
  
  //TH1 *heold = hdatasigtest->createHistogram("heold",testvar);
  //TH1 *heraw = hdata->createHistogram("heraw",*tgtvar,Binning(800,0.,2.));
  TH1 *heraw = hdata->createHistogram("hraw",*rawvar,Binning(800,0.,2.));
  TH1 *hecor = hdata->createHistogram("hecor",*ecorvar);
  
  
  //heold->SetLineColor(kRed);
  hecor->SetLineColor(kBlue);
  heraw->SetLineColor(kMagenta);
  
  hecor->GetXaxis()->SetRangeUser(0.6,1.2);
  //heold->GetXaxis()->SetRangeUser(0.6,1.2);
  
  TCanvas *cresponse = new TCanvas;
  
  hecor->Draw("HIST");
  //heold->Draw("HISTSAME");
  heraw->Draw("HISTSAME");
  cresponse->SaveAs("response.eps");
  cresponse->SetLogy();
  cresponse->SaveAs("responselog.eps");
  
  
  printf("make fine histogram\n");
  TH1 *hecorfine = hdata->createHistogram("hecorfine",*ecorvar,Binning(20e3,0.,2.));

  printf("calc effsigma\n");
  
  double effsigma = effSigma(hecorfine);
  
  printf("effsigma = %5f\n",effsigma);
  
/*  new TCanvas;
  RooPlot *ploteold = testvar.frame(0.6,1.2,100);
  hdatasigtest->plotOn(ploteold);
  ploteold->Draw();    
  
  new TCanvas;
  RooPlot *plotecor = ecorvar->frame(0.6,1.2,100);
  hdatasig->plotOn(plotecor);
  plotecor->Draw(); */   
  
  
}