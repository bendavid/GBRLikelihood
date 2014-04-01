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

void eregtestingExample(bool dobarrel=true, bool doele=true) {
  
  //output dir
  TString dirname = "/data/bendavid/eregexampletest/eregexampletest_test/"; 
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);    
  
  //read workspace from training
  TString fname;
  if (doele && dobarrel) 
    fname = "wereg_ele_eb.root";
  else if (doele && !dobarrel) 
    fname = "wereg_ele_ee.root";
  else if (!doele && dobarrel) 
    fname = "wereg_ph_eb.root";
  else if (!doele && !dobarrel) 
    fname = "wereg_ph_ee.root";
  
  TString infile = TString::Format("/data/bendavid/eregexampletest/%s",fname.Data());
  
  TFile *fws = TFile::Open(infile); 
  RooWorkspace *ws = (RooWorkspace*)fws->Get("wereg");
  
  //read variables from workspace
  RooGBRTargetFlex *meantgt = static_cast<RooGBRTargetFlex*>(ws->arg("sigmeant"));  
  RooRealVar *tgtvar = ws->var("tgtvar");
  
  
  RooArgList vars;
  vars.add(meantgt->FuncVars());
  vars.add(*tgtvar);
   
  //read testing dataset from TTree
  RooRealVar weightvar("weightvar","",1.);

  TTree *dtree;
  
  if (doele) {
    //TFile *fdin = TFile::Open("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-zllm50-v7n_noskim.root");
    TFile *fdin = TFile::Open("/data/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-zllm50-v7n_noskim.root");

    TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("PhotonTreeWriterSingleInvert");
    dtree = (TTree*)ddir->Get("hPhotonTreeSingle");       
  }
  else {
    TFile *fdin = TFile::Open("root://eoscms.cern.ch///eos/cms/store/cmst3/user/bendavid/idTreesAug1/hgg-2013Final8TeV_ID_s12-h124gg-gf-v7n_noskim.root");
    TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("PhotonTreeWriterPreselNoSmear");
    dtree = (TTree*)ddir->Get("hPhotonTreeSingle");       
  }
  
  //selection cuts for testing
  TCut selcut;
  if (dobarrel) 
    selcut = "ph.genpt>25. && ph.isbarrel && ph.ispromptgen"; 
  else
    selcut = "ph.genpt>25. && !ph.isbarrel && ph.ispromptgen"; 
  
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
  
  if (doele) 
    weightvar.SetTitle(prescale100alt*selcut);
  else
    weightvar.SetTitle(selcut);
  
  //make testing dataset
  RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",dtree,vars,weightvar);   

  if (doele) 
    weightvar.SetTitle(prescale1000alt*selcut);
  else
    weightvar.SetTitle(prescale10alt*selcut);
  //make reduced testing dataset for integration over conditional variables
  RooDataSet *hdatasmall = RooTreeConvert::CreateDataSet("hdatasmall",dtree,vars,weightvar);     
    
  //retrieve full pdf from workspace
  RooAbsPdf *sigpdf = ws->pdf("sigpdf");
  
  //input variable corresponding to sceta
  RooRealVar *scetavar = ws->var("var_1");
  
  //regressed output functions
  RooAbsReal *sigmeanlim = ws->function("sigmeanlim");
  RooAbsReal *sigwidthlim = ws->function("sigwidthlim");
  RooAbsReal *signlim = ws->function("signlim");
  RooAbsReal *sign2lim = ws->function("sign2lim");

  //formula for corrected energy/true energy ( 1.0/(etrue/eraw) * regression mean)
  RooFormulaVar ecor("ecor","","1./(@0)*@1",RooArgList(*tgtvar,*sigmeanlim));
  RooRealVar *ecorvar = (RooRealVar*)hdata->addColumn(ecor);
  ecorvar->setRange(0.,2.);
  ecorvar->setBins(800);
  
  //formula for raw energy/true energy (1.0/(etrue/eraw))
  RooFormulaVar raw("raw","","1./@0",RooArgList(*tgtvar));
  RooRealVar *rawvar = (RooRealVar*)hdata->addColumn(raw);
  rawvar->setRange(0.,2.);
  rawvar->setBins(800);

  //clone data and add regression outputs for plotting
  RooDataSet *hdataclone = new RooDataSet(*hdata,"hdataclone");
  RooRealVar *meanvar = (RooRealVar*)hdataclone->addColumn(*sigmeanlim);
  RooRealVar *widthvar = (RooRealVar*)hdataclone->addColumn(*sigwidthlim);
  RooRealVar *nvar = (RooRealVar*)hdataclone->addColumn(*signlim);
  RooRealVar *n2var = (RooRealVar*)hdataclone->addColumn(*sign2lim);
  
  
  //plot target variable and weighted regression prediction (using numerical integration over reduced testing dataset)
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
  
  //plot distribution of regressed functions over testing dataset
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
  
  TCanvas *ceta = new TCanvas;
  RooPlot *ploteta = scetavar->frame(-2.6,2.6,200);
  hdataclone->plotOn(ploteta);
  ploteta->Draw();      
  ceta->SaveAs("eta.eps");  
  

  //create histograms for eraw/etrue and ecor/etrue to quantify regression performance
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