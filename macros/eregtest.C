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

void eregtest() {
  
  TFile *fws = TFile::Open("/afs/cern.ch/work/b/bendavid/bare/eregtesteleJul8_sig5_01_2000/wereg.root");
  RooWorkspace *ws = (RooWorkspace*)fws->Get("wereg");
  
  RooGBRFunction *func = static_cast<RooGBRFunction*>(ws->arg("func"));
  RooRealVar *tgtvar = ws->var("tgtvar");
  
  
  RooArgList vars;
  vars.add(func->Vars());
  vars.add(*tgtvar);
  
  
  RooRealVar weightvar("weightvar","",1.);

  //TFile *fdin = TFile::Open("/home/mingyang/cms/hist/hgg-2013Moriond/merged/hgg-2013Moriond_s12-diphoj-3-v7a_noskim.root");
  TFile *fdin = TFile::Open("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/trainingtreesJul1/hgg-2013Final8TeV_s12-zllm50-v7n_noskim.root");
  TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("PhotonTreeWriterSingleInvert");
  TTree *dtree = (TTree*)ddir->Get("hPhotonTreeSingle");    
  
/*  TFile *fdinsig = TFile::Open("/home/mingyang/cms/hist/hgg-2013Moriond/merged/hgg-2013Moriond_s12-h125gg-gf-v7a_noskim.root");
  TDirectory *ddirsig = (TDirectory*)fdinsig->FindObjectAny("PhotonTreeWriterPreselNoSmear");
  TTree *dtreesig = (TTree*)ddirsig->Get("hPhotonTreeSingle"); */     
  
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
  TCut prescale100alt = "(evt%100==1)";
  TCut prescale1000alt = "(evt%1000==1)";
  //TCut oddevents = prescale100;
  
  weightvar.SetTitle(prescale100alt*selcut);
  RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",dtree,vars,weightvar);   

  weightvar.SetTitle(prescale1000alt*selcut);
  RooDataSet *hdatasmall = RooTreeConvert::CreateDataSet("hdatasmall",dtree,vars,weightvar);     
  
  RooAbsPdf *sigpdf = ws->pdf("sigpdf");
  
  RooAbsReal *sigmeanlim = ws->function("sigmeanlim");
  RooAbsReal *sigwidthlim = ws->function("sigwidthlim");

  

  RooFormulaVar ecor("ecor","","@0*@1",RooArgList(*tgtvar,*sigmeanlim));
  //RooFormulaVar ecor("ecor","","@0",RooArgList(*tgtvar));
  RooRealVar *ecorvar = (RooRealVar*)hdata->addColumn(ecor);
  ecorvar->setRange(0.,2.);
  ecorvar->setBins(800);

/*  RooFormulaVar eraw("eraw","","@0",RooArgList(*tgtvar));
  RooRealVar *erawvar = (RooRealVar*)hdatasig->addColumn(eraw);
  erawvar->setRange(0.,2.);
  erawvar->setBins(400); */ 
  
  RooDataSet *hdataclone = new RooDataSet(*hdata,"hdataclone");
  RooRealVar *meanvar = (RooRealVar*)hdataclone->addColumn(*sigmeanlim);
  RooRealVar *widthvar = (RooRealVar*)hdataclone->addColumn(*sigwidthlim);
  
  new TCanvas;
  RooPlot *plot = tgtvar->frame(0.6,1.2,100);
  hdata->plotOn(plot);
  sigpdf->plotOn(plot,ProjWData(*hdatasmall));
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
  
  //TH1 *heold = hdatasigtest->createHistogram("heold",testvar);
  TH1 *heraw = hdata->createHistogram("heraw",*tgtvar,Binning(800,0.,2.));
  TH1 *hecor = hdata->createHistogram("hecor",*ecorvar);
  
  
  //heold->SetLineColor(kRed);
  hecor->SetLineColor(kBlue);
  heraw->SetLineColor(kMagenta);
  
  hecor->GetXaxis()->SetRangeUser(0.6,1.2);
  //heold->GetXaxis()->SetRangeUser(0.6,1.2);
  
  new TCanvas;
  
  hecor->Draw("HIST");
  //heold->Draw("HISTSAME");
  heraw->Draw("HISTSAME");
  
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