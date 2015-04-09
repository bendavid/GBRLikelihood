
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
#include "RooDataWeightedAverage.h"
#include "AsymPow.h"

#ifndef __CINT__
#include "WMassFitter.h"
#else
class WMassFitter;
#endif /* __CINT __ */



void wfitvvfast() {
  
  TFile *fin = TFile::Open("walttightmerged.root");
  RooWorkspace *win = (RooWorkspace*)fin->Get("waltmerged");
 
  RooDataSet *dataplus = (RooDataSet*)win->data("data_wplus");
  RooDataSet *refdataplus = (RooDataSet*)win->data("refdata_wplus");

  RooDataSet *dataminus = (RooDataSet*)win->data("data_wminus");
  RooDataSet *refdataminus = (RooDataSet*)win->data("refdata_wminus");
  
  std::vector<RooDataSet*> datav;
  datav.push_back(refdataplus);
  datav.push_back(refdataplus);
  datav.push_back(refdataminus);
  datav.push_back(refdataminus);
  
  WMassFitter fitter(datav,"expfalt");
  
// //   TFile *fin = TFile::Open("walttightmergedtemplate.root");
//   TFile *fin = TFile::Open("waltloosemergedtemplate.root");
//   RooWorkspace *win = (RooWorkspace*)fin->Get("waltmerged");
//  
//   RooDataSet *dataplus = (RooDataSet*)win->data("data_wplus");
//   RooDataSet *refdataplus = (RooDataSet*)win->data("refdata_wplus");
// 
//   RooDataSet *dataminus = (RooDataSet*)win->data("data_wminus");
//   RooDataSet *refdataminus = (RooDataSet*)win->data("refdata_wminus");
//   
//   std::vector<RooDataSet*> datav;
//   datav.push_back(dataplus);
//   datav.push_back(refdataplus);
//   datav.push_back(dataminus);
//   datav.push_back(refdataminus);
// 
// //   RooRealVar *mt = static_cast<RooRealVar*>(dataplus->get()->first());
// //   TH1 *hplus = dataplus->createHistogram("hplus",*mt);
// //   new TCanvas;
// //   hplus->Draw("HIST");
// //   return;
//   
//   WMassFitter fitter(datav,"hfuncalt");  
  
  

//   std::vector<RooDataSet*> datav;
//   datav.push_back(dataplus);
//   datav.push_back(refdataplus);
//   datav.push_back(dataminus);
//   datav.push_back(refdataminus);
  
//   WMassFitter fitter(dataplus,refdataplus,dataminus,refdataminus);
//    WMassFitter fitter(refdataplus,refdataplus,refdataminus,refdataminus);
  
//   WMassFitter fitter(datav,"hfuncalt");
  
  fitter.Fit();
  
}