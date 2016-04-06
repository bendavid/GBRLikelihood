//Example macro for training a simple BDT classifier using the GBRLikelihood machinery
//J.Bendavid

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
#include "RooGenericPdf.h"
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
#include "MCGBRIntegrator.h"
#include "TF1.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TStyle.h"
  
void mcinttest() {
     
  TString dirname = "plotstest";
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);

  MCGBRIntegrator mcint("mcint","",1e3,2e4);
//   mcint.SetNEventsBagged(5e5);
//   mcint.SetMinCutSignificance(std::numeric_limits<double>::min());
  mcint.SetMinCutSignificance(5.);
// //   mcint.SetMinCutSignificance(1e-9);
//     mcint.SetMinCutSignificance(10.);
//   mcint.SetMaxDepth(9);
  mcint.SetMinEvents(10);
//   mcint.SetMinEvents(20);
//   mcint.SetShrinkage(0.1);
  mcint.SetShrinkage(0.05);
//   mcint.SetDoEnvelope(true);
//   mcint.SetStagedGeneration(true);
//   mcint.SetNEventsInitial(1e6);
//   mcint.SetMaxNodes(1e3);    
  mcint.TrainForest(250);
  
  
//   const MCGBRTreeD *gentree = mcint.GenTree();
  const MCGBRForest *gentree = mcint.Forest();
  const MCGBRTreeD *lasttree = &mcint.Forest()->Trees().back();
  
  int nbins = 501;
  double low = 0.;
  double high = 1.;
  double step = (high-low)/double(nbins);
  
  double rlow = 0.;
  double rhigh = 2.0;
  double rstep = (rhigh-rlow)/double(nbins);
  
//   gStyle->SetOptStat(0);
  
  TGraph *hfunc = new TGraph(nbins);
  TGraph *hfunc2 = new TGraph(nbins);
  TGraph *hfunc3 = new TGraph(nbins);
  TGraph *hfunc4 = new TGraph(nbins);
//   TGraph2D *hfunc2d = new TGraph2D(nbins*nbins);
  TH2D *hfunc2d = new TH2D("hfunc2d","",nbins,-5.,5.,nbins,-10.,log(2.));
  TH2D *hfunc2d2 = new TH2D("hfunc2d2","",nbins,-5.,5.,nbins,0.,2.);
  TH2D *hfunc2d3 = new TH2D("hfunc2d3","",nbins,-5.,5.,nbins,0.,2.);
  TH2D *hfunc2d4 = new TH2D("hfunc2d4","",nbins,-5.,5.,nbins,-10.,log(2.));
  
  
  TGraph *hfunce = new TGraph(nbins);
  TGraph *hfunce2 = new TGraph(nbins);
  TGraph *hfunce3 = new TGraph(nbins);  
  TGraph *hfunce4 = new TGraph(nbins);  

  std::vector<float> vals(10,0.);
  std::vector<double> valsd(10,0.);
  
  const unsigned int nvars = 4;
  
//   int ipoint2d=0;
 for (int ibin=0; ibin<nbins; ++ibin) {
    float x = low + ibin*step;
    
    for (unsigned int ivar=0; ivar<nvars; ++ivar) {
      vals[ivar] = x;
      valsd[ivar] = x;
//       printf("ivar = %i, x = %5f\n",ivar, x);
    }
    
    double targetmin = mcint.Forest()->GetResponse(vals.data());      
    double target = mcint.ForestGen()->GetResponse(vals.data());      
    
//     printf("targetmin = %5e\n",targetmin);

    double yval = target*exp(-targetmin);
//     double yval = exp(targetmin);
    double yvaltrue = mcint.Camel(nvars,valsd.data());
    double yvalprod = yval*target;
    
    hfunc->SetPoint(ibin,x,yval);
    hfunc2->SetPoint(ibin,x,yvaltrue);
    hfunc3->SetPoint(ibin,x,target);
    hfunc4->SetPoint(ibin,x,yvalprod);

    

  }

 for (int ibin=0; ibin<nbins; ++ibin) {
    float x = low + ibin*step;
    
    vals[0] = x;
    valsd[0] = x;
    
    for (unsigned int ivar=1; ivar<nvars; ++ivar) {
      vals[ivar] = 0.;
      valsd[ivar] = 0.;
//       printf("ivar = %i, x = %5f\n",ivar, x);
    }
    
    double targetmin = mcint.Forest()->GetResponse(vals.data());      
    double target = mcint.ForestGen()->GetResponse(vals.data());      
    
    double yval = target*exp(-targetmin);
//     double yval = exp(targetmin);
    double yvaltrue = mcint.Camel(nvars,valsd.data());
    double yvalprod = yval*target;

    
    
    hfunce->SetPoint(ibin,x,yval);
    hfunce2->SetPoint(ibin,x,yvaltrue);
    hfunce3->SetPoint(ibin,x,target);
    hfunce4->SetPoint(ibin,x,yvalprod);

    

  }  
  
  {
    hfunc2->SetLineColor(kRed);
    hfunc2->SetMarkerColor(kRed);
    
    hfunc3->SetLineColor(kBlue);
    hfunc3->SetMarkerColor(kBlue);
    
    hfunc4->SetLineColor(kMagenta);
    hfunc4->SetMarkerColor(kMagenta);
    
    TCanvas *cfunc = new TCanvas;
    hfunc2->Draw("APL");
    hfunc->Draw("LPSAME");
    hfunc3->Draw("LPSAME");
//     hfunc4->Draw("LPSAME");
    
    hfunc2->GetXaxis()->SetTitle("d, distance along multidimensional diagonal (a.u.)");
    hfunc2->GetYaxis()->SetTitle("function value (a.u.)");
    
    TLegend *leg = new TLegend(0.62,0.6, 0.88,0.9);
    leg->AddEntry(hfunc2,"f(#bar{x}) (Camel)","L");
    leg->AddEntry(hfunc,"e^{h(#bar{x})} (Primary BDT)","L");
    leg->AddEntry(hfunc3,"g(#bar{x}) (Secondary BDT)","L");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();
    
    cfunc->SaveAs("func.pdf");

    TCanvas *cfunclog = new TCanvas;
    hfunc2->Draw("APL");
    hfunc->Draw("LPSAME");
    hfunc3->Draw("LPSAME");
//     hfunc4->Draw("LPSAME");
    
    cfunclog->SetLogy();
    
    TLegend *leglog = new TLegend(0.25,0.15, 0.65,0.45);
    leglog->AddEntry(hfunc2,"f(#bar{x}) (Camel)","L");
    leglog->AddEntry(hfunc,"e^{h(#bar{x})} (Primary BDT)","L");
    leglog->AddEntry(hfunc3,"g(#bar{x}) (Secondary BDT)","L");
    leglog->SetBorderSize(0);
    leglog->SetFillStyle(0);
    leglog->Draw();
    
    cfunclog->SaveAs("funclog.pdf");
  }

  {
    hfunce2->SetLineColor(kRed);
    hfunce2->SetMarkerColor(kRed);
    
    hfunce3->SetLineColor(kBlue);
    hfunce3->SetMarkerColor(kBlue);
        
    hfunce4->SetLineColor(kMagenta);
    hfunce4->SetMarkerColor(kMagenta);
    
    TCanvas *cfunc = new TCanvas;
    hfunce3->Draw("APL");
    hfunce->Draw("LPSAME");
    hfunce2->Draw("LPSAME");
//     hfunce4->Draw("LPSAME");
    
    hfunce2->GetXaxis()->SetTitle("d, distance along multidimensional diagonal (a.u.)");
    hfunce2->GetYaxis()->SetTitle("function value (a.u.)");
    
    TLegend *leg = new TLegend(0.62,0.6, 0.88,0.9);
    leg->AddEntry(hfunce2,"f(#bar{x}) (Camel)","L");
    leg->AddEntry(hfunce,"e^{h(#bar{x})} (Primary BDT)","L");
    leg->AddEntry(hfunce3,"g(#bar{x}) (Secondary BDT)","L");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();
    
    cfunc->SaveAs("func.pdf");

    TCanvas *cfunclog = new TCanvas;
    hfunce3->Draw("APL");
    hfunce->Draw("LPSAME");
    hfunce2->Draw("LPSAME");
//     hfunce4->Draw("LPSAME");
    
    cfunclog->SetLogy();
    
    TLegend *leglog = new TLegend(0.25,0.15, 0.65,0.45);
    leglog->AddEntry(hfunce2,"f(#bar{x}) (Camel)","L");
    leglog->AddEntry(hfunce,"e^{h(#bar{x})} (Primary BDT)","L");
    leglog->AddEntry(hfunce3,"g(#bar{x}) (Secondary BDT)","L");
    leglog->SetBorderSize(0);
    leglog->SetFillStyle(0);
    leglog->Draw();
    
    cfunclog->SaveAs("funclog.pdf");
  }
  
//   new TCanvas;
//   hfunc->Draw("APL");
//   
//   new TCanvas;
//   hfunce->Draw("APL");  
// 
//   new TCanvas;
//   hfunc4->Draw("APL");
//   
//   new TCanvas;
//   hfunce4->Draw("APL");  
  
}
