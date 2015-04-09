
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

void wmergepdf(bool dotight=false) {

  TString fnameprefix;
  if (dotight) {
    fnameprefix = "walttight";
  }
  else {
    fnameprefix = "waltloose";
  }
  
  RooWorkspace *waltmerged = new RooWorkspace("waltmerged");

  for (unsigned int icharge=0; icharge<2; ++icharge) {
    TString sample;
    
    RooArgList nllvars;
    
    if (icharge == 0) {
      sample = "wminus";
    }
    else if (icharge==1) {
      sample = "wplus";
    }
    
    double scale = 1.;
    
    TString pdfname = "expfalt";
    
    
    TString suffixup = TString::Format("%s_126",sample.Data());
    TString suffixdown = TString::Format("%s_76",sample.Data());
    TString suffixdowndown = TString::Format("%s_26",sample.Data());
    TString suffixupup = TString::Format("%s_176",sample.Data());
    
    TFile *finup = TFile::Open(TString::Format("%s_%s.root",fnameprefix.Data(),suffixup.Data()));
    RooWorkspace *wup = (RooWorkspace*)finup->Get(TString::Format("walt_%s",suffixup.Data()));
    RooAbsReal *expfup = (RooAbsReal*)wup->arg(TString::Format("%s_%s",pdfname.Data(),suffixup.Data()));
    nllvars.add(*expfup);
    
    TFile *findown = TFile::Open(TString::Format("%s_%s.root",fnameprefix.Data(),suffixdown.Data()));
    RooWorkspace *wdown = (RooWorkspace*)findown->Get(TString::Format("walt_%s",suffixdown.Data()));
    RooAbsReal *expfdown = (RooAbsReal*)wdown->arg(TString::Format("%s_%s",pdfname.Data(),suffixdown.Data()));
    nllvars.add(*expfdown);
    
    TFile *finupup = TFile::Open(TString::Format("%s_%s.root",fnameprefix.Data(),suffixupup.Data()));
    RooWorkspace *wupup = (RooWorkspace*)finupup->Get(TString::Format("walt_%s",suffixupup.Data()));
    RooAbsReal *expfupup = (RooAbsReal*)wupup->arg(TString::Format("%s_%s",pdfname.Data(),suffixupup.Data()));
    nllvars.add(*expfupup);    

    TFile *findowndown = TFile::Open(TString::Format("%s_%s.root",fnameprefix.Data(),suffixdowndown.Data()));
    RooWorkspace *wdowndown = (RooWorkspace*)findowndown->Get(TString::Format("walt_%s",suffixdowndown.Data()));
    RooAbsReal *expfdowndown = (RooAbsReal*)wdowndown->arg(TString::Format("%s_%s",pdfname.Data(),suffixdowndown.Data()));
    nllvars.add(*expfdowndown);
    
//     RooDataSet *data = (RooDataSet*)wup->data("dsdata");
    RooDataSet *refdata = (RooDataSet*)wup->data("dsref");
    
//     data->merge((RooDataSet*)wdown->data("dsdata"));
    refdata->merge((RooDataSet*)wdown->data("dsref"));

//     data->merge((RooDataSet*)wdowndown->data("dsdata"));
    refdata->merge((RooDataSet*)wdowndown->data("dsref"));

//     data->merge((RooDataSet*)wupup->data("dsdata"));
    refdata->merge((RooDataSet*)wupup->data("dsref"));
    
//     printf("sample %i, scale*sumEntries = %5f, ref sumEntries = %5f\n",icharge,scale*data->sumEntries(), refdata->sumEntries());
    printf("sample %i, ref sumEntries = %5f\n",icharge, refdata->sumEntries());
       
    
    
    TString suffixpdfnom = TString::Format("%s_309",sample.Data());
    
    TFile *finpdfnom = TFile::Open(TString::Format("%s_%s.root",fnameprefix.Data(),suffixpdfnom.Data()));
    RooWorkspace *wpdfnom = (RooWorkspace*)finpdfnom->Get(TString::Format("walt_%s",suffixpdfnom.Data()));
    RooAbsReal *expfpdfnom = (RooAbsReal*)wpdfnom->arg(TString::Format("%s_%s",pdfname.Data(),suffixpdfnom.Data()));
    nllvars.add(*expfpdfnom);
    
//     data->merge((RooDataSet*)wpdfnom->data("dsdata"));
    refdata->merge((RooDataSet*)wpdfnom->data("dsref"));    
    
//     RooArgList pdfall;
//     RooArgList pdfnormall;
    
    const unsigned int npdferr = 26;
    const unsigned int nompdfidx = 309;
//     const unsigned int nompdfidx = 367;
    for (unsigned int ipdferr=0; ipdferr<npdferr; ++ipdferr) {
      int pdfidxup = nompdfidx + 2*ipdferr + 1;
      int pdfidxdown = nompdfidx + 2*ipdferr + 2;
      
      TString suffixpdfup = TString::Format("%s_%i",sample.Data(),pdfidxup);
      TString suffixpdfdown = TString::Format("%s_%i",sample.Data(),pdfidxdown);
      
      TFile *finpdfup = TFile::Open(TString::Format("%s_%s.root",fnameprefix.Data(),suffixpdfup.Data()));
      RooWorkspace *wpdfup = (RooWorkspace*)finpdfup->Get(TString::Format("walt_%s",suffixpdfup.Data()));
      RooAbsReal *expfpdfup = (RooAbsReal*)wpdfup->arg(TString::Format("%s_%s",pdfname.Data(),suffixpdfup.Data()));
      nllvars.add(*expfpdfup);
      
//       data->merge((RooDataSet*)wpdfup->data("dsdata"));
      refdata->merge((RooDataSet*)wpdfup->data("dsref"));         
      
      TFile *finpdfdown = TFile::Open(TString::Format("%s_%s.root",fnameprefix.Data(),suffixpdfdown.Data()));
      RooWorkspace *wpdfdown = (RooWorkspace*)finpdfdown->Get(TString::Format("walt_%s",suffixpdfdown.Data()));
      RooAbsReal *expfpdfdown = (RooAbsReal*)wpdfdown->arg(TString::Format("%s_%s",pdfname.Data(),suffixpdfdown.Data()));
      nllvars.add(*expfpdfdown);
      
//       data->merge((RooDataSet*)wpdfdown->data("dsdata"));
      refdata->merge((RooDataSet*)wpdfdown->data("dsref"));      

      
    }
    
//     data->SetName(TString::Format("data_%s",sample.Data()));
    refdata->SetName(TString::Format("refdata_%s",sample.Data()));
    
//     waltmerged->import(*data);
    waltmerged->import(*refdata);
  }
  
  TString outname = TString::Format("%smerged.root",fnameprefix.Data());
  waltmerged->writeToFile(outname);
}