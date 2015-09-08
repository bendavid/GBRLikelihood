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
#include "RooRevCBFast.h"
#include "RooGausDoubleExp.h"
#include "TStyle.h"
#include "Riostream.h"
 
using namespace RooFit;

 
void fitDataSet(bool dobarrel=true, bool doele=false) {
       

  	//create RooRealVar for target
  	RooRealVar *tgtvar1 = new RooRealVar("tgtvar1","log(etrue/scRawEnergy)",1.);
       if (!dobarrel) tgtvar1->SetTitle("log(etrue/(scRawEnergy + scPreshowerEnergy))"); 

       RooRealVar *tgtvar2 = new RooRealVar("tgtvar2","log(etrue/scRawEnergy)",1.);
  	if (!dobarrel) tgtvar2->SetTitle("log(etrue/(scRawEnergy + scPreshowerEnergy))");  

       RooRealVar *tgtvar3 = new RooRealVar("tgtvar3","log(etrue/scRawEnergy)",1.);
  	if (!dobarrel) tgtvar3->SetTitle("log(etrue/(scRawEnergy + scPreshowerEnergy))");  

       RooRealVar *tgtvar4 = new RooRealVar("tgtvar4","log(etrue/scRawEnergy)",1.);
  	if (!dobarrel) tgtvar4->SetTitle("log(etrue/(scRawEnergy + scPreshowerEnergy))"); 
       
    /*   RooRealVar *tgtvar1 = new RooRealVar("tgtvar1","etrue/scRawEnergy",1.);
  	if (!dobarrel) tgtvar1->SetTitle("etrue/(scRawEnergy + scPreshowerEnergy)"); 

       RooRealVar *tgtvar2 = new RooRealVar("tgtvar2","etrue/scRawEnergy",1.);
  	if (!dobarrel) tgtvar2->SetTitle("etrue/(scRawEnergy + scPreshowerEnergy)");  

       RooRealVar *tgtvar3 = new RooRealVar("tgtvar3","etrue/scRawEnergy",1.);
  	if (!dobarrel) tgtvar3->SetTitle("etrue/(scRawEnergy + scPreshowerEnergy)");  

       RooRealVar *tgtvar4 = new RooRealVar("tgtvar4","etrue/scRawEnergy",1.);
  	if (!dobarrel) tgtvar4->SetTitle("etrue/(scRawEnergy + scPreshowerEnergy)"); 
  	*/
       RooArgList vars1, vars2, vars3, vars4;

  	//add target to full list
  	vars1.addOwned(*tgtvar1);
       vars2.addOwned(*tgtvar2);
       vars3.addOwned(*tgtvar3);
       vars4.addOwned(*tgtvar4);
  	  
  	//RooRealVar for event weight 
  	RooRealVar weightvar_lelr("weightvar","",1.);
       RooRealVar weightvar_lehr("weightvar","",1.);
       RooRealVar weightvar_hehr("weightvar","",1.);
       RooRealVar weightvar_helr("weightvar","",1.);
	
	TChain *tree;
    	tree = new TChain("promptTree");
    	tree->Add("/afs/cern.ch/user/m/musella/public/forKenza/gam_gam_phys14_v5_regtraining_v3.root");  
	

  	TCut selcutLowEtaLowR9, selcutLowEtaHighR9, selcutHighEtaHighR9,selcutHighEtaLowR9;
  	if (dobarrel) 
	{
    		selcutLowEtaLowR9 = "pt>200. && abs(scEta)<1.5 && kSaturated[12]!=1 && r9 <0.94" ;
              selcutLowEtaHighR9 = "pt>200. && abs(scEta)<1.5 && kSaturated[12]!=1  && r9 >=0.94" ; 
              selcutHighEtaHighR9 = "pt>200. && abs(scEta)<1.5 && kSaturated[12]!=1 && r9 >=0.94" ;
              selcutHighEtaLowR9 = "pt>200. && abs(scEta)<1.5 && kSaturated[12]!=1 && r9 <0.94" ;
  	}
  	if(!dobarrel)
	{
    		selcutLowEtaLowR9 = "pt>200. && abs(scEta)>1.5 && kSaturated[12]!=1 && abs(scEta)<=1.8 && abs(scEta)>=1.6 && r9 <0.94" ;
              selcutLowEtaHighR9 = "pt>200. && abs(scEta)>1.5 && kSaturated[12]!=1 && abs(scEta)<=1.8 && abs(scEta)>=1.6 && r9 >=0.94" ; 
              selcutHighEtaHighR9 = "pt>200. && abs(scEta)>1.5 && kSaturated[12]!=1 && abs(scEta)<=2.3 && abs(scEta)>=2.1 && r9 >=0.94" ;
              selcutHighEtaLowR9 = "pt>200. && abs(scEta)>1.5 && kSaturated[12]!=1";// && abs(scEta)<=2.3 && abs(scEta)>=2.1 && r9 <0.94" ;
              
              //selcutHighEtaLowR9 = "pt>200. && abs(scEta)>1.5 && kSaturated[12]!=1 &&  r9 <0.94" ;
              //selcutHighEtaHighR9 = "pt>200. && abs(scEta)>1.5 && kSaturated[12]!=1 && r9 >=0.94" ;   
  	}

       TCut selweight = "1."; 
  	TCut prescale20 = "(event%20==0)";
  	TCut prescale100 = "(event%100==0)";  
  
	//TCut condition = "(etrue/(scRawEnergy + scPreshowerEnergy) >=0.7)";
	
  	//weightvar title used for per-event weights and selection cuts
  	
  	weightvar_lelr.SetTitle(selweight*selcutLowEtaLowR9); ///////////////////////////
       weightvar_lehr.SetTitle(selweight*selcutLowEtaHighR9);
  	weightvar_hehr.SetTitle(selweight*selcutHighEtaHighR9);
       weightvar_helr.SetTitle(selweight*selcutHighEtaLowR9);

  	//create RooDataSet from TChain

       
  	RooDataSet *hdata_lelr = RooTreeConvert::CreateDataSet("hdata1",tree,vars1,weightvar_lelr); 
       RooDataSet *hdata_lehr = RooTreeConvert::CreateDataSet("hdata2",tree,vars2,weightvar_lehr);
       RooDataSet *hdata_hehr = RooTreeConvert::CreateDataSet("hdata3",tree,vars3,weightvar_hehr);
       RooDataSet *hdata_helr = RooTreeConvert::CreateDataSet("hdata4",tree,vars4,weightvar_helr); 
  
       /*TCanvas *c=new TCanvas("c","",1300,800);

       RooPlot *frame1 = tgtvar1->frame(0.6,2.);
       hdata->plotOn(frame1,Name("hdata"));
       hdata2->plotOn(frame1,MarkerColor(kRed), LineColor(kRed),Name("hdata2"));
       frame1->Draw();
       TLegend *leg = new TLegend(0.6,0.8,0.9,0.9);
       leg->AddEntry("hdata","low eta","p");
       leg->AddEntry("hdata2","high eta","p");
       leg->Draw("SAME");*/
       Double_t tgtmin=-0.5,tgtmax=0.6;
       Double_t n1min=1.01,n1max=100, n1val=3., n2min=1.01,n2max=100, n2val=3.;
	
       tgtvar1->setRange(tgtmin,tgtmax);
       tgtvar1->setBins(800.);

       RooRealVar mean1("mean1","mean1",0.,-0.5,0.6);
       RooRealVar sigma1("sigma1","sigma1",0.01,0.0002,0.5);
       RooRealVar alpha11("alpha11","alpha11",3.,0.05,10);
       RooRealVar alpha21("alpha21","alpha21",3.,0.05,10);

       RooRealVar n11("n11","n11",n1val,n1min,n1max);
       RooRealVar n21("n21","n21",n2val,n2min,n2max);
       mean1.setConstant(false);
 	
       tgtvar2->setRange(tgtmin,tgtmax);
       tgtvar2->setBins(800.);
       RooRealVar mean2("mean2","mean2",0.,tgtmin,tgtmax);
       RooRealVar sigma2("sigma2","sigma2",0.01,0.0002,0.5);
       RooRealVar alpha12("alpha12","alpha12",3,0.5,100);
       RooRealVar alpha22("alpha22","alpha22",3,0.5,10);

       RooRealVar n12("n12","n12",n1val,n1min,n1max);
       RooRealVar n22("n22","n22",n2val,n2min,n2max);
       mean2.setConstant(false);

       tgtvar3->setRange(tgtmin,tgtmax);
       tgtvar3->setBins(800.);
       RooRealVar mean3("mean3","mean3",0.,-0.04,0.04);
       RooRealVar sigma3("sigma3","sigma3",0.002,0.0002,0.5);
       RooRealVar alpha13("alpha13","alpha13",3.,0.05,10);
       RooRealVar alpha23("alpha23","alpha23",3.,0.05,10);

       RooRealVar n13("n13","n13",n1val,n1min,n1max);
       RooRealVar n23("n23","n22",n2val,n2min,n2max);
       mean3.setConstant(false);

       tgtvar4->setRange(tgtmin,tgtmax);
       tgtvar4->setBins(800.);
       RooRealVar mean4("mean4","mean4",0.,-0.04,0.04);
       RooRealVar sigma4("sigma4","sigma4",0.002,0.0002,0.5);
       RooRealVar alpha14("alpha14","alpha14",3.,0.05,10);
       RooRealVar alpha24("alpha24","alpha24",3.,0.05,10);

       RooRealVar n14("n14","n14",n1val,n1min,n1max);
       RooRealVar n24("n24","n24",n2val,n1min,n1max);
       mean4.setConstant(false);
/*
       RooRealVar mean_2("mean_2","mean_2",1.,0.2,2);
       RooRealVar sigma_2("sigma_2","sigma_2",0.002,0.0002,0.5);
       RooRealVar alpha1_2("alpha1_2","alpha1_2",1,0.5,100);
       RooRealVar n1_2("n1_2","n1_2",1,0.5,100);
       RooRealVar alpha2_2("alpha2_2","alpha2_2",1,0.5,100);
       RooRealVar n2_2("n2_2","n2_2",1,0.5,100);
       mean_2.setConstant(false);
*/
       RooDoubleCBFast pdfCB1("pdfCB1","CB pdf1",*tgtvar1,mean1,sigma1,alpha11,n11,alpha21,n21);
       
       RooDoubleCBFast pdfCB2("pdfCB2","CB pdf2",*tgtvar2,mean2,sigma2,alpha12,n12,alpha22,n22);
       RooDoubleCBFast pdfCB3("pdfCB3","CB pdf3",*tgtvar3,mean3,sigma3,alpha13,n13,alpha23,n23);
       RooDoubleCBFast pdfCB4("pdfCB4","CB pdf4",*tgtvar4,mean4,sigma4,alpha14,n14,alpha24,n24);
       //RooDoubleCBFast pdfCB_2("pdfCB_2","CB pdf",*tgtvar,mean_2,sigma_2,alpha1_2,n1_2,alpha2_2,n2_2);
       

       //Canvas
       TCanvas *c = new TCanvas("tgtVar","",1300,800);
       c->Divide(2,2);
       

      /* RooDoubleCBFast *pdfCB;
       RooPlot *xframe1 = tgtvar1->frame(Title("fit of data"));
       hdata_lelr->plotOn(xframe1);

       Double_t min(-0.015), max(0.01);
       for(Int_t i=0;i<10;i++) {
              pdfCB = new RooDoubleCBFast("pdfCB","CB pdf",*tgtvar1,mean1,sigma1,alpha11,n11,alpha21,n21);
              pdfCB->fitTo(*hdata_lelr,Range(min,max));  
              min-=0.01;
              max+=0.01;
              //
              if(i+1!=10) {
                     pdfCB->plotOn(xframe1,LineColor(i+1));
              }
              else {
                     pdfCB->plotOn(xframe1,LineColor(28));
              }
              
       }
       xframe1->Draw();*/
       
      // RooRealVar coeff("coeff","coeff",0.5,0.,1.);
       
      // RooAddPdf sumCB("sumCB","sum of CB pdfs",RooArgList(pdfCB,pdfCB_2),coeff);
      // sumCB.fitTo(*hdata);
       

       //sumCB.plotOn(xframe);
       //sumCB.plotOn(xframe,Components(pdfCB),LineStyle(kDashed),LineColor(kRed));
       //sumCB.plotOn(xframe,Components(pdfCB_2),LineStyle(kDashed),LineColor(kRed));


//Canvas
       //TCanvas *canvas_tgtvar2 = new TCanvas("tgtVar2","",800,800);
       c->cd(1);
       gPad->SetGridx();
       gPad->SetGridy();

       RooPlot *xframe1 = tgtvar1->frame(Title("fit of data_lelr"));
       hdata_lelr->plotOn(xframe1);
       pdfCB1.fitTo(*hdata_lelr);  
       pdfCB1.plotOn(xframe1);
       xframe1->Draw();

       c->cd(2);
       gPad->SetGridx();
       gPad->SetGridy();

       RooPlot *xframe2 = tgtvar2->frame(Title("fit of data_lehr"));
       hdata_lehr->plotOn(xframe2);


       //     
       pdfCB2.fitTo(*hdata_lehr);  
       pdfCB2.plotOn(xframe2);
       xframe2->Draw();

       c->cd(3);
       gPad->SetGridx();
       gPad->SetGridy();

       RooPlot *xframe3 = tgtvar3->frame(Title("fit of data_hehr"));
       hdata_hehr->plotOn(xframe3);


       //     
       pdfCB3.fitTo(*hdata_hehr);  
       pdfCB3.plotOn(xframe3);
       xframe3->Draw();

       c->cd(4);
       gPad->SetGridx();
       gPad->SetGridy();

       RooPlot *xframe4 = tgtvar4->frame(Title("fit of data_helr"));
       hdata_helr->plotOn(xframe4);


       //     
       pdfCB4.fitTo(*hdata_helr);  
       pdfCB4.plotOn(xframe4);
       xframe4->Draw();
	

       delete hdata_lelr;  
       delete hdata_lehr;
       delete hdata_hehr;
       delete hdata_helr;
  	return;
  
}			
