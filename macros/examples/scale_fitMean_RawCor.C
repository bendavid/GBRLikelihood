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
#include "TH2.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TGraphErrors.h"

using namespace RooFit;

//function that fit the mean with a gaussian and plot the result for both raw/true and cor/true distributions
void fit_and_plot_mean(Bool_t dobarrel,RooRealVar* pt,int number, TString arg_variable, Double_t xmin, Double_t xmax, Int_t nbins, Double_t* xbins, Double_t step, Double_t init, TString nameXaxis, TString add_to_pngFileName, RooDataSet* hdata, RooRealVar* rawvar, RooRealVar* ecorvar, TCanvas* c_mean);

void scale_fitMean_RawCor(bool dobarrel=true, bool weighted=false) {
  
  	//output dir
  	TString dirname = "Barrel_Log_sig3_alpha2-3_evts15/";
  	gSystem->mkdir(dirname,true);
  	gSystem->cd(dirname);    
  	
	TCut cutUsed = "(event%15==1)";
       
  	//read workspace from training
  	TString fname; //name of the training output file
  	TString add_to_pngFileName, add_to_rootFileName; //extension of root and png files 
  	
  	if (dobarrel) 
 	{	
		fname = "wereg_ph_eb.root";
		add_to_pngFileName = "_eb.png";
              add_to_rootFileName = "_eb.root";
 	}
  	else 
 	{
 		fname = "wereg_ph_ee.root";
		add_to_pngFileName = "_ee.png";
              add_to_rootFileName = "_ee.root";
 	}
       add_to_rootFileName = add_to_rootFileName; //add particular info in file name
       add_to_pngFileName = add_to_pngFileName; 
	
	//Workspace 
  	TString infile = TString::Format("./%s",fname.Data());
  	TFile *fws = TFile::Open(infile); 
  	RooWorkspace *ws = (RooWorkspace*)fws->Get("wereg");
  
  	//read variables from workspace
  	RooGBRTargetFlex *meantgt = static_cast<RooGBRTargetFlex*>(ws->arg("sigmeant"));  
  	RooRealVar *tgtvar = ws->var("tgtvar");
	RooRealVar *etrue = new RooRealVar("etrue","etrue",0.);
	RooRealVar *scRawEnergy = new RooRealVar("scRawEnergy","scRawEnergy",0.);
	RooRealVar *r9 = new RooRealVar("nvtx","nvtx",0.);
	//RooRealVar *maxEnergyXtal = new RooRealVar("maxEnergyXtal","maxEnergyXtal",0.);
  
  	RooArgList vars;
  	vars.add(meantgt->FuncVars());
  	vars.add(*tgtvar);
	vars.add(*etrue);
	vars.add(*scRawEnergy);
	vars.add(*r9);
	//vars.add(*maxEnergyXtal);
       
   
  	//read testing dataset from TTree
  	RooRealVar weightvar("weightvar","",1.);

	TChain *dtree;
    dtree = new TChain("promptTree");
    dtree->Add("/afs/cern.ch/work/k/khoumani/CMSSW_7_2_3/src/HiggsAnalysis/GBRLikelihood/macros/examples/promptTree_FlatPt-300To3000.root");
	dtree->Add("/afs/cern.ch/work/k/khoumani/CMSSW_7_2_3/src/HiggsAnalysis/GBRLikelihood/macros/examples/promptTree_FlatPt-5To300_lessstat.root");
   
       /*if(true_energy) {
	       if(!weighted) {
                     dtree->Draw("etrue>>h1");
                     dtree->Draw("scEta>>h2","etrue<=2500");
                     dtree->Draw("scEta>>h3","etrue>2500 && etrue<=3000");
                     dtree->Draw("scEta>>h4","etrue>3000 && etrue<=5000");
              }
              else {
                     dtree->Draw("etrue>>h1","weight");
                     dtree->Draw("scEta>>h2","weight*(etrue<=2500)");
                     dtree->Draw("scEta>>h3","weight*(etrue>2500 && etrue<=3000)");
                     dtree->Draw("scEta>>h4","weight*(etrue>3000 && etrue<=5000)");
              }
       }
       else {
              if(!weighted) {
                     dtree->Draw("energy>>h1","pt>250. && kSaturated[12]==1 && abs(scEta)<=1.5 && abs(ieta[12])!=1 && abs(ieta[12])!=25 && abs(ieta[12])!=26 && abs(ieta[12])!=45 && abs(ieta[12])!=46 && abs(ieta[12])!=65 && abs(ieta[12])!=66 && abs(ieta[12])!=85");
                     //dtree->Draw("scEta>>h2","pt<=2000");
                     //dtree->Draw("scEta>>h3","pt>2000 && pt<=3000");
                     //dtree->Draw("scEta>>h4","pt>3000 && pt<=5000");
              }
              else {
                     dtree->Draw("pt>>h1","weight");
                     dtree->Draw("scEta>>h2","weight*(pt<=2500)");
                     dtree->Draw("scEta>>h3","weight*(pt>2500 && pt<=3000)");
                     dtree->Draw("scEta>>h4","weight*(pt>3000 && pt<=5000)");
              }
       }*/
       //dtree->Draw("etrue>>h5","weight*(pt<=2000)");
       //dtree->Draw("etrue>>h6","weight*(pt>2000 && pt<=3500)");
       //dtree->Draw("etrue>>h7","weight*(pt>3500 && pt<=5000)");
  	//selection cuts for testing
  	TCut selcut;
       
  	if (dobarrel)
		selcut = "etrue*pt/energy >=250. && etrue*pt/energy<=3000 && abs(scEta)<1.5 && kSaturated[12]!=1 && (eMax/e2nd-1)<0.3";
  	else
		selcut = "etrue*pt/energy >=250. && etrue*pt/energy<=3000 && abs(scEta)>1.5 && kSaturated[12]!=1 ";
  	
       
	TCut selweight;
	if(weighted)
		selweight= "(weight)";
	else
		selweight= "(1.)";
  	
    
  	//weightvar.SetTitle(cutUsed*selcut);
	weightvar.SetTitle(cutUsed*selcut*selweight);
	
  	//create the testing dataset
  	RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",dtree,vars,weightvar);
    
	//retrieve full pdf from workspace
	//RooAbsPdf *sigpdf = ws->pdf("sigpdf");
  	
	//regressed output functions
	RooAbsReal *sigmeanlim = ws->function("sigmeanlim");
	//RooAbsReal *sigwidthlim = ws->function("sigwidthlim");
    
    r9->setRange(0.,40.);
	RooRealVar *pt;
	RooRealVar *scetavar = ws->var("var_1");
	// RooRealVar *r9=ws->var("var_3");
	//r9->setRange(0.75,1.);
	if(dobarrel)
		scetavar->setRange(-1.5,1.5);
	else
		scetavar->setRange(-3,3);
	
	if(!dobarrel) {
		RooRealVar *preshowerEnergy = ws->var("var_28");
		//formula for Pt true_energy/Cosh(eta) = tgtvar*rawvar/cosh(eta) 
		RooFormulaVar ptFormula("pt","","exp(@0)*@1*(1+@2)*1./(cosh(@3))",RooArgList(*tgtvar,*scRawEnergy,*preshowerEnergy,*scetavar));
		pt = (RooRealVar*)hdata->addColumn(ptFormula);
		pt->setRange(250.,5000); //
	}
	else {
		RooFormulaVar ptFormula("pt","","exp(@0)*@1*1./(cosh(@2))",RooArgList(*tgtvar,*scRawEnergy,*scetavar));
		pt = (RooRealVar*)hdata->addColumn(ptFormula);
		pt->setRange(250.,5000); //
	}
	
    
	//create xbins that contains the fit slices limits 
	Int_t nbins=10;

	Double_t *xbins=0, step, init=0;
    
	step = 500.;
    
	xbins = (Double_t*)malloc(sizeof(Double_t)*(nbins+1));
	for(int ii=0;ii<=nbins; ii++) {
		xbins[ii]=init+step*ii;
	}
	
    
	//formula for corrected energy/true energy ( 1.0/(etrue/eraw) * regression mean)
	//RooFormulaVar ecor("ecor","","1./exp(@0)*exp(@1)",RooArgList(*tgtvar,*sigmeanlim));
	RooFormulaVar ecor("ecor","","1./exp(@0)*exp(@1)",RooArgList(*tgtvar,*sigmeanlim));
	RooRealVar *ecorvar = (RooRealVar*)hdata->addColumn(ecor);
	ecorvar->setRange(0.,2);
	ecorvar->setBins(800);
    
	//formula for raw energy/true energy
	//RooFormulaVar raw("raw","","1./exp(@0)",RooArgList(*tgtvar));
	RooFormulaVar raw("raw","","1./exp(@0)",RooArgList(*tgtvar));
	RooRealVar *rawvar = (RooRealVar*)hdata->addColumn(raw);
	rawvar->setRange(0.,2.);
	rawvar->setBins(800);
       
	TCanvas *c_mean = new TCanvas("c_mean","c_mean",800,800);

	fit_and_plot_mean(dobarrel,pt, 1,"pt",0.,5000.,nbins, xbins, step, init, "pT (GeV)", add_to_pngFileName, hdata, rawvar,ecorvar, c_mean);
	c_mean->SaveAs("scale_vs_pt_CMS2.png");
	c_mean->SaveAs("scale_vs_pt_CMS2.root");
	
	//  TCanvas *c_nvtx = new TCanvas("c_nvtx","c_nvtx",800,800);
	//fit_and_plot_mean(dobarrel,nvtx, 1,"nvtx",0.,40.,nbins, xbins, step, init, "nvtx", add_to_pngFileName, hdata, rawvar,ecorvar, c_mean);
	/*	nbins=10;
	init=0;
	step = 40;

	Double_t *xbins2=0;       
    
	xbins2 = (Double_t*)malloc(sizeof(Double_t)*(nbins+1));
	for(int ii=0;ii<=nbins; ii++) {
		xbins2[ii]=init+step*ii;
	}

	//TCanvas *c_r9 = new TCanvas("c_r9","c_r9",800,800);
	fit_and_plot_mean(dobarrel,r9,0,"nvtx",1,40,nbins, xbins2, step, init, "R_{9}", add_to_pngFileName, hdata, rawvar,ecorvar, c_mean);
    
   	c_mean->SaveAs("scale_vs_r9.png");*/
}
void fit_and_plot_mean(Bool_t dobarrel, RooRealVar* pt, int number, TString arg_variable, Double_t xmin, Double_t xmax, Int_t nbins, Double_t* xbins, Double_t step, Double_t init, TString nameXaxis, TString add_to_pngFileName, RooDataSet* hdata, RooRealVar* rawvar, RooRealVar* ecorvar, TCanvas* c_mean)
{
	
	TString str1, str2, str3, str4, variable, str;
	TH1F* h_gaussian_raw = new TH1F("h_gaussian_raw","",nbins,xbins);
	TH1F* h_gaussian_cor = new TH1F("h_gaussian_cor","",nbins,xbins);
	
	TH2 *h2D_average_cor, *h2D_average_raw;
	h2D_average_cor = hdata->createHistogram(*pt,*ecorvar);	
	h2D_average_raw = hdata->createHistogram(*pt,*rawvar);

 	TProfile *h1D_average_cor = h2D_average_cor->ProfileX();
	TProfile *h1D_average_raw = h2D_average_raw->ProfileX();
	
	//h1D_average_cor->SetErrorOption("s");
	h1D_average_cor->Rebin(nbins,"hprof_cor",xbins);
	TProfile *hprof_cor=(TProfile*)gDirectory->Get("hprof_cor");  
	//h1D_average_raw->SetErrorOption("s");
	h1D_average_raw->Rebin(nbins,"hprof_raw",xbins);

	TProfile *hprof_raw=(TProfile*)gDirectory->Get("hprof_raw");       
	
	TCanvas *craw = new TCanvas("craw","craw"); //------------------------------------------------------------CANVAS
	craw->Divide(5,2);
	TCanvas *ccor = new TCanvas("ccor","ccor"); //------------------------------------------------------------CANVAS
	ccor->Divide(5,2);
    
	for(Int_t i=0; i<nbins; i++) {
		
        
		variable = arg_variable;
        
		str1.Form("%d",i);
		str2.Form("%d",i+1);
		str3.Form("%f",step);
		str4.Form("%f",init);
		RooDataSet *hdata_reduced;
        
		str=variable+"> ("+str4+"+("+str3+"*"+str1+")) && "+variable+"< ("+str4+"+("+str3+"*"+str2+")) ";
        
		hdata_reduced = (RooDataSet*)hdata->reduce(str);
        
		RooDataSet *hdata2 = (RooDataSet*)hdata_reduced->reduce(*ecorvar);
		RooRealVar mean_cor("mean_cor","mean_cor",1,0.9,1.1);
		RooRealVar sig_cor("sig_cor","sig_cor",0.01,0.0002,0.8);
		RooRealVar a1_cor("a1_cor","a1_cor",3,0.05,10);
		RooRealVar a2_cor("a2_cor","a2_cor",3,0.05,10);
		RooRealVar n1_cor("n1_cor","n1_cor",3,1.01,500);
		RooRealVar n2_cor("n2_cor","n2_cor",3,1.01,500);
        
		RooDoubleCBFast pdfCB_cor("pdfCB_cor","pdfCB_cor",*ecorvar,mean_cor,sig_cor,a1_cor,n1_cor,a2_cor,n2_cor);
		cout << "dcb corrected" << endl;
		pdfCB_cor.fitTo(*hdata2,Range(0.9,1.1));
		h_gaussian_cor->AddBinContent(i+1,mean_cor.getVal());
		h_gaussian_cor->SetBinError(i+1,mean_cor.getError());
        cout << "histogram corrected filled" << endl;
		RooDataSet *hdata3 = (RooDataSet*)hdata_reduced->reduce(*rawvar);
		RooRealVar mean_raw("mean_raw","mean_raw",1.,0.9,1.1);
		RooRealVar sig_raw("sig_raw","sig_raw",0.01,0.0002,0.8);
		RooRealVar a1_raw("a1_raw","a1_raw",3,0.05,10);
		RooRealVar a2_raw("a2_raw","a2_raw",3,0.05,10);
		RooRealVar n1_raw("n1_raw","n1_raw",3,1.01,500);
		RooRealVar n2_raw("n2_raw","n2_raw",3,1.01,500);
		RooDoubleCBFast pdfCB_raw("pdfCB_raw","pdfCB_raw",*rawvar,mean_raw,sig_raw,a1_raw,n1_raw,a2_raw,n2_raw);
		cout << "dcb raw" << endl;
		pdfCB_raw.fitTo(*hdata3,Range(0.9,1.1));
        
		h_gaussian_raw->AddBinContent(i+1,mean_raw.getVal());
		h_gaussian_raw->SetBinError(i+1,mean_raw.getError());
        cout << "histogram raw filled" << endl;
        
		//-----------------------------------------------------------------------------------
		cout << "start plotting 1" << endl;
		ccor->cd(i+1);
        
		RooPlot *plot = ecorvar->frame(0.,2, 250);
		plot->SetTitle("Distribution of cor/true "+str);
		plot->GetXaxis()->SetTitle("cor/true");
		gPad->SetGridx();
		gPad->SetGridy();
		hdata_reduced->plotOn(plot,Name("MC_sample_cor"));
		pdfCB_cor.plotOn(plot,Name("corrected"),LineColor(kRed));
		plot->Draw();
		TLegend *legend = new TLegend(0.8,0.8,0.9,0.9);
		legend->SetFillColor(kWhite);
		legend->SetLineColor(kWhite);
		legend->AddEntry("MC_sample_cor","MC Sample cor","P");
		legend->AddEntry("corrected","fit","l");
		legend->Draw();
        
		craw->cd(i+1);
        cout << "start plotting 2" << endl;
		//RooPlot *plot3 = rawvar->frame(0.15,0.4, 250);
		RooPlot *plot3 = rawvar->frame(0.,2, 250);
		plot3->SetTitle("Distribution of raw/true "+str);
		plot3->GetXaxis()->SetTitle("raw/true");
		gPad->SetGridx();
		gPad->SetGridy();
		hdata_reduced->plotOn(plot3,Name("MC_sample_raw"));
		pdfCB_raw.plotOn(plot3,Name("raw"),LineColor(kRed));
		plot3->Draw();
		TLegend *legend2 = new TLegend(0.8,0.8,0.9,0.9);
		legend2->SetFillColor(kWhite);
		legend2->SetLineColor(kWhite);
		legend2->AddEntry("MC_sample_raw","MC Sample raw","P");
		legend2->AddEntry("cfit_raw","fit","l");
		legend2->Draw();
        
	}  
	cout << "Save plots of slices fits" << endl;
	craw->SaveAs("scale_vs_"+arg_variable+"_fit_raw"+add_to_pngFileName);
	ccor->SaveAs("scale_vs_"+arg_variable+"_fit_cor"+add_to_pngFileName);
    
	c_mean->cd();
	gPad->SetGridx();
	gPad->SetGridy();
	gStyle->SetOptStat(0);
    
	hprof_raw->SetMarkerStyle(26);
	hprof_raw->SetMarkerColor(kRed);
	hprof_raw->SetLineColor(kRed);
	hprof_raw->SetName("hprof_raw");
	hprof_raw->GetYaxis()->SetRangeUser(0.96,1.04);
	hprof_raw->GetYaxis()->SetTitle("mean");
	if(arg_variable == "pt") {
		if(dobarrel)
			hprof_raw->GetXaxis()->SetRangeUser(250,5000);
		else
			hprof_raw->GetXaxis()->SetRangeUser(250,2500);
		hprof_raw->GetXaxis()->SetTitle("pt (GeV)");
	}      
	hprof_raw->GetYaxis()->SetTitleOffset(1.4);
	hprof_raw->SetName("hprof_raw");
	hprof_raw->SetFillColor(kMagenta);
	hprof_raw->SetFillStyle(3003);
	hprof_raw->Draw("E3");
	hprof_raw->Draw("hist same p");
    
	hprof_cor->SetMarkerStyle(24);
	hprof_cor->SetMarkerColor(kBlue);
	hprof_cor->SetLineColor(kBlue);
	hprof_cor->SetName("hprof_cor");
	//hprof_cor->GetYaxis()->SetRangeUser(0.96,1.04);
	hprof_cor->SetName("hprof_cor");
	hprof_cor->SetFillColor(kBlue);
	hprof_cor->SetFillStyle(3003);
	hprof_cor->Draw("E3 same");
	hprof_cor->Draw("hist same p");
    
	//h_gaussian_cor->GetYaxis()->SetRangeUser(0.98,1.02);
	//if(arg_variable == "r9")
	//       h_gaussian_cor->GetYaxis()->SetRangeUser(0.96,1.02);
	//h_gaussian_cor->GetXaxis()->SetRangeUser(xmin,xmax);
	//h_gaussian_cor->GetYaxis()->SetTitle("mean");
	//h_gaussian_cor->GetYaxis()->SetTitleOffset(1.55);
	//h_gaussian_cor->GetXaxis()->SetTitle(nameXaxis);
    
	h_gaussian_cor->SetMarkerStyle(8);
	h_gaussian_cor->SetMarkerColor(kBlue);
	h_gaussian_cor->SetLineColor(kBlue);
	h_gaussian_cor->SetName("h_gaussian_cor");
	h_gaussian_cor->Draw("e p same");
    
	h_gaussian_raw->SetMarkerStyle(22);
	h_gaussian_raw->SetMarkerColor(kRed);
	h_gaussian_raw->SetLineColor(kRed);
	h_gaussian_raw->SetName("h_gaussian_raw");
	h_gaussian_raw->Draw("e p same");
    
    
	TLegend *legend = new TLegend(0.11,0.7,0.45,0.9);
	legend->SetFillColor(kWhite);
	legend->SetLineColor(kWhite);  
	legend->AddEntry("h_gaussian_raw","raw/true CB mean","p"); 
	legend->AddEntry("h_gaussian_cor","cor/true CB mean","p");
	legend->AddEntry("hprof_raw","raw/true average","p"); 
	legend->AddEntry("hprof_cor","cor/true average","p");
    
	legend->Draw("same");
	
}

