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
#include "TF1.h"
#include "TStyle.h"
using namespace RooFit;
 
//effsigma function from Chris
Double_t effSigma(TH1 * hist);

void plot_Profile(RooDataSet *hdata, bool dobarrel, RooRealVar *xvar,RooRealVar *ecorvar, RooRealVar *rawvar,TString xLabel, TString yLabel,Double_t plotEMin, Double_t plotEMax,Double_t plotXMin_barrel, Double_t plotXMax_barrel);

void plot_Slices(RooDataSet* hdata,bool dobarrel,TObjArray *Arr,RooRealVar *xvar,RooRealVar *ecorvar, RooRealVar *rawvar,TString name_of_hist_mean,TString name_of_hist_chi2, TString xLabel, TString yLabel,Double_t plotEMin, Double_t plotEMax,Double_t plotXMin, Double_t plotXMax, Double_t plotXMin_endcap=0, Double_t plotXMax_endcap=0);

void plot_Mean(RooDataSet* hdata,bool dobarrel,RooRealVar *xvar,RooRealVar *sigmean, TString xLabel, TString yLabel,Double_t plotMeanMin, Double_t plotMeanMax,Double_t plotXMin_barrel, Double_t plotXMax_barrel);

void eregtestingLog(TString dirname = "Barrel_Log_sig3_alpha2-3_evts15_pileup", bool dobarrel=true, bool doele=false) {
       
       
  	//output dir
  	//TString dirname = "Barrel_Log_sig3_alpha3-1_evts_n100";
  	gSystem->mkdir(dirname,true);
  	gSystem->cd(dirname);    
       //1 = testing, 0 = training
	TCut cutUsed = "(event%15==1)";
	Double_t nmax=500.;
	Double_t mintgtvar = -0.5, nBins1 = 100, nBins2=250,nBins3 =800;
  	
	//read workspace from training
  	TString fname;
  	TString add1, addroot; //extension of png files - first part
  	if (doele && dobarrel) 
  	{
  	   	fname = "wereg_ele_eb.root";
		add1 = "_ele_eb.png";
  	}
  	else if (doele && !dobarrel)
  	{
  	  	fname = "wereg_ele_ee.root";
		add1 = "_ele_ee.png";
  	}
  	else if (!doele && dobarrel) 
 	{	
		fname = "wereg_ph_eb.root";
		add1 = "_ph_eb.png";
              addroot = "_ph_eb.root";
 	}
  	else if (!doele && !dobarrel) 
 	{
 		fname = "wereg_ph_ee.root";
		add1 = "_ph_ee.png";
		addroot = "_ph_ee.root";
 	}
       addroot = ""+addroot;
       add1 = ""+add1; 
	TString add2(add1); //extension of png files - second part
       
	//Workspace 
  	TString infile = TString::Format("./%s",fname.Data());
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

  	TChain *tree;
    tree = new TChain("promptTree");
    tree->Add("/afs/cern.ch/work/k/khoumani/CMSSW_7_2_3/src/HiggsAnalysis/GBRLikelihood/macros/examples/promptTree_FlatPt-300To3000.root");
	tree->Add("/afs/cern.ch/work/k/khoumani/CMSSW_7_2_3/src/HiggsAnalysis/GBRLikelihood/macros/examples/promptTree_FlatPt-5To300_lessstat.root");
 
  	//selection cuts for testing
  	TCut selcut;

  	if (dobarrel) 
	{
    		selcut = "etrue*pt/energy >=250. && etrue*pt/energy<=3000 && abs(scEta)<1.5 && kSaturated[12]!=1 && nvtx>=1";
  	}
  	else
	{
    		selcut = "etrue*pt/energy >=250. && etrue*pt/energy<=3000 && abs(scEta)>1.5 && kSaturated[12]!=1 && nvtx >=1" ;   
  	}       
	TCut selweight = "(1)";
  	TCut prescale10alt = "(event%10==1)";
	TCut prescale20alt = "(event%20==1)";
  	TCut prescale25 = "(event%25==0)";
  	TCut oddevents = "(event%2==1)";
  	TCut prescale100alt = "(event%100==1)";
  	TCut prescale1000alt = "(event%1000==1)";
  	TCut prescale50alt = "(event%50==1)";
  
  	if (doele) 
  	  	weightvar.SetTitle(prescale100alt*selcut);
	else
  	  	weightvar.SetTitle(cutUsed*selweight*selcut);
  
  	//make testing dataset
  	RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",tree,vars,weightvar);
  
  	//retrieve full pdf from workspace
  	RooAbsPdf *sigpdf = ws->pdf("sigpdf");
       //sigpdf->Print();
  	//input variable corresponding to sceta
  	RooRealVar *scetavar = ws->var("var_1");

       //R9
       RooRealVar *r9 = ws->var("var_3");
  	
       //regressed output functions
  	RooAbsReal *sigmeanlim = ws->function("sigmeanlim");
  	RooAbsReal *sigwidthlim = ws->function("sigwidthlim");
  	RooAbsReal *sign1lim = ws->function("sign1lim");
  	RooAbsReal *sign2lim = ws->function("sign2lim");

  	//formula for corrected energy/true energy ( 1.0/(etrue/eraw) * regression mean)
  	RooFormulaVar ecor("ecor","","1./(exp(@0))*exp(@1)",RooArgList(*tgtvar,*sigmeanlim));
  	RooRealVar *ecorvar = (RooRealVar*)hdata->addColumn(ecor);
  	ecorvar->setRange(0.,2.);
  	ecorvar->setBins(800);
  	
  	//formula for raw energy/true energy (1.0/(etrue/eraw))
  	RooFormulaVar raw("raw","","1./exp(@0)",RooArgList(*tgtvar));
  	RooRealVar *rawvar = (RooRealVar*)hdata->addColumn(raw);
  	rawvar->setRange(0.,2.);
  	rawvar->setBins(800);
	
  	//clone data and add regression outputs for plotting
  	RooDataSet *hdataclone = new RooDataSet(*hdata,"hdataclone");
  	RooRealVar *meanvar = (RooRealVar*)hdataclone->addColumn(*sigmeanlim);
  	RooRealVar *widthvar = (RooRealVar*)hdataclone->addColumn(*sigwidthlim);
  	RooRealVar *n1var = (RooRealVar*)hdataclone->addColumn(*sign1lim);
  	RooRealVar *n2var = (RooRealVar*)hdataclone->addColumn(*sign2lim);
       
  
  	//plot target variable and weighted regression prediction (using numerical integration over reduced testing dataset)
    TCanvas *craw = new TCanvas("craw","craw",800,600); //------------------------------------------------------------CANVAS
  	//RooPlot *plot = tgtvar->frame(0.6,1.2,100);
  	RooPlot *plot = tgtvar->frame(mintgtvar,.6, nBins3);
  	//plot->SetTitle("Distribution of log(generated_energy/raw_energy)");
  	plot->SetTitle("");
  	plot->GetXaxis()->SetTitle("log(generated_enery / raw_energy)");
	gPad->SetGridx();
	gPad->SetGridy();
  	hdata->plotOn(plot,Name("MC_sample"));
  	sigpdf->plotOn(plot,ProjWData(*hdata),Name("Signal_PDF"));
  	plot->Draw();
	TLegend *legend = new TLegend(0.11,0.8,0.3,0.89);
	legend->SetFillColor(kWhite);
	legend->SetLineColor(kWhite);
	legend->AddEntry("MC_sample","MC simulation","P");
	legend->AddEntry("Signal_PDF","regression PDF","l");
  	legend->Draw();
	craw->SaveAs("RawE"+add1);
  	craw->SetLogy();
  	plot->SetMinimum(0.1);
  	
	TPad *pad_zoom = new TPad("pad_zoom","pad with a zoom on the peak",0.6,0.6,0.9,0.89);
    
	pad_zoom->Draw();
	pad_zoom->cd();
	RooPlot *plot_zoom;
	if(dobarrel)
		plot_zoom = tgtvar->frame(-0.02,0.02,75);
	else
		plot_zoom = tgtvar->frame(-0.04,0.04,75);
	pad_zoom->SetGridx();
	pad_zoom->SetGridy();
	hdata->plotOn(plot_zoom,Name("MC_sample_zoom"));
  	sigpdf->plotOn(plot_zoom,ProjWData(*hdata),Name("Signal_PDF_zoom"));
  	plot_zoom->SetTitle("");
  	plot_zoom->Draw();
  	craw->Update();
  	craw->SaveAs("RawElog"+add1);
  	craw->SaveAs("RawElog"+addroot); 
/*       TCanvas *craw2 = new TCanvas("craw2","craw2"); //------------------------------------------------------------CANVAS
  	//RooPlot *plot = tgtvar->frame(0.6,1.2,100);
  	plot = tgtvar->frame(mintgtvar,.6, nBins2);
  	plot->SetTitle("Distribution of log(generated_energy/raw_energy)");
  	plot->GetXaxis()->SetTitle("log(generated_enery / raw_energy)");
       gPad->SetGridx();
       gPad->SetGridy();
  	hdata->plotOn(plot,Name("MC_sample2"));
  	sigpdf->plotOn(plot,ProjWData(*hdata),Name("Signal_PDF2"));
  	plot->Draw();
       TLegend *legend2 = new TLegend(0.8,0.8,0.9,0.9);
       legend2->SetFillColor(kWhite);
       legend2->SetLineColor(kWhite);
       legend2->AddEntry("MC_sample2","MC Sample","P");
       legend2->AddEntry("Signal_PDF2","Signal PDF","l");
       legend2->Draw();

  	craw2->SaveAs("RawE2"+add1);
  	craw2->SetLogy();
  	plot->SetMinimum(0.1);
  	craw2->SaveAs("RawElog2"+add1);
  	craw2->SaveAs("RawElog2"+addroot);

       TCanvas *craw3 = new TCanvas("craw3","craw3"); //------------------------------------------------------------CANVAS
  	//RooPlot *plot = tgtvar->frame(-0.5,0.6,nBins3);
  	plot = tgtvar->frame(mintgtvar,.6, nBins3);
  	plot->SetTitle("Distribution of log(generated_energy/raw_energy)");
  	plot->GetXaxis()->SetTitle("log(generated_enery / raw_energy)");
       gPad->SetGridx();
       gPad->SetGridy();
  	hdata->plotOn(plot,Name("MC_sample3"));
  	sigpdf->plotOn(plot,ProjWData(*hdata),Name("Signal_PDF3"));
  	plot->Draw();
       TLegend *legend3 = new TLegend(0.8,0.8,0.9,0.9);
       legend3->SetFillColor(kWhite);
       legend3->SetLineColor(kWhite);
       legend3->AddEntry("MC_sample3","MC Sample","P");
       legend3->AddEntry("Signal_PDF3","Signal PDF","l");
       legend3->Draw();

  	craw3->SaveAs("RawE3"+add1);
  	craw3->SetLogy();
  	plot->SetMinimum(0.1);
  	craw3->SaveAs("RawElog3"+add1);
  	craw3->SaveAs("RawElog3"+addroot);
       

 */ 	//plot distribution of regressed functions over testing dataset
 	TCanvas *cmean = new TCanvas("cmean","cmean"); //---------------------------------------------------------CANVAS
  	RooPlot *plotmean = meanvar->frame(-0.5,0.6,100);
  	plotmean->SetTitle("signal mean distribution");
  	plotmean->GetXaxis()->SetTitle("signal_mean");
  	hdataclone->plotOn(plotmean);
  	plotmean->Draw();
  	cmean->SaveAs("mean"+add1);
  	
  	TCanvas *cwidth = new TCanvas("cwidth","cwidth"); //------------------------------------------------------CANVAS
  	//RooPlot *plotwidth = widthvar->frame(Title("hdata vs width"),0.,0.05,100);
  	RooPlot *plotwidth = widthvar->frame(0.,0.05,100);
  	plotwidth->SetTitle("signal width distribution");
  	plotwidth->GetXaxis()->SetTitle("signal_width");
  	hdataclone->plotOn(plotwidth);
  	plotwidth->Draw();
  	cwidth->SaveAs("width"+add1);
  	
  	TCanvas *cn1 = new TCanvas("cn1","cn1"); //------------------------------------------------------------------CANVAS
  	//RooPlot *plotn1 = nvar->frame(Title("hdata vs nvar"),0.,111.,200);
  	//RooPlot *plotn1 = nvar->frame(0.,111.,200);
       RooPlot *plotn1 = n1var->frame(1.,nmax,1000);
  	plotn1->SetTitle("n1 distribution");
  	plotn1->GetXaxis()->SetTitle("n1");
  	hdataclone->plotOn(plotn1);
  	plotn1->Draw();
  	cn1->SaveAs("n1"+add1);
	
  	TCanvas *cn2 = new TCanvas("cn2","cn2"); //---------------------------------------------------------------CANVAS
  	//RooPlot *plotn2 = n2var->frame(Title("hdata vs n2var"),0.,111.,100);
  	RooPlot *plotn2 = n2var->frame(1.,nmax,1000);
  	plotn2->SetTitle("n2 distribution");
  	plotn2->GetXaxis()->SetTitle("n2");
  	hdataclone->plotOn(plotn2);
  	plotn2->Draw();
  	cn2->SaveAs("n2"+add1);
  	
  	TCanvas *ceta = new TCanvas("ceta","ceta"); //------------------------------------------------------------CANVAS
  	RooPlot *ploteta = scetavar->frame(-2.6,2.6,200);
  	ploteta->SetTitle("eta distribution");
  	ploteta->GetXaxis()->SetTitle("eta");
  	hdataclone->plotOn(ploteta);
  	ploteta->Draw();      
  	ceta->SaveAs("eta"+add1);  
  	
	
  	//create histograms for eraw/etrue and ecor/etrue to quantify regression performance
  	TH1 *heraw = hdata->createHistogram("heraw",*rawvar,Binning(800,0.,2.));
  	TH1 *hecor = hdata->createHistogram("hecor",*ecorvar);
	
  	//heold->SetLineColor(kRed);
  	hecor->SetLineColor(kBlue);
  	heraw->SetLineColor(kMagenta);
  	
  	hecor->GetXaxis()->SetRangeUser(0.6,1.2);
  	//heold->GetXaxis()->SetRangeUser(0.6,1.2);
	hecor->SetTitle("magenta: raw/true - blue: mean*raw/true");  	
  	TCanvas *cresponse = new TCanvas("cresponse","cresponse"); //---------------------------------------------CANVAS
  	gPad->SetGridx();
	gPad->SetGridy();
  	hecor->Draw("HIST");
  	//heold->Draw("HISTSAME");
  	heraw->Draw("HISTSAME");
  	cresponse->SaveAs("response"+add1);
  	cresponse->SetLogy();
  	cresponse->SaveAs("responselog"+add1);
       //cresponse->SetLogy(kFALSE);
	cresponse->SaveAs("responselog"+addroot);
  	
  	ofstream file_effsigma;
	file_effsigma.open ("../outputs/effectiveSigma.txt", ios::out | ios::app);
	if (cutUsed == "(event%2==1)")
		file_effsigma << dirname << "_testing,";
	else if (cutUsed == "(event%2==0)")
		file_effsigma << dirname << "_training,";
	
  	printf("make fine histogram\n");
  	TH1 *hecorfine = hdata->createHistogram("hecorfine",*ecorvar,Binning(20e3,0.,2.));
	TH1 *hrawfine = hdata->createHistogram("hrawfine",*rawvar,Binning(20e3,0.,2.));
  	//printf("calc effsigma\n");
  	double effsigmaCor = effSigma(hecorfine);
  	double effsigmaRaw = effSigma(hrawfine);	
  	//printf("effsigmaCor = %5f - effsigmaRaw = %5f\n",effsigmaCor,effsigmaRaw);
  	file_effsigma << effsigmaRaw << "," << effsigmaCor << "\n";
	file_effsigma.close();
    
	/*    Double_t xx=0.2;
       TH1 *h1 = hdata->createHistogram("tgtvar",*tgtvar,Binning(10,-xx,xx));
       RooAbsReal *signal = (RooAbsReal*)sigpdf;
       TH1 *h2 = signal->createHistogram("tgtvar",*tgtvar,Binning(10,-xx,xx));   
       Int_t n=10;
       Double_t res[n], x[n];
       for (Int_t i=0; i<n; i++) x[i]= (Double_t)i/n-xx;
       Double_t p=h1->Chi2Test(h2,"UU P CHI2/NDF",res);
       std::cout << "chi2/ndf "<< p << std::endl;
       //h2->Draw();
       
       TGraph *resgr = new TGraph(n,x,res);
       resgr->GetXaxis()->SetRangeUser(-xx,xx);
       resgr->GetYaxis()->SetRangeUser(-.9,.9);
       resgr->GetYaxis()->SetTitle("Normalized Residuals");
       resgr->SetMarkerStyle(21);
       resgr->SetMarkerColor(2);
       resgr->SetMarkerSize(.9);
       resgr->SetTitle("Normalized Residuals");
       resgr->Draw("APL");
*/
       //h2->SetMarkerColor(kBlue);
       //h2->SetMarkerStyle(8);
       //h2->Scale(1./h2->Integral());
       //h2->Draw("hist p");
       //h1->Scale(1./h1->Integral());
       //h1->SetMarkerColor(kRed);
       //h1->SetMarkerStyle(8);
       //h1->Draw("same hist p");
  	
	//new TCanvas;
  	//RooPlot *ploteold = testvar.frame(0.6,1.2,100);
  	//hdatasigtest->plotOn(ploteold);
  	//ploteold->Draw();    
  	
  	//new TCanvas;
  	//RooPlot *plotecor = ecorvar->frame(0.6,1.2,100);
  	//hdatasig->plotOn(plotecor);
  	//plotecor->Draw();  
  	
 	//________________________________________________________________________________________________________________________
 	//________________________________________________________________________________________________________________________
 	//________________________________________________________________________________________________________________________
 	//
	RooRealVar *scphi = ws->var("var_2");
	//RooRealVar *r9 = ws->var("var_3");
	RooRealVar *scetawidth = ws->var("var_4");
	RooRealVar *scphiwidth = ws->var("var_5");
	RooRealVar *scnclusters = ws->var("var_6");
	RooRealVar *hoveretower = ws->var("var_7");
 	RooRealVar *rho = ws->var("var_8");
	RooRealVar *nvtx = ws->var("var_9");

  	RooFormulaVar meanFormula("mean","","exp(@0)*@1",RooArgList(*sigmeanlim,RooConst(1.)));
  	RooRealVar *sigmean = (RooRealVar*)hdata->addColumn(meanFormula);
	sigmean->setRange(0.,2.); //

	RooRealVar *pt;
	RooRealVar *rawEnergy = ws->var("var_0");

	if(!dobarrel) {
		RooRealVar *preshowerEnergy = ws->var("var_18");
		//formula for Pt true_energy/Cosh(eta) = tgtvar*rawvar/cosh(eta) 
		RooFormulaVar ptFormula("pt","","exp(@0)*@1*(1+@2)*1./(cosh(@3))",RooArgList(*tgtvar,*rawEnergy,*preshowerEnergy,*scetavar));
		pt = (RooRealVar*)hdata->addColumn(ptFormula);
		pt->setRange(0.,2500); //
	}
	else {
		RooFormulaVar ptFormula("pt","","exp(@0)*@1*1./(cosh(@2))",RooArgList(*tgtvar,*rawEnergy,*scetavar));
		pt = (RooRealVar*)hdata->addColumn(ptFormula);
		pt->setRange(0.,5000); //
	}

	//Set ranges for variables
	scetavar->setRange(-2.5,2.5);
	scphi->setRange(-4.,4.);
	r9->setRange(0.7,1.);
	scetawidth->setRange(0.,0.025);
	scphiwidth->setRange(0.,0.2);
	scnclusters->setRange(0.,20.);
	hoveretower->setRange(0.,0.2);
	rho->setRange(0.,30.);
	nvtx->setRange(0.,50.);
	
	//Set Binning
	hoveretower->setBins(50.);
    
	//Define Canvas properties
	Double_t canvas_width(1300);
	Double_t canvas_height(800);
	
//create canvas for ecor (corrected energy)
    
	Double_t plotEMin(0.9), plotEMax(1.1);
	gStyle->SetOptStat(0);
    
    
	TCanvas *canvas_ecor = new TCanvas("corrected_energy","ecor",canvas_width,canvas_height);
    
	canvas_ecor->Divide(3,3);
    
	canvas_ecor->cd(1);
	plot_Profile(hdata, dobarrel,scetavar,ecorvar, rawvar,"eta", "energy/generated_energy",plotEMin, plotEMax,-1.5, 1.5);
	
	canvas_ecor->cd(2);
	plot_Profile(hdata, dobarrel,scphi,ecorvar, rawvar,"phi", "energy/generated_energy",plotEMin, plotEMax,-2., 2.);
	
	canvas_ecor->cd(3);
	plot_Profile(hdata, dobarrel,r9,ecorvar, rawvar,"R9", "energy/generated_energy",plotEMin, plotEMax,0.7,1.);
	
	canvas_ecor->cd(4);
	plot_Profile(hdata, dobarrel,scetawidth,ecorvar, rawvar,"scEtaWidth", "energy/generated_energy",plotEMin, plotEMax,0.005,0.014);
	
	canvas_ecor->cd(5);
	plot_Profile(hdata, dobarrel,scphiwidth,ecorvar, rawvar,"scPhiWidth", "energy/generated_energy",plotEMin, plotEMax,0.,0.12);
	
	canvas_ecor->cd(6);
	if(dobarrel) 
		plot_Profile(hdata, dobarrel,pt,ecorvar, rawvar,"Pt", "energy/generated_energy",plotEMin, plotEMax,0.,5000);
	else
		plot_Profile(hdata, dobarrel,pt,ecorvar, rawvar,"Pt", "energy/generated_energy",plotEMin, plotEMax,0.,2500);
	canvas_ecor->cd(7);
	plot_Profile(hdata, dobarrel,hoveretower,ecorvar, rawvar,"hoveretower", "energy/generated_energy",0.7, plotEMax,0.,0.16);
	
	canvas_ecor->cd(8);
	plot_Profile(hdata, dobarrel,rho,ecorvar, rawvar,"rho", "energy/generated_energy",plotEMin, plotEMax,0.,30);
	
	canvas_ecor->cd(9);
	plot_Profile(hdata, dobarrel,nvtx,ecorvar, rawvar,"number of vertices", "energy/generated_energy",plotEMin, plotEMax,0.,40);

	//save canvas_ecor
	canvas_ecor->SaveAs("corrected_energy"+add2);
	canvas_ecor->SaveAs("corrected_energy"+addroot);
	/*
	  TCanvas *canvas_ecor2 = new TCanvas("corrected_energy2","ecor2",650,400);
	  plot_Profile(hdata, dobarrel,hoveretower,ecorvar, rawvar,"hoveretower", "energy/generated_energy",0.5, 1.2,0.,0.16);
	  canvas_ecor2->SaveAs("corrected_energy2"+add2);
	  canvas_ecor2->SaveAs("corrected_energy2"+addroot);*/


//create canvas for mean_of_gaussian

/*       TObjArray* Arr = new TObjArray();
	TCanvas *canvas_ecor_fit = new TCanvas("ecor_fit","ecor_fit",canvas_width,canvas_height);
	canvas_ecor_fit->Divide(3,3);

	canvas_ecor_fit->cd(1);       
       plot_Slices(hdata,dobarrel,Arr,scetavar,ecorvar, rawvar,"hdata_hist_00000012_1","hdata_hist_00000012_chi2","eta","energy/generated_energy",plotEMin, plotEMax,-1.5,1.5);
       canvas_ecor_fit->cd(2);   
       plot_Slices(hdata,dobarrel,Arr,scphi,ecorvar, rawvar,"hdata_hist_00000014_1","hdata_hist_00000014_chi2","phi","energy/generated_energy",plotEMin, plotEMax,-2,2);
       canvas_ecor_fit->cd(3);       
       plot_Slices(hdata,dobarrel,Arr,r9,ecorvar, rawvar,"hdata_hist_00000016_1","hdata_hist_00000016_chi2","R9","energy/generated_energy",plotEMin, plotEMax,0.7,1);
       canvas_ecor_fit->cd(4);       
       plot_Slices(hdata,dobarrel,Arr,scetawidth,ecorvar, rawvar,"hdata_hist_00000018_1","hdata_hist_00000018_chi2","scEtaWidth","energy/generated_energy",plotEMin, plotEMax,0.005,0.014,0.004,0.04);
	canvas_ecor_fit->cd(5);       
       plot_Slices(hdata,dobarrel,Arr,scphiwidth,ecorvar, rawvar,"hdata_hist_0000001a_1","hdata_hist_0000001a_chi2","scPhiWidth","energy/generated_energy",plotEMin, plotEMax,0.,0.12,0.,0.12);
	canvas_ecor_fit->cd(6);       
       plot_Slices(hdata,dobarrel,Arr,pt,ecorvar, rawvar,"hdata_hist_0000001c_1","hdata_hist_0000001c_chi2","Pt","energy/generated_energy",plotEMin, plotEMax,0.,5000,0.,5000);
	canvas_ecor_fit->cd(7);       
       plot_Slices(hdata,dobarrel,Arr,hoveretower,ecorvar, rawvar,"hdata_hist_0000001e_1","hdata_hist_0000001e_chi2","hoveretower","energy/generated_energy",0.7, plotEMax,0.,0.16,0.,0.16);
	canvas_ecor_fit->cd(8);       
       plot_Slices(hdata,dobarrel,Arr,rho,ecorvar, rawvar,"hdata_hist_00000020_1","hdata_hist_00000020_chi2","rho","energy/generated_energy",plotEMin, plotEMax,0.,30,0.,30);
	canvas_ecor_fit->cd(9);       
       plot_Slices(hdata,dobarrel,Arr,nvtx,ecorvar, rawvar,"hdata_hist_00000022_1","hdata_hist_00000022_chi2","numberOfVertices","energy/generated_energy",plotEMin, plotEMax,0.,40,0.,40);


//save canvas_rawvar
	canvas_ecor_fit->SaveAs("mean_of_gaussian_energy"+add2);
       canvas_ecor_fit->SaveAs("mean_of_gaussian"+addroot);
*//*
//create canvas for mean (#mu)
       Double_t plotMeanMin(0.9), plotMeanMax(1.1);

	TCanvas *canvas_sigmean = new TCanvas("signal_mean","sigmean",canvas_width,canvas_height);
       gPad->SetGridx();
       gPad->SetGridy();
	canvas_sigmean->Divide(3,3);
	canvas_sigmean->cd(1);
	plot_Mean(hdata,dobarrel,scetavar,sigmean, "eta", "Mu",plotMeanMin, plotMeanMax,-1.5, 1.5);
	canvas_sigmean->cd(2);
	plot_Mean(hdata,dobarrel,scphi,sigmean, "phi", "Mu",plotMeanMin, plotMeanMax,-2,2);
       canvas_sigmean->cd(3);
	plot_Mean(hdata,dobarrel,r9,sigmean, "R9", "Mu",plotMeanMin, plotMeanMax,0.7,1);
	canvas_sigmean->cd(4);
	plot_Mean(hdata,dobarrel,scetawidth,sigmean,"scEtaWidth", "Mu",plotMeanMin, plotMeanMax,0.005,0.014);
	canvas_sigmean->cd(5);
	plot_Mean(hdata,dobarrel,scphiwidth,sigmean, "scPhiWidth", "Mu",plotMeanMin, plotMeanMax,0.,0.12);
	canvas_sigmean->cd(6);
	plot_Mean(hdata,dobarrel,pt,sigmean, "pt", "Mu",plotMeanMin, plotMeanMax,0., 5000);
	canvas_sigmean->cd(7);
	plot_Mean(hdata,dobarrel,hoveretower,sigmean, "hoveretower", "Mu",plotMeanMin, 1.5,0.,0.16);
	canvas_sigmean->cd(8);
	plot_Mean(hdata,dobarrel,rho,sigmean, "rho", "Mu",plotMeanMin, plotMeanMax,0.,30.);
	canvas_sigmean->cd(9);
	plot_Mean(hdata,dobarrel,nvtx,sigmean, "nvtx", "Mu",plotMeanMin, plotMeanMax,0.,40.);


//save canvas_sigmean
	canvas_sigmean->SaveAs("signal_mean"+add2);
       canvas_sigmean->SaveAs("signal_mean"+addroot);*/
}

Double_t effSigma(TH1 * hist)
{

	TAxis *xaxis = hist->GetXaxis();
  	Int_t nb = xaxis->GetNbins();
  	if(nb < 10) 
	{
    		cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    		return 0.;
 	}
  
  	Double_t bwid = xaxis->GetBinWidth(1);
  	if(bwid == 0) 
	{
    		cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    		return 0.;
  	}
  	Double_t xmax = xaxis->GetXmax();
  	Double_t xmin = xaxis->GetXmin();
  	Double_t ave = hist->GetMean();
  	Double_t rms = hist->GetRMS();

  	Double_t total=0.;
  	for(Int_t i=0; i<nb+2; i++) 
	{
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
  	for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) 
	{ // Scan window centre
    		Int_t ibm=(ave-xmin)/bwid+1+iscan;
    		Double_t x=(ibm-0.5)*bwid+xmin;
    		Double_t xj=x;
    		Double_t xk=x;
    		Int_t jbm=ibm;
    		Int_t kbm=ibm;
    		Double_t bin=hist->GetBinContent(ibm);
    		total=bin;
    		for(Int_t j=1;j<nb;j++)
		{
    	  		if(jbm < nb) 
			{
    	  		  	jbm++;
    	    			xj+=bwid;
    	    			bin=hist->GetBinContent(jbm);
    	    			total+=bin;
    	    			if(total > rlim) break;
    	  		}
			else ierr=1;
      			if(kbm > 0) 
			{
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
    		if(wid < widmin) 
		{
    		  	widmin=wid;
    		  	ismin=iscan;
    		}   
  	}
  	if(ismin == nrms || ismin == -nrms) ierr=3;
  	if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
  
  	return widmin;
  
}

void plot_Profile(RooDataSet* hdata,bool dobarrel,RooRealVar *xvar,RooRealVar *ecorvar, RooRealVar *rawvar,TString xLabel, TString yLabel,Double_t plotEMin, Double_t plotEMax,Double_t plotXMin_barrel, Double_t plotXMax_barrel) {
	gPad->SetGridx();
	gPad->SetGridy();
	TH2* histo2D_ecor = hdata->createHistogram(*xvar,*ecorvar);	
	TH2* histo2D_raw = hdata->createHistogram(*xvar,*rawvar);
	TProfile *hprof_ecor = histo2D_ecor->ProfileX();
	TProfile *hprof_raw = histo2D_raw->ProfileX();
	   
	hprof_ecor->SetErrorOption("s");
	hprof_raw->SetErrorOption("s");

	hprof_raw->GetYaxis()->SetRangeUser(plotEMin,plotEMax);
	if(dobarrel) 
		hprof_raw->GetXaxis()->SetRangeUser(plotXMin_barrel,plotXMax_barrel);
	
	hprof_raw->GetXaxis()->SetTitle(xLabel);
	hprof_raw->GetYaxis()->SetTitle(yLabel);
	hprof_ecor->SetMarkerStyle(7);
	hprof_ecor->SetMarkerColor(kBlue);
	hprof_raw->SetMarkerColor(kRed);
	hprof_raw->SetMarkerStyle(7);       
	hprof_raw->SetFillColor(kMagenta);
	hprof_raw->SetFillStyle(3003);
	hprof_raw->Draw("E3");
	hprof_raw->Draw("hist same p");
	hprof_ecor->Draw("same");
}

void plot_Slices(RooDataSet* hdata,bool dobarrel,TObjArray *Arr,RooRealVar *xvar,RooRealVar *ecorvar, RooRealVar *rawvar, TString name_of_hist_mean,TString name_of_hist_chi2,TString xLabel, TString yLabel,Double_t plotEMin, Double_t plotEMax,Double_t plotXMin, Double_t plotXMax, Double_t plotXMin_endcap, Double_t plotXMax_endcap) {
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* h2_ecor = hdata->createHistogram(*xvar,*ecorvar);
	Double_t m1(0.98), m2(1.02);
       TString option_fit="RG5"; 
       TF1 *func_gaus = new TF1("func_gaus","gaus(0)",m1,m2);
       func_gaus->SetParameters(1,1,1);
	func_gaus->SetRange(m1,m2);
       h2_ecor->FitSlicesY(func_gaus,0,-1,0,option_fit,Arr);      
       //h2_ecor->FitSlicesY(0,0,-1,0,"G5Q",Arr);
       //Arr->Print();
       TH1D *h2_ecor_1 = (TH1D*)gDirectory->Get(name_of_hist_mean);
       h2_ecor_1->SetMarkerStyle(7);
       h2_ecor_1->SetMarkerColor(kBlue);
       
       TH1D *h2_ecor_chi2 = (TH1D*)gDirectory->Get(name_of_hist_chi2);
       h2_ecor_chi2->SetMarkerStyle(7);
       h2_ecor_chi2->SetMarkerColor(kGreen);
       
       TH2* histo2D_ecor = hdata->createHistogram(*xvar,*ecorvar);	
	TProfile *hprof_ecor = histo2D_ecor->ProfileX();
	hprof_ecor->GetYaxis()->SetRangeUser(plotEMin,plotEMax);
       if(dobarrel) 
              hprof_ecor->GetXaxis()->SetRangeUser(plotXMin,plotXMax);
       else if(!dobarrel && plotXMin_endcap!=0 && plotXMin_endcap!=0)
              hprof_ecor->GetXaxis()->SetRangeUser(plotXMin_endcap,plotXMax_endcap);
	hprof_ecor->GetXaxis()->SetTitle(xLabel);
	hprof_ecor->GetYaxis()->SetTitle(yLabel);
       hprof_ecor->SetMarkerStyle(7);
       hprof_ecor->SetMarkerColor(kRed);     
       hprof_ecor->SetFillColor(kMagenta);
       hprof_ecor->SetFillStyle(3003);
       hprof_ecor->Draw("E3");
       hprof_ecor->Draw("hist same p");
       h2_ecor_1->Draw("same");     
       //h2_ecor_chi2->Draw("hist same p");  
}

void plot_Mean(RooDataSet* hdata,bool dobarrel,RooRealVar *xvar,RooRealVar *sigmean, TString xLabel, TString yLabel,Double_t plotMeanMin, Double_t plotMeanMax,Double_t plotXMin_barrel, Double_t plotXMax_barrel) {
       gPad->SetGridx();
       gPad->SetGridy();
       
	TH2* histo2_sigmean = hdata->createHistogram(*xvar,*sigmean);	
	TProfile *hprof_sigmean_sceta = histo2_sigmean->ProfileX();
	hprof_sigmean_sceta->SetErrorOption("s");
	hprof_sigmean_sceta->GetYaxis()->SetRangeUser(plotMeanMin,plotMeanMax);
	if(dobarrel) hprof_sigmean_sceta->GetXaxis()->SetRangeUser(plotXMin_barrel,plotXMax_barrel);
	hprof_sigmean_sceta->GetXaxis()->SetTitle(xLabel);
	hprof_sigmean_sceta->GetYaxis()->SetTitle(yLabel);
	hprof_sigmean_sceta->Draw();
}


