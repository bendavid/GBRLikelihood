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

using namespace RooFit;
 
//effsigma function from Chris
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

void eregtesting(bool dobarrel=true, bool doele=false) {
  
  	//output dir
  	TString dirname = "ereg_sig3_EvenOddEvts_alpha2-3/";
  	gSystem->mkdir(dirname,true);
  	gSystem->cd(dirname); 
  	
       /*TCut prescale10alt = "(event%10==1)";
       TCut prescale20alt = "(event%20==1)";
  	TCut oddevents = "(event%2==1)";
  	TCut prescale100alt = "(event%100==1)";
  	TCut prescale1000alt = "(event%1000==1)";
  	TCut prescale50alt = "(event%50==1)";*/
       TCut cutUsed = "(event%2==1)";
       Double_t mintgtvar = 0.6, nBins1 = 100, nBins2=250,nBins3 =800;
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

  	TTree *dtree;
  
  	if (doele) {
              TFile *fdin = TFile::Open("");
    	       TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("");
    	       dtree = (TTree*)ddir->Get("");       
  	}
	else {
		TFile *fTree = TFile::Open("/afs/cern.ch/user/m/musella/public/forKenza/gam_gam_phys14_v5_regtraining_v3.root"); 
		dtree = (TTree*)fTree->Get("promptTree");   
	}
  
  	//selection cuts for testing
  	TCut selcut;
  	if (dobarrel) 
    		selcut = "pt>250. && abs(scEta)<=1.5 && kSaturated[12]!=1"; 
  	else
    		selcut = "pt>250. && abs(scEta)>1.5 && kSaturated[12]!=1"; 
       
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
  	  	weightvar.SetTitle(cutUsed*selcut);
  
  	//make testing dataset
  	RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",dtree,vars,weightvar);
  


/*  	//CreateDataSet hdatasmall
  	if (doele) 
    		weightvar.SetTitle(prescale1000alt*selcut);
  	else
    		weightvar.SetTitle(cutUsed*selcut);
  	//make reduced testing dataset for integration over conditional variables
  	RooDataSet *hdatasmall = RooTreeConvert::CreateDataSet("hdatasmall",dtree,vars,weightvar);     
*/   
  	//retrieve full pdf from workspace
  	RooAbsPdf *sigpdf = ws->pdf("sigpdf");
  
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
  	RooRealVar *n1var = (RooRealVar*)hdataclone->addColumn(*sign1lim);
  	RooRealVar *n2var = (RooRealVar*)hdataclone->addColumn(*sign2lim);
  
  
  	//plot target variable and weighted regression prediction (using numerical integration over reduced testing dataset)
/*TCanvas *craw = new TCanvas("craw","craw"); //------------------------------------------------------------CANVAS
  	//RooPlot *plot = tgtvar->frame(0.6,1.2,100);
  	RooPlot *plot = tgtvar->frame(mintgtvar,2., nBins1);
  	plot->SetTitle("Distribution of generated_energy/raw_energy");
  	plot->GetXaxis()->SetTitle("generated_enery / raw_energy");
       gPad->SetGridx();
       gPad->SetGridy();
  	hdata->plotOn(plot,Name("MC_sample"));
  	sigpdf->plotOn(plot,ProjWData(*hdata),Name("Signal_PDF"));
  	plot->Draw();
       TLegend *legend = new TLegend(0.6,0.6,0.9,0.9);
       legend->SetFillColor(kWhite);
       legend->SetLineColor(kWhite);
       legend->AddEntry("MC_sample","MC Sample","P");
       legend->AddEntry("Signal_PDF","Signal PDF","l");
       legend->Draw();

  	craw->SaveAs("RawE"+add1);
  	craw->SetLogy();
  	plot->SetMinimum(0.1);
  	craw->SaveAs("RawElog"+add1);
  	craw->SaveAs("RawElog"+addroot); 

       TCanvas *craw2 = new TCanvas("craw2","craw2"); //------------------------------------------------------------CANVAS
  	//RooPlot *plot = tgtvar->frame(0.6,1.2,100);
  	plot = tgtvar->frame(mintgtvar,2., nBins2);
  	plot->SetTitle("Distribution of generated_energy/raw_energy");
  	plot->GetXaxis()->SetTitle("generated_enery / raw_energy");
       gPad->SetGridx();
       gPad->SetGridy();
  	hdata->plotOn(plot,Name("MC_sample2"));
  	sigpdf->plotOn(plot,ProjWData(*hdata),Name("Signal_PDF2"));
  	plot->Draw();
       TLegend *legend2 = new TLegend(0.6,0.6,0.9,0.9);
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
  	//RooPlot *plot = tgtvar->frame(0.6,1.2,100);
  	plot = tgtvar->frame(mintgtvar,2., nBins3);
  	plot->SetTitle("Distribution of generated_energy/raw_energy");
  	plot->GetXaxis()->SetTitle("generated_enery / raw_energy");
       gPad->SetGridx();
       gPad->SetGridy();
  	hdata->plotOn(plot,Name("MC_sample3"));
  	sigpdf->plotOn(plot,ProjWData(*hdata),Name("Signal_PDF3"));
  	plot->Draw();
       TLegend *legend3 = new TLegend(0.6,0.6,0.9,0.9);
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
       

  	//plot distribution of regressed functions over testing dataset
  	TCanvas *cmean = new TCanvas("cmean","cmean"); //---------------------------------------------------------CANVAS
  	RooPlot *plotmean = meanvar->frame(0.8,2.0,100);
  	plotmean->SetTitle("Data events vs signal_mean");
  	plotmean->GetXaxis()->SetTitle("signal_mean");
  	hdataclone->plotOn(plotmean);
  	plotmean->Draw();
  	cmean->SaveAs("mean"+add1);
  	
  	TCanvas *cwidth = new TCanvas("cwidth","cwidth"); //------------------------------------------------------CANVAS
  	//RooPlot *plotwidth = widthvar->frame(Title("hdata vs width"),0.,0.05,100);
  	RooPlot *plotwidth = widthvar->frame(0.,0.05,100);
  	plotwidth->SetTitle("Data events vs signal_width");
  	plotwidth->GetXaxis()->SetTitle("signal_width");
  	hdataclone->plotOn(plotwidth);
  	plotwidth->Draw();
  	cwidth->SaveAs("width"+add1);
  	
  	TCanvas *cn1 = new TCanvas("cn1","cn1"); //------------------------------------------------------------------CANVAS
  	//RooPlot *plotn1 = nvar->frame(Title("hdata vs nvar"),0.,111.,200);
  	//RooPlot *plotn1 = nvar->frame(0.,111.,200);
       RooPlot *plotn1 = n1var->frame(0.,100,800);
  	plotn1->SetTitle("Data events vs n1_pow");
  	plotn1->GetXaxis()->SetTitle("n1_pow");
  	hdataclone->plotOn(plotn1);
  	plotn1->Draw();
  	cn1->SaveAs("n1"+add1);
	
  	TCanvas *cn2 = new TCanvas("cn2","cn2"); //---------------------------------------------------------------CANVAS
  	//RooPlot *plotn2 = n2var->frame(Title("hdata vs n2var"),0.,111.,100);
  	RooPlot *plotn2 = n2var->frame(0.,100.,400);
  	plotn2->SetTitle("Data events vs n2_pow");
  	plotn2->GetXaxis()->SetTitle("n2_pow");
  	hdataclone->plotOn(plotn2);
  	plotn2->Draw();
  	cn2->SaveAs("n2"+add1);
  	
  	TCanvas *ceta = new TCanvas("ceta","ceta"); //------------------------------------------------------------CANVAS
  	RooPlot *ploteta = scetavar->frame(-2.6,2.6,200);
  	ploteta->SetTitle("Data events vs eta");
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
       hecor->GetXaxis()->SetTitle("magenta: raw/true - blue: mean*raw/true");  	
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
  	
  	printf("make fine histogram\n");
  	TH1 *hecorfine = hdata->createHistogram("hecorfine",*ecorvar,Binning(20e3,0.,2.));
	TH1 *hrawfine = hdata->createHistogram("hrawfine",*rawvar,Binning(20e3,0.,2.));
  	printf("calc effsigma\n");
  	
  	double effsigmaCor = effSigma(hecorfine);
  	double effsigmaRaw = effSigma(hrawfine);	
  	printf("effsigmaCor = %5f - effsigmaRaw = %5f\n",effsigmaCor,effsigmaRaw);*/
  	
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
	RooRealVar *sigmean = (RooRealVar*)hdata->addColumn(*sigmeanlim);

       RooRealVar *pt;
       RooRealVar *rawEnergy = ws->var("var_0");

       if(!dobarrel) {
              RooRealVar *preshowerEnergy = ws->var("var_28");
	       //formula for Pt corrected_energy/Cosh(eta) = ecorvar*tgtvar*rawvar/cosh(eta) 
  	       RooFormulaVar ptFormula("pt","","@0*@1*(1+@2)*1./(cosh(@3))",RooArgList(*sigmeanlim,*rawEnergy,*preshowerEnergy,*scetavar));
     	       pt = (RooRealVar*)hdata->addColumn(ptFormula);
  	       pt->setRange(0.,2500); //
       }
       else {
              RooFormulaVar ptFormula("pt","","@0*@1*1./(cosh(@2))",RooArgList(*sigmeanlim,*rawEnergy,*scetavar));
     	       pt = (RooRealVar*)hdata->addColumn(ptFormula);
  	       pt->setRange(0.,5000); //
       }

	//Set ranges for variables
	sigmean->setRange(0.,2.);
		//------------------
	scetavar->setRange(-2.5,2.5);
	scphi->setRange(-4.,4.);
	r9->setRange(0.65,1.);
	scetawidth->setRange(0.,0.025);
	scphiwidth->setRange(0.,0.2);
	scnclusters->setRange(0.,20.);
	hoveretower->setRange(0.,0.2);
	rho->setRange(0.,40.);
	nvtx->setRange(0.,50.);
	
       //Set Binning
       hoveretower->setBins(50.);
       
	//Define Canvas properties
	Double_t canvas_width(1300);
	Double_t canvas_height(800);
	
//create canvas for ecor (corrected energy)
       Double_t plotERMin(0.6), plotERMax(1.2);
       Double_t plotECMin(0.9), plotECMax(1.1);

	TCanvas *canvas_ecor = new TCanvas("corrected_energy","ecor",canvas_width,canvas_height);
       
	canvas_ecor->Divide(2,2);
       
	//________________
	//  ecor_sceta
	//________________
	canvas_ecor->cd(1);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_ecor_sceta = hdata->createHistogram(*scetavar,*ecorvar);	
	TProfile *hprof_ecor_sceta = histo2_ecor_sceta->ProfileX();
	hprof_ecor_sceta->SetErrorOption("s");
	hprof_ecor_sceta->GetYaxis()->SetRangeUser(plotECMin,plotECMax);
	hprof_ecor_sceta->GetXaxis()->SetTitle("eta");
	hprof_ecor_sceta->GetYaxis()->SetTitle("ecor");
	hprof_ecor_sceta->Draw();

	/*//________________
	//  ecor_scphi
	//________________
	//TCanvas *canvas_ecor_scphi = new TCanvas("ecor_scphi","ecor_scphi");
	canvas_ecor->cd(2);	
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_ecor_scphi = hdata->createHistogram(*scphi,*ecorvar);	
	TProfile *hprof_ecor_scphi = histo2_ecor_scphi->ProfileX();
	hprof_ecor_scphi->SetErrorOption("s");
	hprof_ecor_scphi->GetYaxis()->SetRangeUser(plotECMin,plotECMax);
	hprof_ecor_scphi->GetXaxis()->SetTitle("phi");
	hprof_ecor_scphi->GetYaxis()->SetTitle("ecor");
	hprof_ecor_scphi->Draw();
	//canvas_ecor_scphi->SaveAs("ecor_scphi.png");*/
	
	//________________
	//  ecor_r9
	//________________
	//TCanvas *canvas_ecor_r9 = new TCanvas("ecor_r9","ecor_r9");	
	canvas_ecor->cd(2);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_ecor_r9 = hdata->createHistogram(*r9,*ecorvar);	
	TProfile *hprof_ecor_r9 = histo2_ecor_r9->ProfileX();
	hprof_ecor_r9->SetErrorOption("s");
	hprof_ecor_r9->GetYaxis()->SetRangeUser(plotECMin,plotECMax);
	hprof_ecor_r9->GetXaxis()->SetTitle("R9");
	hprof_ecor_r9->GetYaxis()->SetTitle("ecor");
	hprof_ecor_r9->Draw();
	//canvas_ecor_r9->SaveAs("ecor_r9.png");

/*	//________________
	//  ecor_scetawidth
	//________________
	//TCanvas *canvas_ecor_scetawidth = new TCanvas("ecor_scetawidth","ecor_scetawidth");	
	canvas_ecor->cd(4);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_ecor_scetawidth = hdata->createHistogram(*scetawidth,*ecorvar);	
	TProfile *hprof_ecor_scetawidth = histo2_ecor_scetawidth->ProfileX();
	hprof_ecor_scetawidth->SetErrorOption("s");
	hprof_ecor_scetawidth->GetYaxis()->SetRangeUser(plotECMin,plotECMax);
	hprof_ecor_scetawidth->GetXaxis()->SetTitle("scetawidth");
	hprof_ecor_scetawidth->GetYaxis()->SetTitle("ecor");
	hprof_ecor_scetawidth->Draw();
	//canvas_ecor_scetawidth->SaveAs("ecor_scetawidth.png");

	//________________
	//  ecor_scphiwidth
	//________________
	//TCanvas *canvas_ecor_scphiwidth = new TCanvas("ecor_scphiwidth","ecor_scphiwidth");	
	canvas_ecor->cd(5);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_ecor_scphiwidth = hdata->createHistogram(*scphiwidth,*ecorvar);	
	TProfile *hprof_ecor_scphiwidth = histo2_ecor_scphiwidth->ProfileX();
	hprof_ecor_scphiwidth->SetErrorOption("s");
	hprof_ecor_scphiwidth->GetYaxis()->SetRangeUser(plotECMin,plotECMax);
	hprof_ecor_scphiwidth->GetXaxis()->SetTitle("scphiwidth");
	hprof_ecor_scphiwidth->GetYaxis()->SetTitle("ecor");
	hprof_ecor_scphiwidth->Draw();
	//canvas_ecor_scphiwidth->SaveAs("ecor_scphiwidth.png");
*/
	//________________
	//  ecor_Pt
	//________________
	//TCanvas *canvas_ecor_pt = new TCanvas("ecor_pt","ecor_pt");	
	canvas_ecor->cd(3);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_ecor_pt = hdata->createHistogram(*pt,*ecorvar);	
	TProfile *hprof_ecor_pt = histo2_ecor_pt->ProfileX();
	hprof_ecor_pt->SetErrorOption("s");
	hprof_ecor_pt->GetYaxis()->SetRangeUser(plotECMin,plotECMax);
	hprof_ecor_pt->GetXaxis()->SetTitle("pt");
	hprof_ecor_pt->GetYaxis()->SetTitle("ecor");
	hprof_ecor_pt->Draw();
	//canvas_ecor_pt->SaveAs("ecor_pt.png");

/*	//________________
	//  ecor_hoveretower
	//________________
	//TCanvas *canvas_ecor_hoveretower = new TCanvas("ecor_hoveretower","ecor_hoveretower");	
	canvas_ecor->cd(7);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_ecor_hoveretower = hdata->createHistogram(*hoveretower,*ecorvar);	
	TProfile *hprof_ecor_hoveretower = histo2_ecor_hoveretower->ProfileX();
       hprof_ecor_hoveretower->SetMarkerColor(kRed);
       hprof_ecor_hoveretower->SetMarkerStyle(7);
	hprof_ecor_hoveretower->SetErrorOption("s");
	hprof_ecor_hoveretower->GetYaxis()->SetRangeUser(plotECMin,plotECMax);
	hprof_ecor_hoveretower->GetXaxis()->SetTitle("hoveretower");
	hprof_ecor_hoveretower->GetYaxis()->SetTitle("ecor");
	hprof_ecor_hoveretower->Draw();
	//canvas_ecor_hoveretower->SaveAs("ecor_hovertower.png");

*/
	//________________
	//  ecor_rho
	//________________
	//TCanvas *canvas_ecor_rho = new TCanvas("ecor_rho","ecor_rho");	
	canvas_ecor->cd(4);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_ecor_rho = hdata->createHistogram(*rho,*ecorvar);	
	TProfile *hprof_ecor_rho = histo2_ecor_rho->ProfileX();
	hprof_ecor_rho->SetErrorOption("s");
	hprof_ecor_rho->GetYaxis()->SetRangeUser(plotECMin,plotECMax);
	hprof_ecor_rho->GetXaxis()->SetTitle("rho");
	hprof_ecor_rho->GetYaxis()->SetTitle("ecor");
	hprof_ecor_rho->Draw();

/*	//________________
	//  ecor_nvtx
	//________________
	//TCanvas *canvas_ecor_nvtx = new TCanvas("ecor_nvtx","ecor_nvtx");	
	canvas_ecor->cd(9);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_ecor_nvtx = hdata->createHistogram(*nvtx,*ecorvar);	
	TProfile *hprof_ecor_nvtx = histo2_ecor_nvtx->ProfileX();
	hprof_ecor_nvtx->SetErrorOption("s");
	hprof_ecor_nvtx->GetYaxis()->SetRangeUser(plotECMin,plotECMax);
	hprof_ecor_nvtx->GetXaxis()->SetTitle("nvtx");
	hprof_ecor_nvtx->GetYaxis()->SetTitle("ecor");
	hprof_ecor_nvtx->Draw();*/

//save canvas_ecor
	canvas_ecor->SaveAs("corrected_energy"+add2);
	canvas_ecor->SaveAs("corrected_energy"+addroot);


//create canvas for rawvar (raw energy)
	TCanvas *canvas_rawvar = new TCanvas("raw_energy","rawvar",canvas_width,canvas_height);
	canvas_rawvar->Divide(2,2);

	//________________
	//  rawvar_sceta
	//________________
	canvas_rawvar->cd(1);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_rawvar_sceta = hdata->createHistogram(*scetavar,*rawvar);	
	TProfile *hprof_rawvar_sceta = histo2_rawvar_sceta->ProfileX();
	hprof_rawvar_sceta->SetErrorOption("s");
	hprof_rawvar_sceta->GetYaxis()->SetRangeUser(plotERMin,plotERMax);
	hprof_rawvar_sceta->GetXaxis()->SetTitle("eta");
	hprof_rawvar_sceta->GetYaxis()->SetTitle("rawvar");
	hprof_rawvar_sceta->Draw();
	//canvas_rawvar_sceta->SaveAs("rawvar_sceta.png");

/*	//________________
	//  rawvar_scphi
	//________________
	//TCanvas *canvas_rawvar_scphi = new TCanvas("rawvar_scphi","rawvar_scphi");
	canvas_rawvar->cd(2);	
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_rawvar_scphi = hdata->createHistogram(*scphi,*rawvar);	
	TProfile *hprof_rawvar_scphi = histo2_rawvar_scphi->ProfileX();
	hprof_rawvar_scphi->SetErrorOption("s");
	hprof_rawvar_scphi->GetYaxis()->SetRangeUser(plotERMin,plotERMax);
	hprof_rawvar_scphi->GetXaxis()->SetTitle("phi");
	hprof_rawvar_scphi->GetYaxis()->SetTitle("rawvar");
	hprof_rawvar_scphi->Draw();
	//canvas_rawvar_scphi->SaveAs("rawvar_scphi.png");
	
	*/ 
	//________________
	//  rawvar_r9
	//________________
	//TCanvas *canvas_rawvar_r9 = new TCanvas("rawvar_r9","rawvar_r9");	
	canvas_rawvar->cd(2);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_rawvar_r9 = hdata->createHistogram(*r9,*rawvar);	
	TProfile *hprof_rawvar_r9 = histo2_rawvar_r9->ProfileX();
	hprof_rawvar_r9->SetErrorOption("s");
	hprof_rawvar_r9->GetYaxis()->SetRangeUser(plotERMin,plotERMax);
	hprof_rawvar_r9->GetXaxis()->SetTitle("R9");
	hprof_rawvar_r9->GetYaxis()->SetTitle("rawvar");
	hprof_rawvar_r9->Draw();
	//canvas_rawvar_r9->SaveAs("rawvar_r9.png");
/*
	//________________
	//  rawvar_scetawidth
	//________________
	//TCanvas *canvas_rawvar_scetawidth = new TCanvas("rawvar_scetawidth","rawvar_scetawidth");	
	canvas_rawvar->cd(4);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_rawvar_scetawidth = hdata->createHistogram(*scetawidth,*rawvar);	
	TProfile *hprof_rawvar_scetawidth = histo2_rawvar_scetawidth->ProfileX();
	hprof_rawvar_scetawidth->SetErrorOption("s");
	hprof_rawvar_scetawidth->GetYaxis()->SetRangeUser(plotERMin,plotERMax);
	hprof_rawvar_scetawidth->GetXaxis()->SetTitle("scetawidth");
	hprof_rawvar_scetawidth->GetYaxis()->SetTitle("rawvar");
	hprof_rawvar_scetawidth->Draw();
	//canvas_rawvar_scetawidth->SaveAs("rawvar_scetawidth.png");

	//________________
	//  rawvar_scphiwidth
	//________________
	//TCanvas *canvas_rawvar_scphiwidth = new TCanvas("rawvar_scphiwidth","rawvar_scphiwidth");	
	canvas_rawvar->cd(5);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_rawvar_scphiwidth = hdata->createHistogram(*scphiwidth,*rawvar);	
	TProfile *hprof_rawvar_scphiwidth = histo2_rawvar_scphiwidth->ProfileX();
	hprof_rawvar_scphiwidth->SetErrorOption("s");
	hprof_rawvar_scphiwidth->GetYaxis()->SetRangeUser(plotERMin,plotERMax);
	hprof_rawvar_scphiwidth->GetXaxis()->SetTitle("scphiwidth");
	hprof_rawvar_scphiwidth->GetYaxis()->SetTitle("rawvar");
	hprof_rawvar_scphiwidth->Draw();
	//canvas_rawvar_scphiwidth->SaveAs("rawvar_scphiwidth.png");

*/	//________________
	//  rawvar_pt
	//________________
	//TCanvas *canvas_rawvar_pt = new TCanvas("rawvar_pt","rawvar_pt");	
	canvas_rawvar->cd(3);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_rawvar_pt = hdata->createHistogram(*pt,*rawvar);	
	TProfile *hprof_rawvar_pt = histo2_rawvar_pt->ProfileX();
	hprof_rawvar_pt->SetErrorOption("s");
	hprof_rawvar_pt->GetYaxis()->SetRangeUser(plotERMin,plotERMax);
	hprof_rawvar_pt->GetXaxis()->SetTitle("pt");
	hprof_rawvar_pt->GetYaxis()->SetTitle("rawvar");
	hprof_rawvar_pt->Draw();
	//canvas_rawvar_pt->SaveAs("rawvar_pt.png");

/*	//________________
	//  rawvar_hoveretower
	//________________
	//TCanvas *canvas_rawvar_hoveretower = new TCanvas("rawvar_hoveretower","rawvar_hoveretower");	
	canvas_rawvar->cd(7);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_rawvar_hoveretower = hdata->createHistogram(*hoveretower,*rawvar);	
	TProfile *hprof_rawvar_hoveretower = histo2_rawvar_hoveretower->ProfileX();
	hprof_rawvar_hoveretower->SetErrorOption("s");
	hprof_rawvar_hoveretower->GetYaxis()->SetRangeUser(plotERMin,plotERMax);
	hprof_rawvar_hoveretower->GetXaxis()->SetTitle("hoveretower");
	hprof_rawvar_hoveretower->GetYaxis()->SetTitle("rawvar");
	hprof_rawvar_hoveretower->Draw();
	//canvas_rawvar_hoveretower->SaveAs("rawvar_hoveretower.png");
*/
	//________________
	//  rawvar_rho
	//________________
	//TCanvas *canvas_rawvar_rho = new TCanvas("rawvar_rho","rawvar_rho");	
	canvas_rawvar->cd(4);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_rawvar_rho = hdata->createHistogram(*rho,*rawvar);	
	TProfile *hprof_rawvar_rho = histo2_rawvar_rho->ProfileX();
	hprof_rawvar_rho->SetErrorOption("s");
	hprof_rawvar_rho->GetYaxis()->SetRangeUser(plotERMin,plotERMax);
	hprof_rawvar_rho->GetXaxis()->SetTitle("rho");
	hprof_rawvar_rho->GetYaxis()->SetTitle("rawvar");
	hprof_rawvar_rho->Draw();

/*	//________________
	//  rawvar_nvtx
	//________________
	//TCanvas *canvas_rawvar_nvtx = new TCanvas("rawvar_nvtx","rawvar_nvtx");	
	canvas_rawvar->cd(9);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_rawvar_nvtx = hdata->createHistogram(*nvtx,*rawvar);	
	TProfile *hprof_rawvar_nvtx = histo2_rawvar_nvtx->ProfileX();
	hprof_rawvar_nvtx->SetErrorOption("s");
	hprof_rawvar_nvtx->GetYaxis()->SetRangeUser(plotERMin,plotERMax);
	hprof_rawvar_nvtx->GetXaxis()->SetTitle("nvtx");
	hprof_rawvar_nvtx->GetYaxis()->SetTitle("rawvar");
	hprof_rawvar_nvtx->Draw();
*/

//save canvas_rawvar
	canvas_rawvar->SaveAs("raw_energy"+add2);
       canvas_rawvar->SaveAs("raw_energy"+addroot);

//create canvas for mean (#mu)
/*       Double_t maxMu(1.5), minMu(0.8);

	TCanvas *canvas_sigmean = new TCanvas("signal_mean","sigmean",canvas_width,canvas_height);
       gPad->SetGridx();
       gPad->SetGridy();
	canvas_sigmean->Divide(3,3);
	//________________
	//  sigmean_sceta
	//________________
	canvas_sigmean->cd(1);
       gPad->SetGridx();
       gPad->SetGridy();
       
	TH2* histo2_sigmean_sceta = hdata->createHistogram(*scetavar,*sigmean);	
	TProfile *hprof_sigmean_sceta = histo2_sigmean_sceta->ProfileX();
	hprof_sigmean_sceta->SetErrorOption("s");
	hprof_sigmean_sceta->GetYaxis()->SetRangeUser(minMu,maxMu);
	hprof_sigmean_sceta->GetXaxis()->SetTitle("eta");
	hprof_sigmean_sceta->GetYaxis()->SetTitle("mu");
	hprof_sigmean_sceta->Draw();
	//canvas_sigmean_sceta->SaveAs("sigmean_sceta.png");

	//________________
	//  sigmean_scphi
	//________________
	//TCanvas *canvas_sigmean_scphi = new TCanvas("sigmean_scphi","sigmean_scphi");
	canvas_sigmean->cd(2);
       gPad->SetGridx();
       gPad->SetGridy();	
	TH2* histo2_sigmean_scphi = hdata->createHistogram(*scphi,*sigmean);	
	TProfile *hprof_sigmean_scphi = histo2_sigmean_scphi->ProfileX();
	hprof_sigmean_scphi->SetErrorOption("s");
	hprof_sigmean_scphi->GetYaxis()->SetRangeUser(minMu,maxMu);
	hprof_sigmean_scphi->GetXaxis()->SetTitle("phi");
	hprof_sigmean_scphi->GetYaxis()->SetTitle("mu");
	hprof_sigmean_scphi->Draw();
	//canvas_sigmean_scphi->SaveAs("sigmean_scphi.png");
	
	 
	//________________
	//  sigmean_r9
	//________________
	//TCanvas *canvas_sigmean_r9 = new TCanvas("sigmean_r9","sigmean_r9");	
	canvas_sigmean->cd(3);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_sigmean_r9 = hdata->createHistogram(*r9,*sigmean);	
	TProfile *hprof_sigmean_r9 = histo2_sigmean_r9->ProfileX();
	hprof_sigmean_r9->SetErrorOption("s");
	hprof_sigmean_r9->GetYaxis()->SetRangeUser(minMu,maxMu);
	hprof_sigmean_r9->GetXaxis()->SetTitle("R9");
	hprof_sigmean_r9->GetYaxis()->SetTitle("mu");
	hprof_sigmean_r9->Draw();
	//canvas_sigmean_r9->SaveAs("sigmean_r9.png");

	//________________
	//  sigmean_scetawidth
	//________________
	//TCanvas *canvas_sigmean_scetawidth = new TCanvas("sigmean_scetawidth","sigmean_scetawidth");	
	canvas_sigmean->cd(4);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_sigmean_scetawidth = hdata->createHistogram(*scetawidth,*sigmean);	
	TProfile *hprof_sigmean_scetawidth = histo2_sigmean_scetawidth->ProfileX();
	hprof_sigmean_scetawidth->SetErrorOption("s");
	hprof_sigmean_scetawidth->GetYaxis()->SetRangeUser(minMu,maxMu);
	hprof_sigmean_scetawidth->GetXaxis()->SetTitle("scetawidth");
	hprof_sigmean_scetawidth->GetYaxis()->SetTitle("mu");
	hprof_sigmean_scetawidth->Draw();
	//canvas_sigmean_scetawidth->SaveAs("sigmean_scetawidth.png");

	//________________
	//  sigmean_scphiwidth
	//________________
	//TCanvas *canvas_sigmean_scphiwidth = new TCanvas("sigmean_scphiwidth","sigmean_scphiwidth");	
	canvas_sigmean->cd(5);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_sigmean_scphiwidth = hdata->createHistogram(*scphiwidth,*sigmean);	
	TProfile *hprof_sigmean_scphiwidth = histo2_sigmean_scphiwidth->ProfileX();
	hprof_sigmean_scphiwidth->SetErrorOption("s");
	hprof_sigmean_scphiwidth->GetYaxis()->SetRangeUser(minMu,maxMu);
	hprof_sigmean_scphiwidth->GetXaxis()->SetTitle("scphiwidth");
	hprof_sigmean_scphiwidth->GetYaxis()->SetTitle("mu");
	hprof_sigmean_scphiwidth->Draw();
	//canvas_sigmean_scphiwidth->SaveAs("sigmean_scphiwidth.png");

	//________________
	//  sigmean_pt
	//________________
	//TCanvas *canvas_sigmean_pt = new TCanvas("sigmean_pt","sigmean_pt");	
	canvas_sigmean->cd(6);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_sigmean_pt = hdata->createHistogram(*pt,*sigmean);	
	TProfile *hprof_sigmean_pt = histo2_sigmean_pt->ProfileX();
	hprof_sigmean_pt->SetErrorOption("s");
	hprof_sigmean_pt->GetYaxis()->SetRangeUser(minMu,maxMu);
	hprof_sigmean_pt->GetXaxis()->SetTitle("pt");
	hprof_sigmean_pt->GetYaxis()->SetTitle("mu");
	hprof_sigmean_pt->Draw();
	//canvas_sigmean_pt->SaveAs("sigmean_pt.png");

	//________________
	//  sigmean_hoveretower
	//________________
	//TCanvas *canvas_sigmean_hoveretower = new TCanvas("sigmean_hoveretower","sigmean_hoveretower");	
	canvas_sigmean->cd(7);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_sigmean_hoveretower = hdata->createHistogram(*hoveretower,*sigmean);	
	TProfile *hprof_sigmean_hoveretower = histo2_sigmean_hoveretower->ProfileX();
	hprof_sigmean_hoveretower->SetErrorOption("s");
	hprof_sigmean_hoveretower->GetYaxis()->SetRangeUser(minMu,maxMu);
	hprof_sigmean_hoveretower->GetXaxis()->SetTitle("hoveretower");
	hprof_sigmean_hoveretower->GetYaxis()->SetTitle("mu");
	hprof_sigmean_hoveretower->Draw();
	//canvas_sigmean_hoveretower->SaveAs("sigmean_hoveretower.png");

	//________________
	//  sigmean_rho
	//________________
	//TCanvas *canvas_sigmean_rho = new TCanvas("sigmean_rho","sigmean_rho");	
	canvas_sigmean->cd(8);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_sigmean_rho = hdata->createHistogram(*rho,*sigmean);	
	TProfile *hprof_sigmean_rho = histo2_sigmean_rho->ProfileX();
	hprof_sigmean_rho->SetErrorOption("s");
	hprof_sigmean_rho->GetYaxis()->SetRangeUser(minMu,maxMu);
	hprof_sigmean_rho->GetXaxis()->SetTitle("rho");
	hprof_sigmean_rho->GetYaxis()->SetTitle("mu");
	hprof_sigmean_rho->Draw();

	//________________
	//  sigmean_nvtx
	//________________
	//TCanvas *canvas_sigmean_nvtx = new TCanvas("sigmean_nvtx","sigmean_nvtx");	
	canvas_sigmean->cd(9);
       gPad->SetGridx();
       gPad->SetGridy();
	TH2* histo2_sigmean_nvtx = hdata->createHistogram(*nvtx,*sigmean);	
	TProfile *hprof_sigmean_nvtx = histo2_sigmean_nvtx->ProfileX();
	hprof_sigmean_nvtx->SetErrorOption("s");
	hprof_sigmean_nvtx->GetYaxis()->SetRangeUser(minMu,maxMu);
	hprof_sigmean_nvtx->GetXaxis()->SetTitle("nvtx");
	hprof_sigmean_nvtx->GetYaxis()->SetTitle("mu");
	hprof_sigmean_nvtx->Draw();

//save canvas_sigmean
	canvas_sigmean->SaveAs("signal_mean"+add2);
       canvas_sigmean->SaveAs("signal_mean"+addroot);*/
}
