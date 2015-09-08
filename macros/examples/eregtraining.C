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

 
using namespace RooFit;

 
void eregtraining(bool dobarrel=true, bool doele=false) {
   
       //output directory 
  	TString dirname = "ereg_sig3_Evts20_alpha36-07/"; 
  	gSystem->mkdir(dirname,true);
  	gSystem->cd(dirname);  
  	
       //values param
       Double_t minCutSig = 3., alpha1_const=3.6, alpha2_const=0.7;
       /*
       TCut prescale10 = "(event%10==0)";
  	TCut prescale20 = "(event%20==0)";
  	TCut prescale25 = "(event%25==0)";
  	TCut prescale50 = "(event%50==0)";
  	TCut prescale100 = "(event%100==0)";  
  	TCut prescale1000 = "(event%1000==0)";  
  	TCut evenevents = "(event%2==0)";
  	TCut oddevents = "(event%2==1)";  
  	*/
       TCut cutUsed = "(event%20==0)";
  	//build vectors with list of input variables 
  	std::vector<std::string> *varsf = new std::vector<std::string>;
  	varsf->push_back("scRawEnergy");
  	varsf->push_back("scEta");
  	varsf->push_back("scPhy");
  	varsf->push_back("r9");  
  	varsf->push_back("etaWidth");
  	varsf->push_back("phiWidth");  
  	varsf->push_back("scClustersSize");
  	varsf->push_back("hadTowOverEm"); 
  	varsf->push_back("rho");
  	varsf->push_back("nvtx");
 	
  	varsf->push_back("eta-scEta"); 
  	varsf->push_back("atan2(sin(phi-scPhy),cos(phi-scPhy))");
  	varsf->push_back("scSeedEnergy/scRawEnergy");
  	
  	varsf->push_back("e3x3/e5x5");
  	varsf->push_back("sigmaIetaIeta");   
  	varsf->push_back("sigmaIphiIphi");   
  	varsf->push_back("covarianceIetaIphi");
  	varsf->push_back("eMax/e5x5");
  	varsf->push_back("e2nd/e5x5");
  	varsf->push_back("eTop/e5x5");
  	varsf->push_back("eBottom/e5x5");
  	varsf->push_back("eLeft/e5x5");
  	varsf->push_back("eRight/e5x5");
  	varsf->push_back("e2x5Max/e5x5");
  	varsf->push_back("e2x5Top/e5x5");
  	varsf->push_back("e2x5Bottom/e5x5");
  	varsf->push_back("e2x5Left/e5x5");
  	varsf->push_back("e2x5Right/e5x5");
  	
  	std::vector<std::string> *varseb = new std::vector<std::string>(*varsf);
  	std::vector<std::string> *varsee = new std::vector<std::string>(*varsf);
  	
  	varseb->push_back("e5x5/scSeedEnergy");
  	
  	varseb->push_back("iEta");
  	varseb->push_back("iPhi");
  	varseb->push_back("(iEta-1*abs(iEta)/iEta)%5");
  	varseb->push_back("(iPhi-1)%2");       
  	varseb->push_back("(abs(iEta)<=25)*((iEta-1*abs(iEta)/iEta)%25) + (abs(iEta)>25)*((iEta-26*abs(iEta)/iEta)%20)");
  	varseb->push_back("(iPhi-1)%20"); 
  	varseb->push_back("cryEta");
  	varseb->push_back("cryPhi");
	
  	varsee->push_back("scPreshowerEnergy/scRawEnergy");
  	  
  	//select appropriate input list for barrel or endcap
  	std::vector<std::string> *varslist;
  	if (dobarrel) varslist = varseb;
  	else varslist = varsee;
  	
  	//create RooRealVars for each input variable
  	RooArgList vars;
  	for (unsigned int ivar=0; ivar<varslist->size(); ++ivar) 
	{
  	  	RooRealVar *var = new RooRealVar(TString::Format("var_%i",ivar),varslist->at(ivar).c_str(),0.);
  	  	vars.addOwned(*var);
  	}
  	
  	//make list of input variable RooRealVars
  	RooArgList condvars(vars);
  	
  	//create RooRealVar for target
  	RooRealVar *tgtvar = new RooRealVar("tgtvar","etrue/scRawEnergy",1.);
  	if (!dobarrel) tgtvar->SetTitle("etrue/(scRawEnergy + scPreshowerEnergy)");  
  	
  	//add target to full list
  	vars.addOwned(*tgtvar);
  	  
  	//RooRealVar for event weight 
  	RooRealVar weightvar("weightvar","",1.);
	
		
  	TChain *tree;

      
  	if (doele) 
	{ 
    		tree = new TChain("");
    		tree->Add("");
  	}
  	else 
	{
    		tree = new TChain("promptTree");
    		tree->Add("/afs/cern.ch/user/m/musella/public/forKenza/gam_gam_phys14_v5_regtraining_v3.root"); 
	}
  
  	//training selection cut 
  	TCut selcut;
  	if (dobarrel) 
	{
    		selcut = "pt>200. && abs(scEta)<=1.5";// && kSaturated[12]!=1" ; // && etrue/scRawEnergy <=2. && etrue/scRawEnergy>=0.6";
  	}
  	else 
	{
    		selcut = "pt>200. && abs(scEta)>1.5";// && kSaturated[12]!=1" ; // && etrue/scRawEnergy <=2. && etrue/scRawEnergy>=0.6" ;   
  	}
  

  
  	//TCut selweight = "(weight)";
       TCut selweight = "(1)";
  	TCut prescale10 = "(event%10==0)";
  	TCut prescale20 = "(event%20==0)";
  	TCut prescale25 = "(event%25==0)";
  	TCut prescale50 = "(event%50==0)";
  	TCut prescale100 = "(event%100==0)";  
  	TCut prescale1000 = "(event%1000==0)";  
  	TCut evenevents = "(event%2==0)";
  	TCut oddevents = "(event%2==1)";  
	
	
  	//weightvar title used for per-event weights and selection cuts
  	if (doele) {
  	  	weightvar.SetTitle(prescale100*evenevents*selcut);
  	}
  	else {
  	  	weightvar.SetTitle(cutUsed*selweight*selcut);
  	}
  	//create RooDataSet from TChain
  	RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",tree,vars,weightvar);   
  	
  	
  	//RooRealVars corresponding to regressed parameters (sigma, mean, left tail parameter, right tail parameter)
  	RooRealVar sigwidthtvar("sigwidthtvar","",0.01);
  	sigwidthtvar.setConstant(false);
  	
  	RooRealVar sigmeantvar("sigmeantvar","",1.);
  	sigmeantvar.setConstant(false); 
  	
  	RooRealVar sign1var("sign1var","",3.);
  	sign1var.setConstant(false);     
  	
  	RooRealVar sign2var("sign2var","",3.);
  	sign2var.setConstant(false);     
	
  	//define non-parametric functions for each regressed parameter
  	RooGBRFunctionFlex *sigwidthtfunc = new RooGBRFunctionFlex("sigwidthtfunc","");
  	RooGBRFunctionFlex *sigmeantfunc = new RooGBRFunctionFlex("sigmeantfunc","");
  	RooGBRFunctionFlex *sign1func = new RooGBRFunctionFlex("sign1func","");
  	RooGBRFunctionFlex *sign2func = new RooGBRFunctionFlex("sign2func","");
	
  	//define mapping of input variables to non-parametric functions (in this case trivial since all 4 functions depend on the same inputs, but this is not a requirement)
  	RooGBRTargetFlex *sigwidtht = new RooGBRTargetFlex("sigwidtht","",*sigwidthtfunc,sigwidthtvar,condvars);  
  	RooGBRTargetFlex *sigmeant = new RooGBRTargetFlex("sigmeant","",*sigmeantfunc,sigmeantvar,condvars);  
  	RooGBRTargetFlex *sign1t = new RooGBRTargetFlex("sign1t","",*sign1func,sign1var,condvars);  
  	RooGBRTargetFlex *sign2t = new RooGBRTargetFlex("sign2t","",*sign2func,sign2var,condvars);  
	
  	//define list of mapped functions to regress
  	RooArgList tgts;
 	tgts.add(*sigwidtht);
  	tgts.add(*sigmeant);
  	tgts.add(*sign1t);
  	tgts.add(*sign2t);  
  
  	//define transformations corresponding to parameter bounds for non-parametric outputs  
  	RooRealConstraint sigwidthlim("sigwidthlim","",*sigwidtht,0.0002,0.5);
  	RooRealConstraint sigmeanlim("sigmeanlim","",*sigmeant,0.2,2.0);
  	RooRealConstraint sign1lim("sign1lim","",*sign1t,0.01,100.); 
  	RooRealConstraint sign2lim("sign2lim","",*sign2t,1.01,100.); 
	
  	//define pdf, which depends on transformed outputs (and is intended to be treated as a conditional pdf over the
  	//regression inputs in this case)
  	//The actual pdf below is a double crystal ball, with crossover points alpha_1 and alpha_2 set constant, but all other
  	//parameters regressed
  	RooDoubleCBFast sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,RooConst(alpha1_const),sign1lim,RooConst(alpha2_const),sign2lim);
  	
  	//dummy variable
  	RooConstVar etermconst("etermconst","",0.);  
  	 
  	//dummy variable
  	RooRealVar r("r","",1.);
  	r.setConstant();
	
  	//define list of pdfs
  	std::vector<RooAbsReal*> vpdf;
  	vpdf.push_back(&sigpdf);  
  	
  	//define list of training datasets
  	std::vector<RooAbsData*> vdata;
  	vdata.push_back(hdata);     
  	
  	//define minimum event weight per tree node
  	double minweight = 200;
  	std::vector<double> minweights;
  	minweights.push_back(minweight);
  	
  	//temp output file
  	TFile *fres = new TFile("fres.root","RECREATE");
	
  	//run training
  	if (1) 
	{
  	  	RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff","",tgts,etermconst,r,vdata,vpdf);
  	  	bdtpdfdiff.SetMinCutSignificance(minCutSig);
  		//bdtpdfdiff.SetPrescaleInit(100);
  	  	bdtpdfdiff.SetShrinkage(0.1);
  	  	bdtpdfdiff.SetMinWeights(minweights);
  	  	bdtpdfdiff.SetMaxNodes(750);
  	  	bdtpdfdiff.TrainForest(1e6);   
  	}
  	   
  	//create workspace and output to file
  	RooWorkspace *wereg = new RooWorkspace("wereg");
  	wereg->import(sigpdf);

  	if (doele && dobarrel)
  	  	wereg->writeToFile("wereg_ele_eb.root");    
  	else if (doele && !dobarrel) 
  	  	wereg->writeToFile("wereg_ele_ee.root");    
  	else if (!doele && dobarrel)
  	  	wereg->writeToFile("wereg_ph_eb.root");    
  	else if (!doele && !dobarrel)
  	  	wereg->writeToFile("wereg_ph_ee.root");    
	
  	return;
  	
  
}			
