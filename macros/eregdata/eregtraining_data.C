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
#include "TRandom3.h"
 
using namespace RooFit;

void swapvars(RooArgList &vars) {
 
  for (int ivar=0; ivar<vars.getSize(); ++ivar) {
    TString title = vars.at(ivar)->GetTitle();
    title.ReplaceAll("ph1.","ph3.");
    title.ReplaceAll("ph2.","ph1.");
    title.ReplaceAll("ph3.","ph2.");
    vars.at(ivar)->SetTitle(title);
  }  
  
  
}


void eregtraining_data() {
  
  TString dirname = "/data/bendavid/regflexmasstestingOct23_closure/"; 
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);  
  
  //fill list of input variable expressions (any valid TTree draw expression will work here)
  std::vector<std::string> *varsf = new std::vector<std::string>;
  varsf->push_back("0.*ph.e");
  varsf->push_back("ph.eerr/(ph.e)");
  varsf->push_back("0.*(ph.scrawe)");
  varsf->push_back("ph.sceta");
  varsf->push_back("ph.scphi");
  varsf->push_back("ph.r9");  
  varsf->push_back("ph.scetawidth");
  varsf->push_back("ph.scphiwidth");  
  varsf->push_back("ph.scnclusters");
  varsf->push_back("ph.hoveretower");
  varsf->push_back("rho");
  varsf->push_back("nVtx");  
 
  varsf->push_back("ph.etaseed-ph.sceta");
  varsf->push_back("atan2(sin(ph.phiseed-ph.scphi),cos(ph.phiseed-ph.scphi))");
  varsf->push_back("ph.eseed/(ph.scrawe)");
  
  varsf->push_back("ph.e3x3seed/ph.e5x5seed");
  varsf->push_back("ph.sigietaietaseed");   
  varsf->push_back("ph.sigiphiphiseed");   
  varsf->push_back("ph.covietaiphiseed");
  varsf->push_back("ph.emaxseed/ph.e5x5seed");
  varsf->push_back("ph.e2ndseed/ph.e5x5seed");
  varsf->push_back("ph.etopseed/ph.e5x5seed");
  varsf->push_back("ph.ebottomseed/ph.e5x5seed");
  varsf->push_back("ph.eleftseed/ph.e5x5seed");
  varsf->push_back("ph.erightseed/ph.e5x5seed");
  varsf->push_back("ph.e2x5maxseed/ph.e5x5seed");
  varsf->push_back("ph.e2x5topseed/ph.e5x5seed");
  varsf->push_back("ph.e2x5bottomseed/ph.e5x5seed");
  varsf->push_back("ph.e2x5leftseed/ph.e5x5seed");
  varsf->push_back("ph.e2x5rightseed/ph.e5x5seed");

  varsf->push_back("(ph.e)/(ph.scrawe)");  
  varsf->push_back("ph.isbarrel*ph.e5x5seed/ph.eseed");
  
  varsf->push_back("ph.isbarrel*ph.ietaseed");
  varsf->push_back("ph.isbarrel*ph.iphiseed");
  varsf->push_back("ph.isbarrel*(ph.ietaseed-1*abs(ph.ietaseed)/ph.ietaseed)%5");
  varsf->push_back("ph.isbarrel*(ph.iphiseed-1)%2");       
  varsf->push_back("ph.isbarrel*((abs(ph.ietaseed)<=25)*((ph.ietaseed-1*abs(ph.ietaseed)/ph.ietaseed)%25) + (abs(ph.ietaseed)>25)*((ph.ietaseed-26*abs(ph.ietaseed)/ph.ietaseed)%20))");
  varsf->push_back("ph.isbarrel*(ph.iphiseed-1)%20"); 
  varsf->push_back("ph.isbarrel*ph.etacryseed");
  varsf->push_back("ph.isbarrel*ph.phicryseed");

  varsf->push_back("!ph.isbarrel*(ph.e)/((ph.scrawe) + ph.scpse)");  
  varsf->push_back("!ph.isbarrel*(ph.scpse)/(ph.scrawe)");
    
  std::vector<std::string> *varslist = varsf;
//   if (dobarrel) varslist = varseb;
//   else varslist = varsee;
  //varslist = varseb;
  
  //construct RooRealVars from variable list
  RooArgList vars_ph1;
  RooArgList vars_ph2;  
  RooArgList vars_common;
  RooArgList vars_all;  
  for (unsigned int ivar=0; ivar<varslist->size(); ++ivar) {
    TString sname = varslist->at(ivar).c_str();
    TString title_ph1(sname);
    title_ph1.ReplaceAll("ph.","ph1.");
    TString title_ph2(sname);
    title_ph2.ReplaceAll("ph.","ph2.");
    
    printf("ivar = %i, sname = %s, t1 = %s, t2 = %s\n",ivar,sname.Data(),title_ph1.Data(),title_ph2.Data());
    if (title_ph1==title_ph2) {
      RooRealVar *var = new RooRealVar(TString::Format("var_common_%i",ivar), sname ,0.);
      //vars_ph1.add(*var);
      //vars_ph2.add(*var);
      vars_common.add(*var);
      vars_all.add(*var);
    }
    else {
      RooRealVar *var_ph1 = new RooRealVar(TString::Format("var_ph1_%i",ivar), title_ph1 ,0.);
      RooRealVar *var_ph2 = new RooRealVar(TString::Format("var_ph2_%i",ivar), title_ph2 ,0.);
      vars_ph1.add(*var_ph1);
      vars_ph2.add(*var_ph2);
      vars_all.add(*var_ph1);
      vars_all.add(*var_ph2);      
    }
  }
  
  vars_ph1.Print("V");
  vars_ph2.Print("V");
  vars_all.Print("V");


  //define some additional RooRealVars
  RooRealVar *tgtvar = new RooRealVar("tgtvar","mass",1.);

  const double pdgzmass = 91.1876;
  const double logpdgzmass = log(pdgzmass);
  
  RooRealVar *mass = new RooRealVar("mass","sqrt(2.0*(ph1.e)*(ph2.e)*(1.0-costheta))",90.,70.,110.);
  
  RooRealVar *logmass = new RooRealVar("logmass","0.5*log(2.0) + 0.5*log((ph1.e)) + 0.5*log((ph2.e)) + 0.5*log(1.0-costheta) -log(91.1876)",0.,log(70.)-logpdgzmass,log(110.)-logpdgzmass);
  RooRealVar *relpt_ph1 = new RooRealVar("relpt_ph1","ph1.pt/mass",1.);
  RooRealVar *relpt_ph2 = new RooRealVar("relpt_ph2","ph2.pt/mass",1.);  
  RooRealVar *eta_ph1 = new RooRealVar("eta_ph1","ph1.eta",1.);
  RooRealVar *eta_ph2 = new RooRealVar("eta_ph2","ph2.eta",1.);
  RooRealVar *cosdphi = new RooRealVar("cosdphi","TMath::Cos(ph1.phi-ph2.phi)",1.);
  
  RooRealVar *evt = new RooRealVar("evt","evt",0.);
  RooRealVar *run = new RooRealVar("run","run",0.);
  RooRealVar *lumi = new RooRealVar("lumi","lumi",0.);
  
  RooRealVar *costheta = new RooRealVar("costheta","costheta",0.);
  RooRealVar *deltaE = new RooRealVar("deltaE","0.5*log((ph1.e)) - 0.5*log((ph2.e))",0.);
  RooRealVar *deltaEswap = new RooRealVar("deltaEswap","0.5*log((ph2.e)) - 0.5*log((ph1.e))",0.);  
  RooRealVar *gene_ph1 = new RooRealVar("gene_ph1","ph1.gene",0.);
  RooRealVar *gene_ph2 = new RooRealVar("gene_ph2","ph2.gene",0.);
  RooRealVar *genmass = new RooRealVar("genmass","genmass",0.);
  
  RooRealVar *e_ph1 = new RooRealVar("e_ph1","(ph1.e)",0.);
  RooRealVar *scrawe_ph1 = new RooRealVar("scrawe_ph1","(ph1.scrawe)",0.);
  
  RooRealVar *e_ph2 = new RooRealVar("e_ph2","ph2.e",0.);
  RooRealVar *scrawe_ph2 = new RooRealVar("scrawe_ph2","(ph2.scrawe)",0.);  

  RooArgList vars_common_data;
  vars_common_data.add(*run);
  vars_common.add(*lumi);  
  
  RooArgList vars_common_mass;
  vars_common_mass.add(*costheta);
  //vars_common_mass.add(*deltaE);
  
  RooArgList varsNoEdep_ph1(vars_ph1);
  RooArgList varsNoEdep_ph2(vars_ph2);
  
//   varsNoEdep_ph1.add(*relpt_ph1);
//   varsNoEdep_ph2.add(*relpt_ph2);
  
  vars_ph1.add(*e_ph1);
  vars_ph1.add(*scrawe_ph1);  
  
  vars_ph2.add(*e_ph2);
  vars_ph2.add(*scrawe_ph2);  
  
  
  RooArgList vars;
  vars.add(vars_ph1);
  vars.add(vars_ph2);
  vars.add(vars_common);
  vars.add(vars_common_data);
  vars.add(vars_common_mass);
  vars.add(*deltaE);
  vars.add(*deltaEswap);
  vars.add(*relpt_ph1);
  vars.add(*relpt_ph2);
  vars.add(*logmass);
  vars.add(*mass);
  vars.add(*evt);
  
  //define lists of input variables used later to define pdfs
  RooArgList corvarsNoEdep_ph1;
  corvarsNoEdep_ph1.add(varsNoEdep_ph1);
  corvarsNoEdep_ph1.add(vars_common);
  corvarsNoEdep_ph1.add(vars_common_data);
  corvarsNoEdep_ph1.add(*costheta);
  corvarsNoEdep_ph1.add(*deltaE);

  RooArgList corvarsNoEdep_ph2;
  corvarsNoEdep_ph2.add(varsNoEdep_ph2);
  corvarsNoEdep_ph2.add(vars_common);
  corvarsNoEdep_ph2.add(vars_common_data);  
  corvarsNoEdep_ph1.add(*costheta);  
  corvarsNoEdep_ph2.add(*deltaEswap);
  
  RooArgList corvars_ph1;
  corvars_ph1.add(vars_ph1);
  corvars_ph1.add(vars_common);
  corvars_ph1.add(vars_common_data);

  RooArgList corvars_ph2;
  corvars_ph2.add(vars_ph2);
  corvars_ph2.add(vars_common);
  corvars_ph2.add(vars_common_data);  
    
   
  
  RooArgList condvarsmassNoEdep;
  condvarsmassNoEdep.add(varsNoEdep_ph1);
  condvarsmassNoEdep.add(varsNoEdep_ph2);
  condvarsmassNoEdep.add(vars_common);
  condvarsmassNoEdep.add(*costheta);  
  condvarsmassNoEdep.add(*deltaE);
  
  RooArgList condvarsmassNoEdepData;
  condvarsmassNoEdepData.add(condvarsmassNoEdep);
  condvarsmassNoEdep.add(vars_common_data);    
  
      
  //define RooRealVar for event weight
  RooRealVar weightvar("weightvar","",1.);


  //load trees
  TString treeloc = "RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/JetPub/JetCorrectionMod/SeparatePileUpMod/ElectronIDMod/MuonIDMod/PhotonPairSelectorPreselInvertEleVetoNoSmear/PhotonTreeWriterPreselInvertEleVetoNoSmear/hPhotonTree";
     
  TChain *tree = new TChain(treeloc);
  tree->Add("/data/bendavid/diphoTrees8TeVOct6/hgg-2013Final8TeV_s12-zllm50-v7n_noskim.root");

  TChain *treedata = new TChain(treeloc);
  treedata->Add("/data/bendavid/diphoTrees8TeVOct6/hgg-2013Final8TeV_r12a-pho-j22-v1_noskim.root");  
  treedata->Add("/data/bendavid/diphoTrees8TeVOct6/hgg-2013Final8TeV_r12b-dph-j22-v1_noskim.root");  
  treedata->Add("/data/bendavid/diphoTrees8TeVOct6/hgg-2013Final8TeV_r12c-dph-j22-v1_noskim.root");  
  treedata->Add("/data/bendavid/diphoTrees8TeVOct6/hgg-2013Final8TeV_r12d-dph-j22-v1_noskim.root");   

  //define cuts  
  
  int prescaleinit = 250;
  
  TCut selcut = "mass>70. && mass<110. && ph1.pt>28. && ph2.pt>20.";
  //TCut selcut = "mass>110. && ph1.pt>28. && ph2.pt>20.";
  TCut tcut = "evt%4==0";
  TCut mccut =  "evt%8==0";
  TCut mcrawcut =  "evt%8==2";
  
  TCut selweight = "xsecweight(procidx)";
  TCut prescale10 = "(evt%10==0)";
  TCut prescale20 = "(evt%20==0)";
  TCut prescale25 = "(evt%25==0)";
  TCut prescale50 = "(evt%50==0)";
  TCut prescale100 = "(evt%100==0)";  
  TCut prescale1000 = "(evt%1000==0)";  
  TCut evenevents = "(evt%2==0)";
  TCut oddevents = "(evt%2==1)";
  
  //TCut halfweight = "(0.25*1.0)";

  
  RooDataSet *hdatanoswap = 0;
  RooDataSet *hdatanoswapD = 0;
  RooDataSet *hdataswap = 0;
  RooDataSet *hdataswapD = 0;  
  {
    //load trees into datasets using list of RooRealVars
    //RooRealVar name is arbitrary, RooRealVar Title is used as TTree draw expression (and the weight variable Title becomes the per event weight)
    weightvar.SetTitle(selcut*tcut);
    hdatanoswap = RooTreeConvert::CreateDataSet("hdatanoswap",tree,vars,weightvar);   
    hdatanoswapD = RooTreeConvert::CreateDataSet("hdatanoswapD",treedata,vars,weightvar);  
  }  
  
  //swap variable titles for electron 1 and electron 2
  swapvars(vars);

  printf("swapped vars:\n");
  for (int ivar=0; ivar<vars.getSize(); ++ivar) {
    printf("%i: %s, %s\n",ivar, vars.at(ivar)->GetName(),vars.at(ivar)->GetTitle());
  }

  
  {
    //load trees with electron 1 and 2 swapped
    weightvar.SetTitle(selcut*tcut);
    hdataswap = RooTreeConvert::CreateDataSet("hdataswap",tree,vars,weightvar);   
    hdataswapD = RooTreeConvert::CreateDataSet("hdataswapD",treedata,vars,weightvar);  
  }  
   
  swapvars(vars);  

 
  
  printf("final vars:\n");
  for (int ivar=0; ivar<vars.getSize(); ++ivar) {
    printf("%i: %s, %s\n",ivar, vars.at(ivar)->GetName(),vars.at(ivar)->GetTitle());
  }
  
  assert(hdatanoswap->numEntries() == hdataswap->numEntries());
  assert(hdatanoswapD->numEntries() == hdataswapD->numEntries());
  
  RooArgList varsw(vars);
  varsw.add(weightvar);
  
  RooDataSet *inset = 0;
  
  
  //fill final dataset with randomized choice for "electron 1" and "electron 2"
  RooDataSet *hdata = new RooDataSet("hdata","",varsw,RooFit::WeightVar(weightvar));  
  for (int iev=0; iev<hdatanoswap->numEntries(); ++iev) {
    int dsel = RooRandom::integer(2);
    
    if (dsel==0) {
      inset = hdatanoswap;
    }
    else {
      inset = hdataswap;
    }
    
    const RooArgSet *dset = inset->get(iev);
    double weight = inset->weight();
    
    hdata->add(*dset,weight);
    
  }
  
  RooDataSet *hdataAll = hdatanoswap;
  hdataAll->append(*hdataswap);
  
  delete hdataswap;
  
  RooDataSet *hdataD = new RooDataSet("hdataD","",varsw,RooFit::WeightVar(weightvar));  
  for (int iev=0; iev<hdatanoswapD->numEntries(); ++iev) {
    int dsel = RooRandom::integer(2);
    
    if (dsel==0) {
      inset = hdatanoswapD;
    }
    else {
      inset = hdataswapD;
    }
    
    const RooArgSet *dset = inset->get(iev);
    double weight = inset->weight();
    
    hdataD->add(*dset,weight);
    
  }  
  
  RooDataSet *hdataAllD = hdatanoswapD;
  hdataAllD->append(*hdataswapD);  
    
  delete hdataswapD;
  
  printf("hdata:  sum = %5f, num = %i\n",hdata->sumEntries(),hdata->numEntries());
  printf("hdataD: sum = %5f, num = %i\n",hdataD->sumEntries(),hdataD->numEntries());
  
  //define  RooRealVars corresponding to non-parametric functions (scaling and smearing factors)
  RooRealVar *scalevar_ph1 = new RooRealVar("scalevar_ph1","",1.);
  RooRealVar *scalevar_ph2 = new RooRealVar("scalevar_ph2","",1.);
  RooRealVar *smearvar_ph1 = new RooRealVar("smearvar_ph1","",pow(0.005,2));
  RooRealVar *smearvar_ph2 = new RooRealVar("smearvar_ph2","",pow(0.005,2));

  RooRealVar *scaleNoEdepvar_ph1 = new RooRealVar("scaleNoEdepvar_ph1","",1.);
  RooRealVar *smearNoEdepvar_ph1 = new RooRealVar("smearNoEdepvar_ph1","",pow(0.005,2));  
  RooRealVar *scaleNoEdepvar_ph2 = new RooRealVar("scaleNoEdepvar_ph2","",1.);
  RooRealVar *smearNoEdepvar_ph2 = new RooRealVar("smearNoEdepvar_ph2","",pow(0.005,2));
  
  
  scalevar_ph1->setConstant(false);
  scalevar_ph2->setConstant(false);
  smearvar_ph1->setConstant(false);
  smearvar_ph2->setConstant(false);
  
  scaleNoEdepvar_ph1->setConstant(false);
  smearNoEdepvar_ph1->setConstant(false);
  scaleNoEdepvar_ph2->setConstant(false);
  smearNoEdepvar_ph2->setConstant(false);  
  
  //define non-parametric functions
  RooGBRFunctionFlex *scalefunc = new RooGBRFunctionFlex("scalefunc","");
  RooGBRFunctionFlex *smearfunc = new RooGBRFunctionFlex("smearfunc","");

  RooGBRFunctionFlex *scalefuncNoEdep = new RooGBRFunctionFlex("scalefuncNoEdep","");
  RooGBRFunctionFlex *smearfuncNoEdep = new RooGBRFunctionFlex("smearfuncNoEdep","");    
  
  //define mapping of input variables to non-parametric functions
  RooGBRTargetFlex *scale_ph1 = new RooGBRTargetFlex("scale_ph1","",*scalefunc,*scalevar_ph1,corvars_ph1);
  RooGBRTargetFlex *scale_ph2 = new RooGBRTargetFlex("scale_ph2","",*scalefunc,*scalevar_ph2,corvars_ph2);
  RooGBRTargetFlex *smear_ph1 = new RooGBRTargetFlex("smear_ph1","",*smearfunc,*smearvar_ph1,corvars_ph1);
  RooGBRTargetFlex *smear_ph2 = new RooGBRTargetFlex("smear_ph2","",*smearfunc,*smearvar_ph2,corvars_ph2);  

  RooGBRTargetFlex *scaleNoEdep_ph1 = new RooGBRTargetFlex("scaleNoEdep_ph1","",*scalefuncNoEdep,*scaleNoEdepvar_ph1,corvarsNoEdep_ph1);
  RooGBRTargetFlex *smearNoEdep_ph1 = new RooGBRTargetFlex("smearNoEdep_ph1","",*smearfuncNoEdep,*smearNoEdepvar_ph1,corvarsNoEdep_ph1);    
  RooGBRTargetFlex *scaleNoEdep_ph2 = new RooGBRTargetFlex("scaleNoEdep_ph2","",*scalefuncNoEdep,*scaleNoEdepvar_ph2,corvarsNoEdep_ph2);
  RooGBRTargetFlex *smearNoEdep_ph2 = new RooGBRTargetFlex("smearNoEdep_ph2","",*smearfuncNoEdep,*smearNoEdepvar_ph2,corvarsNoEdep_ph2);    
  
  const double negsmearlim = 0.0;
  
  //define transformations corresponding to parameter bounds for non-parametric outputs
  RooRealConstraint *scalelim_ph1 = new RooRealConstraint("scalelim_ph1","",*scale_ph1,0.7,1.3);
  RooRealConstraint *scalelim_ph2 = new RooRealConstraint("scalelim_ph2","",*scale_ph2,0.7,1.3);  
  RooRealConstraint *smearlim_ph1 = new RooRealConstraint("smearlim_ph1","",*smear_ph1,pow(1e-7,2),pow(0.2,2));
  RooRealConstraint *smearlim_ph2 = new RooRealConstraint("smearlim_ph2","",*smear_ph2,pow(1e-7,2),pow(0.2,2));    

  RooRealConstraint *scaleNoEdeplim_ph1 = new RooRealConstraint("scaleNoEdeplim_ph1","",*scaleNoEdep_ph1,0.7,1.3);  
  RooRealConstraint *smearNoEdeplim_ph1 = new RooRealConstraint("smearNoEdeplim_ph1","",*smearNoEdep_ph1,-pow(negsmearlim,2),pow(0.2,2));    
  RooRealConstraint *scaleNoEdeplim_ph2 = new RooRealConstraint("scaleNoEdeplim_ph2","",*scaleNoEdep_ph2,0.7,1.3);  
  RooRealConstraint *smearNoEdeplim_ph2 = new RooRealConstraint("smearNoEdeplim_ph2","",*smearNoEdep_ph2,-pow(negsmearlim,2),pow(0.2,2));     
  
  //define some additional non-parametric functions (width of pdfs for variable transformation step)
  RooRealVar *scalewidthvar_ph1 = new RooRealVar("scalewidthvar_ph1","",0.05);
  scalewidthvar_ph1->setConstant(false);
  RooGBRFunctionFlex *scalewidthfunc = new RooGBRFunctionFlex("scalewidthfunc","");
  RooGBRTargetFlex *scalewidth_ph1 = new RooGBRTargetFlex("scalewidth_ph1","",*scalewidthfunc,*scalewidthvar_ph1,corvars_ph1);
  RooRealConstraint *scalewidthlim_ph1 = new RooRealConstraint("scalewidthlim_ph1","",*scalewidth_ph1,pow(1e-7,1),pow(0.2,1));
  
  RooRealVar *smearwidthvar_ph1 = new RooRealVar("smearwidthvar_ph1","",pow(0.05,2));
  smearwidthvar_ph1->setConstant(false);
  RooGBRFunctionFlex *smearwidthfunc = new RooGBRFunctionFlex("smearwidthfunc","");
  RooGBRTargetFlex *smearwidth_ph1 = new RooGBRTargetFlex("smearwidth_ph1","",*smearwidthfunc,*smearwidthvar_ph1,corvars_ph1);
  RooRealConstraint *smearwidthlim_ph1 = new RooRealConstraint("smearwidthlim_ph1","",*smearwidth_ph1,pow(1e-4,2),pow(0.2,2));    
  
  RooRealVar *scalewidthvar_ph2 = new RooRealVar("scalewidthvar_ph2","",0.05);
  scalewidthvar_ph2->setConstant(false);
  RooGBRTargetFlex *scalewidth_ph2 = new RooGBRTargetFlex("scalewidth_ph2","",*scalewidthfunc,*scalewidthvar_ph2,corvars_ph2);
  RooRealConstraint *scalewidthlim_ph2 = new RooRealConstraint("scalewidthlim_ph2","",*scalewidth_ph2,pow(1e-7,1),pow(0.2,1));
  
  RooRealVar *smearwidthvar_ph2 = new RooRealVar("smearwidthvar_ph2","",pow(0.05,2));
  smearwidthvar_ph2->setConstant(false);
  RooGBRTargetFlex *smearwidth_ph2 = new RooGBRTargetFlex("smearwidth_ph2","",*smearwidthfunc,*smearwidthvar_ph2,corvars_ph2);
  RooRealConstraint *smearwidthlim_ph2 = new RooRealConstraint("smearwidthlim_ph2","",*smearwidth_ph2,pow(1e-4,2),pow(0.2,2));      
  
  //define non-parametric functions for ln(m) pdfs on mc and data (sums of gaussians)
  int ngaus = 6;  
  double step = 0.2/double(std::max(1,ngaus-1));
  
  RooArgList tgtsmassNoEdep;
  
  RooArgList gauspdfsNoEdep;
  RooArgList gauspdfsdataNoEdep;
  RooArgList gauspdfsdataNoEdepsingle;  
  RooArgList gauscoeffsNoEdep;
    
  for (int igaus=0; igaus<ngaus; ++igaus) {
    RooRealVar *gmeanvar = new RooRealVar(TString::Format("gNoEdep_meanvar_%i",igaus),"",log(80.)-logpdgzmass + step*igaus);
    RooRealVar *gsigmavar = new RooRealVar(TString::Format("gNoEdep_sigmavar_%i",igaus),"",0.02);
    RooRealVar *gfracvar = new RooRealVar(TString::Format("gNoEdep_fracvar_%i",igaus),"",1.0);
    
    gmeanvar->setConstant(false);
    gsigmavar->setConstant(false);
    gfracvar->setConstant(false);
    
    RooGBRFunctionFlex *gmeanfunc = new RooGBRFunctionFlex(TString::Format("gNoEdep_meanfunc_%i",igaus),"");
    RooGBRFunctionFlex *gsigmafunc = new RooGBRFunctionFlex(TString::Format("gNoEdep_sigmafunc_%i",igaus),"");
    RooGBRFunctionFlex *gfracfunc = new RooGBRFunctionFlex(TString::Format("gNoEdep_fracfunc_%i",igaus),"");
    
    RooGBRTargetFlex *gmean = new RooGBRTargetFlex(TString::Format("gNoEdep_mean_%i",igaus),"",*gmeanfunc,*gmeanvar,condvarsmassNoEdep);
    RooGBRTargetFlex *gsigma = new RooGBRTargetFlex(TString::Format("gNoEdep_sigma_%i",igaus),"",*gsigmafunc,*gsigmavar,condvarsmassNoEdep);
    RooGBRTargetFlex *gfrac = new RooGBRTargetFlex(TString::Format("gNoEdep_frac_%i",igaus),"",*gfracfunc,*gfracvar,condvarsmassNoEdep);

    RooRealConstraint *gmeanlim = new RooRealConstraint(TString::Format("gNoEdep_meanlim_%i",igaus),"",*gmean,log(75.)-logpdgzmass,log(105.)-logpdgzmass);
    RooRealConstraint *gsigmalim = new RooRealConstraint(TString::Format("gNoEdep_sigmalim_%i",igaus),"",*gsigma,sqrt(0.5*negsmearlim*negsmearlim + 1e-5),0.3);
    
    RooAbsReal *gfraclim = new RooProduct(TString::Format("gNoEdep_fraclim_%i",igaus),"",RooArgList(*gfrac,*gfrac));
     
    RooFormulaVar *gmeanscale = new RooFormulaVar(TString::Format("gNoEdep_meanscale_%i",igaus),"","@0 - 0.5*log(@1*@2)",RooArgList(*gmeanlim,*scaleNoEdeplim_ph1,*scaleNoEdeplim_ph2));
    RooFormulaVar *gsigmasmear = new RooFormulaVar(TString::Format("gNoEdep_sigmasmear_%i",igaus),"","sqrt(@0*@0 + 0.25*(@1+@2))",RooArgList(*gsigmalim,*smearNoEdeplim_ph1,*smearNoEdeplim_ph2)); 

    RooFormulaVar *gmeanscalesingle = new RooFormulaVar(TString::Format("gNoEdep_meanscalesingle_%i",igaus),"","@0 - log(@1)",RooArgList(*gmeanlim,*scaleNoEdeplim_ph1));
    RooFormulaVar *gsigmasmearsingle = new RooFormulaVar(TString::Format("gNoEdep_sigmasmearsingle_%i",igaus),"","sqrt(@0*@0 + 0.5*@1)",RooArgList(*gsigmalim,*smearNoEdeplim_ph1));     
    
    
    if (igaus==0) {
      gfraclim = new RooConstVar(TString::Format("gNoEdep_fraclimconst_%i",igaus),"",1.);
    }
    else {
      tgtsmassNoEdep.add(*gfrac);
    }
    
    RooGaussianFast *gpdf = new RooGaussianFast(TString::Format("gNoEdep_pdf_%i",igaus),"",*logmass,*gmeanlim,*gsigmalim);
    RooGaussianFast *gpdfdata = new RooGaussianFast(TString::Format("gNoEdep_pdfdata_%i",igaus),"",*logmass,*gmeanscale,*gsigmasmear);
    RooGaussianFast *gpdfdatasingle = new RooGaussianFast(TString::Format("gNoEdep_pdfdatasingle_%i",igaus),"",*logmass,*gmeanscalesingle,*gsigmasmearsingle);
    
    //RooBifurGauss *gpdf = new RooBifurGauss(TString::Format("gNoEdep_df_%i",igaus),"",idmva,*gmeanlim,*gsigmalim,*gsigmaRlim);
    
    gauspdfsNoEdep.add(*gpdf);
    gauscoeffsNoEdep.add(*gfraclim);
    
    gauspdfsdataNoEdep.add(*gpdfdata);
    gauspdfsdataNoEdepsingle.add(*gpdfdatasingle);    
    
    tgtsmassNoEdep.add(*gmean);
    tgtsmassNoEdep.add(*gsigma);
    
  }  
  RooCondAddPdf masspdfNoEdep("masspdfNoEdep","",gauspdfsNoEdep,gauscoeffsNoEdep);  
  RooCondAddPdf masspdfdataNoEdep("masspdfdataNoEdep","",gauspdfsdataNoEdep,gauscoeffsNoEdep);
  RooCondAddPdf masspdfdatasingleNoEdep("masspdfdatasingleNoEdep","",gauspdfsdataNoEdepsingle,gauscoeffsNoEdep);
  

  //define lists of targets
  RooArgList tgtsNoEdep;
  tgtsNoEdep.add(tgtsmassNoEdep);
  tgtsNoEdep.add(*scaleNoEdep_ph1);
  tgtsNoEdep.add(*scaleNoEdep_ph2);
  tgtsNoEdep.add(*smearNoEdep_ph1);
  tgtsNoEdep.add(*smearNoEdep_ph2);
  
  RooArgList tgtsNoEdepsingle;
  tgtsNoEdepsingle.add(tgtsmassNoEdep);
  tgtsNoEdepsingle.add(*scaleNoEdep_ph1);
  tgtsNoEdepsingle.add(*smearNoEdep_ph1);
   
  
  RooConstVar etermconst("etermconst","",0.);  
   
  RooRealVar r("r","",1.);
  r.setConstant();


  double minweightmc = 1000.;
  double minweightdata = minweightmc*hdataD->sumEntries()/hdata->sumEntries();
  
  std::vector<double> minweights;
  minweights.push_back(100.);
  minweights.push_back(100.);
  
  //ntot.setConstant();

  TFile *fres = new TFile("fres.root","RECREATE");

  
  //run initialization fits to get starting values
  if (1) {
    std::vector<RooAbsData*> vdatainit;
    vdatainit.push_back(hdata);
    
    std::vector<RooAbsReal*> vpdfinit;
    vpdfinit.push_back(&masspdfNoEdep);
    
    RooHybridBDTAutoPdf bdtpdfinit("bdtpdfinit","",tgtsmassNoEdep,etermconst,r,vdatainit,vpdfinit);
    bdtpdfinit.SetPrescaleInit(prescaleinit);    
    bdtpdfinit.TrainForest(0);
    
  }

  if (1) {
    std::vector<RooAbsData*> vdatainit;
    vdatainit.push_back(hdata);
    vdatainit.push_back(hdataD);
    
    std::vector<RooAbsReal*> vpdfinit;
    vpdfinit.push_back(&masspdfNoEdep);
    vpdfinit.push_back(&masspdfdatasingleNoEdep);
    
    RooHybridBDTAutoPdf bdtpdfinit("bdtpdfinit2","",tgtsNoEdepsingle,etermconst,r,vdatainit,vpdfinit);
    bdtpdfinit.SetPrescaleInit(prescaleinit);    
    bdtpdfinit.TrainForest(0);
    
    scaleNoEdepvar_ph2->setVal(scaleNoEdepvar_ph1->getVal());
    scaleNoEdepvar_ph2->setError(scaleNoEdepvar_ph1->getError());
    
    smearNoEdepvar_ph2->setVal(smearNoEdepvar_ph1->getVal());
    smearNoEdepvar_ph2->setError(smearNoEdepvar_ph1->getError());      
    
  }  
  
  //run final simultaneous mass fit
  if (1) {  

    std::vector<RooAbsData*> vdata;
    vdata.push_back(hdata);    
    vdata.push_back(hdataD);        
    
    std::vector<RooAbsReal*> vpdf;
    vpdf.push_back(&masspdfNoEdep);  
    vpdf.push_back(&masspdfdataNoEdep);         
    
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfNoEdep","",tgtsNoEdep,etermconst,r,vdata,vpdf);
    bdtpdfdiff.SetMinCutSignificance(5.);
    bdtpdfdiff.SetDoInitialFit(false);
    bdtpdfdiff.SetShrinkage(0.1);
    bdtpdfdiff.SetMinWeights(minweights);
    bdtpdfdiff.SetMaxNodes(150);
    bdtpdfdiff.TrainForest(1e6);
    
  }      
  

  //fill dataset with corrections (first step corrections which depend on dielectron correlations)
  RooFormulaVar *scaleout_ph1 = new RooFormulaVar("scaleout_ph1","","-log(@0)",RooArgList(*scaleNoEdeplim_ph1));
  
  RooRealVar *scaleoutvar_ph1 = (RooRealVar*)hdataAllD->addColumn(*scaleout_ph1);
  scaleoutvar_ph1->removeRange();
  scaleoutvar_ph1->setConstant(false);

  //define pdfs for variable transformation step
  
  RooFormulaVar *scalemean_ph1 = new RooFormulaVar("scalemean_ph1","","-log(@0)",RooArgList(*scalelim_ph1));
  RooGaussianFast *scalepdf_ph1 = new RooGaussianFast("scalepdf_ph1","",*scaleoutvar_ph1,*scalemean_ph1,*scalewidthlim_ph1);

  RooArgList tgtsscale;
  tgtsscale.add(*scale_ph1);
  tgtsscale.add(*scalewidth_ph1);
     
  
  RooRealVar *smearoutvar_ph1 = (RooRealVar*)hdataAllD->addColumn(*smearNoEdeplim_ph1);
  smearoutvar_ph1->removeRange();
  smearoutvar_ph1->setConstant(false);
  
  RooGaussianFast *smearpdf_ph1 = new RooGaussianFast("smearpdf_ph1","",*smearoutvar_ph1,*smearlim_ph1,*smearwidthlim_ph1);

  RooArgList tgtssmear;
  tgtssmear.add(*smear_ph1);
  tgtssmear.add(*smearwidth_ph1);  

  //run fit for variable transformation on scale output
  if (1) {  
    std::vector<RooAbsData*> vdata;
    vdata.push_back(hdataAllD);    
    
    std::vector<RooAbsReal*> vpdf;
    vpdf.push_back(scalepdf_ph1);  
    
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff2","",tgtsscale,etermconst,r,vdata,vpdf);
    bdtpdfdiff.SetMinCutSignificance(5.);
    bdtpdfdiff.SetPrescaleInit(prescaleinit);    
    bdtpdfdiff.SetDoInitialFit(true);
    bdtpdfdiff.SetShrinkage(0.1);
    bdtpdfdiff.SetMinWeights(minweights);
    bdtpdfdiff.SetMaxNodes(1500);
    bdtpdfdiff.TrainForest(1e6);
  }   
  
  scale_ph2->SetUseFunc(true);
  scalewidth_ph2->SetUseFunc(true);
  

  
  
  RooWorkspace *wereg = new RooWorkspace("weregmass");
  wereg->import(*scalelim_ph1,RecycleConflictNodes());
  wereg->import(*scalelim_ph2,RecycleConflictNodes());
  wereg->import(*scalewidthlim_ph1,RecycleConflictNodes());
  wereg->import(*scalewidthlim_ph2,RecycleConflictNodes());
  wereg->import(*scaleNoEdeplim_ph1,RecycleConflictNodes());
  wereg->import(*smearNoEdeplim_ph1,RecycleConflictNodes());
  wereg->import(*scaleNoEdeplim_ph2,RecycleConflictNodes());
  wereg->import(*smearNoEdeplim_ph2,RecycleConflictNodes());  
  wereg->import(*scalepdf_ph1,RecycleConflictNodes());  
  wereg->import(masspdfNoEdep,RecycleConflictNodes());  
  wereg->import(masspdfdataNoEdep,RecycleConflictNodes());  
  
  wereg->defineSet("vars",vars,true);
  
  wereg->writeToFile("weregmass_scale.root");    
  
  //run fit for variable transformation on smearing output
  if (1) {  
    std::vector<RooAbsData*> vdata;
    vdata.push_back(hdataAllD);    
    
    std::vector<RooAbsReal*> vpdf;
    vpdf.push_back(smearpdf_ph1);  
    
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff3","",tgtssmear,etermconst,r,vdata,vpdf);
    bdtpdfdiff.SetMinCutSignificance(5.);
    bdtpdfdiff.SetPrescaleInit(prescaleinit);    
    bdtpdfdiff.SetDoInitialFit(true);
    bdtpdfdiff.SetShrinkage(0.1);
    bdtpdfdiff.SetMinWeights(minweights);
    bdtpdfdiff.SetMaxNodes(1500);
    bdtpdfdiff.TrainForest(1e6);
  }   
  
  smear_ph2->SetUseFunc(true);
  smearwidth_ph2->SetUseFunc(true);  

  wereg->import(*smearlim_ph1,RecycleConflictNodes());  
  wereg->import(*smearlim_ph2,RecycleConflictNodes());
  wereg->import(*smearwidthlim_ph1,RecycleConflictNodes());
  wereg->import(*smearwidthlim_ph2,RecycleConflictNodes());    
  wereg->import(*smearpdf_ph1,RecycleConflictNodes());  
   
  wereg->writeToFile("weregmass.root");    

  //make some plots
  if (1) {
    TCanvas *cmc = new TCanvas;
    RooPlot *plot = logmass->frame(100);
    hdata->plotOn(plot);
    masspdfNoEdep.plotOn(plot,ProjWData(*hdata));
    plot->Draw();
    cmc->SaveAs("mplotmc.eps");

    TCanvas *cdata = new TCanvas;
    RooPlot *plotD = logmass->frame(100);
    hdataD->plotOn(plotD);
    masspdfdataNoEdep.plotOn(plotD,ProjWData(*hdataD));
    plotD->Draw();
    cdata->SaveAs("mplotdata.eps");
  }
   
  if (1) 
  {    
    TCanvas *cmod = new TCanvas;
    RooPlot *plotmod = scaleoutvar_ph1->frame(-0.02,0.02,200);
    hdataAllD->plotOn(plotmod);
    scalepdf_ph1->plotOn(plotmod,ProjWData(*hdataAllD));
    plotmod->Draw();
    cmod->SaveAs("plotmod1.eps");
  }
  
  if (1)
  {    
    TCanvas *cmod = new TCanvas;
    RooPlot *plotmod = smearoutvar_ph1->frame(0.,pow(-0.02,2),200);
    hdataAllD->plotOn(plotmod);
    smearpdf_ph1->plotOn(plotmod,ProjWData(*hdataAllD));
    plotmod->Draw();
    cmod->SaveAs("plotmod2.eps");
  }  
  
  //RooAbsArg::setDirtyInhibit(false);

  
//   RooFormulaVar *scaledmass = new RooFormulaVar("scaledmass","","sqrt(@0*@1)*@2",RooArgList(*scalelim_ph1,*scalelim_ph2,*mass));
//   RooRealVar *scaledmassvar = (RooRealVar*)hdataD->addColumn(*scaledmass);
  
//   new TCanvas;
//   RooPlot *plotm = mass->frame
  
   
  
//   RooAbsReal *condnll = masspdf.createNLL(*hdata,ConditionalObservables(condvarsmass),NumCPU(16));
//   RooAbsReal *condnllD = masspdfdata.createNLL(*hdataD,ConditionalObservables(condvars),NumCPU(16));
//   
//   printf("condnll = %5f, condnllD = %5f, totalnll = %5f\n",condnll->getVal(),condnllD->getVal(), condnll->getVal()+condnllD->getVal());
  
  
  
  //return;
  
   
   
  
  return;
  
  
}