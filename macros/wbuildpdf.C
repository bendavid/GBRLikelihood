
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
#include "RooHistFunc.h"
#include "RooDataHist.h"


void wbuildpdf(int charge = 1, int weight=0, bool dotight=false) {
 
  printf("dotight = %i\n",int(dotight));
  
  TString suffix;
  
  double scale = 0.01;
  double tgtevents = 0.;
  
  TChain *wtree = new TChain("WTreeProducer");
  if (charge==1) {
    wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_1p1.root");
    
    //effacc = 0.4337
    
    
    
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_1p2.root");
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_1p3.root");
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_1p4.root");
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_2p1.root");    
    
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_2p2.root");
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_2p3.root");
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_2p4.root");
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_3p1.root");    
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_3p2.root");
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_3p3.root");
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_3p4.root");    
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_4p1.root");    
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_4p2.root");
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_4p3.root");
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_4p4.root");    
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_5p1.root");    
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_5p2.root");
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_6p1.root");    
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/WTreeProducer_tree_6p2.root");
    
    tgtevents = 18456.*0.1063/(0.1071+0.1068+0.1138)*0.4337*5000.;
    
    suffix = TString::Format("wplus_%i",weight);
  }
  else if (charge==-1) {
    wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_1p1.root");
    
    //effacc = 0.4105
    
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_1p2.root");
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_1p3.root");
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_1p4.root");
//     wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_1p5.root");
    
    
/*    wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_1p6.root");
    wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_2p1.root");    
    wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_2p2.root");
    wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_2p3.root");
    wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_2p4.root");
    wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_2p5.root");
    wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_3p1.root");    
    wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_3p2.root");
    wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_3p3.root");
    wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_3p4.root");    
    wtree->Add("root://eoscms//store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/WTreeProducer_tree_3p5.root"); */  
    
    tgtevents = 12858.*0.1063/(0.1071+0.1068+0.1138)*0.4105*5000.;

    suffix = TString::Format("wminus_%i",weight);
  }
  
  tgtevents *= 0.2;
  
  int ntotal = wtree->GetEntries();
  
  double fracevents = tgtevents/double(ntotal);
  
  assert(fracevents<=1.);
  
//   printf("nevents = %i\n",int(wtree->GetEntries()));
  
//   return;
  
  RooRealVar *mupt = new RooRealVar("mupt","Mu_pt",0.);
  RooRealVar *mueta = new RooRealVar("mueta","Mu_eta",0.);
  RooRealVar *met = new RooRealVar("met","pfmet",0.);
  RooRealVar *sindphi = new RooRealVar("sindphi","TMath::Sin(pfmet_phi-Mu_phi)",0.);
  //RooRealVar *mt = new RooRealVar("mt","sqrt(2.0*Mu_pt*pfmet*(1.-TMath::Cos(Mu_phi-pfmet_phi)))",0.);
  RooRealVar *mt = new RooRealVar("mt","W_mt",0.,50.,110.);
  mt->setBins(120);
  
  RooArgList vars;
  vars.add(*mupt);
  vars.add(*mueta);
  vars.add(*met);
  vars.add(*sindphi);
  
  RooArgList condvars(vars);
  
//   vars.add(*mt);
  
  TCut loosesel = "(Mu_pt > 25. && abs(Mu_eta)<2.1)";
  TCut tightsel = "(Mu_pt > 30. && W_mt > 50. && W_pt < 25.)";
  
  TCut sel;
  if (dotight) {
    sel = loosesel*tightsel;
  }
  else {
    sel = loosesel;
  }
  
//   TCut selref = "(evt%20==0)";
  
  TCut selref(TString::Format("(evt%1000 < %i)",int(fracevents*1000.)));
  
//   TCut seldata = "(evt%20==1)";
//   TCut sel = "(Mu_pt > 0. && evt%200==0)";
  
  RooRealVar weightvar("weightvar","",1.);
  
//   TFile *fin = TFile::Open(TString::Format("walt_%s.root",suffix.Data()));
//   RooWorkspace *walt = (RooWorkspace*)fin->Get(TString::Format("walt_%s",suffix.Data()));
//   RooDataSet *dsnom = (RooDataSet*)walt->data("dsnom");
//   RooDataSet *dsalt = (RooDataSet*)walt->data("dsalt");

//   TCut seldatafull = sel*seldata;  
//   weightvar.SetTitle(seldatafull);
//   RooDataSet *dsdata = RooTreeConvert::CreateDataSet("dsdata",wtree,vars,weightvar); 
//   
//   TH1D *hdata = new TH1D("hdata","",120,50.,110.);
//   wtree->Draw("W_mt>>hdata",seldatafull,"goff");
  
  
//   if (charge==-1) {
//     double xsecratio = 12858./18456.;
//     int lastevent(xsecratio*dsdata->numEntries());
//     RooAbsData *dsdatareduced = dsdata->reduce(RooFit::EventRange(0,lastevent));
//     delete dsdata;
//     dsdata = (RooDataSet*)dsdatareduced;
//   }
  
  TCut selreffull = sel*selref;
  weightvar.SetTitle(selreffull);
  RooDataSet *dsref = RooTreeConvert::CreateDataSet("dsref",wtree,vars,weightvar); 
  
  TCut selaltfull(selreffull*TCut(TString::Format("LHE_weight[%i]",weight)));
//   selaltfull = selreffull;
  weightvar.SetTitle(selaltfull);
  RooDataSet *dsalt = RooTreeConvert::CreateDataSet("dsalt",wtree,vars,weightvar);
  
  printf("dsrel = %5f, dsalt = %5f\n",dsref->sumEntries(),dsalt->sumEntries());
  
  TH1D *href = new TH1D("href","",120,50.,110.);
  wtree->Draw("W_mt>>href",selreffull,"goff");        
  
  TH1D *halt = new TH1D("halt","",120,50.,110.);
  wtree->Draw("W_mt>>halt",selaltfull,"goff");      
  
  RooRealVar *faltvar = new RooRealVar(TString::Format("faltvar_%s",suffix.Data()),"",0.);
  faltvar->setConstant(false);
  
  RooGBRFunctionFlex *faltfunc = new RooGBRFunctionFlex(TString::Format("faltfunc_%s",suffix.Data()),"");
  RooGBRTargetFlex *falt = new RooGBRTargetFlex(TString::Format("falt_%s",suffix.Data()),"",*faltfunc,*faltvar,condvars);
  
  RooArgList tgts;
  tgts.add(*falt);  

  RooFormulaVar *sigreal = new RooFormulaVar("sigreal","","exp(@0)/(1.+exp(@0))",RooArgList(*falt));
  RooFormulaVar *bkgreal = new RooFormulaVar("bkgreal","","1./(1.+exp(@0))",RooArgList(*falt));
    
  RooConstVar etermconst("etermconst","",0.);  
  RooRealVar r("r","",1.);
  r.setConstant();
  
  //define list of pdfs
  std::vector<RooAbsReal*> vpdf;
  vpdf.push_back(sigreal);  
  vpdf.push_back(bkgreal);  

  //define list of training datasets
  std::vector<RooAbsData*> vdata;
  vdata.push_back(dsalt);
  vdata.push_back(dsref);
  
  double minweight = 2000;
  std::vector<double> minweights;
  minweights.push_back(0);  
  minweights.push_back(minweight); 

  RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff","",tgts,etermconst,r,vdata,vpdf,1);
//   bdtpdfdiff.SetMinCutSignificance(0.1);
  bdtpdfdiff.SetShrinkage(0.1);
  bdtpdfdiff.SetMinWeights(minweights);
  bdtpdfdiff.SetMaxNodes(1000);
  bdtpdfdiff.TrainForest(300);  
  
  RooFormulaVar *expfalt = new RooFormulaVar(TString::Format("expfalt_%s",suffix.Data()),"","exp(@0)",*falt);
  
//   dsdata->addColumn(*expfalt);
  dsref->addColumn(*expfalt);
  
  RooConstVar *normalt = new RooConstVar(TString::Format("normalt_%s",suffix.Data()),"",dsalt->sumEntries());
  RooDataHist *dmtalt = new RooDataHist(TString::Format("dmtalt_%s",suffix.Data()),"",*mt,halt);
  RooHistFunc *hfuncalt = new RooHistFunc(TString::Format("hfuncalt_%s",suffix.Data()),"",*mt,*dmtalt);

  RooConstVar *normref = new RooConstVar(TString::Format("normref_%s",suffix.Data()),"",dsref->sumEntries());
  RooDataHist *dmtref = new RooDataHist(TString::Format("dmtref_%s",suffix.Data()),"",*mt,href);
  RooHistFunc *hfuncref = new RooHistFunc(TString::Format("hfuncref_%s",suffix.Data()),"",*mt,*dmtref);

//   RooFormulaVar *hfuncratio = new RooFormulaVar(TString::Format("hfuncratio_%s",suffix.Data()),"","@0/@1",RooArgList(*hfuncalt,*hfuncref));
//   dmtref->addColumn(*hfuncratio);  
  
  RooWorkspace *walt = new RooWorkspace(TString::Format("walt_%s",suffix.Data()));
  walt->import(*falt);
  walt->import(*hfuncalt);
  walt->import(*hfuncref);
  walt->import(*normalt);
//   walt->import(*dsdata);
  walt->import(*dsref);
  walt->import(*dmtalt);
  walt->import(*dmtref);
  
  TString filename;
  if (dotight) {
    filename = TString::Format("walttight_%s.root",suffix.Data());
  }
  else {
    filename = TString::Format("waltloose_%s.root",suffix.Data());
  }
  walt->writeToFile(filename);     
  
}