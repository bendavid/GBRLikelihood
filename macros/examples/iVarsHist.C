//*****************************************************************************
// Macro to draw and save histograms of "input variables" for the regression 
// procedure 
//*****************************************************************************
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCut.h"
#include <stdio.h>
#include "TSystem.h"
//Using a TChain class in case of a Tree on multiple files
using namespace std;
void iVarsHist(bool dobarrel = true)
{

	TString dirname;
	if (dobarrel)
		dirname = "Barrel_inputs_saturatedAndNotSaturated/";
	else
		dirname = "Endcap_inputs_saturatedAndNotSaturated/";
  	gSystem->mkdir(dirname,true);
  	gSystem->cd(dirname);  

       //Get the Tree
	//TChain *fTree;
	//fTree = new TChain("regressionAnalyzer/promptTree"); //Name 
       //fTree->Add("/afs/cern.ch/user/m/musella/public/forKenza/gam_gam_phys14_v5_regtraining_v3.root");
	   //fTree->Add("/afs/cern.ch/user/m/mdonega/public/forKenza/output_SinglePhoton_FlatPt-300To3000_numEvent100.root");
	   //fTree->Add("/tmp/mdonega/output_SinglePhoton_FlatPt-300To3000.root");
	TChain* fTree = new TChain("promptTree");
    fTree->Add("/afs/cern.ch/work/k/khoumani/CMSSW_7_2_3/src/HiggsAnalysis/GBRLikelihood/macros/examples/promptTree_FlatPt-300To3000.root");
	fTree->Add("/afs/cern.ch/work/k/khoumani/CMSSW_7_2_3/src/HiggsAnalysis/GBRLikelihood/macros/examples/promptTree_FlatPt-5To300_lessstat.root");
       //______________________________________________________________________
       // Drawing histograms
       //______________________________________________________________________
       	
       //Print number of overflows and underflows in statBox
	gStyle->SetOptStat(1111110); 

	
	//Define Canvas properties
	Double_t canvas_width(1300);
	Double_t canvas_height(800);
	   
       //Get weights
       //Float_t sum_of_weights = 0., value_of_weight;
       //Compute sum of weights
       //Long64_t number_of_entries = fTree->GetEntries();
       //fTree->SetBranchAddress("weight",&value_of_weight);
       //for(Long64_t ii=0; ii<number_of_entries; ii++) {
         //     fTree->GetEntry(ii);
           //   sum_of_weights +=value_of_weight;
       //}
       //fTree->SetAlias("nEntries","This->GetEntries()");
       //TString str = "";
       //str.Form("%f",fTree->GetEntries()*1./sum_of_weights);


       //TString cut_expression = "( abs(scEta)<=1.5) * nEntries*weight* 1./"+str;
       //TString cut_expression = "( abs(scEta)<=1.5) *weight *"+str;
       //TCut selweight(cut_expression);
       
       //int gg = (fTree->GetEntries());
       //printf("sum_of_weights: %f\n", gg *1./ sum_of_weights );
       //printf("number_of_entries: %d\n", gg );

       
       //name of images
       TString add_name;
       TCut selection;
       if (dobarrel)
              add_name = "_eb_gam_gam.png";
       else
              add_name = "_ee_gam_gam.png";

       //Cut on events
       if (dobarrel)
              selection = "(etrue*pt/energy >=200. && etrue*pt/energy<=3000 && abs(scEta)<1.5 )";
       else
              selection = "(etrue*pt/energy >=200. && etrue*pt/energy<=3000 && abs(scEta)>1.5 )";
       

// --------- Energy 
       // Canvas 1: Energy_Histo
       TCanvas *canvas_energy = new TCanvas("Energy_Histo","Histograms of energy variables",canvas_width,canvas_height);
       TPad *padTop1 = new TPad("padTop1","padTop",0.05,0.5,0.95,0.95,256);
       TPad *padBottom1 = new TPad("padBottom1","padBottom",0.05,0.05,0.95,0.45,256);
       padTop1->Divide(3,1);
       padBottom1->Divide(3,1);

       padTop1->Draw();
       padBottom1->Draw();

       padTop1->cd(1);
       fTree->Draw("scRawEnergy",selection);
       padTop1->cd(2);
       fTree->Draw("etrue",selection);
       padTop1->cd(3);
       fTree->Draw("scSeedEnergy",selection);  
       padBottom1->cd(1);
       fTree->Draw("pt",selection);
       padBottom1->cd(2);
       fTree->Draw("scPreshowerEnergy",selection);
       padBottom1->cd(3);
       fTree->Draw("hadTowOverEm",selection);
      //TH1F *histo = (TH1F*)gDirectory->Get("htemp");
       canvas_energy->SaveAs("energy"+add_name);       
      
       // --------- Phi and Eta
       // Canvas 2: Phi_and_Eta_Histo
       TCanvas *canvas_phi_eta = new TCanvas("phi_eta_Histo","Histograms of phi and eta variables",canvas_width,canvas_height);
       TPad *padTop2 = new TPad("padTop2","padTop",0.05,0.5,0.95,0.95,256);
       TPad *padBottom2 = new TPad("padBottom2","padBottom",0.05,0.05,0.95,0.45,256);
       padTop2->Divide(4,1);
       padBottom2->Divide(4,1);

       padTop2->Draw();
       padBottom2->Draw();

       padTop2->cd(1);
       fTree->Draw("scEta",selection);
       padTop2->cd(2);
       fTree->Draw("eta",selection);
       padTop2->cd(3);
       fTree->Draw("etaWidth",selection);  
       padTop2->cd(4);
       fTree->Draw("cryEta",selection);

       padBottom2->cd(1);
       fTree->Draw("scPhy",selection);
       padBottom2->cd(2);
       fTree->Draw("phi",selection);
       padBottom2->cd(3);
       fTree->Draw("phiWidth",selection);
       padBottom2->cd(4);
       fTree->Draw("cryPhi",selection);
       canvas_phi_eta->SaveAs("phi_eta"+add_name);

       // --------- IPhi and IEta
       // Canvas 3: IPhi_and_IEta_Histo
       TCanvas *canvas_iphi_ieta = new TCanvas("iphi_ieta_Histo","Histograms of iphi and ieta variables",canvas_width,canvas_height);
       TPad *padTop3 = new TPad("padTop3","padTop",0.05,0.5,0.95,0.95,256);
       TPad *padBottom3 = new TPad("padBottom3","padBottom",0.05,0.05,0.95,0.45,256);
       padTop3->Divide(2,1);
       padBottom3->Divide(3,1);

       padTop3->Draw();
       padBottom3->Draw();

       padTop3->cd(1);
       fTree->Draw("iEta",selection);
       padTop3->cd(2);
       fTree->Draw("iPhi",selection);
       
       padBottom3->cd(1);
       fTree->Draw("sigmaIphiIphi",selection);
       padBottom3->cd(2);
       fTree->Draw("sigmaIetaIeta",selection);
       padBottom3->cd(3);
       fTree->Draw("covarianceIetaIphi",selection);
       canvas_iphi_ieta->SaveAs("iphi_ieta"+add_name);

       // --------- rho, R9, csize,nvtx : RRCV
       // Canvas 4: rho, R9, Cluster size, and number of vertices
       TCanvas *canvas_RRCV = new TCanvas("CVR_Histo","Histograms of rho, R9, Cluster size, and number of vertices",canvas_width,canvas_height);
       TPad *padTop4 = new TPad("padTop4","padTop",0.05,0.5,0.95,0.95,256);
       TPad *padBottom4 = new TPad("padBottom4","padBottom",0.05,0.05,0.95,0.45,256);
       padTop4->Divide(3,1);
       padBottom4->Divide(2,1);

       padTop4->Draw();
       padBottom4->Draw();

       padTop4->cd(1);
       fTree->Draw("rho");
       padTop4->cd(2);
       fTree->Draw("r9",selection);
       padTop4->cd(3);
       fTree->Draw("full5x5_r9",selection);

       padBottom4->cd(1);
       fTree->Draw("scClustersSize",selection);
       padBottom4->cd(2);
       fTree->Draw("nvtx");
       canvas_RRCV->SaveAs("rho_R9_csize_nvtx"+add_name);

       // --------- Energy2 
       // Canvas 5: eMax, e2nd, eTop, eBottom, eLeft, eRight

       TCanvas *canvas_energy2 = new TCanvas("Energy2_Histo","Histograms of energy variables (Max, right, Left...)",canvas_width,canvas_height);
       TPad *padTop5 = new TPad("padTop5","padTop",0.05,0.5,0.95,0.95,256);
       TPad *padBottom5 = new TPad("padBottom5","padBottom",0.05,0.05,0.95,0.45,256);
       padTop5->Divide(3,1);
       padBottom5->Divide(3,1);

       padTop5->Draw();
       padBottom5->Draw();

       padTop5->cd(1);
       fTree->Draw("eMax",selection);
       padTop5->cd(2);
       fTree->Draw("e2nd",selection);
       padTop5->cd(3);
       fTree->Draw("eTop",selection); 
       padBottom5->cd(1);
       fTree->Draw("eBottom",selection);
       padBottom5->cd(2);
       fTree->Draw("eLeft",selection);
       padBottom5->cd(3);
       fTree->Draw("eRight",selection);
       canvas_energy2->SaveAs("eMax_2nd_Top_Bottom_Left_Right"+add_name);


       // --------- Energy3 
       // Canvas 6: e2x5Max, e2x5Top, e2x5Bottom, e2x5Left, e2x5Right

       TCanvas *canvas_energy3 = new TCanvas("Energy3_Histo","Histograms of e2x5 variables (Max, right, Left...)",canvas_width,canvas_height);
       TPad *padTop6 = new TPad("padTop6","padTop",0.05,0.5,0.95,0.95,256);
       TPad *padBottom6 = new TPad("padBottom6","padBottom",0.05,0.05,0.95,0.45,256);
       padTop6->Divide(3,1);
       padBottom6->Divide(2,1);

       padTop6->Draw();
       padBottom6->Draw();

       padTop6->cd(1);
       fTree->Draw("e2x5Max",selection);
       padTop6->cd(2);
       fTree->Draw("e2x5Top",selection);
       padTop6->cd(3);
       fTree->Draw("e2x5Bottom",selection); 
       padBottom6->cd(1);
       fTree->Draw("e2x5Left",selection);
       padBottom6->cd(2);
       fTree->Draw("e2x5Right",selection);
       canvas_energy3->SaveAs("e2x5Vars"+add_name);

       // --------- Comparison with full5x5  
       // Canvas 7: full5x5_e3x3, full5x5_e5x5, full5x5_sigmaIetaIeta


       TCanvas *canvas_comparison = new TCanvas("Comparison_with_full5x5_Histo","Histograms of full5x5 variables",canvas_width,canvas_height);
       TPad *padTop7 = new TPad("padTop7","padTop",0.05,0.5,0.95,0.95,256);
       TPad *padBottom7 = new TPad("padBottom7","padBottom",0.05,0.05,0.95,0.45,256);
       padTop7->Divide(3,1);
       padBottom7->Divide(3,1);

       padTop7->Draw();
       padBottom7->Draw();

       padTop7->cd(1);
       fTree->Draw("full5x5_e3x3",selection);
       padTop7->cd(2);
       fTree->Draw("full5x5_e5x5",selection);
       padTop7->cd(3);
       fTree->Draw("full5x5_sigmaIetaIeta",selection);  
       padBottom7->cd(1);
       fTree->Draw("e3x3",selection);
       padBottom7->cd(2);
       fTree->Draw("e5x5",selection);
       padBottom7->cd(3);
       fTree->Draw("sigmaIetaIeta",selection);
       canvas_comparison->SaveAs("comparison_full5x5"+add_name);

       // --------- weight and Luminosity  
       // Canvas 8: weight and luminosity (weight_lumi)


  //     TCanvas *canvas_weight_lumi = new TCanvas("weight_lumi_Histo","Histograms of weight and luminosity variables",canvas_width,canvas_height);
  //     canvas_weight_lumi->Divide(2,1);
  //
  //     canvas_weight_lumi->cd(1);
  //     fTree->Draw("weight",selection);
  //     canvas_weight_lumi->cd(2);
  //     fTree->Draw("lumi",selection);
  //     canvas_weight_lumi->SaveAs("weight_lumi"+add_name); 

//_____________________________
 /*      TString add_name("_restricted_range.png");
       // --------- Energy 
       // Canvas 1: Energy_Histo
       TCanvas *canvas_energy = new TCanvas("Energy_Histo","Histograms of energy variables",canvas_width,canvas_height);
       TPad *padTop1 = new TPad("padTop1","padTop",0.05,0.5,0.95,0.95,256);
       TPad *padBottom1 = new TPad("padBottom1","padBottom",0.05,0.05,0.95,0.45,256);
       padTop1->Divide(3,1);
       padBottom1->Divide(3,1);

       padTop1->Draw();
       padBottom1->Draw();

       padTop1->cd(1);
       fTree->Draw("scRawEnergy","weight* (scRawEnergy<1000 && abs(scEta)<=1.5)");
       padTop1->cd(2);
       fTree->Draw("etrue","weight* (etrue<1000 && abs(scEta)<=1.5)");
       padTop1->cd(3);
       fTree->Draw("scSeedEnergy","weight* (scSeedEnergy<1000 && abs(scEta)<=1.5)");  
       padBottom1->cd(1);
       fTree->Draw("pt","weight* (pt<1000 && abs(scEta)<=1.5)");
       padBottom1->cd(2);
       fTree->Draw("scPreshowerEnergy","weight* (scPreshowerEnergy<100 && abs(scEta)<=1.5)");
       padBottom1->cd(3);
       fTree->Draw("hadTowOverEm","weight* (hadTowOverEm<100 && abs(scEta)<=1.5)");
       //TH1F *histo = (TH1F*)gDirectory->Get("htemp");
       canvas_energy->SaveAs("energy"+add_name);       
       
       // --------- Phi and Eta
       // Canvas 2: Phi_and_Eta_Histo
       TCanvas *canvas_phi_eta = new TCanvas("phi_eta_Histo","Histograms of phi and eta variables",canvas_width,canvas_height);
       TPad *padTop2 = new TPad("padTop2","padTop",0.05,0.5,0.95,0.95,256);
       TPad *padBottom2 = new TPad("padBottom2","padBottom",0.05,0.05,0.95,0.45,256);
       padTop2->Divide(4,1);
       padBottom2->Divide(4,1);

       padTop2->Draw();
       padBottom2->Draw();

       padTop2->cd(1);
       fTree->Draw("scEta","weight* (abs(scEta)<=1.5)");
       padTop2->cd(2);
       fTree->Draw("eta","weight* (abs(scEta)<=1.5)");
       padTop2->cd(3);
       fTree->Draw("etaWidth","weight* (etaWidth<0.04 && abs(scEta)<=1.5)");  
       padTop2->cd(4);
       fTree->Draw("cryEta","weight* (abs(scEta)<=1.5)");

       padBottom2->cd(1);
       fTree->Draw("scPhy","weight* (abs(scEta)<=1.5)");
       padBottom2->cd(2);
       fTree->Draw("phi","weight* (abs(scEta)<=1.5)");
       padBottom2->cd(3);
       fTree->Draw("phiWidth","weight* (phiWidth<0.1 && abs(scEta)<=1.5)");
       padBottom2->cd(4);
       fTree->Draw("cryPhi","weight* (abs(scEta)<=1.5)");
       canvas_phi_eta->SaveAs("phi_eta"+add_name);

       // --------- IPhi and IEta
       // Canvas 3: IPhi_and_IEta_Histo
       TCanvas *canvas_iphi_ieta = new TCanvas("iphi_ieta_Histo","Histograms of iphi and ieta variables",canvas_width,canvas_height);
       TPad *padTop3 = new TPad("padTop3","padTop",0.05,0.5,0.95,0.95,256);
       TPad *padBottom3 = new TPad("padBottom3","padBottom",0.05,0.05,0.95,0.45,256);
       padTop3->Divide(2,1);
       padBottom3->Divide(3,1);

       padTop3->Draw();
       padBottom3->Draw();

       padTop3->cd(1);
       fTree->Draw("iEta","weight* (abs(scEta)<=1.5)");
       padTop3->cd(2);
       fTree->Draw("iPhi","weight* (abs(scEta)<=1.5)");
       
       padBottom3->cd(1);
       fTree->Draw("sigmaIphiIphi","weight* (sigmaIphiIphi<0.35 && abs(scEta)<=1.5)");
       padBottom3->cd(2);
       fTree->Draw("sigmaIetaIeta","weight* (abs(scEta)<=1.5)");
       padBottom3->cd(3);
       fTree->Draw("covarianceIetaIphi","weight* (covarianceIetaIphi<5 && covarianceIetaIphi>-10 && abs(scEta)<=1.5)");
       canvas_iphi_ieta->SaveAs("iphi_ieta"+add_name);

       // --------- rho, R9, csize,nvtx : RRCV
       // Canvas 4: rho, R9, Cluster size, and number of vertices
       TCanvas *canvas_RRCV = new TCanvas("CVR_Histo","Histograms of rho, R9, Cluster size, and number of vertices",canvas_width,canvas_height);
       TPad *padTop4 = new TPad("padTop4","padTop",0.05,0.5,0.95,0.95,256);
       TPad *padBottom4 = new TPad("padBottom4","padBottom",0.05,0.05,0.95,0.45,256);
       padTop4->Divide(3,1);
       padBottom4->Divide(2,1);

       padTop4->Draw();
       padBottom4->Draw();

       padTop4->cd(1);
       fTree->Draw("rho");
       padTop4->cd(2);
       fTree->Draw("r9","weight* (abs(scEta)<=1.5 && r9<1.05 && r9>0.65)");
       padTop4->cd(3);
       fTree->Draw("full5x5_r9","weight* (abs(scEta)<=1.5 && full5x5_r9<1.05 && full5x5_r9>0.65)");

       padBottom4->cd(1);
       fTree->Draw("scClustersSize","weight* (abs(scEta)<=1.5 && scClustersSize<15)");
       padBottom4->cd(2);
       fTree->Draw("nvtx");
       canvas_RRCV->SaveAs("rho_R9_csize_nvtx"+add_name);

       // --------- Energy2 
       // Canvas 5: eMax, e2nd, eTop, eBottom, eLeft, eRight

       TCanvas *canvas_energy2 = new TCanvas("Energy2_Histo","Histograms of energy variables (Max, right, Left...)",canvas_width,canvas_height);
       TPad *padTop5 = new TPad("padTop5","padTop",0.05,0.5,0.95,0.95,256);
       TPad *padBottom5 = new TPad("padBottom5","padBottom",0.05,0.05,0.95,0.45,256);
       padTop5->Divide(3,1);
       padBottom5->Divide(3,1);

       padTop5->Draw();
       padBottom5->Draw();

       padTop5->cd(1);
       fTree->Draw("eMax","weight* (abs(scEta)<=1.5 && eMax<1000)");
       padTop5->cd(2);
       fTree->Draw("e2nd","weight* (abs(scEta)<=1.5 && e2nd<1000)");
       padTop5->cd(3);
       fTree->Draw("eTop","weight* (abs(scEta)<=1.5 && eTop<1000)"); 
       padBottom5->cd(1);
       fTree->Draw("eBottom","weight* (abs(scEta)<=1.5 && eBottom<1000)");
       padBottom5->cd(2);
       fTree->Draw("eLeft","weight* (abs(scEta)<=1.5 && eLeft<1000)");
       padBottom5->cd(3);
       fTree->Draw("eRight","weight* (abs(scEta)<=1.5 && eRight<1000)");
       canvas_energy2->SaveAs("eMax_2nd_Top_Bottom_Left_Right"+add_name);


       // --------- Energy3 
       // Canvas 6: e2x5Max, e2x5Top, e2x5Bottom, e2x5Left, e2x5Right

       TCanvas *canvas_energy3 = new TCanvas("Energy3_Histo","Histograms of e2x5 variables (Max, right, Left...)",canvas_width,canvas_height);
       TPad *padTop6 = new TPad("padTop6","padTop",0.05,0.5,0.95,0.95,256);
       TPad *padBottom6 = new TPad("padBottom6","padBottom",0.05,0.05,0.95,0.45,256);
       padTop6->Divide(3,1);
       padBottom6->Divide(2,1);

       padTop6->Draw();
       padBottom6->Draw();

       padTop6->cd(1);
       fTree->Draw("e2x5Max","weight* (abs(scEta)<=1.5 && e2x5Max<1000)");
       padTop6->cd(2);
       fTree->Draw("e2x5Top","weight* (abs(scEta)<=1.5 && e2x5Top<1000)");
       padTop6->cd(3);
       fTree->Draw("e2x5Bottom","weight* (abs(scEta)<=1.5 && e2x5Bottom<1000)"); 
       padBottom6->cd(1);
       fTree->Draw("e2x5Left","weight* (abs(scEta)<=1.5 && e2x5Left<1000)");
       padBottom6->cd(2);
       fTree->Draw("e2x5Right","weight* (abs(scEta)<=1.5 && e2x5Right<1000)");
       canvas_energy3->SaveAs("e2x5Vars"+add_name);

       // --------- Comparison with full5x5  
       // Canvas 7: full5x5_e3x3, full5x5_e5x5, full5x5_sigmaIetaIeta


       TCanvas *canvas_comparison = new TCanvas("Comparison_with_full5x5_Histo","Histograms of full5x5 variables",canvas_width,canvas_height);
       TPad *padTop7 = new TPad("padTop7","padTop",0.05,0.5,0.95,0.95,256);
       TPad *padBottom7 = new TPad("padBottom7","padBottom",0.05,0.05,0.95,0.45,256);
       padTop7->Divide(3,1);
       padBottom7->Divide(3,1);

       padTop7->Draw();
       padBottom7->Draw();

       padTop7->cd(1);
       fTree->Draw("full5x5_e3x3","weight* (abs(scEta)<=1.5 && full5x5_e3x3<1000)");
       padTop7->cd(2);
       fTree->Draw("full5x5_e5x5","weight* (abs(scEta)<=1.5 && full5x5_e5x5<1000)");
       padTop7->cd(3);
       fTree->Draw("full5x5_sigmaIetaIeta","weight* (abs(scEta)<=1.5 && full5x5_sigmaIetaIeta<0.35)");  
       padBottom7->cd(1);
       fTree->Draw("e3x3","weight* (abs(scEta)<=1.5 && e3x3<1000)");
       padBottom7->cd(2);
       fTree->Draw("e5x5","weight* (abs(scEta)<=1.5 && e5x5<1000)");
       padBottom7->cd(3);
       fTree->Draw("sigmaIetaIeta","weight* (abs(scEta)<=1.5 && sigmaIetaIeta<0.35)");
       canvas_comparison->SaveAs("comparison_full5x5"+add_name);

       // --------- weight and Luminosity  
       // Canvas 8: weight and luminosity (weight_lumi)


       TCanvas *canvas_weight_lumi = new TCanvas("weight_lumi_Histo","Histograms of weight and luminosity variables",canvas_width,canvas_height);
       canvas_weight_lumi->Divide(2,1);

       canvas_weight_lumi->cd(1);
       fTree->Draw("weight","weight* (abs(scEta)<=1.5)");
       canvas_weight_lumi->cd(2);
       fTree->Draw("lumi","weight* (abs(scEta)<=1.5)");
       canvas_weight_lumi->SaveAs("weight_lumi"+add_name); */

}
