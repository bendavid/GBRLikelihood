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
#include "TH2.h"
//Using a TChain class in case of a Tree on multiple files
using namespace std;
void hist_iEta_iPhi(bool dobarrel = true)
{

       TString dirname;
       if (dobarrel)
              dirname = "barrel/";
       else
              dirname = "endcap/";
  	gSystem->mkdir(dirname,true);
  	gSystem->cd(dirname);  

       //Get the Tree
       TChain *tree;
       tree = new TChain("promptTree"); //Name 
       tree->Add("/afs/cern.ch/user/m/musella/public/forKenza/gam_gam_phys14_v5_regtraining_v3.root");


       //Print number of overflows and underflows in statBox
       gStyle->SetOptStat(1111110); 


       //Define Canvas properties
	Double_t canvas_width(1300);
	Double_t canvas_height(800);

       TString add_name;
       TCut selection1, selection2;
       if (dobarrel)
              add_name = "_eb_gam_gam.png";
       else
              add_name = "_ee_gam_gam.png";

       //Cut on events
       if (dobarrel) {
              //selection = "(abs(scEta)<=1.5 && kSaturated[12]!=1 )";
              selection1 = "(abs(scEta)<=1.5 && kSaturated[12]!=1 && r9 >= 0.94)";
              selection2 = "(abs(scEta)<=1.5 && kSaturated[12]!=1 && r9 < 0.94)";
       }
       else {
              //selection = "(abs(scEta)>1.5 && kSaturated[12]!=1 )";
       }

       TCanvas *canvas_iphi_ieta = new TCanvas("iphi_ieta_Histo","",canvas_width,canvas_height);
       canvas_iphi_ieta->Divide(1,2);
       canvas_iphi_ieta->cd(1);
       gPad->SetGridx();
       gPad->SetGridy();
       //TH2 *hist = new TH2();
       tree->Draw("iPhi:iEta >> hist1(171,-85.5,85.5,361,-0.5,360.5)",selection1,"colz");
       canvas_iphi_ieta->cd(2);
       gPad->SetGridx();
       gPad->SetGridy();
       //TH2 *hist = new TH2();
       tree->Draw("iPhi:iEta >> hist2(171,-85.5,85.5,361,-0.5,360.5)",selection2,"colz");

}




