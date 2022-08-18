#define monopho17MVA_cxx
#include "monopho17MVA.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include <fstream>


/*
**** 11182021 cutbased medium ID ****
sieie < 0.0105
H/E < 0.05
HCAL iso < 6.40
Track iso < 0.89
ECAL iso < 4.78
*/

using namespace std;


Double_t EcalEA(Double_t eta){
    
    Float_t EffectiveArea = 0.0;

    if(fabs(eta) >= 0   && fabs(eta) < 0.5 ) EffectiveArea = 0.108;
    if(fabs(eta) >= 0.5   && fabs(eta) < 1.0 ) EffectiveArea = 0.110;
    if(fabs(eta) >= 1.0   && fabs(eta) < 1.4442 ) EffectiveArea = 0.096;
        
    return EffectiveArea;
}


Double_t EcalEA_ptscale(Double_t eta, Double_t pt){
    Float_t EffectiveArea = 0.0;

    if(fabs(eta) >= 0 && fabs(eta) < 1.4442 ) EffectiveArea = 0.002948 * pt;

    return EffectiveArea;
}


Double_t HcalEA(Double_t eta){
    Float_t EffectiveArea = 0.0;
    
    if(fabs(eta) >= 0   && fabs(eta) < 0.5 ) EffectiveArea = 0.069;
    if(fabs(eta) >= 0.5   && fabs(eta) < 1.0 ) EffectiveArea = 0.107;
    if(fabs(eta) >= 1.0   && fabs(eta) < 1.4442 ) EffectiveArea = 0.142;

    return EffectiveArea;
}

Double_t HcalEA_ptscale(Double_t eta, Double_t pt){
    Float_t EffectiveArea = 0.0;

    if(fabs(eta) >=0 && fabs(eta) < 1.4442  ) EffectiveArea =  0.0121 * pt + 0.00001335 * pt * pt;
    
    return EffectiveArea;

}

Double_t TkrEffAreas(Double_t eta){
    Float_t EffectiveArea = 0.0;
    
    if(fabs(eta) >= 0   && fabs(eta) < 0.5 ) EffectiveArea = 0.043;
    if(fabs(eta) >= 0.5   && fabs(eta) < 1.0 ) EffectiveArea = 0.040;
    if(fabs(eta) >= 1.0   && fabs(eta) < 1.4442 ) EffectiveArea = 0.038;
    
    return EffectiveArea;
}

//if there are only oot || if they don't match
Float_t newuMET(Float_t met, Float_t metphi, Float_t phi, Float_t et)
{
	Float_t ophoEt = et;
	Float_t ophoPhi = phi;
	Float_t pfMET = met;
	Float_t pfMETPhi = metphi;

	Float_t pfMETX = met * cos(pfMETPhi);
	Float_t pfMETY = met * sin(pfMETPhi);
	Float_t OOTEtX = et * cos(ophoPhi);
	Float_t OOTEtY = et * sin(ophoPhi);

	pfMETX -= OOTEtX;
	pfMETY -= OOTEtY;

	Float_t newuMET = sqrt(pfMETX * pfMETX + pfMETY * pfMETY);

	return newuMET;
}

//if they are both and they do match
Float_t newmMET(Float_t met, Float_t metphi, Float_t phi, Float_t et)
{
	Float_t phoEt;
	Float_t phoPhi;
	Float_t pfMET = met;
	Float_t pfMETPhi = metphi;
	Float_t ophoEt = et;
	Float_t ophoPhi = phi;

	Float_t pfMETX = met * cos(pfMETPhi);
	Float_t pfMETY = met * sin(pfMETPhi);
	Float_t ITEtX = phoEt * cos(phoPhi);
	Float_t ITEtY = phoEt * sin(phoPhi);

	pfMETX += ITEtX;
	pfMETY += ITEtY;

	Float_t OOTEtX = et * cos(ophoPhi);
	Float_t OOTEtY = et * sin(ophoPhi);

	pfMETX -= OOTEtX;
	pfMETY -= OOTEtY;

	Float_t newmMET = sqrt(pfMETX * pfMETX + pfMETY * pfMETY);

	return newmMET;
}


//for checking matching

Float_t DeltaR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2)
{

	Float_t ophoPhi = phi1;
	Float_t phoPhi = phi2;
	Float_t dphi = fabs(phoPhi - ophoPhi);
	Float_t tdphi = dphi;
    //phi wrap
	if (dphi > TMath::Pi()) tdphi = TMath::Pi()*2.0 - dphi;
	dphi = tdphi;

	Float_t ophoEta = eta1;
	Float_t phoEta = eta2;
	Float_t deta = fabs(phoEta - ophoEta);

	Float_t deltaR = sqrt(deta*deta + dphi * dphi);
	return deltaR;
}


//calulate invariant mass Z

Float_t InvariMass(Float_t Et1, Float_t Et2, Float_t Phi1, Float_t Phi2, Float_t Eta1, Float_t Eta2)
{
	Float_t Theta1 = 2 * atan(exp(-1.*Eta1));
	Float_t Theta2 = 2 * atan(exp(-1.*Eta2));
	Float_t phoPhi1 = Phi1;
	Float_t phoPhi2 = Phi2;
	Float_t Etot1 = Et1 / sin(Theta1); //E_tot1
	Float_t Etot2 = Et2 / sin(Theta2); //E_tot2
	
	//reconstruct the vectors for x, y and z

	
	Float_t phoX1 = Etot1 * cos(Phi1) * sin(Theta1); 
	Float_t phoY1 = Etot1 * sin(Phi1) * sin(Theta1);
	Float_t phoZ1 = Etot1 * cos(Theta1);

	Float_t phoX2 = Etot2 * cos(Phi2) * sin(Theta2);
	Float_t phoY2 = Etot2 * sin(Phi2) * sin(Theta2);
	Float_t phoZ2 = Etot2 * cos(Theta2);

	Float_t EX1 = Et1*sin(Phi1);
	Float_t EY1 = Et1*cos(Phi1);
	Float_t EZ1 = Etot1*cos(Theta1);

	Float_t EX2 = Et2*sin(Phi2);
	Float_t EY2 = Et2*cos(Phi2);
	Float_t EZ2 = Etot2*cos(Theta2);

	Float_t E1 = sqrt(Etot1*Etot1);
	Float_t E2 = sqrt(Etot2*Etot2);

	//Float_t InvariMassSq = 2 * E1*E2 - phoMag1 * phoMag1 - phoMag2 * phoMag2 - 2 * DotProd12;

	//1	
	Float_t InvariMassSq = (Etot1 + Etot2)*(Etot1 + Etot2) - (EX1 + EX2)*(EX1 + EX2) - (EY1 + EY2)*(EY1 + EY2) - (EZ1 + EZ2)*(EZ1 + EZ2);

	Float_t InvMass = sqrt(InvariMassSq);

	
	
	return InvMass;
	
}


void monopho17MVA::Loop()
{
//   In a ROOT session, you can do:
//      root> .L monopho17MVA.C
//      root> monopho17MVA t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
    TFile *monophoseedtime = new TFile("monopho17UL_cutbased_1118B.root", "RECREATE");


    //**********************CANDIDATE**********************
    TH1F *ewing_candidate = new TH1F("ewing_candidate", "ewing_candidate", 2000, 0, 1.0);
	ewing_candidate->GetXaxis()->SetTitle("ewing");
	ewing_candidate->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_candidate_iso = new TH1F("ewing_candidate_iso", "ewing_candidate_iso", 2000, 0, 1.0);
	ewing_candidate_iso->GetXaxis()->SetTitle("ewing");
	ewing_candidate_iso->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_candidate_3ns_iso = new TH1F("ewing_candidate_3ns_iso", "ewing_candidate_3ns_iso", 2000, 0, 1.0);
	ewing_candidate_3ns_iso->GetXaxis()->SetTitle("ewing");
	ewing_candidate_3ns_iso->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_candidate_ss = new TH1F("ewing_candidate_ss", "ewing_candidate_ss", 2000, 0, 1.0);
	ewing_candidate_ss->GetXaxis()->SetTitle("ewing");
	ewing_candidate_ss->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_candidate_ss_iso = new TH1F("ewing_candidate_ss_iso", "ewing_candidate_ss_iso", 2000, 0, 1.0);
	ewing_candidate_ss_iso->GetXaxis()->SetTitle("ewing");
	ewing_candidate_ss_iso->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_candidate_ss_3ns_iso = new TH1F("ewing_candidate_ss_3ns_iso", "ewing_candidate_ss_3ns_iso", 2000, 0, 1.0);
	ewing_candidate_ss_3ns_iso->GetXaxis()->SetTitle("ewing");
	ewing_candidate_ss_3ns_iso->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_candidate_HCAL = new TH1F("ewing_candidate_HCAL", "ewing_candidate_HCAL", 2000, 0, 1.0);
	ewing_candidate_HCAL->GetXaxis()->SetTitle("ewing");
	ewing_candidate_HCAL->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_candidate_ECAL = new TH1F("ewing_candidate_ECAL", "ewing_candidate_ECAL", 2000, 0, 1.0);
	ewing_candidate_ECAL->GetXaxis()->SetTitle("ewing");
	ewing_candidate_ECAL->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_candidate_track = new TH1F("ewing_candidate_track", "ewing_candidate_track", 2000, 0, 1.0);
	ewing_candidate_track->GetXaxis()->SetTitle("ewing");
	ewing_candidate_track->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_candidate_ss_HCAL = new TH1F("ewing_candidate_ss_HCAL", "ewing_candidate_ss_HCAL", 2000, 0, 1.0);
	ewing_candidate_ss_HCAL->GetXaxis()->SetTitle("ewing");
	ewing_candidate_ss_HCAL->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_candidate_ss_ECAL = new TH1F("ewing_candidate_ss_ECAL", "ewing_candidate_ss_ECAL", 2000, 0, 1.0);
	ewing_candidate_ss_ECAL->GetXaxis()->SetTitle("ewing");
	ewing_candidate_ss_ECAL->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_candidate_ss_track = new TH1F("ewing_candidate_ss_track", "ewing_candidate_ss_track", 2000, 0, 1.0);
	ewing_candidate_ss_track->GetXaxis()->SetTitle("ewing");
	ewing_candidate_ss_track->GetYaxis()->SetTitle("Entries");





    TH1F *timing_candidate_iso = new TH1F("timing_candidate_iso", "timing_candidate_iso", 200, -25, 25);
	timing_candidate_iso->GetXaxis()->SetTitle("Timing");
	timing_candidate_iso->GetYaxis()->SetTitle("Entries");

    TH1F *timing_candidate_iso_ewing1 = new TH1F("timing_candidate_iso_ewing1", "timing_candidate_iso_ewing1", 200, -25, 25);
	timing_candidate_iso_ewing1->GetXaxis()->SetTitle("Timing");
	timing_candidate_iso_ewing1->GetYaxis()->SetTitle("Entries");

    TH1F *timing_candidate_3ns_iso = new TH1F("timing_candidate_3ns_iso", "timing_candidate_3ns_iso", 200, -25, 25);
	timing_candidate_3ns_iso->GetXaxis()->SetTitle("Timing");
	timing_candidate_3ns_iso->GetYaxis()->SetTitle("Entries");

    TH1F *timing_candidate_3ns_iso_ewing1 = new TH1F("timing_candidate_3ns_iso_ewing1", "timing_candidate_3ns_iso_ewing1", 200, -25, 25);
	timing_candidate_3ns_iso_ewing1->GetXaxis()->SetTitle("Timing");
	timing_candidate_3ns_iso_ewing1->GetYaxis()->SetTitle("Entries");

    TH1F *timing_candidate_3ns_iso_ewing2 = new TH1F("timing_candidate_3ns_iso_ewing2", "timing_candidate_3ns_iso_ewing2", 200, -25, 25);
	timing_candidate_3ns_iso_ewing2->GetXaxis()->SetTitle("Timing");
	timing_candidate_3ns_iso_ewing2->GetYaxis()->SetTitle("Entries");

    TH1F *timing_candidate_3ns_iso_ewing3 = new TH1F("timing_candidate_3ns_iso_ewing3", "timing_candidate_3ns_iso_ewing3", 200, -25, 25);
	timing_candidate_3ns_iso_ewing3->GetXaxis()->SetTitle("Timing");
	timing_candidate_3ns_iso_ewing3->GetYaxis()->SetTitle("Entries");

    TH1F *timing_candidate_ss_iso = new TH1F("timing_candidate_ss_iso", "timing_candidate_ss_iso", 200, -25, 25);
	timing_candidate_ss_iso->GetXaxis()->SetTitle("Timing");
	timing_candidate_ss_iso->GetYaxis()->SetTitle("Entries");

    TH1F *timing_candidate_ss_iso_ewing1 = new TH1F("timing_candidate_ss_iso_ewing1", "timing_candidate_ss_iso_ewing1", 200, -25, 25);
	timing_candidate_ss_iso_ewing1->GetXaxis()->SetTitle("Timing");
	timing_candidate_ss_iso_ewing1->GetYaxis()->SetTitle("Entries");

    TH1F *timing_candidate_ss_3ns_iso = new TH1F("timing_candidate_ss_3ns_iso", "timing_candidate_ss_3ns_iso", 200, -25, 25);
	timing_candidate_ss_3ns_iso->GetXaxis()->SetTitle("Timing");
	timing_candidate_ss_3ns_iso->GetYaxis()->SetTitle("Entries");

    TH1F *timing_candidate_ss_3ns_iso_ewing1 = new TH1F("timing_candidate_ss_3ns_iso_ewing1", "timing_candidate_ss_3ns_iso_ewing1", 200, -25, 25);
	timing_candidate_ss_3ns_iso_ewing1->GetXaxis()->SetTitle("Timing");
	timing_candidate_ss_3ns_iso_ewing1->GetYaxis()->SetTitle("Entries");

    TH1F *timing_candidate_ss_3ns_iso_ewing2 = new TH1F("timing_candidate_ss_3ns_iso_ewing2", "timing_candidate_ss_3ns_iso_ewing2", 200, -25, 25);
	timing_candidate_ss_3ns_iso_ewing2->GetXaxis()->SetTitle("Timing");
	timing_candidate_ss_3ns_iso_ewing2->GetYaxis()->SetTitle("Entries");

    TH1F *timing_candidate_ss_3ns_iso_ewing3 = new TH1F("timing_candidate_ss_3ns_iso_ewing3", "timing_candidate_ss_3ns_iso_ewing3", 200, -25, 25);
	timing_candidate_ss_3ns_iso_ewing3->GetXaxis()->SetTitle("Timing");
	timing_candidate_ss_3ns_iso_ewing3->GetYaxis()->SetTitle("Entries");


    TH2F *timing_vs_LICTD_candidate_iso = new TH2F("timing_vs_LICTD_candidate_iso", "timing_vs_LICTD_candidate_iso", 200, -25.0, 25.0, 200, -25, 25);
	timing_vs_LICTD_candidate_iso->GetXaxis()->SetTitle("Seedtime");
	timing_vs_LICTD_candidate_iso->GetYaxis()->SetTitle("LICTD");

    TH2F *timing_vs_LICTD_candidate_iso_ewing1 = new TH2F("timing_vs_LICTD_candidate_iso_ewing1", "timing_vs_LICTD_candidate_iso_ewing1", 200, -25.0, 25.0, 200, -25, 25);
	timing_vs_LICTD_candidate_iso_ewing1->GetXaxis()->SetTitle("Seedtime");
	timing_vs_LICTD_candidate_iso_ewing1->GetYaxis()->SetTitle("LICTD");

    TH2F *timing_vs_LICTD_candidate_ss_iso = new TH2F("timing_vs_LICTD_candidate_ss_iso", "timing_vs_LICTD_candidate_ss_iso", 200, -25.0, 25.0, 200, -25, 25);
	timing_vs_LICTD_candidate_ss_iso->GetXaxis()->SetTitle("Seedtime");
	timing_vs_LICTD_candidate_ss_iso->GetYaxis()->SetTitle("LICTD");

    TH2F *timing_vs_LICTD_candidate_ss_iso_ewing1 = new TH2F("timing_vs_LICTD_candidate_ss_iso_ewing1", "timing_vs_LICTD_candidate_ss_iso_ewing1", 200, -25.0, 25.0, 200, -25, 25);
	timing_vs_LICTD_candidate_ss_iso_ewing1->GetXaxis()->SetTitle("Seedtime");
	timing_vs_LICTD_candidate_ss_iso_ewing1->GetYaxis()->SetTitle("LICTD");


    



    //**************SPIKES**********************
    
    TH1F *ewing_spike = new TH1F("ewing_spike", "ewing_spike", 2000, 0, 1.0);
	ewing_spike->GetXaxis()->SetTitle("ewing");
	ewing_spike->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_spike_iso = new TH1F("ewing_spike_iso", "ewing_spike_iso", 2000, 0, 1.0);
	ewing_spike_iso->GetXaxis()->SetTitle("ewing");
	ewing_spike_iso->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_spike_early_iso = new TH1F("ewing_spike_early_iso", "ewing_spike_early_iso", 2000, 0, 1.0);
	ewing_spike_early_iso->GetXaxis()->SetTitle("ewing");
	ewing_spike_early_iso->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_spike_ss = new TH1F("ewing_spike_ss", "ewing_spike_ss", 2000, 0, 1.0);
	ewing_spike_ss->GetXaxis()->SetTitle("ewing");
	ewing_spike_ss->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_spike_ss_iso = new TH1F("ewing_spike_ss_iso", "ewing_spike_ss_iso", 2000, 0, 1.0);
	ewing_spike_ss_iso->GetXaxis()->SetTitle("ewing");
	ewing_spike_ss_iso->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_spike_ss_early_iso = new TH1F("ewing_spike_ss_early_iso", "ewing_spike_ss_early_iso", 2000, 0, 1.0);
	ewing_spike_ss_early_iso->GetXaxis()->SetTitle("ewing");
	ewing_spike_ss_early_iso->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_spike_HCAL = new TH1F("ewing_spike_HCAL", "ewing_spike_HCAL", 2000, 0, 1.0);
	ewing_spike_HCAL->GetXaxis()->SetTitle("ewing");
	ewing_spike_HCAL->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_spike_ECAL = new TH1F("ewing_spike_ECAL", "ewing_spike_ECAL", 2000, 0, 1.0);
	ewing_spike_ECAL->GetXaxis()->SetTitle("ewing");
	ewing_spike_ECAL->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_spike_track = new TH1F("ewing_spike_track", "ewing_spike_track", 2000, 0, 1.0);
	ewing_spike_track->GetXaxis()->SetTitle("ewing");
	ewing_spike_track->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_spike_ss_HCAL = new TH1F("ewing_spike_ss_HCAL", "ewing_spike_ss_HCAL", 2000, 0, 1.0);
	ewing_spike_ss_HCAL->GetXaxis()->SetTitle("ewing");
	ewing_spike_ss_HCAL->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_spike_ss_ECAL = new TH1F("ewing_spike_ss_ECAL", "ewing_spike_ss_ECAL", 2000, 0, 1.0);
	ewing_spike_ss_ECAL->GetXaxis()->SetTitle("ewing");
	ewing_spike_ss_ECAL->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_spike_ss_track = new TH1F("ewing_spike_ss_track", "ewing_spike_ss_track", 2000, 0, 1.0);
	ewing_spike_ss_track->GetXaxis()->SetTitle("ewing");
	ewing_spike_ss_track->GetYaxis()->SetTitle("Entries");





    TH1F *timing_spike_iso = new TH1F("timing_spike_iso", "timing_spike_iso", 200, -25, 25);
	timing_spike_iso->GetXaxis()->SetTitle("Timing");
	timing_spike_iso->GetYaxis()->SetTitle("Entries");

    TH1F *timing_spike_iso_ewing1 = new TH1F("timing_spike_iso_ewing1", "timing_spike_iso_ewing1", 200, -25, 25);
	timing_spike_iso_ewing1->GetXaxis()->SetTitle("Timing");
	timing_spike_iso_ewing1->GetYaxis()->SetTitle("Entries");

    TH1F *timing_spike_early_iso = new TH1F("timing_spike_early_iso", "timing_spike_early_iso", 200, -25, 25);
	timing_spike_early_iso->GetXaxis()->SetTitle("Timing");
	timing_spike_early_iso->GetYaxis()->SetTitle("Entries");

    TH1F *timing_spike_early_iso_ewing1 = new TH1F("timing_spike_early_iso_ewing1", "timing_spike_early_iso_ewing1", 200, -25, 25);
	timing_spike_early_iso_ewing1->GetXaxis()->SetTitle("Timing");
	timing_spike_early_iso_ewing1->GetYaxis()->SetTitle("Entries");

    TH1F *timing_spike_early_iso_ewing2 = new TH1F("timing_spike_early_iso_ewing2", "timing_spike_early_iso_ewing2", 200, -25, 25);
	timing_spike_early_iso_ewing2->GetXaxis()->SetTitle("Timing");
	timing_spike_early_iso_ewing2->GetYaxis()->SetTitle("Entries");

    TH1F *timing_spike_early_iso_ewing3 = new TH1F("timing_spike_early_iso_ewing3", "timing_spike_early_iso_ewing3", 200, -25, 25);
	timing_spike_early_iso_ewing3->GetXaxis()->SetTitle("Timing");
	timing_spike_early_iso_ewing3->GetYaxis()->SetTitle("Entries");

    TH1F *timing_spike_ss_iso = new TH1F("timing_spike_ss_iso", "timing_spike_ss_iso", 200, -25, 25);
	timing_spike_ss_iso->GetXaxis()->SetTitle("Timing");
	timing_spike_ss_iso->GetYaxis()->SetTitle("Entries");

    TH1F *timing_spike_ss_iso_ewing1 = new TH1F("timing_spike_ss_iso_ewing1", "timing_spike_ss_iso_ewing1", 200, -25, 25);
	timing_spike_ss_iso_ewing1->GetXaxis()->SetTitle("Timing");
	timing_spike_ss_iso_ewing1->GetYaxis()->SetTitle("Entries");

    TH1F *timing_spike_ss_early_iso = new TH1F("timing_spike_ss_early_iso", "timing_spike_ss_early_iso", 200, -25, 25);
	timing_spike_ss_early_iso->GetXaxis()->SetTitle("Timing");
	timing_spike_ss_early_iso->GetYaxis()->SetTitle("Entries");

    TH1F *timing_spike_ss_early_iso_ewing1 = new TH1F("timing_spike_ss_early_iso_ewing1", "timing_spike_ss_early_iso_ewing1", 200, -25, 25);
	timing_spike_ss_early_iso_ewing1->GetXaxis()->SetTitle("Timing");
	timing_spike_ss_early_iso_ewing1->GetYaxis()->SetTitle("Entries");

    TH1F *timing_spike_ss_early_iso_ewing2 = new TH1F("timing_spike_ss_early_iso_ewing2", "timing_spike_ss_early_iso_ewing2", 200, -25, 25);
	timing_spike_ss_early_iso_ewing2->GetXaxis()->SetTitle("Timing");
	timing_spike_ss_early_iso_ewing2->GetYaxis()->SetTitle("Entries");

    TH1F *timing_spike_ss_early_iso_ewing3 = new TH1F("timing_spike_ss_early_iso_ewing3", "timing_spike_ss_early_iso_ewing3", 200, -25, 25);
	timing_spike_ss_early_iso_ewing3->GetXaxis()->SetTitle("Timing");
	timing_spike_ss_early_iso_ewing3->GetYaxis()->SetTitle("Entries");


    TH2F *timing_vs_LICTD_spike_iso = new TH2F("timing_vs_LICTD_spike_iso", "timing_vs_LICTD_spike_iso", 200, -25.0, 25.0, 200, -25, 25);
	timing_vs_LICTD_spike_iso->GetXaxis()->SetTitle("Seedtime");
	timing_vs_LICTD_spike_iso->GetYaxis()->SetTitle("LICTD");

    TH2F *timing_vs_LICTD_spike_iso_ewing1 = new TH2F("timing_vs_LICTD_spike_iso_ewing1", "timing_vs_LICTD_spike_iso_ewing1", 200, -25.0, 25.0, 200, -25, 25);
	timing_vs_LICTD_spike_iso_ewing1->GetXaxis()->SetTitle("Seedtime");
	timing_vs_LICTD_spike_iso_ewing1->GetYaxis()->SetTitle("LICTD");

    TH2F *timing_vs_LICTD_spike_ss_iso = new TH2F("timing_vs_LICTD_spike_ss_iso", "timing_vs_LICTD_spike_ss_iso", 200, -25.0, 25.0, 200, -25, 25);
	timing_vs_LICTD_spike_ss_iso->GetXaxis()->SetTitle("Seedtime");
	timing_vs_LICTD_spike_ss_iso->GetYaxis()->SetTitle("LICTD");

    TH2F *timing_vs_LICTD_spike_ss_iso_ewing1 = new TH2F("timing_vs_LICTD_spike_ss_iso_ewing1", "timing_vs_LICTD_spike_ss_iso_ewing1", 200, -25.0, 25.0, 200, -25, 25);
	timing_vs_LICTD_spike_ss_iso_ewing1->GetXaxis()->SetTitle("Seedtime");
	timing_vs_LICTD_spike_ss_iso_ewing1->GetYaxis()->SetTitle("LICTD");


    //**********************PromptZ**********************
    TH1F *InvMass_Z = new TH1F("InvMass_Z", "InvMass_Z", 2000, 25, 25);
	InvMass_Z->GetXaxis()->SetTitle("Invaraint mass of Z (GeV)");
	InvMass_Z->GetYaxis()->SetTitle("Entries");
    
    TH1F *ewing_promptZ = new TH1F("ewing_promptZ", "ewing_promptZ", 2000, 0, 1.0);
	ewing_promptZ->GetXaxis()->SetTitle("ewing");
	ewing_promptZ->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_promptZ_iso = new TH1F("ewing_promptZ_iso", "ewing_promptZ_iso", 2000, 0, 1.0);
	ewing_promptZ_iso->GetXaxis()->SetTitle("ewing");
	ewing_promptZ_iso->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_promptZ_3ns_iso = new TH1F("ewing_promptZ_3ns_iso", "ewing_promptZ_3ns_iso", 2000, 0, 1.0);
	ewing_promptZ_3ns_iso->GetXaxis()->SetTitle("ewing");
	ewing_promptZ_3ns_iso->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_promptZ_HCAL = new TH1F("ewing_promptZ_HCAL", "ewing_promptZ_HCAL", 2000, 0, 1.0);
	ewing_promptZ_HCAL->GetXaxis()->SetTitle("ewing");
	ewing_promptZ_HCAL->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_promptZ_ECAL = new TH1F("ewing_promptZ_ECAL", "ewing_promptZ_ECAL", 2000, 0, 1.0);
	ewing_promptZ_ECAL->GetXaxis()->SetTitle("ewing");
	ewing_promptZ_ECAL->GetYaxis()->SetTitle("Entries");

    TH1F *ewing_promptZ_track = new TH1F("ewing_promptZ_track", "ewing_promptZ_track", 2000, 0, 1.0);
	ewing_promptZ_track->GetXaxis()->SetTitle("ewing");
	ewing_promptZ_track->GetYaxis()->SetTitle("Entries");




    TH1F *timing_promptZ_iso = new TH1F("timing_promptZ_iso", "timing_promptZ_iso", 200, -25, 25);
	timing_promptZ_iso->GetXaxis()->SetTitle("Timing");
	timing_promptZ_iso->GetYaxis()->SetTitle("Entries");

    TH1F *timing_promptZ_iso_ewing1 = new TH1F("timing_promptZ_iso_ewing1", "timing_promptZ_iso_ewing1", 200, -25, 25);
	timing_promptZ_iso_ewing1->GetXaxis()->SetTitle("Timing");
	timing_promptZ_iso_ewing1->GetYaxis()->SetTitle("Entries");

    TH1F *timing_promptZ_3ns_iso = new TH1F("timing_promptZ_3ns_iso", "timing_promptZ_3ns_iso", 200, -25, 25);
	timing_promptZ_3ns_iso->GetXaxis()->SetTitle("Timing");
	timing_promptZ_3ns_iso->GetYaxis()->SetTitle("Entries");

    TH1F *timing_promptZ_3ns_iso_ewing1 = new TH1F("timing_promptZ_3ns_iso_ewing1", "timing_promptZ_3ns_iso_ewing1", 200, -25, 25);
	timing_promptZ_3ns_iso_ewing1->GetXaxis()->SetTitle("Timing");
	timing_promptZ_3ns_iso_ewing1->GetYaxis()->SetTitle("Entries");

    TH1F *timing_promptZ_3ns_iso_ewing2 = new TH1F("timing_promptZ_3ns_iso_ewing2", "timing_promptZ_3ns_iso_ewing2", 200, -25, 25);
	timing_promptZ_3ns_iso_ewing2->GetXaxis()->SetTitle("Timing");
	timing_promptZ_3ns_iso_ewing2->GetYaxis()->SetTitle("Entries");

    TH1F *timing_promptZ_3ns_iso_ewing3 = new TH1F("timing_promptZ_3ns_iso_ewing3", "timing_promptZ_3ns_iso_ewing3", 200, -25, 25);
	timing_promptZ_3ns_iso_ewing3->GetXaxis()->SetTitle("Timing");
	timing_promptZ_3ns_iso_ewing3->GetYaxis()->SetTitle("Entries");



    //**********************Beamhalo**********************
    TH1F *timing_beamhalo = new TH1F("timing_beamhalo", "timing_beamhalo", 200, -25, 25);
	timing_beamhalo->GetXaxis()->SetTitle("Timing");
	timing_beamhalo->GetYaxis()->SetTitle("Entries");

    TH1F *timing_beamhalo_iso_ewing1 = new TH1F("timing_beamhalo_iso_ewing1", "timing_beamhalo_iso_ewing1", 200, -25, 25);
	timing_beamhalo_iso_ewing1->GetXaxis()->SetTitle("Timing");
	timing_beamhalo_iso_ewing1->GetYaxis()->SetTitle("Entries");


    TH1F *promptz_sieie = new TH1F("promptz_sieie", "promptz_sieie", 200, 0, 0.03);
	promptz_sieie->GetXaxis()->SetTitle("SigmaIEtaIEta");
	promptz_sieie->GetYaxis()->SetTitle("Entries");

    TH1F *beamhalo_sieie = new TH1F("beamhalo_sieie", "beamhalo_sieie", 200, 0, 0.03);
	beamhalo_sieie->GetXaxis()->SetTitle("SigmaIEtaIEta");
	beamhalo_sieie->GetYaxis()->SetTitle("Entries");

    TH2F *promptz_sieie_vs_MIPE = new TH2F("promptz_sieie_vs_MIPE", "promptz_sieie_vs_MIPE", 200, 0, 0.03, 200, 0, 10);
	promptz_sieie_vs_MIPE->GetXaxis()->SetTitle("SigmaIEtaIEta");
	promptz_sieie_vs_MIPE->GetYaxis()->SetTitle("MIPtotalEnergy");

    auto promptz_sieie_vs_MIPE_profile = new TProfile("promptz_sieie_vs_MIPE_profile", "promptz_sieie_vs_MIPE_profile", 2000, 0, 0.03, 0, 10);
	promptz_sieie_vs_MIPE_profile->GetXaxis()->SetTitle("SigmaIEtaIEta");
	promptz_sieie_vs_MIPE_profile->GetYaxis()->SetTitle("MIPTotalEnergy (GeV)");



    
    //ofstream file1; //, file2, file3, file4, file5, file6, file7, file8, file9;
    //file1.open("pickevents.txt", ios::out);
    /*
    file2.open("ootmetfiltersinfo.txt", ios::out);
    file3.open("spike_info_after.txt", ios::out);
    file4.open("spike_info_first.txt", ios::out);
    file5.open("spike_info_sameclus.txt", ios::out);
    file6.open("err_etawing_cal.txt", ios::out);
    file7.open("err_etawing_cal_new.txt", ios::out);
    file8.open("info_crystal_seed_lictdcal.txt", ios::out);
    file9.open("info_crystal_seed_lictdcal_around0.txt", ios::out);
    */

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        if ((jentry % 10000) == 0)
            // to print the number of processed entries
            std::cout << "Processed: " << jentry << std::endl;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        


        
        //============Selections==============

        /*
        ================================================================
        ================================================================
        ================================================================
        ================================================================
        ================================================================
        */

        //declare variables
        Int_t nphotons = 0;
        Int_t nophotons = 0;
        
        Int_t mitIDXSEL = -1; //set index for photon, initial value could be something non-physical, negative something.
		Int_t mootIDXSEL = -1;
        Int_t mootZIDXSEL = -1;
        Int_t mitZIDXSEL = -1;
        
		Bool_t misOOT = kFALSE; //if it's OOT photon, select the OOT, can be true or false but false makes more sense.
		Int_t mIDXSEL = -1;
        Int_t mZIDXSEL = -1;

        vector<Int_t> cellsEB_right;
        vector<Int_t> cellsEB_left;

        Int_t rightcellsIDX_only = -1;
        Int_t leftcellsIDX_only = -1;

        Int_t rightcellsIDX = -1;
        Int_t leftcellsIDX = -1;

        Float_t maxCellsE = 0;

        Float_t it_initimediff = 0;
        Float_t lictdit = 0;

        Float_t oot_initimediff = 0;
        Float_t lictdoot = 0;

        Float_t case3_initimediff = 0;
        Float_t lictdcase3 = 0;

        Float_t case4_initimediff = 0;
        Float_t lictdcase4 = 0;

        Float_t case4_initimediff_ivrsp = 0;
        Float_t lictdcase4_ivrsp = 0;

        Float_t case5_initimediff = 0;
        Float_t lictdcase5 = 0;





        
        /*
        ===========----------=====--------=========----========--=========----------========
        ==============---=========--==============--===--=====--=--=======---=====--========
        ==============---=========--------=======--=====--===--===--======----------========
        ==============---=========--=============--======--=--====--======---===============
        ==============---=========--------=======--======----=====--======---===============
        */
        //=========================preselection============================

        //adding these tempaltes info for doing timing fit with LICTD info

        vector<Int_t> itpho;
		vector<Int_t> itootmatched;
		vector<Int_t> ootpho;
		vector<Int_t> ootitmatched;

        vector<Int_t> iiZpho;
		vector<Int_t> itootZmatched;
		vector<Int_t> ooZpho;
		vector<Int_t> ootitZmatched;

        //in-time collection
		for (int ipho = 0; ipho < nPho; ipho++)
        {
            if ((*phoEt)[ipho] > 225 && fabs((*phoEta)[ipho]) < 1.442 && (*phoR9)[ipho] > 0.8 
            && (*phoHoverE)[ipho] < 0.05 && ((HLTPho >> 10 & 1) == 1))//&& HLTPho < 2048
            {
            
				itpho.push_back(ipho);
				itootmatched.push_back(-1);

                //prompt Z in-time collection
                float PhoSep = 999;
				for (int iipho = ipho + 1; iipho < nPho; iipho++)
				{
					Float_t deltaR = DeltaR((*phoEta)[ipho], (*phoPhi)[ipho], (*phoEta)[iipho], (*phoPhi)[iipho]);
					if (deltaR < PhoSep)
					{
						PhoSep = deltaR;
					}
                
                    //collect the second photon
                    if ((*phoEt)[iipho] > 10.0 && fabs((*phoEta)[iipho]) < 1.442 && (*phoR9)[iipho] > 0.8 && (*phoHoverE)[iipho] < 0.05 && PhoSep > 0.2)
                    {
                        iiZpho.push_back(iipho);
                        itootZmatched.push_back(-1);
                    }
                }
                
                nphotons++;
			}
		}

        //out-of-time collection
		for (int opho = 0; opho < onPho; opho++)
		{
            if ((*ophoEt)[opho] > 225 && fabs((*ophoEta)[opho]) < 1.442 
            && (*ophoR9)[opho] > 0.8 && (*ophoHoverE)[opho] < 0.05 && ((HLTPho >> 10 & 1) == 1))//&& HLTPho < 2048
            {
            
				ootpho.push_back(opho);
				ootitmatched.push_back(-1);


                //prompt Z out-of-time collection
                float PhoSep = 999;
				for (int oopho = opho + 1; oopho < nPho; oopho++)
				{
					Float_t deltaR = DeltaR((*ophoEta)[opho], (*ophoPhi)[opho], (*ophoEta)[oopho], (*ophoPhi)[oopho]);
					if (deltaR < PhoSep)
					{
						PhoSep = deltaR;
					}

                    //collect the second photon
                    if ((*ophoEt)[oopho] > 10.0 && fabs((*ophoEta)[oopho]) < 1.442 && (*ophoR9)[oopho] > 0.8 && (*ophoHoverE)[oopho] < 0.05)
                    {
                        ooZpho.push_back(oopho);
                        ootitZmatched.push_back(-1);
                    }
                }

                nophotons++;

			}
		}

        //for getting the minimum deltaR
		//for (unsigned i = 0; i < itpho.size(); i++)
		for (unsigned i = 0; i < itpho.size(); ++i)
		{
			float PhoSep = 999;
			Int_t itootmatchedidx = -1; //defualt that they are not matched
			//Int_t itootseparateidx = 1;

			//for (unsigned j = 0; j < ootpho.size(); j++)
			for (unsigned j = 0; j < ootpho.size(); ++j)
			{
				Float_t deltaR = DeltaR((*phoEta)[i], (*phoPhi)[i], (*ophoEta)[j], (*ophoPhi)[j]);
				if (deltaR < PhoSep)
				{
					PhoSep = deltaR;
				}
				
				if (PhoSep < 0.2)
				{
					itootmatchedidx = j;
					itootmatched[i] = itootmatchedidx;
					ootitmatched[itootmatchedidx] = i;
				}
			}
		}

        //second photon matched for Z prompt
        for (unsigned i = 0; i < iiZpho.size(); ++i)
		{
			float PhoSep = 999;
			Int_t itootZmatchedidx = -1; //defualt that they are not matched
			//Int_t itootseparateidx = 1;

			//for (unsigned j = 0; j < ootpho.size(); j++)
			for (unsigned j = 0; j < ooZpho.size(); ++j)
			{
				Float_t deltaR = DeltaR((*phoEta)[i], (*phoPhi)[i], (*ophoEta)[j], (*ophoPhi)[j]);
				if (deltaR < PhoSep)
				{
					PhoSep = deltaR;
				}
				
				if (PhoSep < 0.2)
				{
					itootZmatchedidx = j;
					itootZmatched[i] = itootZmatchedidx;
					ootitZmatched[itootZmatchedidx] = i;
				}
			}
		}

        //===========================================================================


        //===========================analysis=====================================


        //---------------------------case 1, only it, no oot---------------------------
		if (itpho.size() != 0 && ootpho.size() == 0 && itootmatched.size() != 0)
		{
			mIDXSEL = itpho[0];
        
            //candidate event
			if (pfMET > 210 && (*phohasPixelSeed)[mIDXSEL] == 0 && metFilters == 0 &&
            (*phoMIPTotEnergy)[mIDXSEL] < 4.9 &&
            //(*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 &&x
            //(*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001 &&
            (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.0105)
			{      
                Float_t uncorrectedPhoEt = ((*phoCalibEt)[mIDXSEL]);
                uncorrectedPhoEt = uncorrectedPhoEt - 0.015 * uncorrectedPhoEt;
                Float_t phoPFECALClusIsoCorr = (*phoPFClusEcalIso)[mIDXSEL] - rho*EcalEA((*phoSCEta)[mIDXSEL]) - EcalEA_ptscale((*phoSCEta)[mIDXSEL], uncorrectedPhoEt); 
                Float_t phoPFHCALClusIsoCorr = (*phoPFClusHcalIso)[mIDXSEL] - rho*HcalEA((*phoSCEta)[mIDXSEL]) - HcalEA_ptscale((*phoSCEta)[mIDXSEL], uncorrectedPhoEt); 
                Float_t phoTkrIsoCorr        = (*phoTrkSumPtHollowConeDR03)[mIDXSEL] - rho * TkrEffAreas((*phoSCEta)[mIDXSEL]);


                //cand1 searching adjencent cells
                for (int icell = 0; icell < nAllCellsEB; icell++)
                {                    
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] - 1)
                    {                    
                        cellsEB_left.push_back(icell);
                    }            

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 
                    }        
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == itpho[0])
                            {
                                //cal aboslute timing diff
                                float it_abstimediff = abs((*phoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (it_abstimediff > it_initimediff)
                                {
                                    it_initimediff = it_abstimediff;
                                    lictdit = ((*phoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    if (nphotons == 1) //on candidate only
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (phoPFECALClusIsoCorr < 4.78 &&
                            phoPFHCALClusIsoCorr < 6.40 &&
                            phoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*phoSeedTime)[mIDXSEL], lictdit);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL], lictdit);
                                }

                                //matrix method
                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(phoPFHCALClusIsoCorr);
                            }

                            if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(phoTkrIsoCorr);
                            }
                        }

                        if ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*phoSeedTime)[mIDXSEL], lictdit);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL], lictdit);
                                    }

                                    //matrix method
                                    if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(phoPFHCALClusIsoCorr);
                                }

                                if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(phoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (phoPFECALClusIsoCorr < 4.78 &&
                        phoPFHCALClusIsoCorr < 6.40 &&
                        phoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*phoSeedTime)[mIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*phoSeedTime)[mIDXSEL], lictdit);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL], lictdit);
                            }

                            //matrix method
                            if ((*phoSeedTime)[mIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*phoSeedTime)[mIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(phoPFHCALClusIsoCorr);
                            }

                            if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(phoTkrIsoCorr);
                            }
                        }

                        if ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*phoSeedTime)[mIDXSEL], lictdit);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL], lictdit);
                                    }

                                    //matrix method
                                    if ((*phoSeedTime)[mIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*phoSeedTime)[mIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(phoPFHCALClusIsoCorr);
                                }

                                if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(phoTkrIsoCorr);
                                }
                            }
                        }
                
                }


                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == itpho[0])
                            {
                                //cal aboslute timing diff
                                float it_abstimediff = abs((*phoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (it_abstimediff > it_initimediff)
                                {
                                    it_initimediff = it_abstimediff;
                                    lictdit = ((*phoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    if (nphotons == 1) //on candidate only
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (phoPFECALClusIsoCorr < 4.78 &&
                            phoPFHCALClusIsoCorr < 6.40 &&
                            phoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*phoSeedTime)[mIDXSEL], lictdit);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL], lictdit);
                                }

                                //matrix method
                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(phoPFHCALClusIsoCorr);
                            }

                            if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(phoTkrIsoCorr);
                            }
                        }

                        if ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*phoSeedTime)[mIDXSEL], lictdit);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL], lictdit);
                                    }

                                    //matrix method
                                    if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*phoSeedTime)[mIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(phoPFHCALClusIsoCorr);
                                }

                                if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(phoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (phoPFECALClusIsoCorr < 4.78 &&
                        phoPFHCALClusIsoCorr < 6.40 &&
                        phoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*phoSeedTime)[mIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*phoSeedTime)[mIDXSEL], lictdit);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL], lictdit);
                            }

                            //matrix method
                            if ((*phoSeedTime)[mIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*phoSeedTime)[mIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(phoPFHCALClusIsoCorr);
                            }

                            if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(phoTkrIsoCorr);
                            }
                        }

                        if ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*phoSeedTime)[mIDXSEL], lictdit);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL], lictdit);
                                    }

                                    //matrix method
                                    if ((*phoSeedTime)[mIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(phoPFHCALClusIsoCorr);
                                }

                                if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(phoTkrIsoCorr);
                                }
                            }
                        }
                }


                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

            

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == itpho[0])
                            {
                                //cal aboslute timing diff
                                float it_abstimediff = abs((*phoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (it_abstimediff > it_initimediff)
                                {
                                    it_initimediff = it_abstimediff;
                                    lictdit = ((*phoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    if (nphotons == 1) //on candidate only
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (phoPFECALClusIsoCorr < 4.78 &&
                            phoPFHCALClusIsoCorr < 6.40 &&
                            phoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*phoSeedTime)[mIDXSEL], lictdit);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL], lictdit);
                                }

                                //matrix method
                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(phoPFHCALClusIsoCorr);
                            }

                            if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(phoTkrIsoCorr);
                            }
                        }

                        if ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*phoSeedTime)[mIDXSEL], lictdit);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL], lictdit);
                                    }

                                    //matrix method
                                    if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*phoSeedTime)[mIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(phoPFHCALClusIsoCorr);
                                }

                                if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(phoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (phoPFECALClusIsoCorr < 4.78 &&
                        phoPFHCALClusIsoCorr < 6.40 &&
                        phoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*phoSeedTime)[mIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*phoSeedTime)[mIDXSEL], lictdit);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL], lictdit);
                            }

                            //matrix method
                            if ((*phoSeedTime)[mIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*phoSeedTime)[mIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(phoPFHCALClusIsoCorr);
                            }

                            if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(phoTkrIsoCorr);
                            }
                        }

                        if ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*phoSeedTime)[mIDXSEL], lictdit);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL], lictdit);
                                    }

                                    //matrix method
                                    if ((*phoSeedTime)[mIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*phoSeedTime)[mIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(phoPFHCALClusIsoCorr);
                                }

                                if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(phoTkrIsoCorr);
                                }
                            }
                        }
                }

            }
        
        

    
            //promptZ selection
            if (iiZpho.size() != 0 && ooZpho.size() == 0 && itootZmatched.size() != 0)
            {
                mZIDXSEL = iiZpho[0];
                //Prompt Z template:
                if ((*phohasPixelSeed)[mIDXSEL] == 1 &&
                    //(*phoMIPTotEnergy)[mIDXSEL] < 4.9 &&
                    //(*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 &&
                    //(*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.0105 &&
                    metFilters == 0 && //106X 256 for passing everything&&
                    (*phohasPixelSeed)[mZIDXSEL] == 1 &&
                    //(*phoMIPTotEnergy)[mZIDXSEL] < 4.9 &&
                    //(*phoSigmaIEtaIEtaFull5x5)[mZIDXSEL] > 0.001 &&
                    //(*phoSigmaIPhiIPhiFull5x5)[mZIDXSEL] > 0.001 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mZIDXSEL] < 0.0105)
                {
                   
                    Float_t uncorrectedPhoEt = ((*phoCalibEt)[mIDXSEL]);
                    uncorrectedPhoEt = uncorrectedPhoEt - 0.015 * uncorrectedPhoEt;
                    Float_t phoPFECALClusIsoCorr = (*phoPFClusEcalIso)[mIDXSEL] - rho*EcalEA((*phoSCEta)[mIDXSEL]) - EcalEA_ptscale((*phoSCEta)[mIDXSEL], uncorrectedPhoEt); 
                    Float_t phoPFHCALClusIsoCorr = (*phoPFClusHcalIso)[mIDXSEL] - rho*HcalEA((*phoSCEta)[mIDXSEL]) - HcalEA_ptscale((*phoSCEta)[mIDXSEL], uncorrectedPhoEt); 
                    Float_t phoTkrIsoCorr        = (*phoTrkSumPtHollowConeDR03)[mIDXSEL] - rho * TkrEffAreas((*phoSCEta)[mIDXSEL]);

                    Float_t InvM = InvariMass((*phoEt)[mIDXSEL], (*phoEt)[mZIDXSEL], (*phoPhi)[mIDXSEL], (*phoPhi)[mZIDXSEL], (*phoEta)[mIDXSEL], (*phoEta)[mZIDXSEL]);
                    
                    InvMass_Z->Fill(InvM);

                    if (InvM > 85 && InvM < 100)
                    {   
                        //promptZ1 searching adjencent cells
                        for (int icell = 0; icell < nAllCellsEB; icell++)
                        {       
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] + 1)
                            {
                                cellsEB_right.push_back(icell); 
                            }

                            //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                            {
                                cellsEB_right.push_back(icell); 
                            }

                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] - 1)
                            {               
                                cellsEB_left.push_back(icell);
                            }       

                            //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                            {
                                cellsEB_left.push_back(icell); 
                            }  
                        }

                        //Now I have to take if there's only one and both
                        //If there's only right
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                        {
                            rightcellsIDX_only = cellsEB_right[0];
                            maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);
                
                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*phoMIPTotEnergy)[mIDXSEL] < 4.9 && (*phoMIPTotEnergy)[mZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                                promptz_sieie->Fill((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL], (*phoMIPTotEnergy)[mIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL], (*phoMIPTotEnergy)[mIDXSEL]);
                                    }
                                }
                        }
                        
                        //if there is only left
                        if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                        {
                            leftcellsIDX_only = cellsEB_left[0];
                            maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*phoMIPTotEnergy)[mIDXSEL] < 4.9 && (*phoMIPTotEnergy)[mZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                                promptz_sieie->Fill((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL], (*phoMIPTotEnergy)[mIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL], (*phoMIPTotEnergy)[mIDXSEL]);
                                    }
                                }
                        }
                        

                        //if there are both
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                        {
                            rightcellsIDX = cellsEB_right[0];
                            leftcellsIDX = cellsEB_left[0];

                            if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[rightcellsIDX];
                            }
                            else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[leftcellsIDX];
                            }

                            Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);


                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*phoMIPTotEnergy)[mIDXSEL] < 4.9 && (*phoMIPTotEnergy)[mZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*phoSeedTime)[mIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                                promptz_sieie->Fill((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL], (*phoMIPTotEnergy)[mIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL], (*phoMIPTotEnergy)[mIDXSEL]);
                                    }
                                }
                        }
                    }
                }
            }


            //BeamHalo template
            if (pfMET > 210 &&
                (*phohasPixelSeed)[mIDXSEL] == 0 &&
                (*phoMIPTotEnergy)[mIDXSEL] > 4.9 &&
                //(*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 &&
                //(*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
                metFilters == 8)//106X beamhalo flag
            {
                
                Float_t uncorrectedPhoEt = ((*phoCalibEt)[mIDXSEL]);
                uncorrectedPhoEt = uncorrectedPhoEt - 0.015 * uncorrectedPhoEt;
                Float_t phoPFECALClusIsoCorr = (*phoPFClusEcalIso)[mIDXSEL] - rho*EcalEA((*phoSCEta)[mIDXSEL]) - EcalEA_ptscale((*phoSCEta)[mIDXSEL], uncorrectedPhoEt); 
                Float_t phoPFHCALClusIsoCorr = (*phoPFClusHcalIso)[mIDXSEL] - rho*HcalEA((*phoSCEta)[mIDXSEL]) - HcalEA_ptscale((*phoSCEta)[mIDXSEL], uncorrectedPhoEt); 
                Float_t phoTkrIsoCorr        = (*phoTrkSumPtHollowConeDR03)[mIDXSEL] - rho * TkrEffAreas((*phoSCEta)[mIDXSEL]);

                timing_beamhalo->Fill((*phoSeedTime)[mIDXSEL]);


                for (int icell = 0; icell < nAllCellsEB; icell++)
                {       
                   if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] + 1)
                   {
                       cellsEB_right.push_back(icell); 
                   }

                   //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                   if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                   {
                       cellsEB_right.push_back(icell); 
                   }

                   if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] - 1)
                   {               
                       cellsEB_left.push_back(icell);
                   }       

                   //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                   if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                   {
                       cellsEB_left.push_back(icell); 
                   }  
               }

               //Now I have to take if there's only one and both
               //If there's only right
               if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
               {
                   rightcellsIDX_only = cellsEB_right[0];
                   maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                   Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);
                  
                   if (phoPFECALClusIsoCorr < 4.78 &&
                   phoPFHCALClusIsoCorr < 6.40 &&
                   phoTkrIsoCorr < 0.89
                   )
                   {
                       if (Ewing > 0.01)
                       {
                           timing_beamhalo_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                           beamhalo_sieie->Fill((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL]); 
                       }
                   }
               }

               //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                    if (phoPFECALClusIsoCorr < 4.78 &&
                    phoPFHCALClusIsoCorr < 6.40 &&
                    phoTkrIsoCorr < 0.89
                    )
                    {
                        if (Ewing > 0.01)
                        {
                            timing_beamhalo_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                            beamhalo_sieie->Fill((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL]); 
                        }
                    }

                }

                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                
                    if (phoPFECALClusIsoCorr < 4.78 &&
                    phoPFHCALClusIsoCorr < 6.40 &&
                    phoTkrIsoCorr < 0.89
                    )
                    {
                        if (Ewing > 0.01)
                        {
                            timing_beamhalo_iso_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                            beamhalo_sieie->Fill((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL]); 
                        }
                    }
                }

            }
        }
    
            

        //---------------------------case 2, only oot, no it---------------------------
		if (itpho.size() == 0 && ootpho.size() != 0 && ootitmatched.size() != 0)
		{
			mIDXSEL = ootpho[0];
			misOOT = kTRUE;
			Float_t NewMET = newuMET(pfMET, pfMETPhi, (*ophoPhi)[mIDXSEL], (*ophoEt)[mIDXSEL]);

			//Candidate Events
			if (NewMET > 210 && (*ophohasPixelSeed)[mIDXSEL] == 0 && metFilters == 0 &&
            (*ophoMIPTotEnergy)[mIDXSEL] < 4.9 &&
            //(*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 &&
            //(*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001 &&
            (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.0105)//106X 256 for passing everything
			{
                
                Float_t uncorrectedPhoEt = ((*ophoEt)[mIDXSEL]);
                uncorrectedPhoEt = uncorrectedPhoEt - 0.015 * uncorrectedPhoEt;
                Float_t ophoPFECALClusIsoCorr = (*ophoPFClusEcalIso)[mIDXSEL] - rho*EcalEA((*ophoSCEta)[mIDXSEL]) - EcalEA_ptscale((*ophoSCEta)[mIDXSEL], uncorrectedPhoEt); 
                Float_t ophoPFHCALClusIsoCorr = (*ophoPFClusHcalIso)[mIDXSEL] - rho*HcalEA((*ophoSCEta)[mIDXSEL]) - HcalEA_ptscale((*ophoSCEta)[mIDXSEL], uncorrectedPhoEt); 
                Float_t ophoTkrIsoCorr        = (*ophoTrkSumPtHollowConeDR03)[mIDXSEL] - rho * TkrEffAreas((*ophoSCEta)[mIDXSEL]);

                        
                //cand1 searching adjencent cells
                for (int icell = 0; icell < nAllCellsEB; icell++)
                {                    
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] - 1)
                    {                    
                        cellsEB_left.push_back(icell);
                    }            

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 
                    }        
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == ootpho[0])
                            {
                                //cal aboslute timing diff
                                float oot_abstimediff = abs((*ophoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (oot_abstimediff > oot_initimediff)
                                {
                                    oot_initimediff = oot_abstimediff;
                                    lictdoot = ((*ophoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                 
                        }                
                    }

                    if (nophotons == 1) //on candidate only
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (ophoPFECALClusIsoCorr < 4.78 &&
                            ophoPFHCALClusIsoCorr < 6.40 &&
                            ophoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);
                                }

                                //matrix method
                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);
                                    }

                                    //matrix method
                                    if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (ophoPFECALClusIsoCorr < 4.78 &&
                        ophoPFHCALClusIsoCorr < 6.40 &&
                        ophoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);
                            }

                            //matrix method
                            if ((*ophoSeedTime)[mIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*ophoSeedTime)[mIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);
                                    }

                                    //matrix method
                                    if ((*ophoSeedTime)[mIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*ophoSeedTime)[mIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                
                }

                            
                            
                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

            

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == ootpho[0])
                            {
                                //cal aboslute timing diff
                                float oot_abstimediff = abs((*ophoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (oot_abstimediff > oot_initimediff)
                                {
                                    oot_initimediff = oot_abstimediff;
                                    lictdoot = ((*ophoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    if (nophotons == 1) //on candidate only
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (ophoPFECALClusIsoCorr < 4.78 &&
                            ophoPFHCALClusIsoCorr < 6.40 &&
                            ophoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);
                                }

                                //matrix method
                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);
                                    }

                                    //matrix method
                                    if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (ophoPFECALClusIsoCorr < 4.78 &&
                        ophoPFHCALClusIsoCorr < 6.40 &&
                        ophoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);
                            }

                            //matrix method
                            if ((*ophoSeedTime)[mIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*ophoSeedTime)[mIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);
                                    }

                                    //matrix method
                                    if ((*ophoSeedTime)[mIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*ophoSeedTime)[mIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                
                }

                            
                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

            

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == ootpho[0])
                            {
                                //cal aboslute timing diff
                                float oot_abstimediff = abs((*ophoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (oot_abstimediff > oot_initimediff)
                                {
                                    oot_initimediff = oot_abstimediff;
                                    lictdoot = ((*ophoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    if (nophotons == 1) //on candidate only
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (ophoPFECALClusIsoCorr < 4.78 &&
                            ophoPFHCALClusIsoCorr < 6.40 &&
                            ophoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);
                                }

                                //matrix method
                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);
                                    }

                                    //matrix method
                                    if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (ophoPFECALClusIsoCorr < 4.78 &&
                        ophoPFHCALClusIsoCorr < 6.40 &&
                        ophoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);
                            }

                            //matrix method
                            if ((*ophoSeedTime)[mIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*ophoSeedTime)[mIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL], lictdoot);
                                    }

                                    //matrix method
                                    if ((*ophoSeedTime)[mIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*ophoSeedTime)[mIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                
                }

            }
        
        
        
            //promptZ
            if (iiZpho.size() == 0 && ooZpho.size() != 0 && ootitZmatched.size() != 0)
            {
                mZIDXSEL = ooZpho[0];
                //Prompt Z template:
                if ((*ophohasPixelSeed)[mIDXSEL] == 1 &&
                    //(*ophoMIPTotEnergy)[mIDXSEL] < 4.9 &&
                    //(*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 &&
                    //(*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.0105 &&
                    metFilters == 0 && //106X 256 for passing everything&&
                    (*ophohasPixelSeed)[mZIDXSEL] == 1 &&
                    //(*ophoMIPTotEnergy)[mZIDXSEL] < 4.9 &&
                    //(*ophoSigmaIEtaIEtaFull5x5)[mZIDXSEL] > 0.001 &&
                    //(*ophoSigmaIPhiIPhiFull5x5)[mZIDXSEL] > 0.001 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mZIDXSEL] < 0.0105)
                {
                    Float_t uncorrectedPhoEt = ((*ophoEt)[mIDXSEL]);
                    uncorrectedPhoEt = uncorrectedPhoEt - 0.015 * uncorrectedPhoEt;
                    Float_t ophoPFECALClusIsoCorr = (*ophoPFClusEcalIso)[mIDXSEL] - rho*EcalEA((*ophoSCEta)[mIDXSEL]) - EcalEA_ptscale((*ophoSCEta)[mIDXSEL], uncorrectedPhoEt); 
                    Float_t ophoPFHCALClusIsoCorr = (*ophoPFClusHcalIso)[mIDXSEL] - rho*HcalEA((*ophoSCEta)[mIDXSEL]) - HcalEA_ptscale((*ophoSCEta)[mIDXSEL], uncorrectedPhoEt); 
                    Float_t ophoTkrIsoCorr        = (*ophoTrkSumPtHollowConeDR03)[mIDXSEL] - rho * TkrEffAreas((*ophoSCEta)[mIDXSEL]);

                    Float_t InvM = InvariMass((*ophoEt)[mIDXSEL], (*ophoEt)[mZIDXSEL], (*ophoPhi)[mIDXSEL], (*ophoPhi)[mZIDXSEL], (*ophoEta)[mIDXSEL], (*ophoEta)[mZIDXSEL]);
                    
                    InvMass_Z->Fill(InvM);

                    //promptZ for timing fit
                    if (InvM > 85 && InvM < 100)
                    {
                        //promptZ1 searching adjencent cells
                        for (int icell = 0; icell < nAllCellsEB; icell++)
                        {       
                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] + 1)
                            {
                                cellsEB_right.push_back(icell); 
                            }

                            //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                            {
                                cellsEB_right.push_back(icell); 
                            }

                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] - 1)
                            {               
                                cellsEB_left.push_back(icell);
                            }       

                            //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                            {
                                cellsEB_left.push_back(icell); 
                            }  
                        }

                        //Now I have to take if there's only one and both
                        //If there's only right
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                        {
                            rightcellsIDX_only = cellsEB_right[0];
                            maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*ophoMIPTotEnergy)[mIDXSEL] < 4.9 && (*ophoMIPTotEnergy)[mZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                                promptz_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL], (*ophoMIPTotEnergy)[mIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL], (*ophoMIPTotEnergy)[mIDXSEL]);
                                    }
                                }
                        }

                        //if there is only left
                        if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                        {
                            leftcellsIDX_only = cellsEB_left[0];
                            maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                        
                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*ophoMIPTotEnergy)[mIDXSEL] < 4.9 && (*ophoMIPTotEnergy)[mZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                                promptz_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL], (*ophoMIPTotEnergy)[mIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL], (*ophoMIPTotEnergy)[mIDXSEL]);
                                    }
                                }
                        }

                        //if there are both
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                        {
                            rightcellsIDX = cellsEB_right[0];
                            leftcellsIDX = cellsEB_left[0];

                            if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[rightcellsIDX];
                            }
                            else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[leftcellsIDX];
                            }

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*ophoMIPTotEnergy)[mIDXSEL] < 4.9 && (*ophoMIPTotEnergy)[mZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*ophoSeedTime)[mIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                                promptz_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL], (*ophoMIPTotEnergy)[mIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL], (*ophoMIPTotEnergy)[mIDXSEL]);
                                    }
                                }
                        }
                    }
                }
            }

            //BeamHalo template
            if (NewMET > 210 &&
                (*ophohasPixelSeed)[mIDXSEL] == 0 &&
                (*ophoMIPTotEnergy)[mIDXSEL] > 4.9 &&
                //(*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 &&
                //(*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001
                metFilters == 8)//106X beamhalo flag
            {
                Float_t uncorrectedoPhoEt = ((*ophoEt)[mIDXSEL]);
                uncorrectedoPhoEt = uncorrectedoPhoEt - 0.015 * uncorrectedoPhoEt;
                Float_t ophoPFECALClusIsoCorr = (*ophoPFClusEcalIso)[mIDXSEL] - rho*EcalEA((*ophoSCEta)[mIDXSEL]) - EcalEA_ptscale((*ophoSCEta)[mIDXSEL], uncorrectedoPhoEt); 
                Float_t ophoPFHCALClusIsoCorr = (*ophoPFClusHcalIso)[mIDXSEL] - rho*HcalEA((*ophoSCEta)[mIDXSEL]) - HcalEA_ptscale((*ophoSCEta)[mIDXSEL], uncorrectedoPhoEt); 
                Float_t ophoTkrIsoCorr        = (*ophoTrkSumPtHollowConeDR03)[mIDXSEL] - rho * TkrEffAreas((*ophoSCEta)[mIDXSEL]);

                timing_beamhalo->Fill((*ophoSeedTime)[mIDXSEL]);

                
                for (int icell = 0; icell < nAllCellsEB; icell++)
                {       
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] - 1)
                    {               
                        cellsEB_left.push_back(icell);
                    }       

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 
                    }  
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                    if (ophoPFECALClusIsoCorr < 4.78 &&
                    ophoPFHCALClusIsoCorr < 6.40 &&
                    ophoTkrIsoCorr < 0.89)
                    {
                        if (Ewing > 0.01)
                        {
                            beamhalo_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]); 
                            timing_beamhalo_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }
                }

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                    if (ophoPFECALClusIsoCorr < 4.78 &&
                    ophoPFHCALClusIsoCorr < 6.40 &&
                    ophoTkrIsoCorr < 0.89)
                    {
                        if (Ewing > 0.01)
                        {
                            beamhalo_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]); 
                            timing_beamhalo_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                
                }

                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                    if (ophoPFECALClusIsoCorr < 4.78 &&
                    ophoPFHCALClusIsoCorr < 6.40 &&
                    ophoTkrIsoCorr < 0.89)
                    {
                        if (Ewing > 0.01)
                        {
                            beamhalo_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]); 
                            timing_beamhalo_iso_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }
                }
            }
        }

        //============cases 3 when there are both=============
        if (itpho.size() != 0 && ootpho.size() != 0)
		{
			//case 3, only it pass trigger
			mootIDXSEL = ootitmatched[0];
			mitIDXSEL = itootmatched[0];

			//Candidate Events
			if (pfMET > 210 && (*phohasPixelSeed)[mitIDXSEL] == 0 && metFilters == 0 &&
            (*phoMIPTotEnergy)[mitIDXSEL] < 4.9 &&
            //(*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 &&
            //(*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001 &&
            (*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.0105)
			{
                
                Float_t uncorrectedPhoEt = ((*phoCalibEt)[mitIDXSEL]);
                uncorrectedPhoEt = uncorrectedPhoEt - 0.015 * uncorrectedPhoEt;
                Float_t phoPFECALClusIsoCorr = (*phoPFClusEcalIso)[mitIDXSEL] - rho*EcalEA((*phoSCEta)[mitIDXSEL]) - EcalEA_ptscale((*phoSCEta)[mitIDXSEL], uncorrectedPhoEt); 
                Float_t phoPFHCALClusIsoCorr = (*phoPFClusHcalIso)[mitIDXSEL] - rho*HcalEA((*phoSCEta)[mitIDXSEL]) - HcalEA_ptscale((*phoSCEta)[mitIDXSEL], uncorrectedPhoEt); 
                Float_t phoTkrIsoCorr        = (*phoTrkSumPtHollowConeDR03)[mitIDXSEL] - rho * TkrEffAreas((*phoSCEta)[mitIDXSEL]);

                //cand1 searching adjencent cells
                for (int icell = 0; icell < nAllCellsEB; icell++)
                {                    
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mitIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && (*phoSeedIEta)[mitIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mitIDXSEL] - 1)
                    {                    
                        cellsEB_left.push_back(icell);
                    }            

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && (*phoSeedIEta)[mitIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 
                    }        
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == itootmatched[0])
                            {
                                //cal aboslute timing diff
                                float case3_abstimediff = abs((*phoSeedTime)[mitIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (case3_abstimediff > case3_initimediff)
                                {
                                    case3_initimediff = case3_abstimediff;
                                    lictdcase3 = ((*phoSeedTime)[mitIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    if (nphotons == 1 && nophotons == 1)
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (phoPFECALClusIsoCorr < 4.78 &&
                            phoPFHCALClusIsoCorr < 6.40 &&
                            phoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);
                                }

                                //matrix method
                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(phoPFHCALClusIsoCorr);
                            }

                            if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(phoTkrIsoCorr);
                            }
                        }

                        if ((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);
                                    }

                                    //matrix method
                                    if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(phoPFHCALClusIsoCorr);
                                }

                                if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(phoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (phoPFECALClusIsoCorr < 4.78 &&
                        phoPFHCALClusIsoCorr < 6.40 &&
                        phoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);
                            }

                            //matrix method
                            if ((*phoSeedTime)[mitIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*phoSeedTime)[mitIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(phoPFHCALClusIsoCorr);
                            }

                            if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(phoTkrIsoCorr);
                            }
                        }

                        if ((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);
                                    }

                                    //matrix method
                                    if ((*phoSeedTime)[mitIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*phoSeedTime)[mitIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(phoPFHCALClusIsoCorr);
                                }

                                if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(phoTkrIsoCorr);
                                }
                            }
                        }
                
                }


                            

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

            

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == itootmatched[0])
                            {
                                //cal aboslute timing diff
                                float case3_abstimediff = abs((*phoSeedTime)[mitIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (case3_abstimediff > case3_initimediff)
                                {
                                    case3_initimediff = case3_abstimediff;
                                    lictdcase3 = ((*phoSeedTime)[mitIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    if (nphotons == 1 && nophotons == 1)
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (phoPFECALClusIsoCorr < 4.78 &&
                            phoPFHCALClusIsoCorr < 6.40 &&
                            phoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);
                                }

                                //matrix method
                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(phoPFHCALClusIsoCorr);
                            }

                            if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(phoTkrIsoCorr);
                            }
                        }

                        if ((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);
                                    }

                                    //matrix method
                                    if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(phoPFHCALClusIsoCorr);
                                }

                                if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(phoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (phoPFECALClusIsoCorr < 4.78 &&
                        phoPFHCALClusIsoCorr < 6.40 &&
                        phoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);
                            }

                            //matrix method
                            if ((*phoSeedTime)[mitIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*phoSeedTime)[mitIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(phoPFHCALClusIsoCorr);
                            }

                            if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(phoTkrIsoCorr);
                            }
                        }

                        if ((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);
                                    }

                                    //matrix method
                                    if ((*phoSeedTime)[mitIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*phoSeedTime)[mitIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(phoPFHCALClusIsoCorr);
                                }

                                if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(phoTkrIsoCorr);
                                }
                            }
                        }
                
                }


                            

                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

            

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == itootmatched[0])
                            {
                                //cal aboslute timing diff
                                float case3_abstimediff = abs((*phoSeedTime)[mitIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (case3_abstimediff > case3_initimediff)
                                {
                                    case3_initimediff = case3_abstimediff;
                                    lictdcase3 = ((*phoSeedTime)[mitIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    if (nphotons == 1 && nophotons == 1)
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (phoPFECALClusIsoCorr < 4.78 &&
                            phoPFHCALClusIsoCorr < 6.40 &&
                            phoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);
                                }

                                //matrix method
                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(phoPFHCALClusIsoCorr);
                            }

                            if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(phoTkrIsoCorr);
                            }
                        }

                        if ((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);
                                    }

                                    //matrix method
                                    if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(phoPFHCALClusIsoCorr);
                                }

                                if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(phoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (phoPFECALClusIsoCorr < 4.78 &&
                        phoPFHCALClusIsoCorr < 6.40 &&
                        phoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);
                            }

                            //matrix method
                            if ((*phoSeedTime)[mitIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*phoSeedTime)[mitIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(phoPFHCALClusIsoCorr);
                            }

                            if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(phoTkrIsoCorr);
                            }
                        }

                        if ((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL], lictdcase3);
                                    }

                                    //matrix method
                                    if ((*phoSeedTime)[mitIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*phoSeedTime)[mitIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (phoPFECALClusIsoCorr < 4.78 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(phoPFHCALClusIsoCorr);
                                }

                                if (phoPFHCALClusIsoCorr < 6.40 && phoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 4.78 && phoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(phoTkrIsoCorr);
                                }
                            }
                        }
                
                }

            }


            //-----------------promptZ cases 3, 4 ,5-----------------
            if (iiZpho.size() != 0 && ooZpho.size() != 0)
            {
                mootZIDXSEL = ootitZmatched[0];
		    	mitZIDXSEL = itootZmatched[0];

                //Prompt Z template case 3
                if ((*phohasPixelSeed)[mitIDXSEL] == 1 &&
                    //(*phoMIPTotEnergy)[mitIDXSEL] < 4.9 &&
                    //(*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 &&
                    //(*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.0105 &&
                    metFilters == 0 && //106X 256 for passing everything&&
                    (*phohasPixelSeed)[mitZIDXSEL] == 1 &&
                    //(*phoMIPTotEnergy)[mitZIDXSEL] < 4.9 &&
                    //(*phoSigmaIEtaIEtaFull5x5)[mitZIDXSEL] > 0.001 &&
                    //(*phoSigmaIPhiIPhiFull5x5)[mitZIDXSEL] > 0.001 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mitZIDXSEL] < 0.0105)
                {
                    
                    Float_t uncorrectedPhoEt = ((*phoCalibEt)[mitIDXSEL]);
                    uncorrectedPhoEt = uncorrectedPhoEt - 0.015 * uncorrectedPhoEt;
                    Float_t phoPFECALClusIsoCorr = (*phoPFClusEcalIso)[mitIDXSEL] - rho*EcalEA((*phoSCEta)[mitIDXSEL]) - EcalEA_ptscale((*phoSCEta)[mitIDXSEL], uncorrectedPhoEt); 
                    Float_t phoPFHCALClusIsoCorr = (*phoPFClusHcalIso)[mitIDXSEL] - rho*HcalEA((*phoSCEta)[mitIDXSEL]) - HcalEA_ptscale((*phoSCEta)[mitIDXSEL], uncorrectedPhoEt); 
                    Float_t phoTkrIsoCorr        = (*phoTrkSumPtHollowConeDR03)[mitIDXSEL] - rho * TkrEffAreas((*phoSCEta)[mitIDXSEL]);

                    Float_t InvM = InvariMass((*phoEt)[mitIDXSEL], (*phoEt)[mitZIDXSEL], (*phoPhi)[mitIDXSEL], (*phoPhi)[mitZIDXSEL], (*phoEta)[mitIDXSEL], (*phoEta)[mitZIDXSEL]);
                    
                    InvMass_Z->Fill(InvM);

                    if (InvM > 85 && InvM < 100)
                    {
                    
                        //promptZ1 searching adjencent cells
                        for (int icell = 0; icell < nAllCellsEB; icell++)
                        {       
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mitIDXSEL] + 1)
                            {
                                cellsEB_right.push_back(icell); 
                            }

                            //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && (*phoSeedIEta)[mitIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                            {
                                cellsEB_right.push_back(icell); 
                            }

                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mitIDXSEL] - 1)
                            {               
                                cellsEB_left.push_back(icell);
                            }       

                            //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && (*phoSeedIEta)[mitIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                            {
                                cellsEB_left.push_back(icell); 
                            }  
                        }

                        //Now I have to take if there's only one and both
                        //If there's only right
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                        {
                            rightcellsIDX_only = cellsEB_right[0];
                            maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*phoMIPTotEnergy)[mitIDXSEL] < 4.9 && (*phoMIPTotEnergy)[mitZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                                promptz_sieie->Fill((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL], (*phoMIPTotEnergy)[mitIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL], (*phoMIPTotEnergy)[mitIDXSEL]);
                                    }
                                }
                        }

                        //if there is only left
                        if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                        {
                            leftcellsIDX_only = cellsEB_left[0];
                            maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                           
                        
                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*phoMIPTotEnergy)[mitIDXSEL] < 4.9 && (*phoMIPTotEnergy)[mitZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                                promptz_sieie->Fill((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL], (*phoMIPTotEnergy)[mitIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL], (*phoMIPTotEnergy)[mitIDXSEL]);
                                    }
                                }
                        }

                        //if there are both
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                        {
                            rightcellsIDX = cellsEB_right[0];
                            leftcellsIDX = cellsEB_left[0];

                            if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[rightcellsIDX];
                            }
                            else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[leftcellsIDX];
                            }

                            Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);


                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (phoPFECALClusIsoCorr < 4.78 &&
                                phoPFHCALClusIsoCorr < 6.40 &&
                                phoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*phoMIPTotEnergy)[mitIDXSEL] < 4.9 && (*phoMIPTotEnergy)[mitZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*phoSeedTime)[mitIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                                promptz_sieie->Fill((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL], (*phoMIPTotEnergy)[mitIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL], (*phoMIPTotEnergy)[mitIDXSEL]);
                                    }
                                }
                        }
                    }

                }

                Float_t NewmpMET = newmMET(pfMET, pfMETPhi, (*ophoPhi)[mootIDXSEL], (*ophoEt)[mootIDXSEL]);

                //Prompt Z template case 4
                if ((*ophohasPixelSeed)[mootIDXSEL] == 1 &&
                    //(*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 &&
                    //(*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 &&
                    //(*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] < 0.0105 &&
                    metFilters == 0 && //106X 256 for passing everything&&
                    (*ophohasPixelSeed)[mootZIDXSEL] == 1 &&
                    //(*ophoMIPTotEnergy)[mootZIDXSEL] < 4.9 &&
                    //(*ophoSigmaIEtaIEtaFull5x5)[mootZIDXSEL] > 0.001 &&
                    //(*ophoSigmaIPhiIPhiFull5x5)[mootZIDXSEL] > 0.001 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mootZIDXSEL] < 0.0105)
                {
                    
                    Float_t uncorrectedoPhoEt = ((*ophoEt)[mootIDXSEL]);
                    uncorrectedoPhoEt = uncorrectedoPhoEt - 0.015 * uncorrectedoPhoEt;
                    Float_t ophoPFECALClusIsoCorr = (*ophoPFClusEcalIso)[mootIDXSEL] - rho*EcalEA((*ophoSCEta)[mootIDXSEL]) - EcalEA_ptscale((*ophoSCEta)[mootIDXSEL], uncorrectedoPhoEt); 
                    Float_t ophoPFHCALClusIsoCorr = (*ophoPFClusHcalIso)[mootIDXSEL] - rho*HcalEA((*ophoSCEta)[mootIDXSEL]) - HcalEA_ptscale((*ophoSCEta)[mootIDXSEL], uncorrectedoPhoEt); 
                    Float_t ophoTkrIsoCorr        = (*ophoTrkSumPtHollowConeDR03)[mootIDXSEL] - rho * TkrEffAreas((*ophoSCEta)[mootIDXSEL]);


                    Float_t InvM = InvariMass((*ophoEt)[mootIDXSEL], (*ophoEt)[mootZIDXSEL], (*ophoPhi)[mootIDXSEL], (*ophoPhi)[mootZIDXSEL], (*ophoEta)[mootIDXSEL], (*ophoEta)[mootZIDXSEL]);
                    
                    InvMass_Z->Fill(InvM);


                    if (InvM > 85 && InvM < 100)
                    {
                        //promptZ1 searching adjencent cells
                        for (int icell = 0; icell < nAllCellsEB; icell++)
                        {       
                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] + 1)
                            {
                                cellsEB_right.push_back(icell); 
                            }

                            //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                            {
                                cellsEB_right.push_back(icell); 
                            }

                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] - 1)
                            {               
                                cellsEB_left.push_back(icell);
                            }       

                            //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                            {
                                cellsEB_left.push_back(icell); 
                            }  
                        }

                        //Now I have to take if there's only one and both
                        //If there's only right
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                        {
                            rightcellsIDX_only = cellsEB_right[0];
                            maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);


                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 && (*ophoMIPTotEnergy)[mootZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                                promptz_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL], (*ophoMIPTotEnergy)[mootIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL], (*ophoMIPTotEnergy)[mootIDXSEL]);
                                    }
                                }
                        }
                
                

                        //if there is only left
                        if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                        {
                            leftcellsIDX_only = cellsEB_left[0];
                            maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 && (*ophoMIPTotEnergy)[mootZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                                promptz_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL], (*ophoMIPTotEnergy)[mootIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL], (*ophoMIPTotEnergy)[mootIDXSEL]);
                                    }
                                }
                        }
                
                
                        //if there are both
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                        {
                            rightcellsIDX = cellsEB_right[0];
                            leftcellsIDX = cellsEB_left[0];

                            if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[rightcellsIDX];
                            }
                            else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[leftcellsIDX];
                            }

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);


                            
                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 && (*ophoMIPTotEnergy)[mootZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                                promptz_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL], (*ophoMIPTotEnergy)[mootIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL], (*ophoMIPTotEnergy)[mootIDXSEL]);
                                    }
                                }
                        }
                
                    }
                }

                //Prompt Z template case 5
                if ((*phohasPixelSeed)[mitIDXSEL] == 1 &&
                    //(*phoMIPTotEnergy)[mitIDXSEL] < 4.9 &&
                    //(*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 &&
                    //(*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.0105 &&
                    metFilters == 0 && //106X 256 for passing everything&&
                    (*phohasPixelSeed)[mitZIDXSEL] == 1 &&
                    //(*phoMIPTotEnergy)[mitZIDXSEL] < 4.9 &&
                    //(*phoSigmaIEtaIEtaFull5x5)[mitZIDXSEL] > 0.001 &&
                    //(*phoSigmaIPhiIPhiFull5x5)[mitZIDXSEL] > 0.001 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mitZIDXSEL] < 0.0105 &&

                    (*ophohasPixelSeed)[mootIDXSEL] == 1 &&
                    //(*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 &&
                    //(*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 &&
                    //(*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] < 0.0105 &&

                    (*ophohasPixelSeed)[mootZIDXSEL] == 1 &&
                    //(*ophoMIPTotEnergy)[mootZIDXSEL] < 4.9 &&
                    //(*ophoSigmaIEtaIEtaFull5x5)[mootZIDXSEL] > 0.001 &&
                    //(*ophoSigmaIPhiIPhiFull5x5)[mootZIDXSEL] > 0.001 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mootZIDXSEL] < 0.0105)
                {
                    Float_t uncorrectedoPhoEt = ((*ophoEt)[mootIDXSEL]);
                    uncorrectedoPhoEt = uncorrectedoPhoEt - 0.015 * uncorrectedoPhoEt;
                    Float_t ophoPFECALClusIsoCorr = (*ophoPFClusEcalIso)[mootIDXSEL] - rho*EcalEA((*ophoSCEta)[mootIDXSEL]) - EcalEA_ptscale((*ophoSCEta)[mootIDXSEL], uncorrectedoPhoEt); 
                    Float_t ophoPFHCALClusIsoCorr = (*ophoPFClusHcalIso)[mootIDXSEL] - rho*HcalEA((*ophoSCEta)[mootIDXSEL]) - HcalEA_ptscale((*ophoSCEta)[mootIDXSEL], uncorrectedoPhoEt); 
                    Float_t ophoTkrIsoCorr        = (*ophoTrkSumPtHollowConeDR03)[mootIDXSEL] - rho * TkrEffAreas((*ophoSCEta)[mootIDXSEL]);


                    Float_t InvM = InvariMass((*ophoEt)[mootIDXSEL], (*ophoEt)[mootZIDXSEL], (*ophoPhi)[mootIDXSEL], (*ophoPhi)[mootZIDXSEL], (*ophoEta)[mootIDXSEL], (*ophoEta)[mootZIDXSEL]);
                    
                    InvMass_Z->Fill(InvM);


                    if (InvM > 85 && InvM < 100)
                    {
                        //promptZ1 searching adjencent cells
                        for (int icell = 0; icell < nAllCellsEB; icell++)
                        {       
                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] + 1)
                            {
                                cellsEB_right.push_back(icell); 
                            }

                            //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                            {
                                cellsEB_right.push_back(icell); 
                            }

                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] - 1)
                            {               
                                cellsEB_left.push_back(icell);
                            }       

                            //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                            {
                                cellsEB_left.push_back(icell); 
                            }  
                        }

                        //Now I have to take if there's only one and both
                        //If there's only right
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                        {
                            rightcellsIDX_only = cellsEB_right[0];
                            maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*phoMIPTotEnergy)[mitIDXSEL] < 4.9 && (*phoMIPTotEnergy)[mitZIDXSEL] && 
                                    (*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 && (*ophoMIPTotEnergy)[mootZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                                promptz_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL], (*ophoMIPTotEnergy)[mootIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL], (*ophoMIPTotEnergy)[mootIDXSEL]);
                                    }
                                }
                        }

                        //if there is only left
                        if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                        {
                            leftcellsIDX_only = cellsEB_left[0];
                            maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*phoMIPTotEnergy)[mitIDXSEL] < 4.9 && (*phoMIPTotEnergy)[mitZIDXSEL] && 
                                    (*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 && (*ophoMIPTotEnergy)[mootZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                                promptz_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL], (*ophoMIPTotEnergy)[mootIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL], (*ophoMIPTotEnergy)[mootIDXSEL]);
                                    }
                                }
                        }

                        //if there are both
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                        {
                            rightcellsIDX = cellsEB_right[0];
                            leftcellsIDX = cellsEB_left[0];

                            if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[rightcellsIDX];
                            }
                            else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[leftcellsIDX];
                            }

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);


                            ewing_promptZ->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    if ((*phoMIPTotEnergy)[mitIDXSEL] < 4.9 && (*phoMIPTotEnergy)[mitZIDXSEL] && 
                                    (*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 && (*ophoMIPTotEnergy)[mootZIDXSEL] < 4.9)
                                    {
                                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                        {
                                            ewing_promptZ_3ns_iso->Fill(Ewing);

                                            timing_promptZ_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                            if (Ewing > 0.01)
                                            {
                                                timing_promptZ_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                                promptz_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL]);
                                            }
                                        }
                                    }

                                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0 && Ewing > 0.01)
                                    {
                                        promptz_sieie_vs_MIPE->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL], (*ophoMIPTotEnergy)[mootIDXSEL]);
                                        promptz_sieie_vs_MIPE_profile->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL], (*ophoMIPTotEnergy)[mootIDXSEL]);
                                    }
                                }
                        }
                    }
                }
            }


            //BeamHalo template
			if (pfMET > 210 &&
				(*phohasPixelSeed)[mitIDXSEL] == 0 &&
				(*phoMIPTotEnergy)[mitIDXSEL] > 4.9 &&
				//(*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 &&
				//(*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001
				metFilters == 8)//106X beamhalo flag
			{
                
                
                Float_t uncorrectedPhoEt = ((*phoCalibEt)[mitIDXSEL]);  
                uncorrectedPhoEt = uncorrectedPhoEt - 0.015 * uncorrectedPhoEt;
                Float_t phoPFECALClusIsoCorr = (*phoPFClusEcalIso)[mitIDXSEL] - rho*EcalEA((*phoSCEta)[mitIDXSEL]) - EcalEA_ptscale((*phoSCEta)[mitIDXSEL], uncorrectedPhoEt); 
                Float_t phoPFHCALClusIsoCorr = (*phoPFClusHcalIso)[mitIDXSEL] - rho*HcalEA((*phoSCEta)[mitIDXSEL]) - HcalEA_ptscale((*phoSCEta)[mitIDXSEL], uncorrectedPhoEt); 
                Float_t phoTkrIsoCorr        = (*phoTrkSumPtHollowConeDR03)[mitIDXSEL] - rho * TkrEffAreas((*phoSCEta)[mitIDXSEL]);

                timing_beamhalo->Fill((*phoSeedTime)[mitIDXSEL]);

                for (int icell = 0; icell < nAllCellsEB; icell++)
                {       
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mitIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && (*phoSeedIEta)[mitIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mitIDXSEL] - 1)
                    {               
                        cellsEB_left.push_back(icell);
                    }       

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && (*phoSeedIEta)[mitIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 
                    }  
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);
                

                    if (phoPFECALClusIsoCorr < 4.78 &&
                    phoPFHCALClusIsoCorr < 6.40 &&
                    phoTkrIsoCorr < 0.89)
                    {
                        if (Ewing > 0.01)
                        {
                            timing_beamhalo_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                            beamhalo_sieie->Fill((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL]); 
                        }
                    }

                    
                }

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);


                    if (phoPFECALClusIsoCorr < 4.78 &&
                    phoPFHCALClusIsoCorr < 6.40 &&
                    phoTkrIsoCorr < 0.89)
                    {
                        if (Ewing > 0.01)
                        {
                            timing_beamhalo_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                            beamhalo_sieie->Fill((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL]); 
                        }
                    }

                
                }

                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);


                    if (phoPFECALClusIsoCorr < 4.78 &&
                    phoPFHCALClusIsoCorr < 6.40 &&
                    phoTkrIsoCorr < 0.89)
                    {
                        if (Ewing > 0.01)
                        {
                            timing_beamhalo_iso_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                            beamhalo_sieie->Fill((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL]); 
                        }
                    }

                }
            }



            //Case 4, only OOT
            Float_t NewmMET = newmMET(pfMET, pfMETPhi, (*ophoPhi)[mootIDXSEL], (*ophoEt)[mootIDXSEL]);

            //Candidate Events
            if (NewmMET > 210 && (*ophohasPixelSeed)[mootIDXSEL] == 0 && metFilters == 0 &&
            (*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 &&
            //(*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 &&
            //(*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001 &&
            (*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] < 0.0105)//106X 256 for passing everything
            {
              
                Float_t uncorrectedPhoEt = ((*ophoEt)[mootIDXSEL]);
                uncorrectedPhoEt = uncorrectedPhoEt - 0.015 * uncorrectedPhoEt;
                Float_t ophoPFECALClusIsoCorr = (*ophoPFClusEcalIso)[mootIDXSEL] - rho*EcalEA((*ophoSCEta)[mootIDXSEL]) - EcalEA_ptscale((*ophoSCEta)[mootIDXSEL], uncorrectedPhoEt); 
                Float_t ophoPFHCALClusIsoCorr = (*ophoPFClusHcalIso)[mootIDXSEL] - rho*HcalEA((*ophoSCEta)[mootIDXSEL]) - HcalEA_ptscale((*ophoSCEta)[mootIDXSEL], uncorrectedPhoEt); 
                Float_t ophoTkrIsoCorr        = (*ophoTrkSumPtHollowConeDR03)[mootIDXSEL] - rho * TkrEffAreas((*ophoSCEta)[mootIDXSEL]);


                //cand1 searching adjencent cells
                for (int icell = 0; icell < nAllCellsEB; icell++)
                {                    
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] - 1)
                    {                    
                        cellsEB_left.push_back(icell);
                    }            

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 
                    }        
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

            

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == ootitmatched[0])
                            {
                                //cal aboslute timing diff
                                float case4_abstimediff = abs((*ophoSeedTime)[mootIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (case4_abstimediff > case4_initimediff)
                                {
                                    case4_initimediff = case4_abstimediff;
                                    lictdcase4 = ((*ophoSeedTime)[mootIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    if (nphotons == 1 && nophotons == 1)
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (ophoPFECALClusIsoCorr < 4.78 &&
                            ophoPFHCALClusIsoCorr < 6.40 &&
                            ophoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);
                                }

                                //matrix method
                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);
                                    }

                                    //matrix method
                                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (ophoPFECALClusIsoCorr < 4.78 &&
                        ophoPFHCALClusIsoCorr < 6.40 &&
                        ophoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);
                            }

                            //matrix method
                            if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*ophoSeedTime)[mootIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);
                                    }

                                    //matrix method
                                    if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*ophoSeedTime)[mootIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                
                }

                            

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

            

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == ootitmatched[0])
                            {
                                //cal aboslute timing diff
                                float case4_abstimediff = abs((*ophoSeedTime)[mootIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (case4_abstimediff > case4_initimediff)
                                {
                                    case4_initimediff = case4_abstimediff;
                                    lictdcase4 = ((*ophoSeedTime)[mootIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    if (nphotons == 1 && nophotons == 1)
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (ophoPFECALClusIsoCorr < 4.78 &&
                            ophoPFHCALClusIsoCorr < 6.40 &&
                            ophoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);
                                }

                                //matrix method
                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);
                                    }

                                    //matrix method
                                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (ophoPFECALClusIsoCorr < 4.78 &&
                        ophoPFHCALClusIsoCorr < 6.40 &&
                        ophoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);
                            }

                            //matrix method
                            if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*ophoSeedTime)[mootIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);
                                    }

                                    //matrix method
                                    if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*ophoSeedTime)[mootIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                
                }


                            
                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

            

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == ootitmatched[0])
                            {
                                //cal aboslute timing diff
                                float case4_abstimediff = abs((*ophoSeedTime)[mootIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (case4_abstimediff > case4_initimediff)
                                {
                                    case4_initimediff = case4_abstimediff;
                                    lictdcase4 = ((*ophoSeedTime)[mootIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    if (nphotons == 1 && nophotons == 1)
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (ophoPFECALClusIsoCorr < 4.78 &&
                            ophoPFHCALClusIsoCorr < 6.40 &&
                            ophoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);
                                }

                                //matrix method
                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);
                                    }

                                    //matrix method
                                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (ophoPFECALClusIsoCorr < 4.78 &&
                        ophoPFHCALClusIsoCorr < 6.40 &&
                        ophoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);
                            }

                            //matrix method
                            if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*ophoSeedTime)[mootIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase4);
                                    }

                                    //matrix method
                                    if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*ophoSeedTime)[mootIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                
                }
            }


            //BeamHalo template
			if (NewmMET > 210 &&
				(*ophohasPixelSeed)[mootIDXSEL] == 0 &&
				(*ophoMIPTotEnergy)[mootIDXSEL] > 4.9 &&
				//(*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 &&
				//(*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001
				metFilters == 8)//106X beamhalo flag
			{
                
                
                Float_t uncorrectedoPhoEt = ((*ophoEt)[mootIDXSEL]);
                uncorrectedoPhoEt = uncorrectedoPhoEt - 0.015 * uncorrectedoPhoEt;
                Float_t ophoPFECALClusIsoCorr = (*ophoPFClusEcalIso)[mootIDXSEL] - rho*EcalEA((*ophoSCEta)[mootIDXSEL]) - EcalEA_ptscale((*ophoSCEta)[mootIDXSEL], uncorrectedoPhoEt); 
                Float_t ophoPFHCALClusIsoCorr = (*ophoPFClusHcalIso)[mootIDXSEL] - rho*HcalEA((*ophoSCEta)[mootIDXSEL]) - HcalEA_ptscale((*ophoSCEta)[mootIDXSEL], uncorrectedoPhoEt); 
                Float_t ophoTkrIsoCorr        = (*ophoTrkSumPtHollowConeDR03)[mootIDXSEL] - rho * TkrEffAreas((*ophoSCEta)[mootIDXSEL]);

                timing_beamhalo->Fill((*ophoSeedTime)[mootIDXSEL]);

                for (int icell = 0; icell < nAllCellsEB; icell++)
                {       
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] - 1)
                    {               
                        cellsEB_left.push_back(icell);
                    }       

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 
                    }  
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);
                
                    if (ophoPFECALClusIsoCorr < 4.78 &&
                    ophoPFHCALClusIsoCorr < 6.40 &&
                    ophoTkrIsoCorr < 0.89)
                    {
                        if (Ewing > 0.01)
                        {
                            timing_beamhalo_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                            beamhalo_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL]); 
                        }
                    }

                    
                }

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                    if (ophoPFECALClusIsoCorr < 4.78 &&
                    ophoPFHCALClusIsoCorr < 6.40 &&
                    ophoTkrIsoCorr < 0.89)
                    {
                        if (Ewing > 0.01)
                        {
                            timing_beamhalo_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                            beamhalo_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL]); 
                        }
                    }

             
                }

                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);


                    if (ophoPFECALClusIsoCorr < 4.78 &&
                    ophoPFHCALClusIsoCorr < 6.40 &&
                    ophoTkrIsoCorr < 0.89)
                    {
                        if (Ewing > 0.01)
                        {
                            timing_beamhalo_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                            beamhalo_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL]); 
                        }
                    }

                    
                }
            }
            


            //Case 5, it + oot but take oot
            //Candidate Events
            if (pfMET > 210 &&
                (*phohasPixelSeed)[mitIDXSEL] == 0 &&
                (*phoMIPTotEnergy)[mitIDXSEL] < 4.9 && (*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 &&
                metFilters == 0 && //106X 256 for passing everything&&
                NewmMET > 210 &&
                //(*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 &&
                //(*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001 &&
                (*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.0105 &&
                ///(*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 &&
                //(*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001 &&
                (*ophoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.0105 &&
                (*ophohasPixelSeed)[mootIDXSEL] == 0)
            {
                
                
                Float_t uncorrectedPhoEt = ((*ophoEt)[mootIDXSEL]);
                uncorrectedPhoEt = uncorrectedPhoEt - 0.015 * uncorrectedPhoEt;
                Float_t ophoPFECALClusIsoCorr = (*ophoPFClusEcalIso)[mootIDXSEL] - rho*EcalEA((*ophoSCEta)[mootIDXSEL]) - EcalEA_ptscale((*ophoSCEta)[mootIDXSEL], uncorrectedPhoEt); 
                Float_t ophoPFHCALClusIsoCorr = (*ophoPFClusHcalIso)[mootIDXSEL] - rho*HcalEA((*ophoSCEta)[mootIDXSEL]) - HcalEA_ptscale((*ophoSCEta)[mootIDXSEL], uncorrectedPhoEt); 
                Float_t ophoTkrIsoCorr        = (*ophoTrkSumPtHollowConeDR03)[mootIDXSEL] - rho * TkrEffAreas((*ophoSCEta)[mootIDXSEL]);


                        
                //cand1 searching adjencent cells
                for (int icell = 0; icell < nAllCellsEB; icell++)
                {                    
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] - 1)
                    {                    
                        cellsEB_left.push_back(icell);
                    }            

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 
                    }        
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

            

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == ootitmatched[0])
                            {
                                //cal aboslute timing diff
                                float case5_abstimediff = abs((*ophoSeedTime)[mootIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (case5_abstimediff > case5_initimediff)
                                {
                                    case5_initimediff = case5_abstimediff;
                                    lictdcase5 = ((*ophoSeedTime)[mootIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    if (nphotons == 1 && nophotons == 1)
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (ophoPFECALClusIsoCorr < 4.78 &&
                            ophoPFHCALClusIsoCorr < 6.40 &&
                            ophoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);
                                }

                                //matrix method
                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);
                                    }

                                    //matrix method
                                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (ophoPFECALClusIsoCorr < 4.78 &&
                        ophoPFHCALClusIsoCorr < 6.40 &&
                        ophoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);
                            }

                            //matrix method
                            if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*ophoSeedTime)[mootIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);
                                    }

                                    //matrix method
                                    if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*ophoSeedTime)[mootIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                
                }


                         
                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

            

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == ootitmatched[0])
                            {
                                //cal aboslute timing diff
                                float case5_abstimediff = abs((*ophoSeedTime)[mootIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (case5_abstimediff > case5_initimediff)
                                {
                                    case5_initimediff = case5_abstimediff;
                                    lictdcase5 = ((*ophoSeedTime)[mootIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    if (nphotons == 1 && nophotons == 1)
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (ophoPFECALClusIsoCorr < 4.78 &&
                            ophoPFHCALClusIsoCorr < 6.40 &&
                            ophoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);
                                }

                                //matrix method
                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);
                                    }

                                    //matrix method
                                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (ophoPFECALClusIsoCorr < 4.78 &&
                        ophoPFHCALClusIsoCorr < 6.40 &&
                        ophoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);
                            }

                            //matrix method
                            if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*ophoSeedTime)[mootIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);
                                    }

                                    //matrix method
                                    if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*ophoSeedTime)[mootIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                
                }

                           

                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

            

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == ootitmatched[0])
                            {
                                //cal aboslute timing diff
                                float case5_abstimediff = abs((*ophoSeedTime)[mootIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (case5_abstimediff > case5_initimediff)
                                {
                                    case5_initimediff = case5_abstimediff;
                                    lictdcase5 = ((*ophoSeedTime)[mootIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    if (nphotons == 1 && nophotons == 1)
                    {
                        ewing_candidate->Fill(Ewing);

                        //timing
                        if (ophoPFECALClusIsoCorr < 4.78 &&
                            ophoPFHCALClusIsoCorr < 6.40 &&
                            ophoTkrIsoCorr < 0.89
                            )
                            {
                                ewing_candidate_iso->Fill(Ewing);

                                timing_candidate_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                timing_vs_LICTD_candidate_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);

                                if (Ewing > 0.01)
                                {
                                    timing_candidate_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_candidate_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);
                                }

                                //matrix method
                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    ewing_candidate_3ns_iso->Fill(Ewing);

                                    timing_candidate_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }

                                    if (Ewing > 0.02)
                                    {
                                        timing_candidate_3ns_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }

                                    if (Ewing > 0.03)
                                    {
                                        timing_candidate_3ns_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                    }
                                }
                            }

                        if (Ewing > 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_candidate_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_candidate_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001)
                        {
                            ewing_candidate_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_candidate_ss_iso->Fill(Ewing);

                                    timing_candidate_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_candidate_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);

                                    if (Ewing > 0.01)
                                    {
                                        timing_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        timing_vs_LICTD_candidate_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);
                                    }

                                    //matrix method
                                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                    {
                                        ewing_candidate_ss_3ns_iso->Fill(Ewing);

                                        timing_candidate_ss_3ns_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        if (Ewing > 0.01)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_candidate_ss_3ns_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }
                                    }
                                }

                            if (Ewing > 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_candidate_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_candidate_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                    }

                    //spike timing
                    ewing_spike->Fill(Ewing);

                    if (ophoPFECALClusIsoCorr < 4.78 &&
                        ophoPFHCALClusIsoCorr < 6.40 &&
                        ophoTkrIsoCorr < 0.89
                        )
                        {
                            ewing_spike_iso->Fill(Ewing);

                            timing_spike_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                            timing_vs_LICTD_spike_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);

                            if (Ewing < 0.01)
                            {
                                timing_spike_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                timing_vs_LICTD_spike_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);
                            }

                            //matrix method
                            if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                            {
                                ewing_spike_early_iso->Fill(Ewing);

                                timing_spike_early_iso->Fill((*ophoSeedTime)[mootIDXSEL]);


                                if (Ewing > 0.01)
                                {
                                    timing_spike_early_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }

                                if (Ewing > 0.02)
                                {
                                    timing_spike_early_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }

                                if (Ewing > 0.03)
                                {
                                    timing_spike_early_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                        }

                        if (Ewing < 0.01)
                        {
                            if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_HCAL->Fill(ophoPFHCALClusIsoCorr);
                            }

                            if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                            {
                                ewing_spike_ECAL->Fill(ophoPFECALClusIsoCorr);
                            }

                            if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                            {
                                ewing_spike_track->Fill(ophoTkrIsoCorr);
                            }
                        }

                        if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001)
                        {
                            ewing_spike_ss->Fill(Ewing);

                            //timing
                            if (ophoPFECALClusIsoCorr < 4.78 &&
                                ophoPFHCALClusIsoCorr < 6.40 &&
                                ophoTkrIsoCorr < 0.89
                                )
                                {
                                    ewing_spike_ss_iso->Fill(Ewing);

                                    timing_spike_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL]);

                                    timing_vs_LICTD_spike_ss_iso->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);

                                    if (Ewing < 0.01)
                                    {
                                        timing_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                        timing_vs_LICTD_spike_ss_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL], lictdcase5);
                                    }

                                    //matrix method
                                    if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                                    {
                                        ewing_spike_ss_early_iso->Fill(Ewing);

                                        timing_spike_ss_early_iso->Fill((*ophoSeedTime)[mootIDXSEL]);


                                        if (Ewing > 0.01)
                                        {
                                            timing_spike_ss_early_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.02)
                                        {
                                            timing_spike_ss_early_iso_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }

                                        if (Ewing > 0.03)
                                        {
                                            timing_spike_ss_early_iso_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                                        }
                                    }
                                }


                            if (Ewing < 0.01)
                            {
                                if (ophoPFECALClusIsoCorr < 4.78 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_HCAL->Fill(ophoPFHCALClusIsoCorr);
                                }

                                if (ophoPFHCALClusIsoCorr < 6.40 && ophoTkrIsoCorr < 0.89)
                                {
                                    ewing_spike_ss_ECAL->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 4.78 && ophoPFHCALClusIsoCorr < 6.40)
                                {
                                    ewing_spike_ss_track->Fill(ophoTkrIsoCorr);
                                }
                            }
                        }
                
                }
            }

            //BeamHalo template
			if (pfMET > 210 &&
				(*phohasPixelSeed)[mitIDXSEL] == 0 &&
				(*phoMIPTotEnergy)[mitIDXSEL] > 4.9 &&
				//(*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 &&
				//(*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001 &&
                //(*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 &&
				//(*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001
				metFilters == 8 && //106X beamhalo flag&&
				NewmMET > 210 &&
				(*ophohasPixelSeed)[mootIDXSEL] == 0 &&
				(*ophoMIPTotEnergy)[mootIDXSEL] > 4.9
				)
			{
                
                
                Float_t uncorrectedoPhoEt = ((*ophoEt)[mootIDXSEL]);
                uncorrectedoPhoEt = uncorrectedoPhoEt - 0.015 * uncorrectedoPhoEt;
                Float_t ophoPFECALClusIsoCorr = (*ophoPFClusEcalIso)[mootIDXSEL] - rho*EcalEA((*ophoSCEta)[mootIDXSEL]) - EcalEA_ptscale((*ophoSCEta)[mootIDXSEL], uncorrectedoPhoEt); 
                Float_t ophoPFHCALClusIsoCorr = (*ophoPFClusHcalIso)[mootIDXSEL] - rho*HcalEA((*ophoSCEta)[mootIDXSEL]) - HcalEA_ptscale((*ophoSCEta)[mootIDXSEL], uncorrectedoPhoEt); 
                Float_t ophoTkrIsoCorr        = (*ophoTrkSumPtHollowConeDR03)[mootIDXSEL] - rho * TkrEffAreas((*ophoSCEta)[mootIDXSEL]);

                
                timing_beamhalo->Fill((*ophoSeedTime)[mootIDXSEL]);

                for (int icell = 0; icell < nAllCellsEB; icell++)
                {       
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] - 1)
                    {               
                        cellsEB_left.push_back(icell);
                    }       

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 
                    }  
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);
                

                    if (ophoPFECALClusIsoCorr < 4.78 &&
                    ophoPFHCALClusIsoCorr < 6.40 &&
                    ophoTkrIsoCorr < 0.89)
                    {
                        if (Ewing > 0.01)
                        {
                            timing_beamhalo_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                            beamhalo_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL]); 
                        }
                    }

                    
                }

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);


                    if (ophoPFECALClusIsoCorr < 4.78 &&
                    ophoPFHCALClusIsoCorr < 6.40 &&
                    ophoTkrIsoCorr < 0.89)
                    {
                        if (Ewing > 0.01)
                        {
                            timing_beamhalo_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                            beamhalo_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL]); 
                        }
                    }

                }

                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);


                    if (ophoPFECALClusIsoCorr < 4.78 &&
                    ophoPFHCALClusIsoCorr < 6.40 &&
                    ophoTkrIsoCorr < 0.89)
                    {
                        if (Ewing > 0.01)
                        {
                            timing_beamhalo_iso_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                            beamhalo_sieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL]); 
                        }
                    }
                }
            
            }
        
        }
    
        
    }
        


    /*
    =================================================================================
    =================================================================================
    =================================================================================
    =================================================================================
    =================================================================================
    */
    
    ewing_candidate->Write();
    ewing_candidate_iso->Write();
    timing_candidate_iso->Write();
    timing_candidate_iso_ewing1->Write();

    timing_vs_LICTD_candidate_iso->SetOption("COLZ");
    timing_vs_LICTD_candidate_iso->Write();

    timing_vs_LICTD_candidate_iso_ewing1->SetOption("COLZ");
    timing_vs_LICTD_candidate_iso_ewing1->Write();

    ewing_candidate_3ns_iso->Write();
    timing_candidate_3ns_iso->Write();

    timing_candidate_3ns_iso_ewing1->Write();
    timing_candidate_3ns_iso_ewing2->Write();
    timing_candidate_3ns_iso_ewing3->Write();

    ewing_candidate_HCAL->Write();
    ewing_candidate_ECAL->Write();
    ewing_candidate_track->Write();

    ewing_candidate_ss->Write();
    ewing_candidate_ss_iso->Write();
    timing_candidate_ss_iso->Write();
    timing_candidate_ss_iso_ewing1->Write();

    timing_vs_LICTD_candidate_ss_iso->SetOption("COLZ");
    timing_vs_LICTD_candidate_ss_iso->Write();

    timing_vs_LICTD_candidate_ss_iso_ewing1->SetOption("COLZ");
    timing_vs_LICTD_candidate_ss_iso_ewing1->Write();

    ewing_candidate_ss_3ns_iso->Write();
    timing_candidate_ss_3ns_iso->Write();

    timing_candidate_ss_3ns_iso_ewing1->Write();
    timing_candidate_ss_3ns_iso_ewing2->Write();
    timing_candidate_ss_3ns_iso_ewing3->Write();

    ewing_candidate_ss_HCAL->Write();
    ewing_candidate_ss_ECAL->Write();
    ewing_candidate_ss_track->Write();



    ewing_spike->Write();
    ewing_spike_iso->Write();
    timing_spike_iso->Write();
    timing_spike_iso_ewing1->Write();


    timing_vs_LICTD_spike_iso->SetOption("COLZ");
    timing_vs_LICTD_spike_iso->Write();
    
    timing_vs_LICTD_spike_iso_ewing1->SetOption("COLZ");
    timing_vs_LICTD_spike_iso_ewing1->Write();
    
    ewing_spike_early_iso->Write();
    timing_spike_early_iso->Write();

    timing_spike_early_iso_ewing1->Write();
    timing_spike_early_iso_ewing2->Write();
    timing_spike_early_iso_ewing3->Write();

    ewing_spike_HCAL->Write();
    ewing_spike_ECAL->Write();
    ewing_spike_track->Write();

    ewing_spike_ss->Write();
    ewing_spike_ss_iso->Write();

    timing_spike_ss_iso->Write();
    timing_spike_ss_iso_ewing1->Write();

    timing_vs_LICTD_spike_ss_iso->SetOption("COLZ");
    timing_vs_LICTD_spike_ss_iso->Write();

    timing_vs_LICTD_spike_ss_iso_ewing1->SetOption("COLZ");
    timing_vs_LICTD_spike_ss_iso_ewing1->Write();

    ewing_spike_ss_early_iso->Write();
    timing_spike_ss_early_iso->Write();

    timing_spike_ss_early_iso_ewing1->Write();
    timing_spike_ss_early_iso_ewing2->Write();
    timing_spike_ss_early_iso_ewing3->Write();

    ewing_spike_ss_HCAL->Write();
    ewing_spike_ss_ECAL->Write();
    ewing_spike_ss_track->Write();

    InvMass_Z->SetLineColor(1);
	InvMass_Z->Write();

    
    ewing_promptZ->Write();
    ewing_promptZ_iso->Write();
    timing_promptZ_iso->Write();
    timing_promptZ_iso_ewing1->Write();

    ewing_promptZ_3ns_iso->Write();
    timing_promptZ_3ns_iso->Write();

    timing_promptZ_3ns_iso_ewing1->Write();
    timing_promptZ_3ns_iso_ewing2->Write();
    timing_promptZ_3ns_iso_ewing3->Write();

    ewing_promptZ_HCAL->Write();
    ewing_promptZ_ECAL->Write();
    ewing_promptZ_track->Write();

    timing_beamhalo->Write();
    timing_beamhalo_iso_ewing1->Write();


    beamhalo_sieie->Write();

    promptz_sieie->Write();

    promptz_sieie_vs_MIPE->SetOption("COLZ");
    promptz_sieie_vs_MIPE->Write();

    promptz_sieie_vs_MIPE_profile->SetLineColor(1);
	promptz_sieie_vs_MIPE_profile->Write();

    cout << "monopho17UL_cutbased_1118B.root" << endl;



}

/*


========For timing fit========
candidate_events
promptZ_events
spike_events --> OLD spike template for timing fit
BeamHalo_template
timing_spikeincandnotimingcut_ewing1 
--> NEW spike template, candidate selection with ewing < 0.01, contains single and double spikes


========For ewing method========
candidate_eta_wing
candidate_eta_wing_prompt
timing_candidate_ewing1_allcases
timing_candidate_ewing3_allcases


promptZ_eta_wing
promptZ_eta_wing_prompt
timing_promptZ_ewing1_allcases
timing_promptZ_ewing3_allcases


spikeincandidate_eta_wing 
--> OLD spike template, requires in the candidate sample with timing < -12.5 ns
timing_spikeincand_ewing1_allcases
timing_spikeincand_ewing3_allcases


ewing_singlespike_allcases 
--> NEW spike template, requires single spike (lictd < -10) in candidate event and timing < -12.5 ns
timing_singlespike_ewing1_allcases
timing_singlespike_ewing3_allcases


ewing_doublespike_allcases
--> also check this double spike case (lictd > -10) with timing < -12.5 ns in candidate sample
timing_doublespike_ewing3_allcases
timing_doublespike_ewing1_allcases



12.04.2020
to do:
clean the code for timing fit only.
check the ewing distribution for new spike template and candidate event(without |3|ns selection)
candidate event with ewing > 0.01 cut. Compare it with new spike template and do timing fit again.
ewing distribution from single and double spike in can and spike tempalte?



12.08.2020
added beamhalo cross check (p.39)


12.14.2020
fixed the ewing cuts for new spike



01292021
adding one photon selection test

02022021
looking at the number of photon with pt > 225
===============================================
04302021
this is for the original template (with BDT cuts)
===============================================


05032021
adding BDT score distribution on every template


05062021
applying BDT score to candidate only, not on templates

05072021
adding ECAL HCAL TRACK ISO for timing only, to see how it affects


05252021
candidate event selecting ewing > 0.01 for timing fit

apply muon and electron veto -> Gulia's slides (https://indico.cern.ch/event/1016471/contributions/4266383/attachments/2204447/3729507/Monophoton_09_03_2021.pdf)

maybe next time -> tau veto

05182021
phoet 230 GeV -> 225 GeV

05212021
eta wing study
candidate, spike, prompt ewing > 0.01
plot ewing distribution in each cut
plot timing distribution with ewing > 0.01 & in each cut

05242021
fixed the ewing filled problem
fixed the ewing syntax
added a timing for spike with single spike selected but without ewing

05252021
H/E cut from 0.02197 -> 0.042

05282021
added InvMass of Z for different photon ID

06062021
removed sieie sipip > 0.001 cuts for candidate and templates

06072021
spike estimate using timing fit and ewing

06142021
added spike pt distribution
removed (*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.01055 since we have new photon ID for this
-> need to be discussed the result, keep it for now

06152021
the late time (5 - 10 ns) contribution without the sieie < 0.01055 looks like BH.
To further study it, make:

1. eta and phi distribution of those guys
2. seedtime vs eta 

other idea is to run the analysis reuqires only one photon -> this runs separately


06162021
fixed the bug which I should not have the cuts on 0.01015


06172021
late time halo studies:
Discussion in the meeting:

Why do these halos pass the BDT (but others do not)?
What can we check:

1. BDT score for those halo events and the template Halo region, compare the BDT score between them

2. Those guys vs the close to in-time halo	


==== 0622 ====
veto the photons for candidate only
studying the late time stuff:
1. n - 1 plot for each BDT scores
2. n - 1 vs sieie

in the region of
a. late time 5ns - 10ns in candidate events
b. late time 5ns - 10ns in beamhalo template
c. close to prompt window -5ns - -3ns in beamhalo template
d. early time -15ns - -10ns in beamhalo template

==== 07262021 ====
fixed the single spike selected for ewing, just the timing


==== 07272021 ====
add the old photon id selection plots to make comparison
this version is for 09Aug2019_UL2017-v1

the selection that we removed the sieie && sipip for templates and candidate sample
we use ewing cuts > 0.01 for candidate, prompt and beamhalo while ewing < 0.01 for spike

for timing fit, take the candidate, prompt and beamhalo with ewing > 0.01 and spike < 0.01
for ewing matrix method, take the candidate, prompt, spike with ewing > 0.01


histogram for studies:
==== timing fit ====
---- old ----
candidate - timing_ewing_candidate_OLDIDs
spike - timing_spike_OLDIDs
promptZ - timing_ewing_promptZ_OLDIDs
beamhalo - timing_beamhalo_iso_ewing1

---- new ----
candidate - timing_ewing_candidate_BDT_ECAL_HCAL_TRACK
spike - timing_spike_BDT_ECALISO_HCALISO_TRACKISO / timing_spike
promptZ - timing_ewing_promptZ_BDT_ECAL_HCAL_TRACK / timing_ewing_promptZ
beamhalo - timing_beamhalo_iso_BDT_ewing1_ECALISO_HCALISO_TRACKISO / timing_beamhalo_iso_BDT_ewing1


==== ewing matrix ====
---- old ----
candidate - ewing_candidate_3ns_OLDIDs, timing_ewing_candidate_3ns_OLDIDs
spike - ewing_spike_OLDIDs, timing_ewing_spike_OLDIDs
promptZ - ewing_promptZ_OLDIDs, timing_ewing_promptZ_OLDIDs


---- new ----
candidate - ewing_candidate_3ns_BDT_ECAL_HCAL_TRACK, timing_ewing_candidate_BDT_ECAL_HCAL_TRACK
spike - ewing_spike, timing_ewing_spike
promptZ - ewing_promptZ, timing_ewing_promptZ



==== 072820201 ====
added the ewing of candidate in prompt window


==== 07302021 ====
apply isolation cuts everywhere, the only difference between old and new selections are sieie vs BDT
add ewing > 0.02, 0.03 for ewing matrix method

==== 08022021 ====
added the beamhalo isolation cuts
added the prompt and spike before BDT cuts timing distribution

separated spike and candidate which one photon requirement is only for candidate
fixed case5 to select both in-time and out-of-time
added the old selection sieie/sipip > / < 0.001 cuts, to check the old results

==== 08032021 ====
added the LICTD plots for candidate > 0.01 with old selections and spike < 0.01

fixed the ewing method for spike is the same selection to candidate but in the early time


=== 08092021 ===
added single and double spike template
double spike: -17 < time < -9, -6 < LICTD < 0
single spike: -16 < time < -10, -20 < LICTD < -12 

=== 081020201 ===
removed the double/single spikes timing requirements.

=== 08122021 ===
added a sieie && sipip > 0.001 for candidate event
keep the templates the same but just to add sieie && sipip cut to candidate
and do the fitting


=== 08162021 ===
added sieie && sipip > 0.001 for spike template
specifically, this is for the ewing matrix method, this will change the efficiency of the spike
I applied it for both timing and ewing just to make sure, ideally, this cuts do not affect timing


=== 08172021 ===
added mip cut < 4.9, CSC tag on beamhalo and compare (suppress the peak just before the prompt?)


=== 09142021 ===
keep the old ID

make a comparison between old and new
with the phoCalibET for photon and phoEt for oot in the BDT score calculation


11182021
R9 distirbution from the spike template (early)
photon ET distribution from the spike template (early)

*/



