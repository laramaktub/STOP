// C / C++ headers
#include <cmath>
#include <iostream>
using namespace std;

#include <algorithm>
#include "TColor.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPolyLine.h"
#include "TROOT.h"
#include "TString.h"
#include "TPRegexp.h"

// ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TBranch.h>
#include <set>
#include <string>
#include <TLorentzVector.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
// Define format and input file
#include "Reader_prefinal1024.h" 


using namespace TMVA;

float triggereff(float, float,float,float,float);

// ###################
// #  Main function  #
// ###################

int main (int argc, char *argv[])
{

  // ################################
  // ##       Open the tree        ##
  // ################################

     

   TH1D* BDT1 = new TH1D("BDT1",",BDT1",80,-1,1) ;
   TH1D* BDT2 = new TH1D("BDT2",",BDT2",80,-1,1) ;
   TH1D* BDT3 = new TH1D("BDT3",",BDT3",80,-1,1) ;
   TH1D* BDT4 = new TH1D("BDT4",",BDT4",80,-1,1) ;
   TH1D* BDT5 = new TH1D("BDT5",",BDT5",80,-1,1) ;
   TH1D* BDT6 = new TH1D("BDT6",",BDT6",80,-1,1) ;

   TFile *fin = TFile::Open(argv[1]);
   TTree* theInputTree = (TTree*) fin->Get("babyTuple"); 
   TFile *fout = new TFile(argv[2],"recreate");

   Float_t met_phi;
   Float_t met_signif;
   Float_t phi_lepton_met;
   Float_t lepton_pT;
   Float_t lepton_phi;
   Float_t nbjets;
   Float_t mlpjpa;
   Float_t dphiljpa;
   Float_t dphiljb;
   Float_t dphimetl;
   Float_t missmass;
   Float_t mjetlepmax;
   Float_t m3;
   Float_t m3_pt;
   Float_t mi34;
   Float_t mi34b;
   Float_t mi34b_m2;
   Float_t m3topW ;
   Float_t jet1_pT ;
   Float_t jet2_pT ;
   Float_t jet3_pT ;
   Float_t jet4_pT ;
   Float_t jet1_phi ;
   Float_t jet2_phi ;
   Float_t jet3_phi ;
   Float_t jet4_phi ;
    Float_t jet1_eta ;
   Float_t jet2_eta ;
   Float_t jet3_eta ;
   Float_t jet4_eta ;
   
  
 Float_t njets_JESup;
 Float_t njets_JESdown;
 Float_t puweight;
 Float_t xsweight;
 Float_t nvertex;
 Float_t leptonflav;
 Float_t trweight;
 Float_t lepton_eta;

 
   Float_t b2_pt;
   Float_t b1_eta;
   Float_t b1plusb2_pt;
   Float_t dR_b1b2;
   Float_t dPhi_b1b2;
 
   Float_t b1_pt_2;
   Float_t b2_pt_2;
   Float_t b1_eta_2;
   Float_t b1plusb2_pt_2;
   Float_t dR_b1b2_2;
   Float_t dPhi_b1b2_2;
   Float_t mlb1_2;
 
   Float_t HT_av;
   Float_t HT_MET_lep_pt;
   Float_t HT;
 
 
   Float_t mlb_hemi;
   Float_t HTfrac;
   Float_t HTfrac_FNAL;
   Float_t chargino_mixing;
   Float_t mT2;
   Float_t mT2bl;
   Float_t mT2W;
   Float_t dR_LepB;
   Float_t dR_LepJet;
   Float_t dPhi_JetMet;
   Float_t METoverSqrtHT;
   Float_t Chi2SNT;
   
   Float_t met;
   Float_t mT;
   Float_t b1_pt;
   Float_t m3b;
   Float_t njets;
   Float_t stop_mass;
   Float_t lsp_mass;
   
   
   Float_t bdt_R1;
   Float_t bdt_R2;
   Float_t bdt_R3;
   Float_t bdt_R4;
   Float_t bdt_R5;
   Float_t bdt_R6;
   
    TTree *BDTtree     = new TTree("BDTtree","Tree of variables");
  BDTtree->Branch("met", &met, "met/F");
  BDTtree->Branch("mT", &mT, "mT/F");
  BDTtree->Branch("mT2W", &mT2W, "mT2W/F");
  BDTtree->Branch("HTfrac", &HTfrac, "HTfrac/F");
  BDTtree->Branch("b1_pt", &b1_pt, "b1_pt/F");
  BDTtree->Branch("lepton_pT", &lepton_pT, "lepton_pT/F");
  BDTtree->Branch("dPhi_JetMet", &dPhi_JetMet, "dPhi_JetMet/F");
  BDTtree->Branch("dR_LepB", &dR_LepB, "dR_LepB/F");
  BDTtree->Branch("m3b", &m3b, "m3b/F");
  BDTtree->Branch("mlb_hemi", &mlb_hemi, "mlb_hemi/F");
  BDTtree->Branch("jet1_pT", &jet1_pT, "jet1_pT/F");
  BDTtree->Branch("jet2_pT", &jet2_pT, "jet2_pT/F");
  BDTtree->Branch("jet3_pT", &jet3_pT, "jet3_pT/F");
  BDTtree->Branch("jet1_eta", &jet1_eta, "jet1_eta/F");
  BDTtree->Branch("jet2_eta", &jet2_eta, "jet2_eta/F");
  BDTtree->Branch("jet3_eta", &jet3_eta, "jet3_eta/F");
  BDTtree->Branch("jet1_phi", &jet1_phi, "jet1_phi/F");
  BDTtree->Branch("jet2_phi", &jet2_phi, "jet2_phi/F");
  BDTtree->Branch("jet3_phi", &jet3_phi, "jet3_phi/F");
  BDTtree->Branch("njets", &njets, "njets/F");
  BDTtree->Branch("METoverSqrtHT", &METoverSqrtHT, "METoverSqrtHT/F");
  BDTtree->Branch("dR_LepJet", &dR_LepJet, "dR_LepJet/F");
  BDTtree->Branch("HT_MET_lep_pt", &HT_MET_lep_pt, "HT_MET_lep_pt/F");
  BDTtree->Branch("Chi2SNT", &Chi2SNT, "Chi2SNT/F");
  BDTtree->Branch("stop_mass", &stop_mass, "stop_mass/F");
  BDTtree->Branch("lsp_mass", &lsp_mass, "lsp_mass/F");
  BDTtree->Branch("njets_JESup", &njets_JESup, "njets_JESup/F");
  BDTtree->Branch("njets_JESdown", &njets_JESdown, "njets_JESdown/F");
  BDTtree->Branch("puweight", &puweight, "puweight/F");
  BDTtree->Branch("xsweight", &xsweight, "xsweight/F");
  BDTtree->Branch("nvertex", &nvertex, "nvertex/F");
  BDTtree->Branch("leptonflav", &leptonflav, "leptonflav/F");
  BDTtree->Branch("trweight", &trweight, "trweight/F");
  BDTtree->Branch("lepton_eta", &lepton_eta, "lepton_eta/F");
  BDTtree->Branch("lepton_phi", &lepton_phi, "lepton_phi/F");
  
  BDTtree->Branch("bdt_R1", &bdt_R1, "bdt_R1/F");
  BDTtree->Branch("bdt_R2", &bdt_R2, "bdt_R2/F");
  BDTtree->Branch("bdt_R3", &bdt_R3, "bdt_R3/F");
  BDTtree->Branch("bdt_R4", &bdt_R4, "bdt_R4/F");
  BDTtree->Branch("bdt_R5", &bdt_R5, "bdt_R5/F");
  BDTtree->Branch("bdt_R6", &bdt_R6, "bdt_R6/F");
  
  
   babyEvent myEvent;
   intermediatePointers pointers;
   InitializeBranches(theInputTree,&myEvent,&pointers);
  


 Reader* readerReg1;
   readerReg1= new Reader("V");

    string NN_varsReg1 = string(argv[3]);
    
    
    TString str_NNvariablesReg1 = TString(NN_varsReg1);
    TPMERegexp _variablesReg1(",");

    set<string> NNvariablesReg1;
    int nvariablesReg1 = 0;
 
    nvariablesReg1 = _variablesReg1.Split(str_NNvariablesReg1);
 
    for(int i=0; i<nvariablesReg1; i++){
    NNvariablesReg1.insert(_variablesReg1[i].Data());
    }
    
    
    
 
   Reader* reader; 
   reader = new Reader("V");

    string NN_vars = string(argv[4]);
    TString str_NNvariables = TString(NN_vars);
    TPMERegexp _variables(",");

    set<string> NNvariables;
    int nvariables = 0;
 
    nvariables = _variables.Split(str_NNvariables);
 
    for(int i=0; i<nvariables; i++){
    NNvariables.insert(_variables[i].Data());
    }
 

   cout << "" << endl;
   cout << "**********************************************************" << endl;
   cout << "*   Will use these variables in BDT  * "                    << NN_vars << endl;
   cout << "**********************************************************" << endl;
   cout << "" << endl;


    string setup_directory = string(argv[5]);


   TH1D *Stop = new TH1D("Stop","",1000,0,1000);
   TH1D *Neut = new TH1D("Neut","",1000,0,1000);


  if(NNvariables.find("met")!=NNvariables.end())
    reader->AddVariable("met",&met);
  if(NNvariables.find("mT")!=NNvariables.end())
    reader->AddVariable("mT",&mT);
  if(NNvariables.find("lepton_pT")!=NNvariables.end())
    reader->AddVariable("lepton_pT",&lepton_pT);
  if(NNvariables.find("lepton_phi")!=NNvariables.end())
    reader->AddVariable("lepton_phi",&lepton_phi);
  if(NNvariables.find("njets")!=NNvariables.end())
    reader->AddVariable("njets",&njets);
  if(NNvariables.find("nbjets")!=NNvariables.end())
    reader->AddVariable("nbjets",&nbjets);
  if(NNvariables.find("m3b")!=NNvariables.end())
    reader->AddVariable("m3b",&m3b);
  if(NNvariables.find("phi_lepton_met")!=NNvariables.end())
    reader->AddVariable("phi_lepton_met",&phi_lepton_met);
     if(NNvariables.find("jet1_pT")!=NNvariables.end())
    reader->AddVariable("jet1_pT",&jet1_pT);
    if(NNvariables.find("jet2_pT")!=NNvariables.end())
    reader->AddVariable("jet2_pT",&jet1_pT);
      if(NNvariables.find("jet3_pT")!=NNvariables.end())
    reader->AddVariable("jet3_pT",&jet1_pT);
  if(NNvariables.find("b1_pt")!=NNvariables.end())
    reader->AddVariable("b1_pt",&b1_pt);
 
    if(NNvariables.find("dR_b1b2")!=NNvariables.end())
    reader->AddVariable("dR_b1b2",&dR_b1b2);
  if(NNvariables.find("dR_b1b2_2")!=NNvariables.end())
    reader->AddVariable("dR_b1b2_2",&dR_b1b2_2);
  if(NNvariables.find("dPhi_b1b2")!=NNvariables.end())
    reader->AddVariable("dPhi_b1b2",&dPhi_b1b2);
  if(NNvariables.find("dPhi_b1b2_2")!=NNvariables.end())
    reader->AddVariable("dPhi_b1b2_2",&dPhi_b1b2_2);
  if(NNvariables.find("HT_MET_lep_pt")!=NNvariables.end())
    reader->AddVariable("HT_MET_lep_pt",&HT_MET_lep_pt);
  if(NNvariables.find("mlb_hemi")!=NNvariables.end())
    reader->AddVariable("mlb_hemi",&mlb_hemi);
  if(NNvariables.find("HTfrac")!=NNvariables.end())
    reader->AddVariable("HTfrac",&HTfrac);
  if(NNvariables.find("mT2W")!=NNvariables.end())
    reader->AddVariable("mT2W",&mT2W);
  if(NNvariables.find("dR_LepB")!=NNvariables.end())
    reader->AddVariable("dR_LepB",&dR_LepB);
  if(NNvariables.find("dR_LepJet")!=NNvariables.end())
    reader->AddVariable("dR_LepJet",&dR_LepJet);
  if(NNvariables.find("dPhi_JetMet")!=NNvariables.end())
    reader->AddVariable("dPhi_JetMet",&dPhi_JetMet);
  if(NNvariables.find("METoverSqrtHT")!=NNvariables.end())
    reader->AddVariable("METoverSqrtHT",&METoverSqrtHT);
  if(NNvariables.find("HT")!=NNvariables.end())
    reader->AddVariable("HT",&HT);
    
    if(NNvariables.find("Chi2SNT")!=NNvariables.end())
    reader->AddVariable("Chi2SNT",&Chi2SNT);
    
    
     
      if(NNvariablesReg1.find("met")!=NNvariablesReg1.end())
    readerReg1->AddVariable("met",&met);
    if(NNvariablesReg1.find("mT")!=NNvariablesReg1.end())
    readerReg1->AddVariable("mT",&mT); 
      if(NNvariablesReg1.find("lepton_pT")!=NNvariablesReg1.end())
    readerReg1->AddVariable("lepton_pT",&lepton_pT);
    if(NNvariablesReg1.find("njets")!=NNvariablesReg1.end())
    readerReg1->AddVariable("njets",&njets);
     if(NNvariablesReg1.find("jet1_pT")!=NNvariablesReg1.end())
    readerReg1->AddVariable("jet1_pT",&jet1_pT);
if(NNvariablesReg1.find("jet2_pT")!=NNvariablesReg1.end())
    readerReg1->AddVariable("jet2_pT",&jet2_pT);
if(NNvariablesReg1.find("jet3_pT")!=NNvariablesReg1.end())
    readerReg1->AddVariable("jet3_pT",&jet3_pT);


      if(NNvariablesReg1.find("b1_pt")!=NNvariablesReg1.end())
    readerReg1->AddVariable("b1_pt",&b1_pt);
       if(NNvariablesReg1.find("mlb_hemi")!=NNvariablesReg1.end())
    readerReg1->AddVariable("mlb_hemi",&mlb_hemi);
   if(NNvariablesReg1.find("HTfrac")!=NNvariablesReg1.end())
    readerReg1->AddVariable("HTfrac",&HTfrac);
    if(NNvariablesReg1.find("mT2W")!=NNvariablesReg1.end())
    readerReg1->AddVariable("mT2W",&mT2W);  
    
   
     
  
  if(NNvariablesReg1.find("lepton_phi")!=NNvariablesReg1.end())
    readerReg1->AddVariable("lepton_phi",&lepton_phi);
 
  
     if(NNvariablesReg1.find("dPhi_JetMet")!=NNvariablesReg1.end())
    readerReg1->AddVariable("dPhi_JetMet",&dPhi_JetMet);
     

   
  if(NNvariablesReg1.find("nbjets")!=NNvariablesReg1.end())
    readerReg1->AddVariable("nbjets",&nbjets);
  if(NNvariablesReg1.find("m3b")!=NNvariablesReg1.end())
    readerReg1->AddVariable("m3b",&m3b);
  if(NNvariablesReg1.find("phi_lepton_met")!=NNvariablesReg1.end())
    readerReg1->AddVariable("phi_lepton_met",&phi_lepton_met);
 
    if(NNvariablesReg1.find("dR_b1b2")!=NNvariablesReg1.end())
    readerReg1->AddVariable("dR_b1b2",&dR_b1b2);
  if(NNvariablesReg1.find("dR_b1b2_2")!=NNvariablesReg1.end())
    readerReg1->AddVariable("dR_b1b2_2",&dR_b1b2_2);
  if(NNvariablesReg1.find("dPhi_b1b2")!=NNvariablesReg1.end())
    readerReg1->AddVariable("dPhi_b1b2",&dPhi_b1b2);
  if(NNvariablesReg1.find("dPhi_b1b2_2")!=NNvariablesReg1.end())
    readerReg1->AddVariable("dPhi_b1b2_2",&dPhi_b1b2_2);
  if(NNvariablesReg1.find("HT_MET_lep_pt")!=NNvariablesReg1.end())
    readerReg1->AddVariable("HT_MET_lep_pt",&HT_MET_lep_pt);
 

  
  if(NNvariablesReg1.find("dR_LepB")!=NNvariablesReg1.end())
    readerReg1->AddVariable("dR_LepB",&dR_LepB);
  if(NNvariablesReg1.find("dR_LepJet")!=NNvariablesReg1.end())
    readerReg1->AddVariable("dR_LepJet",&dR_LepJet);
  
  if(NNvariablesReg1.find("METoverSqrtHT")!=NNvariablesReg1.end())
    readerReg1->AddVariable("METoverSqrtHT",&METoverSqrtHT);
  if(NNvariablesReg1.find("HT")!=NNvariablesReg1.end())
    readerReg1->AddVariable("HT",&HT);
    if(NNvariablesReg1.find("njets_JESdown")!=NNvariablesReg1.end())
    readerReg1->AddVariable("njets_JESdown",&njets_JESdown); 
    if(NNvariablesReg1.find("njets_JESup")!=NNvariablesReg1.end())
    readerReg1->AddVariable("njets_JESup",&njets_JESup);
    
  
if(NNvariablesReg1.find("Chi2SNT")!=NNvariablesReg1.end())
    readerReg1->AddVariable("Chi2SNT",&Chi2SNT);



   string Decay_Mode = string(argv[6]);

 
   TString BDT_dir = "/afs/cern.ch/work/l/lara/CMSSW_5_3_11/src/testbabyReaderSTOPS/runBDT/";

     readerReg1->BookMVA("bdt_R1",  BDT_dir + "Reg1_"+Decay_Mode+"/OUTPUT/" + setup_directory + "/weights/BDT_BDT.weights.xml");
   reader->BookMVA("bdt_R2",  BDT_dir + "Reg2_"+Decay_Mode+"/OUTPUT/" + setup_directory + "/weights/BDT_BDT.weights.xml");
   reader->BookMVA("bdt_R3",  BDT_dir + "Reg3_"+Decay_Mode+"/OUTPUT/" + setup_directory + "/weights/BDT_BDT.weights.xml");
   reader->BookMVA("bdt_R4",  BDT_dir + "Reg4_"+Decay_Mode+"/OUTPUT/" + setup_directory + "/weights/BDT_BDT.weights.xml");
   reader->BookMVA("bdt_R5",  BDT_dir + "Reg5_"+Decay_Mode+"/OUTPUT/" + setup_directory + "/weights/BDT_BDT.weights.xml");
   reader->BookMVA("bdt_R6",  BDT_dir + "Reg6_"+Decay_Mode+"/OUTPUT/" + setup_directory + "/weights/BDT_BDT.weights.xml");

  // ########################################
  // ##        Run over the events         ##
  // ########################################
  
  cout << theInputTree->GetEntries() <<endl;


  for (int i = 0 ; i < theInputTree->GetEntries() ; i++){

     
	
	
	  ReadEvent(theInputTree,i,&pointers,&myEvent);


      
        bool isSignal = false;

        if (myEvent.mStop == -1) { isSignal = false;}
                else isSignal = true;


         if (myEvent.nJets < 4) continue;
        if (myEvent.numberOfLepton != 1) continue;
	
	if (myEvent.MET < 80) continue;
        if (myEvent.MT < 100) continue;
        if (myEvent.nBTag < 1) continue;
	
	
	
	// Electron triggers
	if ((abs(myEvent.leadingLeptonPDGId) == 11) && (myEvent.triggerElec ==false)) continue;
	// Muon triggers

	if (abs(myEvent.leadingLeptonPDGId) == 13)
	{
		if ((myEvent.leadingLepton.Pt() < 25) && (myEvent.xtriggerMuon ==false)) continue;
		else if ((myEvent.leadingLepton.Pt() > 25) && (myEvent.triggerMuon ==false)) continue;
	}

	if (myEvent.isolatedTrackVeto == false ) continue;
	if (myEvent.tauVeto == false) continue;

        
        if (isSignal && (myEvent.event%2)!=0) continue;


        if ( isSignal && !((myEvent.mStop == atof(argv[6])) && (myEvent.mNeutralino == atof(argv[7]))) )continue;



      
        met = myEvent.MET;
        mT = myEvent.MT;
        mT2W = myEvent.MT2W;
        HTfrac = myEvent.HTRatio;
        b1_pt = myEvent.leadingBPt;
        lepton_pT = myEvent.leadingLeptonPt;
        dPhi_JetMet = myEvent.deltaPhiMETJets;
        dR_LepB = myEvent.deltaRLeptonLeadingB;
        m3b = myEvent.M3b;
        mlb_hemi = myEvent.Mlb_hemi;
        jet1_pT = myEvent.leadingJetPt;
	jet2_pT = myEvent.jets.at(1).Pt();
	jet3_pT = myEvent.jets.at(2).Pt();
	jet1_eta = myEvent.jets.at(0).Eta();
	jet2_eta = myEvent.jets.at(1).Eta();
	jet3_eta = myEvent.jets.at(2).Eta();
	jet1_phi = myEvent.jets.at(0).Phi();
	jet2_phi = myEvent.jets.at(1).Phi();
	jet3_phi = myEvent.jets.at(2).Phi();
        njets = myEvent.nJets;
        METoverSqrtHT = myEvent.METoverSqrtHT;
        HT_MET_lep_pt = myEvent.HTPlusLeptonPtPlusMET;
        Chi2SNT = myEvent.hadronicChi2; 
	njets_JESdown=myEvent.nJets_JESdown;
	njets_JESup= myEvent.nJets_JESup;
	puweight= myEvent.weightPileUp;
	xsweight=myEvent.weightCrossSection;
	nvertex=myEvent.numberOfPrimaryVertices;
	leptonflav=myEvent.leadingLeptonPDGId;
	trweight=triggereff(myEvent.leadingLeptonPt, myEvent.leadingLepton.Eta(), myEvent.jets.at(3).Pt(),leptonflav, myEvent.run );
	lepton_eta=myEvent.leadingLepton.Eta();
	lepton_phi=myEvent.leadingLepton.Phi();

        stop_mass = myEvent.mStop;
        lsp_mass = myEvent.mNeutralino;
	
      
        
 

 
      bdt_R1 = readerReg1->EvaluateMVA( "bdt_R1" );
      bdt_R2 = reader->EvaluateMVA( "bdt_R2" );
      bdt_R3 = reader->EvaluateMVA( "bdt_R3" );
      bdt_R4 = reader->EvaluateMVA( "bdt_R4" );
      bdt_R5 = reader->EvaluateMVA( "bdt_R5" );
      bdt_R6 = reader->EvaluateMVA( "bdt_R6" );
      
      


      double weight = myEvent.weightCrossSection;
     // double lumiweight = (20000. * weight*2);
      double lumiweight = (20000. * weight);

      if (isSignal) lumiweight=lumiweight*2;

      Stop->Fill(myEvent.mStop);
      Neut->Fill(myEvent.mNeutralino);
      
      BDTtree->Fill();



	cout << "bdt1: "<< bdt_R1 << endl; 
	cout << "bdt2: "<< bdt_R2 << endl; 
   
      BDT1->Fill(bdt_R1, lumiweight);
      BDT2->Fill(bdt_R2, lumiweight);
      BDT3->Fill(bdt_R3, lumiweight);
      BDT4->Fill(bdt_R4, lumiweight);
      BDT5->Fill(bdt_R5, lumiweight);
      BDT6->Fill(bdt_R6, lumiweight);
	

  } 
    
      fout->cd();
      Stop->Write(); 
      Neut->Write(); 
      BDT1->Write(); 
      BDT2->Write(); 
      BDT3->Write(); 
      BDT4->Write(); 
      BDT5->Write(); 
      BDT6->Write(); 
      fout->Write();
      fout->Close();

  return (0);

}


float triggereff(float pt, float eta, float ptjet4, float flavour,float runnumber){

float triggeff=1.0;



if (fabs(flavour)==13 && pt <= 24){

float muoneff=1.0;
float jeteff=1.0;

  if ( fabs(eta)<0.8 ) {
    if (pt>=20 && pt<22)  muoneff=  0.94;     
    if (pt>=22 && pt<24)  muoneff=  0.97;      
    if (pt>=24 && pt<26)  muoneff=  0.96;
  } 
  else if ( fabs(eta)>=0.8 && fabs(eta)<1.5 ) {
    if (pt>=20 && pt<22)  muoneff=  0.87;
    if (pt>=22 && pt<24)  muoneff=  0.90;
    if (pt>=24 && pt<26)  muoneff=  0.91;
  }
  else if ( fabs(eta)>=1.5 && fabs(eta)<2.1 ) {
    if (pt>=20 && pt<22)  muoneff=  0.76;
    if (pt>=22 && pt<24)  muoneff=  0.83;
    if (pt>=24 && pt<26)  muoneff=  0.87;
  }
  
 /// Jet Turn ons PFJet30 Runs < 198049 
if (runnumber <= 198049){
  jeteff=(1/(1 + exp(10.95-ptjet4*0.4503)));  }
 
if (runnumber > 198049 && runnumber < 199608){
  jeteff=(1/(1 + exp(1.393-ptjet4*0.2067)));}
  
  if (runnumber >= 199608){ 
  
  jeteff=(0.6412/(0.6414 + exp(9.173-ptjet4*0.3877)));}
  
  
  triggeff=muoneff*jeteff;
  
  
  } //close the if for the x trigger
  
  else if (pt > 24 || fabs(flavour==11)) {
    if (fabs(flavour) == 11) 
  {
    if ( fabs(eta)<1.5) 
    {
      if ( pt>=20 && pt<22 )   triggeff= 0.00;
      if ( pt>=22 && pt<24 )   triggeff= 0.00;
      if ( pt>=24 && pt<26 )   triggeff= 0.00;
      if ( pt>=26 && pt<28 )   triggeff= 0.08;
      if ( pt>=28 && pt<30 )   triggeff= 0.61;
      if ( pt>=30 && pt<32 )   triggeff= 0.86;
      if ( pt>=32 && pt<34 )   triggeff= 0.88;
      if ( pt>=34 && pt<36 )   triggeff= 0.90;
      if ( pt>=36 && pt<38 )   triggeff= 0.91;
      if ( pt>=38 && pt<40 )   triggeff= 0.92;
      if ( pt>=40 && pt<50 )   triggeff= 0.94;
      if ( pt>=50 && pt<60 )   triggeff= 0.95;
      if ( pt>=60 && pt<80 )   triggeff= 0.96;
      if ( pt>=80 && pt<100 )  triggeff= 0.96;
      if ( pt>=100 && pt<150 ) triggeff= 0.96;
      if ( pt>=150 && pt<200 ) triggeff= 0.97;
      if ( pt>=200 )           triggeff= 0.97;
    } 
    else if ( fabs(eta)>=1.5 && fabs(eta)<2.1) 
    {
      if ( pt>=20 && pt<22 )   triggeff= 0.00;
      if ( pt>=22 && pt<24 )   triggeff= 0.00;
      if ( pt>=24 && pt<26 )   triggeff= 0.02;
      if ( pt>=26 && pt<28 )   triggeff= 0.18;
      if ( pt>=28 && pt<30 )   triggeff= 0.50;
      if ( pt>=30 && pt<32 )   triggeff= 0.63;
      if ( pt>=32 && pt<34 )   triggeff= 0.68;
      if ( pt>=34 && pt<36 )   triggeff= 0.70;
      if ( pt>=36 && pt<38 )   triggeff= 0.72;
      if ( pt>=38 && pt<40 )   triggeff= 0.74;
      if ( pt>=40 && pt<50 )   triggeff= 0.76;
      if ( pt>=50 && pt<60 )   triggeff= 0.77;
      if ( pt>=60 && pt<80 )   triggeff= 0.78;
      if ( pt>=80 && pt<100 )  triggeff= 0.80;
      if ( pt>=100 && pt<150 ) triggeff= 0.79;
      if ( pt>=150 && pt<200 ) triggeff= 0.76;
      if ( pt>=200 )           triggeff= 0.81;
    }
  } 
  //muon efficiencies
  else if (fabs(flavour) == 13) 
  {
    if ( fabs(eta)<0.8 ) 
    {
      if (pt>=20 && pt<22)   triggeff= 0.00;     
      if (pt>=22 && pt<24)   triggeff= 0.03;      
      if (pt>=24 && pt<26)   triggeff= 0.87;
      if (pt>=26 && pt<28)   triggeff= 0.90;
      if (pt>=28 && pt<30)   triggeff= 0.91;
      if (pt>=30 && pt<32)   triggeff= 0.91;
      if (pt>=32 && pt<34)   triggeff= 0.92;
      if (pt>=34 && pt<36)   triggeff= 0.93;
      if (pt>=36 && pt<38)   triggeff= 0.93;
      if (pt>=38 && pt<40)   triggeff= 0.93;
      if (pt>=40 && pt<50)   triggeff= 0.94;
      if (pt>=50 && pt<60)   triggeff= 0.95;
      if (pt>=60 && pt<80)   triggeff= 0.95;
      if (pt>=80 && pt<100)  triggeff= 0.94;
      if (pt>=100 && pt<150) triggeff= 0.94;
      if (pt>=150 && pt<200) triggeff= 0.93;
      if (pt>=200)           triggeff= 0.92;
    } 
    else if ( fabs(eta)>=0.8 && fabs(eta)<1.5 ) 
    {
      if (pt>=20 && pt<22)   triggeff= 0.00;
      if (pt>=22 && pt<24)   triggeff= 0.05;
      if (pt>=24 && pt<26)   triggeff= 0.78;
      if (pt>=26 && pt<28)   triggeff= 0.81;
      if (pt>=28 && pt<30)   triggeff= 0.81;
      if (pt>=30 && pt<32)   triggeff= 0.81;
      if (pt>=32 && pt<34)   triggeff= 0.82;
      if (pt>=34 && pt<36)   triggeff= 0.82;
      if (pt>=36 && pt<38)   triggeff= 0.83;
      if (pt>=38 && pt<40)   triggeff= 0.83;
      if (pt>=40 && pt<50)   triggeff= 0.84;
      if (pt>=50 && pt<60)   triggeff= 0.84;
      if (pt>=60 && pt<80)   triggeff= 0.84;
      if (pt>=80 && pt<100)  triggeff= 0.84;
      if (pt>=100 && pt<150) triggeff= 0.84;
      if (pt>=150 && pt<200) triggeff= 0.84;
      if (pt>=200)           triggeff= 0.82;
    } 
    else if ( fabs(eta)>=1.5 && fabs(eta)<2.1 ) 
    {
      if (pt>=20 && pt<22)   triggeff= 0.00;
      if (pt>=22 && pt<24)   triggeff= 0.11;
      if (pt>=24 && pt<26)   triggeff= 0.76;
      if (pt>=26 && pt<28)   triggeff= 0.78;
      if (pt>=28 && pt<30)   triggeff= 0.79;
      if (pt>=30 && pt<32)   triggeff= 0.80;
      if (pt>=32 && pt<34)   triggeff= 0.80;
      if (pt>=34 && pt<36)   triggeff= 0.81;
      if (pt>=36 && pt<38)   triggeff= 0.81;
      if (pt>=38 && pt<40)   triggeff= 0.82;
      if (pt>=40 && pt<50)   triggeff= 0.82;
      if (pt>=50 && pt<60)   triggeff= 0.83;
      if (pt>=60 && pt<80)   triggeff= 0.83;
      if (pt>=80 && pt<100)  triggeff= 0.83;
      if (pt>=100 && pt<150) triggeff= 0.83;
      if (pt>=150 && pt<200) triggeff= 0.82;
      if (pt>=200)           triggeff= 0.82;
    }
  }
  
  
  }
  
  return triggeff; 
  
  }
  
  
 
  
  


