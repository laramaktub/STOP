#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <set>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPRegexp.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

//#include "/afs/cern.ch/work/s/sigamani/public/CMSSW_5_3_7_STOP/src/xsections.h"

using namespace TMVA;
using namespace std;

string strip(string line)
{
  int l = line.size();
  if ( l == 0 ) return std::string("");
  int n = 0;
  while (((line[n] == 0)    ||
	  (line[n] == ' ' ) ||
	  (line[n] == '\n') ||
	  (line[n] == '\t')) && n < l) n++;

  int m = l-1;
  while (((line[m] == 0)    ||
	  (line[m] == ' ')  ||
	  (line[m] == '\n') ||
	  (line[m] == '\t')) && m > 0) m--;
  return line.substr(n,m-n+1);
}

vector<string> getFilenames(string filelist)
{
  std::ifstream stream(filelist.c_str());
  if ( !stream.good() )
    cout << "unable to open file: " << filelist << endl;

  // Get list of ntuple files to be processed

  std::vector<std::string> v;
  std::string filename;
  while ( stream >> filename )
	if ( strip(filename) != "" ) v.push_back(filename);
  return v;
}


void NN(TString outputfilename, string sigfilename, string bkgfilelistname, TString NN_vars = "met,mT,lepton_pT", 
	TString NN_book_options = "!H:VarTransform=Norm:NeuronType=tanh:NCycles=600:HiddenLayers=N+5,N+3:NormMode=EqualNumEvents", 
	TString train_test = "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:SplitSeed=100"){
  TFile* outputFile = new TFile(outputfilename,"RECREATE");
  
  // create TMVA factory
  Factory* factory = new Factory("BDT",outputFile);

  // open files with tree for signal and background
  vector<string> sigfiles = getFilenames(sigfilename);
  vector<TFile*> sigSrc;
  for(unsigned int i=0; i<sigfiles.size(); i++){
    TFile* sigSrc_i = TFile::Open(sigfiles[i].c_str());
    sigSrc.push_back(sigSrc_i);
  }
  vector<string> bkgfiles = getFilenames(bkgfilelistname);
  vector<TFile*> bkgSrc;
  for(unsigned int i=0; i<bkgfiles.size(); i++){
    TFile* bkgSrc_i = TFile::Open(bkgfiles[i].c_str());
    bkgSrc.push_back(bkgSrc_i);
  }
  cout << "including " << bkgSrc.size() << " background sources" << endl;

  // Get the signal and background trees from TFile source(s);
  // multiple trees can be registered with the Factory
  vector<TTree*> sigTree;
  for(unsigned int i=0; i<sigSrc.size(); i++){
    TTree* sigTree_i = (TTree*)sigSrc[i]->Get( "BDTtree" );
    sigTree.push_back(sigTree_i);
  }
  vector<TTree*> bkgTree;
  for(unsigned int i=0; i<bkgSrc.size(); i++){
    TTree* bkgTree_i = (TTree*)bkgSrc[i]->Get( "BDTtree" );
    bkgTree.push_back(bkgTree_i);
  }

  // Set the event weights per tree (these weights are applied in
  // addition to individual event weights that can be specified)
  double weight=1.;

  Double_t sigWeight = weight;

  vector<Double_t> bkgWeight;
  for(unsigned int i=0;i<bkgTree.size(); i++){
    //first get name of bkg
    string name = bkgfiles[i];
    //give weight according to the name

    if(name.find("TTJetsSemiLep") != string::npos)
      weight = 53.2*1.7e-02; // cross-section * QC+4J selection efficiency
    else if(name.find("TTJetsFullLep") != string::npos)
      weight = 13.43*8.4e-03;
    else if(name.find("TTJetsHad") != string::npos)
      weight = 53.4*2.9e-05;
/*    else if(name.find("topbar_") != string::npos)
      weight = 1.44*393/137980;
    else if(name.find("topbarT_") != string::npos)
      weight = 22.65*3820/1944822;
    else if(name.find("topTW_") != string::npos)
      weight = 2*7.87*20023/814390;
    else if(name.find("zjets_") != string::npos)
      weight = 3048.*8549/36059026;
    else if(name.find("ww_") != string::npos)
      weight = 43.*7789/4225857;
    else if(name.find("wz_") != string::npos)
      weight = 18.2*6248/4265171;
    else if(name.find("zz_") != string::npos)
      weight = 5.9*2628/4190973;
    else if(name.find("qcdmu_") != string::npos)
      weight = 2.966e8*0.0002855*80/25066511;
    else if(name.find("qcd30_") != string::npos)
      weight = 5.944e7*0.061*0;
    else if(name.find("qcd80_") != string::npos)
      weight = 898200.*0.159*10/2654051;
    else if(name.find("qcd170_") != string::npos)
      weight = 22140.*0.1474*89/2720054;
    else if(name.find("qcd250_") != string::npos)
      weight = 2900.*0.1269*182/2761474;
    else if(name.find("qcd350_") != string::npos)
      weight = 520*0.1058*350/2769821;
*/
    bkgWeight.push_back(weight);
    cout << "weight of bg " << i << " " << name << " = " << weight << endl;
    }


  outputFile->cd();

  // Register the trees
  for(unsigned int i=0;i<sigTree.size(); i++){
   // factory->AddSignalTree(sigTree[i],sigWeight);
   factory->AddTree(sigTree[i], "Signal", sigWeight, "(event%2)==1", "train");
   factory->AddTree(sigTree[i], "Signal", sigWeight, "(event%2)==0", "test");
  }
  for(unsigned int i=0;i<bkgTree.size(); i++){
   factory->AddBackgroundTree(bkgTree[i],bkgWeight[i]);
  // factory->AddTree(bkgTree[i], "Background", bkgWeight[i], "(event%2)==1", "train");
   //factory->AddTree(bkgTree[i], "Background", bkgWeight[i], "(event%2)==0", "test");
  }




  // decode NN_vars
  set<TString> variables;
  int nvariables = 0;
  TPMERegexp _variables(",");
  nvariables = _variables.Split(NN_vars);
  for(int i=0; i<nvariables; i++){
    variables.insert(_variables[i]);
  }

  if(variables.find("event")!=variables.end())
    factory->AddVariable("event",'I');
  if(variables.find("met")!=variables.end())
    factory->AddVariable("met",'F');
  if(variables.find("met_phi")!=variables.end())
    factory->AddVariable("met_phi",'F');
  if(variables.find("met_signif")!=variables.end())
    factory->AddVariable("met_signif",'F');
  if(variables.find("mT")!=variables.end())
    factory->AddVariable("mT",'F');
  if(variables.find("lepton_pT")!=variables.end())
    factory->AddVariable("lepton_pT",'F');
  if(variables.find("lepton_phi")!=variables.end())
    factory->AddVariable("lepton_phi",'F');
  if(variables.find("njets")!=variables.end())
    factory->AddVariable("njets",'I');
  if(variables.find("nbjets")!=variables.end())
    factory->AddVariable("nbjets",'I');
  if(variables.find("m3")!=variables.end())
    factory->AddVariable("m3",'F');
  if(variables.find("m3b")!=variables.end())
    factory->AddVariable("m3b",'F');
  if(variables.find("mjetlepmax")!=variables.end())
    factory->AddVariable("mjetlepmax",'F');
  if(variables.find("mi34b_m2")!=variables.end())
    factory->AddVariable("mi34b_m2",'F');
  if(variables.find("mi34")!=variables.end())
    factory->AddVariable("mi34",'F');
  if(variables.find("m3_pt")!=variables.end())
    factory->AddVariable("m3_pt",'F');
  if(variables.find("mi34b")!=variables.end())
    factory->AddVariable("mi34b",'F');
  if(variables.find("phi_lepton_met")!=variables.end())
    factory->AddVariable("phi_lepton_met",'F');
  if(variables.find("m2")!=variables.end())
    factory->AddVariable("m2",'F');
  if(variables.find("mi34m")!=variables.end())
    factory->AddVariable("mi34m",'F');
  if(variables.find("jet1_pT")!=variables.end())
    factory->AddVariable("jet1_pT",'F');
  if(variables.find("jet2_pT")!=variables.end())
    factory->AddVariable("jet2_pT",'F');
  if(variables.find("jet3_pT")!=variables.end())
    factory->AddVariable("jet3_pT",'F');
  if(variables.find("jet4_pT")!=variables.end())
    factory->AddVariable("jet4_pT",'F');
  if(variables.find("jet1_phi")!=variables.end())
    factory->AddVariable("jet1_phi",'F');
  if(variables.find("jet2_phi")!=variables.end())
    factory->AddVariable("jet2_phi",'F');
  if(variables.find("jet3_phi")!=variables.end())
    factory->AddVariable("jet3_phi",'F');
  if(variables.find("jet4_phi")!=variables.end())
    factory->AddVariable("jet4_phi",'F');
  if(variables.find("missmass")!=variables.end())
    factory->AddVariable("missmass",'F');
  if(variables.find("mlpjpa")!=variables.end())
    factory->AddVariable("mlpjpa",'F');
  if(variables.find("dphiljpa")!=variables.end())
    factory->AddVariable("dphiljpa",'F');
  if(variables.find("dphiljb")!=variables.end())
    factory->AddVariable("dphiljb",'F');
  if(variables.find("m3top")!=variables.end())
    factory->AddVariable("m3top",'F');
  if(variables.find("m3topW")!=variables.end())
    factory->AddVariable("m3topW",'F');
  if(variables.find("b1_pt")!=variables.end())
    factory->AddVariable("b1_pt",'F');
  if(variables.find("b1_pt_2")!=variables.end())
    factory->AddVariable("b1_pt_2",'F');
  if(variables.find("b2_pt")!=variables.end())
    factory->AddVariable("b2_pt",'F');
  if(variables.find("b2_pt_2")!=variables.end())
    factory->AddVariable("b2_pt_2",'F');
  if(variables.find("b1plusb2_pt")!=variables.end())
    factory->AddVariable("b1plusb2_pt",'F');
  if(variables.find("b1plusb2_pt_2")!=variables.end())
    factory->AddVariable("b1plusb2_pt_2",'F');
  if(variables.find("b1_eta")!=variables.end())
    factory->AddVariable("b1_eta",'F');
  if(variables.find("b1_eta_2")!=variables.end())
    factory->AddVariable("b1_eta_2",'F');
  if(variables.find("dR_b1b2")!=variables.end())
    factory->AddVariable("dR_b1b2",'F');
  if(variables.find("dR_b1b2_2")!=variables.end())
    factory->AddVariable("dR_b1b2_2",'F');
  if(variables.find("dPhi_b1b2")!=variables.end())
    factory->AddVariable("dPhi_b1b2",'F');
  if(variables.find("dPhi_b1b2_2")!=variables.end())
    factory->AddVariable("dPhi_b1b2_2",'F');
  if(variables.find("mlb1")!=variables.end())
    factory->AddVariable("mlb1",'F');
  if(variables.find("mlb1_2")!=variables.end())
    factory->AddVariable("mlb1_2",'F');
  if(variables.find("HT")!=variables.end())
    factory->AddVariable("HT",'F');
  if(variables.find("HT_av")!=variables.end())
    factory->AddVariable("HT_av",'F');
  if(variables.find("HT_MET_lep_pt")!=variables.end())
    factory->AddVariable("HT_MET_lep_pt",'F');
  if(variables.find("mlb_hemi")!=variables.end())
    factory->AddVariable("mlb_hemi",'F');
  if(variables.find("HTfrac")!=variables.end())
    factory->AddVariable("HTfrac",'F');
  if(variables.find("HTfrac_FNAL")!=variables.end())
    factory->AddVariable("HTfrac_FNAL",'F');
  if(variables.find("mT2W")!=variables.end())
    factory->AddVariable("mT2W",'F');
  if(variables.find("dR_LepB")!=variables.end())
    factory->AddVariable("dR_LepB",'F');
  if(variables.find("dR_LepJet")!=variables.end())
    factory->AddVariable("dR_LepJet",'F');
  if(variables.find("dPhi_JetMet")!=variables.end())
    factory->AddVariable("dPhi_JetMet",'F');
  if(variables.find("METoverSqrtHT")!=variables.end())
    factory->AddVariable("METoverSqrtHT",'F');
  if(variables.find("Chi2SNT")!=variables.end())
    factory->AddVariable("Chi2SNT",'F');


/*
  // Register the trees
  for(unsigned int i=0;i<sigTree.size(); i++){
   // factory->AddSignalTree(sigTree[i],sigWeight);
   factory->AddTree(sigTree[i], "Signal", sigWeight, "(event%2)==1", "train");
   factory->AddTree(sigTree[i], "Signal", sigWeight, "(event%2)==0", "test");
  }
  for(unsigned int i=0;i<bkgTree.size(); i++){
  //  factory->AddBackgroundTree(bkgTree[i],bkgWeight[i]);
   factory->AddTree(bkgTree[i], "Background", bkgWeight[i], "(event%2)==1", "train");
   factory->AddTree(bkgTree[i], "Background", bkgWeight[i], "(event%2)==0", "test");
  }
*/
  cout << "Added variables" << endl;

  factory->PrepareTrainingAndTestTree("",train_test);
  cout << "Prepared training and testing" << endl;

  factory->BookMethod(TMVA::Types::kBDT,"BDT",NN_book_options);
  cout << "booked method" << endl;

  factory->TrainAllMethods();
  cout << "trained method" << endl;

  factory->TestAllMethods();
  cout << "tested method" << endl;

  factory->EvaluateAllMethods();
  cout << "evaluated method" << endl;

}
