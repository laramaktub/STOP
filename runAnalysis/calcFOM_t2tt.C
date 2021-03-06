#include <iostream>   
#include <algorithm>  
#include <iomanip>
#include <string>
#include <fstream>
#include <map>
#include <sstream>
#include "TColor.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPolyLine.h"
#include "TROOT.h"
#include <stdio.h>
#include <stdlib.h>




using namespace std;

double f_syst = 0.15;

TString signalregiont2tt(int,int);


double calcSoSqrtB(double S, double B){

if (B < 1) B = 1;

      if (S >= 3)
         return S / sqrt(B + f_syst*f_syst * B*B);
      else

    return 0.;
}



int maxbin(TH1D* histo){
  int bi = 0;
  double max = 0;
  for( int i = 1; i <= histo->GetNbinsX(); i++){

    //if( histo->GetBinContent(i) > max && histo->GetBinError(i) < histo->GetBinContent(i) ){
    if( histo->GetBinContent(i) > max ){
      max = histo->GetBinContent(i);
      bi = i;
    }
  }
  return bi;
}



void outputBDTplots(TString setup, TString MVA, TString dir, int MSTOP, int MLSP,bool lineacut,bool logbool){


	  gStyle->SetOptStat(0);
	  gStyle->SetCanvasColor(0);
	  gStyle->SetPadColor(0);
	  gStyle->SetMarkerStyle(15);
	  gStyle->SetMarkerSize(0.25);
	  gStyle->SetTextFont(42);
	  gStyle->SetMarkerColor(37);


          std::ostringstream ostr1;
          ostr1 << MSTOP; 
          std::string stop = ostr1.str();

          std::ostringstream ostr2;
          ostr2 << MLSP;
          std::string neut = ostr2.str();

          TString ntpdir ="ntuples_bestset_t2tt_good"; 
	  TString dataset_name;
          TString datasetnombre;
	
	
	  if (dir == "T2bw025") {dataset_name = "t2bw_025";}
	  if (dir == "T2bw050") {dataset_name = "t2bw_050";}
	  if (dir == "T2bw075") {dataset_name = "t2bw_075";}
	  if (dir == "T2tt") {dataset_name = "t2tt_all";}
	  if (dir == "T2tt") {datasetnombre = "t2tt_half";}
	  
          TFile ttbar(ntpdir+"/"+setup+"/"+dir+"/ttbar/output/ttbar_all_0.root");
          TFile wjets(ntpdir+"/"+setup+"/"+dir+"/wjets/output/wjets_all_0.root");
          TFile others(ntpdir+"/"+setup+"/"+dir+"/others/output/others_all_0.root");
          TFile sig(ntpdir+"/"+setup+"/"+dir+"/signal/"+TString(stop)+"/"+TString(neut)+"/output/"+dataset_name+"_0.root");



	  TH1D* TTBar= (TH1D*)ttbar.Get(MVA);
	  TH1D* WJets= (TH1D*)wjets.Get(MVA);
	  TH1D* Others= (TH1D*)others.Get(MVA);
	  TH1D* signal= (TH1D*)sig.Get(MVA);

	  TTBar->SetFillColor(40);
	  TTBar->SetFillStyle(3001);
	  TTBar->SetLineWidth(0.);
	  signal->SetLineColor(kRed);
	  signal->SetLineWidth(2.);


    	  WJets->SetFillColor(41);
    	  Others->SetFillColor(43);
	  
		  
	  
          THStack *stack= new THStack("stack", "");
          stack->Add(Others);
          stack->Add(WJets);
          stack->Add(TTBar);

         

          int nbins = TTBar->GetNbinsX();
          double hmin = TTBar->GetXaxis()->GetBinLowEdge(1);
          double hmax = TTBar->GetXaxis()->GetBinUpEdge(nbins);


          TH1D* SoB = new TH1D("","",nbins,hmin,hmax);
          int max_bin;
          double cutvalue; 
	  double nsignal = signal->Integral();




          SoB->SetLineColor(kBlue);
          SoB->SetLineWidth(2.);




	      for(int b=1; b<=nbins; ++b){

	      //cout << "Bkg: " << TTBar->Integral(b,nbins+1) << endl;
	      //cout << "Sig: " << signal->Integral(b,nbins+1) << endl;

	      double sig_ = signal->Integral(b,nbins+1);
	      double bkg_ = TTBar->Integral(b,nbins+1) + WJets->Integral(b,nbins+1) + Others->Integral(b,nbins+1);
	      
	      double SoSqrtB = calcSoSqrtB(sig_, bkg_ );

	      SoB->SetBinContent(b,SoSqrtB);

	      }


	/*      double bkgcut=0.;

 if (MLSP==150 && (MSTOP ==300 || MSTOP == 400 ||  MSTOP == 500 || MSTOP == 600 || MSTOP == 700 || MSTOP == 800)){
	     
	     double binR1=TTBar->FindBin(0.25);
	     double binR2=TTBar->FindBin(0.4);
	     double binR3=TTBar->FindBin(0.5);
	     double binR4=TTBar->FindBin(0.5);
	     double binR5=TTBar->FindBin(0.4);
	     double binR6=TTBar->FindBin(0.1);
	     
	     TH1F *histobackground=(TH1F*)TTBar->Clone();
	     histobackground->Add(WJets);
	     histobackground->Add(Others);
	     
	     //histobackground->SetAxisRange(0.3,1,"X");
	    	      
	      if (MSTOP==300)  cout << TTBar->Integral(binR1,nbins+1)+Others->Integral(binR1,nbins+1)+WJets->Integral(binR1,nbins+1) << " +/- "  << sqrt(histobackground->Integral(binR1,nbins+1)) <<  endl; 
	      if (MSTOP==400)  cout << histobackground->Integral(binR2,nbins+1) << " +/- "  << sqrt(histobackground->Integral(binR2,nbins+1)) <<  endl;  
	      if (MSTOP==500)  cout << histobackground->Integral(binR3,nbins+1) << " +/- "  << sqrt(histobackground->Integral(binR3,nbins+1)) <<  endl;  
	      if (MSTOP==600)  cout << histobackground->Integral(binR4,nbins+1) << " +/- "  << sqrt(histobackground->Integral(binR4,nbins+1)) <<  endl;  
	      if (MSTOP==700)  cout << histobackground->Integral(binR5,nbins+1) << " +/- "  << sqrt(histobackground->Integral(binR5,nbins+1)) <<  endl;  
	      if (MSTOP==800)  cout << histobackground->Integral(binR6,nbins+1) << " +/- "  << sqrt(histobackground->Integral(binR6,nbins+1)) <<  endl;  

	       }*/
 


	   max_bin = maxbin(SoB);
	   cutvalue = TTBar->GetBinLowEdge(max_bin);
           

	   double sob = SoB->GetBinContent(max_bin);

       	   char cutval_[32];
	   snprintf(cutval_, 32, "%.1g", cutvalue);
       	   char fom_[32];
	   snprintf(fom_, 32, "%.2g", sob);
       	   char signal_[32];
	   snprintf(signal_, 32, "%.2g", nsignal);

           TLegendEntry *legge;
           TLegend *leg;
           leg = new TLegend(0.6,0.55,0.9,0.85);
           leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.043);
           legge = leg->AddEntry(TTBar,   "t#bar{t} + Jets", "F");
           legge = leg->AddEntry(WJets,   "W + Jets", "F");
           legge = leg->AddEntry(Others,   "Others", "F");
           legge = leg->AddEntry(signal, "sig("+TString(stop)+","+TString(neut)+")", "l");
           leg->SetFillColor(0);

 

           TCanvas c2("c2","c2",800,600);
	   
	 	   
	   stack->Draw("");
	   signal->Draw("same");
	  
	  stack->GetYaxis()->SetTitle("Entries");
	  stack->GetXaxis()->SetTitle("BDT output");
	  
	  signal->GetXaxis()->SetRangeUser(-0.4,0.8);
	  stack->GetXaxis()->SetRangeUser(-0.4,0.8);
	  
	  signal->SetMinimum(0.3);
	  stack->SetMinimum(0.3);
	  
	  //signal->GetYaxis()->SetNdivisions(10);
	  //stack->GetYaxis()->SetNdivisions(10);
	  
	  
	   if(logbool) c2.SetLogy();
           leg->Draw();
	   
	  
	   
	   double maximobin=stack->GetMaximum();
	   double cutvalueplot = atof(TString(cutval_));
	   TLine *linea = new TLine(cutvalueplot,0,cutvalueplot,maximobin);
	   linea->SetLineWidth(2);
	   linea->SetLineColor(4);
	   if (lineacut) linea->Draw();
           TLatex l1;
           l1.SetTextAlign(12);
           l1.SetTextSize(0.04);
           l1.SetNDC();
           l1.DrawLatex(0.155, 0.98, "CMS Simulation, 20 fb^{-1}");
           l1.DrawLatex(0.7, 0.98, "#sqrt{s} = 8 TeV");

           TLatex l2;
           l2.SetTextAlign(12);
           l2.SetTextSize(0.04);
           l2.SetNDC();
           l2.DrawLatex(0.22, 0.3, "#color[2]{Optimal cut: \t "+TString(cutval_)+"}");
           l2.DrawLatex(0.22, 0.25, "#color[2]{FOM@cut: \t "+TString(fom_)+"}");
           l2.DrawLatex(0.22, 0.2, "#color[2]{Tot. N_{Signal}: \t "+TString(signal_)+"}");
	   
          if (!logbool)  c2.Print("~/www/STOP/BDTTraining/8TeV/"+datasetnombre+"/BestSet/"+TString(MVA)+"/histo_"+TString(stop)+"_"+TString(neut)+"_lineal.png");
          if(logbool)   c2.Print("~/www/STOP/BDTTraining/8TeV/"+datasetnombre+"/BestSet/"+TString(MVA)+"/histo_"+TString(stop)+"_"+TString(neut)+"_log.png");
}


void makeAlloutputBDTplots(TString setup, TString BDT, TString dir, bool lineacut,bool logbool){

              for(int x=175; x<=800; x+=25){
                           
                       for(int y=25; y<=700; y+=25){
				if (x - y > 99){
				//if (x==725 && y==25){
				outputBDTplots(setup,BDT,dir,x,y,lineacut,logbool);
				}
			}
		}
}



double FOM(TString setup, TString BDT, TString dir, int MSTOP, int MLSP, int num){


	  gStyle->SetOptStat(0);
	  gStyle->SetCanvasColor(0);
	  gStyle->SetPadColor(0);
	  gStyle->SetMarkerStyle(15);
	  gStyle->SetMarkerSize(0.25);
	  gStyle->SetTextFont(42);
	  gStyle->SetMarkerColor(37);

	  
          std::ostringstream ostr1;
          ostr1 << MSTOP; 
          std::string stop = ostr1.str();

          std::ostringstream ostr2;
          ostr2 << MLSP;
          std::string neut = ostr2.str();


           TString dataset_name;
	     TString datasetnombre;
          if (dir == "T2bw025") {dataset_name = "t2bw_025";}
	  if (dir == "T2bw050") {dataset_name = "t2bw_050";}
	  if (dir == "T2bw075") {dataset_name = "t2bw_075";}
          if (dir == "T2tt") {dataset_name = "t2tt_all";}
	  if (dir == "T2tt") {datasetnombre = "t2tt_half";}

	  TString ntpdir = "ntuples_bestset_t2tt_good";

	  TFile ttbar(ntpdir+"/"+setup+"/"+dir+"/ttbar/output/ttbar_all_0.root");
	  TFile wjets(ntpdir+"/"+setup+"/"+dir+"/wjets/output/wjets_all_0.root");
	  TFile others(ntpdir+"/"+setup+"/"+dir+"/others/output/others_all_0.root");
	  TFile sig(ntpdir+"/"+setup+"/"+dir+"/"+"signal/"+TString(stop)+"/"+TString(neut)+"/output/"+dataset_name+"_0.root"); 
	  
	  
	//  TFile sig(ntpdir+"/"+setup+"/"+dir+"/"+TString(stop)+"/"+TString(neut)+"/output/"+dir+"_0.root"); 	
//	  TFile sig("../Leg3_NN_Strass/"+setup+"/"+TString(stop)+"/"+TString(neut)+"/output/SMS_t2tt_0.root"); 	



	  TH1D* TTBar= (TH1D*)ttbar.Get(BDT);
	  TH1D* WJets= (TH1D*)wjets.Get(BDT);
	  TH1D* Others= (TH1D*)others.Get(BDT);
	  TH1D* signal= (TH1D*)sig.Get(BDT);
	  
	

          int nbins = TTBar->GetNbinsX();
          double hmin = TTBar->GetXaxis()->GetBinLowEdge(1);
          double hmax = TTBar->GetXaxis()->GetBinUpEdge(nbins);
	  
	  



          TH1D* FOM = new TH1D("","",nbins,hmin,hmax);
          int max_bin;
          double cutvalue; 
	  double nsignal = signal->Integral();

 for(int b=1; b<=nbins; ++b){

	      double sig_ = signal->Integral(b,nbins+1);
	      double bkg_ = TTBar->Integral(b,nbins+1) + WJets->Integral(b,nbins+1) + Others->Integral(b,nbins+1) ;
	      
	      double SoSqrtB = calcSoSqrtB(sig_, bkg_ );

	      FOM->SetBinContent(b,SoSqrtB);

	      }


	    //max_bin = maxbin(FOM);
	    //cutvalue = TTBar->GetBinLowEdge(max_bin);
	    cutvalue=0.35;
	    if (BDT=="BDT1") cutvalue=0.35;
	    cout << "cut value " << cutvalue << endl;
	    
	    max_bin=TTBar->FindBin(cutvalue);
	    cout << "max bin " << max_bin << endl;
	    
            double fom = FOM->GetBinContent(max_bin);
	    
	    


   if (num==0) return nsignal;
   if (num==1) return fom;
   if (num==2) return cutvalue;
}




void make2Dplot(TString setup, TString BDT,TString dir, int num){

          gStyle->SetOptStat(0);
          gStyle->SetCanvasColor(0);
          gStyle->SetPadColor(0);
          gStyle->SetMarkerStyle(15);
          gStyle->SetMarkerSize(0.25);
          gStyle->SetTextFont(42);
          gStyle->SetMarkerColor(37);
          if (num <3) {gStyle->SetPaintTextFormat("4.1f");}

                           
              TH2D* TwoDPlot = new TH2D("","",26,157.5, 812.5, 28, 12.5,712.5);

 	 
              for(int x=175; x<=800; x+=25){
	
    	               for(int y=25; y<=700; y+=25){
			       if (x - y > 99){ 
				TwoDPlot->Fill(x,y,FOM(setup,BDT,dir,x,y,num));

			       }	
              }
	   }

           TCanvas c1("c1","c1",800,600);
	   c1.SetLeftMargin(0.1706731);
	   c1.SetRightMargin(0.1983173);
	   c1.SetTopMargin(0.04895105);
	   c1.SetBottomMargin(0.1416084);
	   c1.Range(-289.7381,-191.8196,1334.643,1074.487);
           TwoDPlot->SetMarkerSize(1.);
           TwoDPlot->SetMarkerColor(kWhite);
           TwoDPlot->Draw("COLZ TEXT");
           TwoDPlot->GetYaxis()->SetTitle("LSP mass");
           TwoDPlot->GetXaxis()->SetTitle("Stop mass");
           TwoDPlot->GetZaxis()->SetTitle("Entries");
           TwoDPlot->GetZaxis()->SetTitleOffset(1.1);
           if (num==0) {TwoDPlot->GetZaxis()->SetTitle("Entries"); TwoDPlot->GetZaxis()->SetRangeUser(0,5000);}
           if (num==1) {TwoDPlot->GetZaxis()->SetTitle("FOM"); TwoDPlot->GetZaxis()->SetRangeUser(0,15);}
           if (num==2) {TwoDPlot->GetZaxis()->SetTitle("Cutval"); TwoDPlot->GetZaxis()->SetRangeUser(0,1);}
           TLatex l1;
           l1.SetTextAlign(12);
           l1.SetTextSize(0.04);
           l1.SetNDC();
           l1.DrawLatex(0.155, 0.98, "CMS Simulation, 20 fb^{-1}");
           l1.DrawLatex(0.7, 0.98, "#sqrt{s} = 8 TeV");
	   
	   	    
TString dataset_name;
TString datasetnombre;
          if (dir == "T2bw025") {dataset_name = "t2bw_025";}
          if (dir == "T2bw050") {dataset_name = "t2bw_050";}
          if (dir == "T2bw075") {dataset_name = "t2bw_075";}
	  if (dir == "T2tt") {dataset_name = "t2tt_all";}
	  if (dir == "T2tt") {datasetnombre = "t2tt_half";}
	
          //c1.Print("~/www/STOP/BDTTraining/08TeV/"+dataset_name+"/"+setup+"/"+BDT+"/"+"FOM_Lara.png");
c1.Print("~/www/STOP/BDTTraining/8TeV/"+datasetnombre+"/BestSet/"+BDT+"/"+"FOM_Lara.png");
}



void make2Dplot_MAX(TString setup, TString dir, int num, int MLSP,TString opt){

          gStyle->SetOptStat(0);
          gStyle->SetCanvasColor(0);
          gStyle->SetPadColor(0);
          gStyle->SetMarkerStyle(15);
          gStyle->SetMarkerSize(0.25);
          gStyle->SetTextFont(42);
          gStyle->SetMarkerColor(37);
	  if (num <3) {gStyle->SetPaintTextFormat("4.1f");}
	  
	  TString BDT ="";
	  
	  
	char cadena[128];
   // Crea un fichero de salida
   ofstream fs(setup+".txt"); 

   fs.setf(ios::fixed);
   fs.precision(1);
   

              TH2D* TwoDPlot = new TH2D("","",26,157.5, 812.5, 28, 12.5,712.5);

 	 
              for(int x=175; x<=800; x+=25){
           //   for(int x=175; x<=225; x+=25){

	
    	               for(int y=25; y<=700; y+=25){
		       
			       if (x - y > 99){ 
			       
			       
			               BDT = signalregiont2tt(x,y);
				       
				       cout << BDT << endl;
 
					
					double array[1] = {
							   FOM(setup,BDT,dir,x,y,num) 
							  
							   };
							   
				        double arraycut[1]= {
							   FOM(setup,BDT,dir,x,y,1) 
							  
							   };
					                                  
					


					double temp = 0.;
					int mvaval = 0;
					double tempcut=0.;

						// Get the maximum of each point for all MVAs
						  for(int i=0;i<1;i++){

							if(arraycut[i]>temp){ 

								temp	= arraycut[i]; 
								mvaval 	= i; 
								tempcut=array[i];
								
							}  
						  }
						  
						
						  
			
	
								if (num ==1) {TwoDPlot->Fill(x,y,temp);}
								if (num ==2) {TwoDPlot->Fill(x,y,tempcut);}
								if (num ==3) {
								if (BDT=="BDT1") mvaval=1;
								if (BDT=="BDT3") mvaval=3;
								if (BDT=="BDT4") mvaval=4;
								TwoDPlot->Fill(x,y,mvaval);}
								//	createTableCLs(x,y,temp);
								
			       }	
              	        }
			
			
			
	   }
	   
	                                                     
	   
	   fs.close();


           TCanvas c1("c1","c1",800,600);
	   c1.SetLeftMargin(0.1706731);
	   c1.SetRightMargin(0.1983173);
	   c1.SetTopMargin(0.04895105);
	   c1.SetBottomMargin(0.1416084);
	   c1.Range(-289.7381,-191.8196,1334.643,1074.487);
           TwoDPlot->SetMarkerSize(1.); 
           TwoDPlot->SetMarkerColor(kWhite); 
           TwoDPlot->Draw("COLZ TEXT");

           TwoDPlot->GetYaxis()->SetTitle("LSP mass"); 
           TwoDPlot->GetXaxis()->SetTitle("Stop mass");  
           if (num==0) {TwoDPlot->GetZaxis()->SetTitle("Entries"); TwoDPlot->GetZaxis()->SetRangeUser(0,5000);}
           if (num==1) {TwoDPlot->GetZaxis()->SetTitle("FOM"); TwoDPlot->GetZaxis()->SetRangeUser(0,15);}
           if (num==2) {TwoDPlot->GetZaxis()->SetTitle("Optimal cut point"); TwoDPlot->GetZaxis()->SetRangeUser(0,1);}
           if (num==3) {TwoDPlot->GetZaxis()->SetTitle("Best performing BDT training"); TwoDPlot->GetZaxis()->SetRangeUser(0,6);}
	   
	                                             
                                                                
	                                                        TFile* f = new TFile(dir+"_"+opt+".root","RECREATE");
	                                                        TwoDPlot->Write("twodplot");
	                                                        f->Write();
	                                                        f->Close();

           TLatex l1;
           l1.SetTextAlign(12);
           l1.SetTextSize(0.04);
           l1.SetNDC();
           l1.DrawLatex(0.155, 0.98, "CMS Simulation, 20 fb^{-1}");
           l1.DrawLatex(0.7, 0.98, "#sqrt{s} = 8 TeV");
           l1.SetTextSize(0.03);
	   
	   
	   
	
	 
	 	
	   
	   
	   if (num==3){
	   
	   
       tex = new TLatex(0.7449749,0.8251748,"BDT1");
tex->SetNDC();
   tex->SetTextAlign(12);
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetTextAngle(41.15586);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.7449749,0.7097902,"BDT2");
tex->SetNDC();
   tex->SetTextAlign(12);
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetTextAngle(41.15586);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.7437186,0.5891608,"BDT3");
tex->SetNDC();
   tex->SetTextAlign(12);
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetTextAngle(41.15586);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.7424623,0.4755245,"BDT4");
tex->SetNDC();
   tex->SetTextAlign(12);
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetTextAngle(41.15586);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.7424623,0.3583916,"BDT5");
tex->SetNDC();
   tex->SetTextAlign(12);
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetTextAngle(41.15586);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.7424623,0.2447552,"BDT6");
tex->SetNDC();
   tex->SetTextAlign(12);
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetTextAngle(41.15586);
   tex->SetLineWidth(2); 
     tex->Draw();
   
   }
	    
	    
TString dataset_name;
TString datasetnombre;
          if (dir == "T2bw025") {dataset_name = "t2bw_025";}
          if (dir == "T2bw050") {dataset_name = "t2bw_050";}
          if (dir == "T2bw075") {dataset_name = "t2bw_075";}
	if (dir == "T2tt") {dataset_name = "t2tt_all";}
	  if (dir == "T2tt") {datasetnombre = "t2tt_half";} 
	  
//if (num==0) c1.Print("~/www/STOP/BDTTraining/08TeV/"+dataset_name+"/"+setup+"/Entries_Lara.png");
//if (num==1) c1.Print("~/www/STOP/BDTTraining/08TeV/"+dataset_name+"/"+setup+"/FOM_Lara.png");
//if (num==2) c1.Print("~/www/STOP/BDTTraining/08TeV/"+dataset_name+"/"+setup+"/OptimalCut_Lara.png");
//if (num==3) c1.Print("~/www/STOP/BDTTraining/08TeV/"+dataset_name+"/"+setup+"/BestBDT_Lara.png");

if (num==0) c1.Print("~/www/STOP/BDTTraining/8TeV/"+datasetnombre+"/BestSet/Entries_Lara_SR_fixcut.png");
if (num==1) c1.Print("~/www/STOP/BDTTraining/8TeV/"+datasetnombre+"/BestSet/FOM_Lara_SR_fixcut.png");
if (num==2) c1.Print("~/www/STOP/BDTTraining/8TeV/"+datasetnombre+"/BestSet/OptimalCut_Lara_SR_fixcut.png");
if (num==3) c1.Print("~/www/STOP/BDTTraining/8TeV/"+datasetnombre+"/BestSet/BestBDT_Lara_SR_fixcut.png");


}





void make2d(){


make2Dplot("setup_2", "MVA1", 1);
make2Dplot("setup_2", "MVA2", 1);
make2Dplot("setup_2", "MVA3", 1);
make2Dplot("setup_2", "MVA4", 1);
make2Dplot("setup_2", "MVA5", 1);
make2Dplot("setup_2", "MVA6", 1);


/*make2Dplot("setup_2", "MVA1", 2);
make2Dplot("setup_2", "MVA2", 2);
make2Dplot("setup_2", "MVA3", 2);
make2Dplot("setup_2", "MVA4", 2);
make2Dplot("setup_2", "MVA5", 2);
make2Dplot("setup_2", "MVA6", 2);
*/
}





void createTableCLs(int S, int N, double signal, double ttbar){
 
  char datacardname[100];
  sprintf(datacardname,"CLsS%dN%d.txt",S,N);
 
  ofstream  tablesFile(datacardname);
  tablesFile.setf(ios::fixed);
  tablesFile.precision(1);
 
  tablesFile << "imax 1  number of channels" << endl
             << "jmax 2  number of backgrounds" << endl
             << "kmax 4  number of nuisance parameters (sources of systematical uncertainties)" << endl
             << "------------"<<endl
             << "bin 1"<<endl    
  //           << "observation \t" << fillTable(map_yields_mu,"Data",cutname)+fillTable(map_yields_el,"Data",cutname) << endl
             << "------------" << endl
             << "bin            \t 1              \t 1              \t 1" << endl
             << "process        \t signal         \t TTJets         \t others" << endl
             << "process        \t 0              \t 1              \t 2" << endl
             << "rate           \t " << signal << "  \t \t \t "
  //           << fillTable(map_yields_mu,"TTbar",cutname) + fillTable(map_yields_el,"TTbar",cutname) <<  " \t \t \t "
  //           << fillTable(map_yields_mu,"others",cutname) + fillTable(map_yields_el,"others",cutname)<< endl
             << "------------" << endl
             << "lumi    \t lnN \t 1.022        \t 1.022          \t 1.022         \t  lumi affects both signal and mc. lnN = lognormal" << endl
             << "eff_unc \t lnN \t 1.3          \t -              \t -             \t  stop cross section + signal efficiency + other minor ones." << endl
             << "stat_unc \t lnN \t -             \t 1.045          \t 1.8           \t  4.5\% and 80\% stat uncertainty on TTJets and others" << endl
             << "bg_others \t lnN \t -            \t -              \t 1.30          \t  30\% uncertainty on the rest of the backgrounds" << endl;
 
  tablesFile.close();
 
}

TString signalregiont2tt (int stopmass, int lspmass)
{

TString BDT="";

  int signalregion = 0;

  if (stopmass >= 175 && stopmass <= 325)
    signalregion = 1;

  else if (stopmass >= 350 && stopmass <= 700)
    {
      if (lspmass > stopmass - 225)
	signalregion = 1;
      else
	signalregion = 3;

    }

  else if (stopmass >= 725)
    {
      if (lspmass >= 500)
	signalregion = 1;
      else
	signalregion = 4;
    }

  if (signalregion==1) BDT="BDT1";
  if (signalregion==3) BDT="BDT3";
  if (signalregion==4) BDT="BDT4";
  
  
  return BDT;

}
