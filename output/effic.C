include "recon.h"
#include <string>
#include <ctime>
#include <iostream>
#include <fstream>



using namespace std;


void effic(char config = 'd'){
  
  TString title,answer;
  cout<<"Output ROOT Tree name: ";
  cin>>title;
  cout<<endl;

  TString exist = gSystem->FindFile("./",title+".root");
  if(exist!=""){
    do{
      cout<<"Root file exists. Substitute? (y/n): ";
      cin>>answer;
    }while(answer!="y"&answer!="n");
    if(answer=="n") return;
  }
  
  TFile *ofile = new TFile(title+".root","RECREATE","configurations");
  ofile->Close();

  const Int_t effic_config = 25;
  
  //Will contain three relative efficiency value per configuration: (To be expanded)
  // #1 -> for Confx/ConfA
  // #2 -> for Confxonly/ConfA
  // #3 -> for Confxonly/Confx

  
  
  
  Float_t **efficiencies = 0;
  efficiencies = new Float_t*[3];
  for (Int_t i=0;i<3;i++){
    efficiencies[i] = new Float_t[25];
    for(Int_t j=0; j<25; j++){
      efficiencies[i][j]=999;
    }
  }


  char config_name [25];
  Int_t config_int;
  

  if(config == 'd'){

    //runs all configurations in current folder
    for(config='B';config<='Z';config++){
      config_int = config-66;
      analyse_config(efficiencies,config,config_int, title);
      config_name[config_int]=config;
    }
    
   
  
  }
  else{
    //runs configuration specified configuration, or send error
    //  analyse_config(config,title,efficiencies);
    config_int = config-66;
    analyse_config(efficiencies,config,config_int,title);
    config_name[config_int]=config;
  }


  Int_t n_points = 0;
  Float_t xvals [25]={};
  Float_t y1vals [25]={}; Float_t y2vals[25]={}; Float_t y3vals[25]={};
  
  for(Int_t i = 0; i<25; i++){
  
    if(efficiencies[0][i]!=999){
      xvals[n_points] = i+66;
      y1vals[n_points] = efficiencies[0][i];
      y2vals[n_points] = efficiencies[1][i];
      y3vals[n_points] = efficiencies[2][i];
      n_points++;
    }
  }


  TGraph* comp1 = new TGraph(n_points,xvals,y1vals);
  comp1->SetTitle("Rel efficiency of Cofig X vs Config A");
  comp1->SetLineColor(38);
  comp1->SetLineStyle(1);
  comp1->SetLineWidth(5);
  comp1->SetMarkerStyle(21);
  comp1->SetMarkerColor(9);
  comp1->SetFillStyle(0);
  comp1->SetFillColor(0);
  

  TGraph* comp2 = new TGraph(n_points,xvals,y2vals);
  comp2->SetTitle("Rel efficiency of LaBr3 detectors in Config X vs Config A");
  comp2->SetLineColor(30);
  comp2->SetLineStyle(2);
  comp2->SetLineWidth(4);
  comp2->SetMarkerStyle(22);
  comp2->SetMarkerColor(8);
  comp2->SetFillStyle(0);
  comp2->SetFillColor(0);
  
  TGraph* comp3 = new TGraph(n_points,xvals,y3vals);
  comp3->SetTitle("Rel. efficiency of LaBr3 detectors in Config X vs Config X");
  comp3->SetLineColor(1);
  comp3->SetLineStyle(3);
  comp3->SetLineWidth(2);
  comp3->SetMarkerStyle(2);
  comp3->SetMarkerColor(1);
  comp3->SetFillStyle(0);
  comp3->SetFillColor(0);
  
  TMultiGraph *mg = new TMultiGraph("Relative efficiencies","Relative efficiencies");

  mg->Add(comp1);
  mg->Add(comp2);
  mg->Add(comp3);


  /*
  TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry("comp1","Relative efficiency of Cofiguration X vs Configuration A");
  legend->AddEntry("comp2","Relative efficiency of LaBr3 detectors in Config X vs Configuration A");
  legend->AddEntry("comp3","Relative efficiency of LaBr3 detectors in Config X vs Configuration X");
  */


  TCanvas *comp_canv = new TCanvas("comp","comp",600,500);
  gStyle->SetOptStat(0);
  mg->Draw("APL");
  // legend->Draw();
  TLegend *leg = comp_canv->BuildLegend();
  leg->SetFillStyle(0);


  TFile *ofile = new TFile(title+".root","UPDATE","Configurations");
  comp_canv->Write();
  
  ofile->Close();


}




















void  analyse_config(Float_t **efficiencies,TString config  = "A", Int_t config_int, TString title){
  
  
  TString inputFilePath = "~/geant_ruiz/output/";
 
  TString & inputFileName = "Config_"+config+".root";
  TString inputFile = inputFilePath+inputFileName;
  // Open all_BGO file
  TString exist = gSystem->FindFile(inputFilePath,inputFileName);
  if(exist=="") break;
  
  TFile *fA = new TFile(inputFilePath+"Config_A.root");
  ;  if (fA->IsZombie()) return;

  Float_t e_bgo_A;

  // Open trees and branches
  TTree *h1001A = (TTree*)fA->Get("h1001");
  //  h1001A->SetBranchAddress("e_bgo_first",&e_bgo_A);
  h1001A->SetBranchAddress("e0_conv",&e_bgo_A);

  
  TFile *fX = new TFile(inputFile);
  if(fX->IsZombie())return;

  Int_t n_entries, num_scnt, islabr3, isbgo;
  char configs [30] = {};
  Float_t e_bgo_X;

  TTree *h1001X = (TTree*)fX->Get("h1001");
  h1001X->SetBranchAddress("num_bgo_first",&num_scnt);
  h1001X->SetBranchAddress("islabr3",&islabr3);
  h1001X->SetBranchAddress("isbgo",&isbgo);
  //  h1001X->SetBranchAddress("e_bgo_first",&e_bgo_X);
  h1001X->SetBranchAddress("e0_conv",&e_bgo_X);  

  n_entries = (Int_t)h1001X->GetEntries();

  
  
  for(Int_t i = 0; i<30; i++){
    for(Int_t j=0; j<n_entries; j++){
      h1001X->GetEntry(j);
      
      if(num_scnt==i+1 & (isbgo^islabr3)){
	if(isbgo) configs[i]='B';
	else if(islabr3) configs[i]='L';
	break;
      }
      else configs[i]=999;  //no material data found for this scintillator.
    }

  }
  
  Float_t low_limit_x=0;
  Float_t high_limit_x=0;
  Float_t low_limit_a=0;
  Float_t high_limit_a=0;
  for(Int_t i=0;i<n_entries;i++){
    h1001A->GetEntry(i);
    h1001X->GetEntry(i);
    
    if(e_bgo_A<low_limit_a) low_limit_a=e_bgo_A;
    if(e_bgo_X<low_limit_x) low_limit_x=e_bgo_X;
    if(e_bgo_A>high_limit_a) high_limit_a=e_bgo_A;
    if(e_bgo_X>high_limit_x) high_limit_x=e_bgo_X;
  }


  TH1F *e_gamma_Xhist = new TH1F("e_gamma_firstX","temp title;temp axis;temp axis",100,low_limit_x,high_limit_x);
  TH1F *e_gamma_Ahist = new TH1F("e_gamma_firstA","temp title;temp axis;temp axis",100,low_limit_a,high_limit_a);
  TH1F *e_gamma_Xonly = new TH1F("e_gamma_onlyX","temp title;temp axis;temp axis",100,low_limit_x,high_limit_x);
 

  for(Int_t i=0; i<n_entries; i++){
    h1001A->GetEntry(i);
    h1001X->GetEntry(i);

    
    if(e_bgo_X!=0){
      e_gamma_Xhist->Fill(e_bgo_X);
      if(islabr3) e_gamma_Xonly->Fill(e_bgo_X);
    }
    if(e_bgo_A!=0) e_gamma_Ahist->Fill(e_bgo_A);
  }

  THStack *a = new THStack("a"+config,"Entire configuration "+config);
  THStack *b = new THStack("b"+config,"Only LaBr3 in configuration "+config);
  THStack *c = new THStack("c"+config,"Both");

  //Config_X
  e_gamma_Xhist->SetLineColor(kRed);
  e_gamma_Xhist->SetFillStyle(3001);
  e_gamma_Xhist->SetFillColor(kRed);
  //Config A
  e_gamma_Ahist->SetLineColor(kBlue);
  e_gamma_Ahist->SetFillStyle(3001);
  e_gamma_Ahist->SetFillColor(kBlue);
  //Only LaBr3
  e_gamma_Xonly->SetLineColor(kGreen);
  e_gamma_Xonly->SetFillStyle(3001);
  e_gamma_Xonly->SetFillColor(kGreen);
  
  a->Add(e_gamma_Ahist);
  a->Add(e_gamma_Xhist);
  b->Add(e_gamma_Ahist);
  b->Add(e_gamma_Xonly);
  c->Add(e_gamma_Ahist);
  c->Add(e_gamma_Xhist);
  c->Add(e_gamma_Xonly);
    
  
  TFile *ofile = new TFile(title+".root","UPDATE","Configurations");
  
  TCanvas *c1  = new TCanvas("Configuration"+config,"multipads",40,40,678,800);
  c1->Divide(2,1,0.005,0.005);

  c1->cd(1);
  a->Draw("nostack");


  c1->cd(2);

  c->Draw("nostack");
  /*
  c1->cd(3);
  c->Draw("nostack");
  */
  
  c1->Write();

  ofile->Close();
  
  Float_t IntX = e_gamma_Xhist->Integral();
  Float_t IntA = e_gamma_Ahist->Integral();
  Float_t IntXonly = e_gamma_Xonly->Integral();

  cout<<"<--------- Configuration "<<config<<" --------->"<<endl;
  cout<<"Counts of Configuration A\t\t"<<IntA<<endl;
  cout<<"Counts of Configuration "<<config<<"\t\t"<<IntX<<endl;
  cout<<"Counts of Configuration "<<config<<" (only LaBr3)\t"<<IntXonly<<endl;
  cout<<"Relative efficiency\t\t\t"<<IntX/IntA<<endl;
  cout<<"Relative efficiency (only LaBr3)\t"<<IntXonly/IntA<<endl;
  cout<<"further statistic (specify)\t\t"<<IntXonly/IntX<<endl;
  

  efficiencies[0][config_int]=IntX/IntA;
  efficiencies[1][config_int]=IntXonly/IntA;
  efficiencies[2][config_int]=IntXonly/IntX;
  

  
}
