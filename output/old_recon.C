//CURRENT REACTION SETTING: 23Mg (p,g) 24Al

/*
USAGE OF ARGUMENTS: the argument to the recon() function is the name of the input file with .root extension.

If current ROOT Tree has "binitenerg" branch, use recon.C instead of old_recon.C for more accurate initial beam energy.


NAMING OF INPUT FILE: input file can be called "centred.root" if resonance is centred. 
Otherwise, it has to be called "ff_direction.root" where ff is the distance of resonance position from centre of gas target (in cm) 
and direction is "upstream" or "downstream" from centre of gas target.

eg. if resonance happens at 3.5 cm downstream (or upstream) from centre of gas target, filename is "3.5_downstream.root" (or "3.5_upsteam.root"). 
The initial beam energy will be calculated automatically and the reconstruction adjusted.
*/




#include "recon.h"



void old_recon(TString inputFileName = "centred.root"){

   E_i = filename_to_E(inputFileName);



  cout<<"ESTIMATED INITIAL BEAM ENERGY =\t"<<E_i<<endl;

 
  TString inputFilePath = "./";

  TFile *f = new TFile(inputFilePath+inputFileName);

  if (f->IsZombie()) return;
  
  
  Float_t m = 22.9941237;
  Float_t zint,gammatof,dt,z,zbeam;
  
  TTree *h1000 = (TTree*)f->Get("h1000");
  TTree *h1001 = (TTree*)f->Get("h1001");

  
  h1000->SetBranchAddress("zint",&zint);
  h1001->SetBranchAddress("gammatof",&gammatof);
  h1000->SetBranchAddress("z",&zbeam);
  
  m=m*amu;//Converts mass in MeV/c^2
  

  Int_t nentries = (Int_t)h1000->GetEntries();
  
  TH1F *zbeam_hist = new TH1F("zbeam","filler",10000,-200.,-50.);

  for (Int_t i=0;i<nentries;i++){
    h1000->GetEntry(i);
    zbeam_hist->Fill(zbeam);
  }
  zbeam = zbeam_hist->GetMean();

  

  if(zbeam==0)
    {
      zbeam=-88.*.999; 

      cout<<"Initial beam position estimated since missing on current ROOT Tree"<<endl;
    }


  //Calculating expected initial velocity of incident particles (non-relativistically)
  Float_t v_i = sqrt(2*E_i/m)*c; // cm/ns
  Float_t t1 = ((-zbeam-half_thickness)/v_i); //time to reach beginning of gas target (ns)
  


  TH1F *recon = new TH1F("zrecon","Z-coordinate reconstruction;Distance from centre of gas target [cm];Counts",10000,-7.,7.);
  TH1F *zint_hist = new TH1F("zint_hist","Z-coordinate of resonance position (zint);Distance from centre of gas target [cm];Counts",10000,-7.,7.);
 
  
  nentries = (Int_t)h1001->GetEntries();

  //Data for graph
  const Int_t num_points_max = nentries;
  Float_t nevts [num_points_max] = {};
  Float_t meanz [num_points_max] = {};
  Float_t xerr [num_points_max] = {};
  Float_t yerr [num_points_max] = {};
  Int_t counter = 0;

  
  Float_t rel_diff [num_points_max] = {};
  Float_t nevts2 [num_points_max] = {};
  Float_t yerr2 [num_points_max] ={};

  
  
  for (Int_t i=0;i<nentries;i++){
    h1000->GetEntry(i);
    h1001->GetEntry(i);
    if(gammatof<200. && gammatof>80.){
      dt = gammatof-t1;
      z = (-SPow*pow(dt,2))/(2*m)*pow(c,2) + v_i*dt-half_thickness; // cm from centre of gas target
      recon->Fill(z);
      
      zint_hist->Fill(zint);

      if(i<51||(i%50==0 && i<301)||i%500==0){
	nevts[counter]=i;
	meanz[counter]=recon->GetMean();
	
	yerr[counter]=recon->GetMeanError();
	//yerr[counter]=recon->GetRMS();
	counter++;
	//cout<<counter<<"\tzint\t"<<zint<<endl;
      }
      
      rel_diff[i] = (zint-z)/zint;
      nevts2[i] = i+1;
      


    }
  }

  
  // TFile outputfile(outputFileName.c_str(),"RECREATE");


  c1 = new TCanvas("canvas","multipads",900,700);
  c1->Divide(2,2,0.005,0.005);

  c1->cd(1);
  //  c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
  c1->cd(1).SetLogx();  
  gPad->SetTickx(2);
  gr = new TGraphErrors(counter,nevts,meanz,xerr,yerr);
  gr->SetTitle("Z-reconstruction mean vs number of events;N events;Distance from centre of gas target [cm]");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("AP");
  c1->Update();
  
  c1->cd(2);
  //  c2 = new TCanvas("c2","Afefef Simple Graph with error bars",200,10,700,500);
  recon->Draw();
 

  c1->cd(3);
  zint_hist->Draw();



  c1->cd(4);  
  h1000->Draw("zint","react");
  /*
  c1->cd(4);
  c1->cd(4).SetLogx();
  gr2 = new TGraphErrors(counter,nevts
  */
  cout<<"E_i\t"<<E_i<<"\t\tMeV\nSPow\t"<<SPow<<"\tMeV/cm\nm\t"<<m<<"\t\tMeV/c^2"<<endl;

}
















/*


Double_t pos_to_E(Double_t position){

  Double_t Initial_beam_energy = (position + half_thickness) * SPow + Resonance_energy;

  return Initial_beam_energy;
}


Double_t filename_to_E(TString inputFileName){
  

  if(inputFileName.Index("centred") != string::npos ){
    cout << "The resonance is centred."<<endl;
    
    E_i = pos_to_E(0);
  }
  else if(inputFileName.Index("_upstream") != string::npos ){
    cout << "The resonance is upstream."<<endl;
    
    TObjArray *tx = inputFileName.Tokenize("_");
    
    Double_t value = -((TObjString *)(tx->At(0)))->String().Atof();
    
    E_i = pos_to_E(value);
    
  }
  else if(inputFileName.Index("_downstream") != string::npos ){
    cout << "The resonance is downstream."<<endl;
    
    TObjArray *tx = inputFileName.Tokenize("_");
    
    Double_t value = ((TObjString *)(tx->At(0)))->String().Atof();
    
    E_i = pos_to_E(value);
  }
  else{
    cout << "Please, check file name for correctness." << endl;
  }

  return E_i;
}
*/
