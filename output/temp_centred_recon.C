//CURRENT REACTION SETTING: 23Mg (p,g) 24Al

/*
USAGE OF ARGUMENTS: the argument to the recon() function is the name of the input file with .root extension.

This version of recon.C can only be used with ROOT Trees that contain the "binitenerg" branches, from which the initial beam energy can be extracted. Otherwiese, use old_recon.C.



PHOTON TIME OF FLIGHT ADJUSTMENTS: in this version of recon.C, a simple adjustments is carried out to account for the time of flight of radiated photons in the reconstruction process. The position vectors for the centres of the photodetectors are listed in the array coords[30][3]. The average time of flight of a photon (assuming vacuum) between the centre of the gas target (therefore assuming centred resonance) and the photodetectos is measured and subtracted from gammatof. ISSUES: this does not account for the volumes of the photodetectors and for the different probability in being detected by closer or farther detectors. Also, it is more accurate for centred resonances.

*/


#include "recon.h"



void temp_recon(TString inputFileName = "centred.root"){

  //E_i = filename_to_E(inputFileName);

  Float_t coords[30][3] = {{   0, -4.8,-15.3},
			   {   0,-10.1,-12.2},	//1
			   {   0,  5.0,-12.2},
			   {   0,  9.9, -9.2}, 	//3
			   {   0,  8.0, -3.1},
			   {   0,  8.0,  3.1}, 	//5
			   {   0,  9.9,  9.2},
			   {   0,-10.1, 12.2},	//7
			   {   0,  5.0, 12.2},
			   {   0, -4.8, 15.3}, 	//9
			   {  -4, -2.6, -9.2},
			   {   4, -2.6, -9.2}, 	//11
			   {  -4, -7.9, -6.1},
			   {   4, -7.9, -6.1},	//13
			   {  -4,  2.7, -6.1},
			   {   4,  2.7, -6.1},	//15
			   {  -4, -2.6, -3.1},
			   {   4, -2.6, -3.1}, 	//17
			   {  -4, -7.9,    0},
			   {   4, -7.9,    0}, 	//19
			   {  -4,  2.7,    0},
			   {   4,  2.7,    0}, 	//21
			   {  -4, -2.6,  3.1},
			   {   4, -2.6,  3.1}, 	//23
			   {  -4, -7.9,  6.1},
			   {   4, -7.9,  6.1}, 	//25
			   {  -4,  2.7,  6.1},
			   {   4,  2.7,  6.1},	//27
			   {  -4, -2.6,  9.2},
			   {   4, -2.6,  9.2}};	//29




  Float_t t_gamma=0;

  for(Int_t i=0;i<30;i++){
    t_gamma = t_gamma + sqrt((coords[i][0])**2+(coords[i][1])**2+(coords[i][2])**2);
  }

  t_gamma = t_gamma/(30*c);

  cout<<"t_gamma is\t"<<t_gamma<<" ns"<<endl;



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
  h1000->SetBranchAddress("binitenerg",&E_i);

  h1000->GetEntry(0);

  cout<<"ESTIMATED INITIAL BEAM ENERGY =\t"<<E_i<<endl;



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
  


  TH1F *recon = new TH1F("zrecon","Z-coordinate reconstruction;Distance from centre of gas target [cm];Counts",280,-7.,7.);
  TH1F *zint_hist = new TH1F("zint_hist","Z-coordinate of resonance position (zint);Distance from centre of gas target [cm];Counts",280,-7.,7.);
 
  
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
      //dt = gammatof-t1;
      dt = gammatof - t1 - t_gamma;
      z = (-SPow*pow(dt,2))/(2*m)*pow(c,2) + v_i*dt-half_thickness; // cm from centre of gas target
      recon->Fill(z);
      
      zint_hist->Fill(zint);

      if(i<51||(i%50==0 && i<301)||i%500==0){
	nevts[counter]=i;
	meanz[counter]=recon->GetMean();
	
	//	yerr[counter]=recon->GetMeanError();
	yerr[counter]=recon->GetRMS();
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


