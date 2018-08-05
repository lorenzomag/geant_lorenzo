
Double_t filename_to_E(TString inputFileName);
Double_t pos_to_E(Double_t position);
Float_t indef_integral(Float_t x, Float_t y, Float_t z, Float_t ext);


//Constants

Float_t SPow = 0.029402297; // Gas stopping power MeV/cm
Float_t amu = 931.5; // MeV/c^2
Float_t c = 29979245800.*pow(10,-9); // cm/ns
Float_t half_thickness = 6.14136028; // half-thickness of gas target
Double_t Resonance_energy = 11.41164; //for 23Mg(p,g)24Al
Float_t E_i;





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





Float_t indef_integral(Float_t x, Float_t y, Float_t z, Float_t ext){

  Float_t m = (x**2+y**2)**2;
  
  ext = -ext;
  
  Float_t Int_a = 1./2. * (ext-z)*sqrt(m+(z-ext)**2) - 1./2. * m * log(sqrt(m + (z-ext)**2)+z-ext);
  Float_t Int_c = 1./2. * sqrt(m + (z - ext)**2)* (-z + ext) - 1./2. * m * log(z + sqrt(m + (z - ext)**2) - ext);

  ext = -ext;
  
  Float_t Int_b = 1./2. * (ext-z)*sqrt(m+(z-ext)**2) - 1./2. * m * log(sqrt(m + (z-ext)**2)+z-ext);

  Float_t Integral = Int_b-Int_a;

  cout<<"Int_a\t"<<Int_a<<"\nInt_c\t"<<Int_c<<"\nInt_b\t"<<Int_b<<"\nIntegral\t"<<Integral<<"\n"<<endl;

  return Integral;
}
