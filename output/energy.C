#include "recon.h"


void energy(){

  Float_t distance;
  Double_t energy, current_pos;
  TString input;

  cout<<"Please, insert desired position (in cm) of resonance wrt cetre of gas target or request a 'table' of energies:\t";
  cin>>input;
  cout<<endl;


  if(input == "table"){
    for(Int_t i; i<27; i++){
      current_pos = 6.5.-i*0.5;
      energy = pos_to_E(current_pos);

      cout<<current_pos<<"\tcm\t ==>\t"<<energy<<"\tMeV"<<endl;
    }
  }
  else{
    distance = input.Atof();

   
    energy = pos_to_E(distance);
    cout<<"\nRequired energy:\n"<<energy<<"\tMeV"<<endl;
  }

  cout<<endl;

}
