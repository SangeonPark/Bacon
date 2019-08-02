//Driver to run Bacon

#include "../include/GenLoader.hh"
#include "../include/EvtLoader.hh"
#include "../include/MuonLoader.hh"
#include "../include/ElectronLoader.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <iostream>

//Object Processors
GenLoader       *fGen      = 0;
EvtLoader       *fEvt      = 0;
MuonLoader      *fMuon     = 0;
ElectronLoader  *fElectron = 0;

TTree* load(std::string iName) {
  TFile *lFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  return lTree;
}
int main( int argc, char **argv ) {
  gROOT->ProcessLine("#include <vector>");
  int maxEvents     = atoi(argv[1]);
  std::string lName = argv[2];
  bool        lGen  = atoi(argv[3]);
  //void runBacon(int iNEvents=10,std::string lName="test.root",bool lGen=false) {
  TTree *lTree = load(lName);
  if(lTree->GetEntries() < maxEvents || maxEvents == -1) maxEvents = lTree->GetEntries();
  //Declare Readers
  fEvt      = new EvtLoader     (lTree);
  fMuon     = new MuonLoader    (lTree);
  fElectron = new ElectronLoader(lTree);
  if(lGen) fGen      = new GenLoader     (lTree);

  TFile *lFile = new TFile("Output.root","RECREATE");
  TTree *lOut  = new TTree("Tree","Tree");
  //Setup Tree
  fEvt ->setupTree      (lOut);
  fEvt ->setupRecoilTree();
  fMuon->setupTree      (lOut);
  if(lGen) fGen ->setupTree      (lOut);
  if(lGen) fGen ->setupRecoilTree();
  for(int i0 = 0; i0 < maxEvents; i0++) {
    if(i0 % 1000 == 0) std::cout << "===> Processed " << i0 << " - Done : " << (float(i0)/float(maxEvents)) << std::endl;
    //Load Ntuple & Select 1 good muon
    fMuon    ->load(i0);
    if(!fMuon->selectDiMu()) continue;
    //Load the rest and fill
    fEvt    ->load(i0);
    fEvt    ->fillEvent();
    //Get Muon and comupte recoil variables
    TLorentzVector lMuon = fMuon->muon();
    fEvt    ->fillRecoil(lMuon);
    if(lGen) fGen->load(i0);
    if(lGen) fGen->selectBoson();
    if(lGen) fGen->fillRecoil(lMuon);
    lOut->Fill();
  }
  lFile->cd();
  lOut->Write();
  lFile->Close();
}
