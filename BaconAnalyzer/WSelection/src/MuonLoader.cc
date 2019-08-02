#include "../include/MuonLoader.hh"
#include "TMath.h"

using namespace baconhep;

MuonLoader::MuonLoader(TTree *iTree) {
  fMuons  = new TClonesArray("baconhep::TMuon");
  iTree->SetBranchAddress("Muon",       &fMuons);
  fMuonBr  = iTree->GetBranch("Muon");
}
MuonLoader::~MuonLoader() {
  delete fMuons;
  delete fMuonBr;
}
void MuonLoader::reset() {
  fPt   = 0;
  fEta  = 0;
  fPhi  = 0;
}
void MuonLoader::setupTree(TTree *iTree) {
  reset();
  fTree = iTree;
  fTree->Branch("pt_1"  ,&fPt ,"fPt/F");
  fTree->Branch("eta_1" ,&fEta,"fEta/F");
  fTree->Branch("phi_1" ,&fPhi,"fPhi/F");
}
void MuonLoader::load(int iEvent) {
  fMuons   ->Clear();
  fMuonBr ->GetEntry(iEvent);
}
bool MuonLoader::selectSingleMu() {
  reset();
  TMuon *lMuon = 0;
  int lCount = 0;
  for  (int i0 = 0; i0 < fMuons->GetEntriesFast(); i0++) {
    TMuon *pMuon = (TMuon*)((*fMuons)[i0]);
    if(!passLoose(pMuon)) continue;
    lMuon = pMuon;
    lCount++;
    if(lCount > 1) return false;
  }
  if(lMuon == 0) return false;
  fPt  = lMuon->pt;
  fEta = lMuon->eta;
  fPhi = lMuon->phi;
  return true;
}
bool MuonLoader::selectDiMu() {
  reset();
  TMuon *lMuon1 = 0;
  TMuon *lMuon2 = 0;
  int lCount = 0;
  for  (int i0 = 0; i0 < fMuons->GetEntriesFast(); i0++) {
    TMuon *pMuon = (TMuon*)((*fMuons)[i0]);
    if(!passTight(pMuon)) continue;
    if(lCount == 0) lMuon1 = pMuon;
    if(lCount == 1) lMuon2 = pMuon;
    lCount++;
    if(lCount > 2) return false;
  }
  if(lMuon1 == 0) return false;
  TLorentzVector v1;
  TLorentzVector v2;
  TLorentzVector lDiMuon;
  v1.SetPtEtaPhiM(lMuon1->pt,lMuon1->eta,lMuon1->phi,0.105);
  v2.SetPtEtaPhiM(lMuon2->pt,lMuon2->eta,lMuon2->phi,0.105);
  lDiMuon = v1 + v2;
  if(lDiMuon.Pt()<30 || abs(lDiMuon.M()-91.2)>20) return false;
  fPt = lDiMuon.Pt();
  fEta = lDiMuon.Eta();
  fPhi = lDiMuon.Phi();
  fM = lDiMuon.M();
  return true;
}
bool MuonLoader::vetoMu(TMuon *iMuon) {
  for  (int i0 = 0; i0 < fMuons->GetEntriesFast(); i0++) {
    TMuon *pMuon = (TMuon*)((*fMuons)[i0]);
    if(passLoose(pMuon)) return true;
  }
  return false;
}
//H=>ZZ Mu Id
bool MuonLoader::passLoose(TMuon *muon) {
  if(!(muon->typeBits     & kGlobal || muon->typeBits & kTracker))  return false;
  if(!(muon->selectorBits & kAllArbitrated))                        return false;
  if(!(muon->typeBits     & kPFMuon))                               return false;
  if(fabs(muon->dz)> 1.0)                                           return false;
  double chargedIso = muon->chHadIso04;
  double neutralIso = TMath::Max(muon->gammaIso04 + muon->neuHadIso04 - 0.5 * muon->puIso04, 0.0);
  double totalIso   = chargedIso+neutralIso;
  if(totalIso/muon->pt > 0.4) return false;
  return true;
}
bool MuonLoader::passTight(TMuon *iMuon) {
  if(!(iMuon->typeBits & kGlobal))  return false;
  if(fabs(iMuon->dz)  > 0.2)        return false;
  if(fabs(iMuon->d0)  > 0.045)      return false;
  if(iMuon->muNchi2        > 10)    return false;
  if(iMuon->nValidHits     < 1)     return false;
  if(iMuon->nMatchStn      < 2)     return false;
  if(iMuon->nPixHits       < 1)     return false;
  if(iMuon->nTkLayers      < 6)     return false;
  if(!(iMuon->typeBits & kPFMuon))  return false;


  double chargedIso = iMuon->chHadIso04;
  double neutralIso = TMath::Max(iMuon->gammaIso04 + iMuon->neuHadIso04 - 0.5 * iMuon->puIso04, 0.0);
  double totalIso   = chargedIso+neutralIso;
  if(totalIso/iMuon->pt > 0.15) return false;
  return true;
}
TLorentzVector MuonLoader::muon() {
  TLorentzVector lMuon;
  lMuon.SetPtEtaPhiM(fPt,fEta,fPhi,0.105);
  return lMuon;
}
TLorentzVector MuonLoader::dimuon() {
  TLorentzVector lMuon;
  lMuon.SetPtEtaPhiM(fPt,fEta,fPhi,fM);
  return lMuon;
}
