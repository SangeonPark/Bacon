#ifndef BACONPROD_NTUPLER_FILLERPF_HH
#define BACONPROD_NTUPLER_FILLERPF_HH

#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"

class TClonesArray;


namespace baconhep
{
  class FillerPF
  {
    public:
      FillerPF();
      ~FillerPF();
      
       void fill(TClonesArray       *array,    // output array to be filled
		 TClonesArray       *iVtxCol,
		 const edm::Event   &iEvent);  // event info
    //Useful tools
    float depthDeltaR(const reco::PFCandidate *iPF,const reco::PFRecHitCollection &iPFCol,double iDR=0.08) ;
    float timeDeltaR (const reco::PFCandidate *iPF,const reco::PFRecHitCollection &iPFCol,double iDR=0.08) ;
    float depth(const reco::PFCandidate *iPF);
    float time (const reco::PFCandidate *iPF);
      
      // EDM object collection names
      std::string fPFName;
      std::string fPVName;
      bool        fAddDepthTime;
  };
}
#endif
