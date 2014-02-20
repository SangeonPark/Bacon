#ifndef BACONPROD_NTUPLER_FILLEREVENTINFO_HH
#define BACONPROD_NTUPLER_FILLEREVENTINFO_HH

#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class TEventInfo;  // foward declaration
  class FillerEventInfo
  {
    public:
      FillerEventInfo();
      ~FillerEventInfo();
      
      void fill(TEventInfo         *evtInfo,       // output object to be filled
                const edm::Event   &iEvent,        // EDM event info
		const reco::Vertex &pv,            // event primary vertex
		const bool          hasGoodPV,     // flag for if PV passing cuts is found
		const TriggerBits   triggerBits);  // bits for corresponding fired triggers
	       
      void computeTrackMET(const reco::Vertex &pv, 
                           const reco::PFCandidateCollection *pfCandCol,
                           float &out_met, float &out_metphi);
    
    
      // EDM object collection names
      std::string fPFCandName;
      std::string fPUInfoName;
      std::string fBSName;
      std::string fPFMETName;
      std::string fMVAMETName;
      std::string fMVAMETUName;
      std::string fRhoIsoName;
      std::string fRhoJetName;
      bool        fFillMET;
  };
}
#endif
