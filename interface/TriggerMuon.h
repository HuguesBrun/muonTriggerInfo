// -*- C++ -*-
//
// Package:    TriggerMuon
// Class:      TriggerMuon
// 
/**\class TriggerMuon TriggerMuon.cc hugues/TriggerMuon/src/TriggerMuon.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hugues Brun
//         Created:  Mon Jul 29 11:43:25 CEST 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"


#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

// ROOT stuff !
#include "TH1D.h"
#include <map>
#include "TFile.h"
#include <math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TString.h"
#include "TTree.h"
//
// class declaration
//
typedef std::vector<edm::InputTag> vtag;

class TriggerMuon : public edm::EDAnalyzer {
   public:
      explicit TriggerMuon(const edm::ParameterSet&);
      ~TriggerMuon();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void beginEvent();
    virtual void endEvent();
      virtual float deltaR(float , float , float , float );
      virtual float deltaPhi(float, float);


      // ----------member data ---------------------------
    vtag muonProducers_;
    // input tags

    edm::InputTag               triggerResultsLabel_;
    edm::InputTag               triggerSummaryLabel_;
    
    std::string outputFile_; // output file
    
    std::vector<edm::InputTag> moduleLabels;
    
    std::vector<TString> HLT_name;
    std::vector<int> theBitCorr;
    std::vector<std::string> HLT_triggerObjects;
    
    HLTConfigProvider hltConfig;
    
    
    
    // tree leaf declaration :
    int T_Event_RunNumber;
    int T_Event_EventNumber;
    int T_Event_LuminosityBlock;
    
    int T_Event_HLT_Mu17_Mu8;
    int T_Event_HLT_Mu17_TkMu8;
    int T_Event_HLT_Mu17;
    int T_Event_HLT_Mu8;
    int T_Event_HLT_Mu8_Ele17;
    int T_Event_HLT_Ele8_Mu17;
    int T_Event_HLT_IsoMu24;
    int T_Event_HLT_IsoMu24_2p1;
    
    
    //muons infos
    std::vector<float>*T_Muon_Eta;
 	std::vector<float>*T_Muon_Phi;
	std::vector<float>*T_Muon_Energy;
	std::vector<float>*T_Muon_Et;
	std::vector<float>*T_Muon_Pt;
	std::vector<float>*T_Muon_Px;
	std::vector<float>*T_Muon_Py;
	std::vector<float>*T_Muon_Pz;
	std::vector<float>*T_Muon_Mass;
    
    std::vector<bool> *T_Muon_IsGlobalMuon;
    std::vector<bool> *T_Muon_IsTrackerMuon;
    std::vector<bool> *T_Muon_IsPFMuon;
    std::vector<bool> *T_Muon_IsCaloMuon;
    std::vector<bool> *T_Muon_IsStandAloneMuon;
    std::vector<bool> *T_Muon_IsMuon;
    std::vector<bool> *T_Muon_IsGlobalMuon_PromptTight;
    std::vector<bool> *T_Muon_IsTrackerMuonArbitrated;
    std::vector<int>  *T_Muon_numberOfChambers;
    std::vector<int>  *T_Muon_numberOfChambersRPC;
    std::vector<int>  *T_Muon_numberOfMatches;
    std::vector<int>  *T_Muon_numberOfMatchedStations;
    std::vector<int>  *T_Muon_charge;
    
    
    std::vector<bool> *T_Muon_TMLastStationTight;
    std::vector<float> *T_Muon_globalTrackChi2;
    std::vector<int>  *T_Muon_validMuonHits;
    std::vector<float> *T_Muon_trkKink;
    std::vector<int>  *T_Muon_trkNbOfTrackerLayers;
    std::vector<int>  *T_Muon_trkNbOfValidTrackeHits;
    std::vector<int>  *T_Muon_trkValidPixelHits;
    std::vector<float> *T_Muon_trkError;
    std::vector<float> *T_Muon_dB;
    std::vector<float> *T_Muon_dzPV;
    std::vector<float> *T_Muon_dBstop;
    std::vector<float> *T_Muon_dzstop;
    
    // PF isolation
    std::vector<float> *T_Muon_chargedHadronIsoR04;
    std::vector<float> *T_Muon_neutralHadronIsoR04;
    std::vector<float> *T_Muon_photonIsoR04;
    std::vector<float> *T_Muon_chargedHadronIsoPUR04;
    
    std::vector<float> *T_Muon_chargedHadronIsoR03;
    std::vector<float> *T_Muon_neutralHadronIsoR03;
    std::vector<float> *T_Muon_photonIsoR03;
    std::vector<float> *T_Muon_chargedHadronIsoPUR03;
    
    std::vector<float> *T_Muon_isoR03_emEt;
    std::vector<float> *T_Muon_isoR03_hadEt;
    std::vector<float> *T_Muon_isoR03_hoEt;
    std::vector<float> *T_Muon_isoR03_sumPt;
    std::vector<int> *T_Muon_isoR03_nTracks;
    std::vector<int> *T_Muon_isoR03_nJets;
    std::vector<float> *T_Muon_isoRingsMVA;
    
    //trigger matching
    std::vector<int> *T_Muon_HLT_Mu17_TkMu8_Mu17Leg;
    std::vector<int> *T_Muon_HLT_Mu17_TkMu8_Mu8Leg;
    std::vector<int> *T_Muon_HLT_Mu17_Mu8_Mu17Leg;
    std::vector<int> *T_Muon_HLT_Mu17_Mu8_Mu8Leg;
    std::vector<int> *T_Muon_HLT_Mu17_obj;
    std::vector<int> *T_Muon_HLT_Mu8_obj;
    std::vector<int> *T_Muon_HLT_Mu8_Ele17_Mu8Leg;
    std::vector<int> *T_Muon_HLT_Ele8_Mu17_Mu17Leg;
    std::vector<int> *T_Muon_HLT_IsoMu24;
    std::vector<int> *T_Muon_HLT_IsoMu24_2p1;
    
    
    // root file to store histograms
    TFile*  rootFile_;
    
    //Tree
    TTree* mytree_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

