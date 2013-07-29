#include "../interface/TriggerMuon.h"



TriggerMuon::TriggerMuon(const edm::ParameterSet& iConfig)

{
    
    triggerResultsLabel_    = iConfig.getParameter<edm::InputTag>("TriggerResults");
    triggerSummaryLabel_    = iConfig.getParameter<edm::InputTag>("HLTTriggerSummaryAOD");
    muonProducers_			= iConfig.getParameter<vtag>("muonProducer");
    
    outputFile_   = iConfig.getParameter<std::string>("outputFile");
    rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");
    
    
    HLT_name.push_back("HLT_Mu17_TkMu8_v");//0
    HLT_name.push_back("HLT_Mu17_Mu8_v");//1
    HLT_name.push_back("HLT_Mu17_v");//2
    HLT_name.push_back("HLT_Mu8_v");//3
    HLT_name.push_back("HLT_IsoMu24_v");//4
    HLT_name.push_back("HLT_IsoMu24_eta2p1_v");//5
    
    
    HLT_triggerObjects.push_back("hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17");// 0 -> DoubleMu17Mu8_Mu17
    HLT_triggerObjects.push_back("hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17");// 1 -> DoubleMu17Mu8_Mu17
    HLT_triggerObjects.push_back("hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8");// 2-> DoubleMu17Mu8_Mu8
    HLT_triggerObjects.push_back("hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8");//3-> DoubleMu17Mu8_Mu8
    HLT_triggerObjects.push_back("hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17");// 4-> DoubleMu17TkMu8_Mu17leg
    HLT_triggerObjects.push_back("hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17");// 5-> DoubleMu17TkMu8_Mu17leg
    HLT_triggerObjects.push_back("hltDiMuonGlbFiltered17TrkFiltered8");// 6-> DoubleMu17TkMu8_TkMu8leg
    HLT_triggerObjects.push_back("hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15");// 7 -> IsoMu24
    HLT_triggerObjects.push_back("hltL3fL1sMu12L3Filtered17");// 8 -> Mu17
    HLT_triggerObjects.push_back("hltL3fL1sMu3L3Filtered8");// 9 -> Mu8

    
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
    

}


TriggerMuon::~TriggerMuon()
{
 
delete rootFile_;

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerMuon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    beginEvent();
    //recup the muon collection 
    edm::Handle < std::vector <reco::Muon> > recoMuons;
    edm::InputTag muonProducer = muonProducers_.at(0);
	iEvent.getByLabel(muonProducer, recoMuons);
    
    
    //read the trigger results
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(triggerResultsLabel_, triggerResults);
    
    edm::Handle<trigger::TriggerEvent> triggerSummary;
    iEvent.getByLabel(triggerSummaryLabel_, triggerSummary);
    
    // get the HLT config is changed 
    bool changedConfig = false;
    if (!hltConfig.init(iEvent.getRun(), iSetup, triggerResultsLabel_.process(), changedConfig)) {
        edm::LogError("HLTMatchingFilter") << "Initialization of HLTConfigProvider failed!!";
        return;
    }
    if (changedConfig){
        std::cout << "the present menu is " << hltConfig.tableName() << std::endl;
        moduleLabels.clear();
        for (size_t m = 0 ; m < HLT_name.size() ; m++){
            for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
                if (TString(hltConfig.triggerNames()[j]).Contains(HLT_name[m])){
                    cout << j << " = " << hltConfig.triggerNames()[j] << endl;
                    theBitCorr.push_back(j);
                }
            }
        }
        for (unsigned int j=0; j<HLT_triggerObjects.size(); j++)
            moduleLabels.push_back(edm::InputTag(HLT_triggerObjects[j], "", triggerResultsLabel_.process()));
    }
    
    /// fill the in selected Objet the HLT filter we will use for the matching
    trigger::TriggerObjectCollection allTriggerObjects = triggerSummary->getObjects();
    trigger::TriggerObjectCollection selectedObjects;
    vector<int> theHLTcorr;
    for (size_t t=0; t<moduleLabels.size(); t++) {
        
        size_t filterIndex = (*triggerSummary).filterIndex(moduleLabels[t]);
        if (filterIndex < (*triggerSummary).sizeFilters()) {
            const trigger::Keys &keys = (*triggerSummary).filterKeys(filterIndex);
            
            for (size_t j = 0; j < keys.size(); j++) {
                trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
                selectedObjects.push_back(foundObject);
                theHLTcorr.push_back(t);
            }
        }
    }
    
    
    //fill event variables : 
    T_Event_RunNumber = iEvent.id().run();
    T_Event_EventNumber = iEvent.id().event();
    T_Event_LuminosityBlock = iEvent.id().luminosityBlock();
    
    T_Event_HLT_Mu17_Mu8 =         triggerResults->accept(theBitCorr[1]);
    T_Event_HLT_Mu17_TkMu8 =         triggerResults->accept(theBitCorr[0]);
    T_Event_HLT_Mu17 =         triggerResults->accept(theBitCorr[2]);
    T_Event_HLT_Mu8 =         triggerResults->accept(theBitCorr[3]);
    T_Event_HLT_IsoMu24 =         triggerResults->accept(theBitCorr[4]);
    T_Event_HLT_IsoMu24_2p1 =         triggerResults->accept(theBitCorr[5]);

    int nbMuons = recoMuons->size();
    for (int k = 0 ; k < nbMuons ; k++){// loop on the muons in the event
        const reco::Muon* muon = &((*recoMuons)[k]);
        
        T_Muon_Eta->push_back(muon->eta());
        T_Muon_Phi->push_back(muon->phi());
        T_Muon_IsGlobalMuon->push_back(muon->isGlobalMuon());
        T_Muon_IsPFMuon->push_back(muon->isPFMuon());
        T_Muon_IsTrackerMuon->push_back(muon->isTrackerMuon());
        T_Muon_IsCaloMuon->push_back(muon->isCaloMuon());
        T_Muon_IsStandAloneMuon->push_back(muon->isStandAloneMuon());
        T_Muon_IsMuon->push_back(muon->isMuon());
        T_Muon_Energy->push_back(muon->energy());
        T_Muon_Et->push_back(muon->et());
        T_Muon_Pt->push_back(muon->pt());
        T_Muon_Px->push_back(muon->px());
        T_Muon_Py->push_back(muon->py());
        T_Muon_Pz->push_back(muon->pz());
        T_Muon_Mass->push_back(muon->mass());
        T_Muon_charge->push_back(muon->charge());
        
        T_Muon_numberOfChambers->push_back(muon->numberOfChambers());
        T_Muon_numberOfChambersRPC->push_back(muon->numberOfChambersNoRPC());
        T_Muon_numberOfMatches->push_back(muon->numberOfMatches());
        T_Muon_numberOfMatchedStations->push_back(muon->numberOfMatchedStations());
        bool isMatchTheStation = muon::isGoodMuon(*muon, muon::TMOneStationTight);
        bool isGlobalMuonPT = muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight);
        bool isGlobalMuonArbitrated = muon::isGoodMuon(*muon, muon::TrackerMuonArbitrated);
        T_Muon_TMLastStationTight->push_back(isMatchTheStation);
        T_Muon_IsGlobalMuon_PromptTight->push_back(isGlobalMuonPT);
        T_Muon_IsTrackerMuonArbitrated->push_back(isGlobalMuonArbitrated);
        
        if (muon->globalTrack().isNull()) T_Muon_globalTrackChi2->push_back(-1); else T_Muon_globalTrackChi2->push_back(muon->globalTrack()->normalizedChi2());
        if (muon->globalTrack().isNull()) T_Muon_validMuonHits->push_back(-1); else T_Muon_validMuonHits->push_back(muon->globalTrack()->hitPattern().numberOfValidMuonHits());
        T_Muon_trkKink->push_back(muon->combinedQuality().trkKink);
        if (muon->muonBestTrack().isNull()) {
            T_Muon_trkNbOfTrackerLayers->push_back(-1);
            T_Muon_trkError->push_back(-1);
            T_Muon_trkValidPixelHits->push_back(-1);
            T_Muon_trkNbOfValidTrackeHits->push_back(-1);
        }
        else {
            T_Muon_trkNbOfTrackerLayers->push_back(muon->muonBestTrack()->hitPattern().trackerLayersWithMeasurement());
            T_Muon_trkError->push_back(muon->muonBestTrack()->ptError());
            T_Muon_trkValidPixelHits->push_back(muon->muonBestTrack()->hitPattern().numberOfValidPixelHits());
            T_Muon_trkNbOfValidTrackeHits->push_back(muon->muonBestTrack()->hitPattern().numberOfValidTrackerHits());
        }
        T_Muon_isoR03_emEt->push_back(muon->isolationR03().emEt);
        T_Muon_isoR03_hadEt->push_back(muon->isolationR03().hadEt);
        T_Muon_isoR03_hoEt->push_back(muon->isolationR03().hoEt);
        T_Muon_isoR03_sumPt->push_back(muon->isolationR03().sumPt);
        T_Muon_isoR03_nTracks->push_back(muon->isolationR03().nTracks);
        T_Muon_isoR03_nJets->push_back(muon->isolationR03().nJets);
        T_Muon_chargedHadronIsoR04->push_back(muon->pfIsolationR04().sumChargedHadronPt);
        T_Muon_neutralHadronIsoR04->push_back(muon->pfIsolationR04().sumNeutralHadronEt);
        T_Muon_photonIsoR04->push_back(muon->pfIsolationR04().sumPhotonEt);
        T_Muon_chargedHadronIsoPUR04->push_back(muon->pfIsolationR04().sumPUPt);
        T_Muon_chargedHadronIsoR03->push_back(muon->pfIsolationR03().sumChargedHadronPt);
        T_Muon_neutralHadronIsoR03->push_back(muon->pfIsolationR03().sumNeutralHadronEt);
        T_Muon_photonIsoR03->push_back(muon->pfIsolationR03().sumPhotonEt);
        T_Muon_chargedHadronIsoPUR03->push_back(muon->pfIsolationR03().sumPUPt);
        
        int pass_HLT_Mu17_TkMu8_Mu17Leg = 0;
        int pass_HLT_Mu17_TkMu8_Mu8Leg = 0;
        int pass_HLT_Mu17_Mu8_Mu17Leg = 0;
        int pass_HLT_Mu17_Mu8_Mu8Leg = 0;
        int pass_HLT_Mu17_Mu17_obj = 0;
        int pass_HLT_Mu17_Mu8_obj = 0;
        int pass_HLT_IsoMu24 = 0;
        int pass_HLT_IsoMu24_2p1 = 0;

        
        for (size_t t = 0 ; t < selectedObjects.size() ; t++){
            // cout << "eta = " << selectedObjects[t].eta() << " phi = " << selectedObjects[t].phi() << "muon Pt =" << muon->pt()<< "filter = " << HLT_triggerObjects[theHLTcorr[t]] << endl;
            // cout << "pt trigger" << selectedObjects[t].pt() << endl;
            float HLTdeltaR = deltaR(muon->phi(), selectedObjects[t].phi(), muon->eta(), selectedObjects[t].eta());
            float relatDeltaPt = fabs(selectedObjects[t].pt()-muon->pt())/muon->pt();
            if (HLTdeltaR < 0.1){
                if ((theHLTcorr[t] == 0)||(theHLTcorr[t] == 1))  pass_HLT_Mu17_Mu8_Mu17Leg = 1;
                if ((theHLTcorr[t] == 2)||(theHLTcorr[t] == 3))  pass_HLT_Mu17_Mu8_Mu8Leg = 1;
                if ((theHLTcorr[t] == 4)||(theHLTcorr[t] == 5))  pass_HLT_Mu17_TkMu8_Mu17Leg = 1;
                if (theHLTcorr[t] == 6)                          pass_HLT_Mu17_TkMu8_Mu8Leg = 1;
                if (theHLTcorr[t] == 7)                          pass_HLT_IsoMu24 = 1;
                if (theHLTcorr[t] == 8)                          pass_HLT_Mu17_Mu17_obj = 1;
                if (theHLTcorr[t] == 9)                          pass_HLT_Mu17_Mu8_obj = 1;
            }
        }
        T_Muon_HLT_Mu17_TkMu8_Mu17Leg->push_back(pass_HLT_Mu17_TkMu8_Mu17Leg);
        T_Muon_HLT_Mu17_TkMu8_Mu8Leg->push_back(pass_HLT_Mu17_TkMu8_Mu8Leg);
        T_Muon_HLT_Mu17_Mu8_Mu17Leg->push_back(pass_HLT_Mu17_Mu8_Mu17Leg);
        T_Muon_HLT_Mu17_Mu8_Mu8Leg->push_back(pass_HLT_Mu17_Mu8_Mu8Leg);
        T_Muon_HLT_Mu17_obj->push_back(pass_HLT_Mu17_Mu17_obj);
        T_Muon_HLT_Mu8_obj->push_back(pass_HLT_Mu17_Mu8_obj);
        T_Muon_HLT_IsoMu24->push_back(pass_HLT_IsoMu24);
    }
    
    mytree_->Fill();
    endEvent();
}


// ------------ method called once each job just before starting event loop  ------------
void 
TriggerMuon::beginJob()
{
    mytree_ = new TTree("eventsTree","");
    
    mytree_->Branch("T_Event_RunNumber", &T_Event_RunNumber, "T_Event_RunNumber/I");
    mytree_->Branch("T_Event_EventNumber", &T_Event_EventNumber, "T_Event_EventNumber/I");
    mytree_->Branch("T_Event_LuminosityBlock", &T_Event_LuminosityBlock, "T_Event_LuminosityBlock/I");
    
    mytree_->Branch("T_Event_HLT_Mu17_Mu8",&T_Event_HLT_Mu17_Mu8,"T_Event_HLT_Mu17_Mu8/I");
    mytree_->Branch("T_Event_HLT_Mu17_TkMu8",&T_Event_HLT_Mu17_TkMu8,"T_Event_HLT_Mu17_TkMu8/I");
    mytree_->Branch("T_Event_HLT_Mu17",&T_Event_HLT_Mu17,"T_Event_HLT_Mu17/I");
    mytree_->Branch("T_Event_HLT_Mu8",&T_Event_HLT_Mu8,"T_Event_HLT_Mu8/I");
    mytree_->Branch("T_Event_HLT_IsoMu24",&T_Event_HLT_IsoMu24,"T_Event_HLT_IsoMu24/I");
    mytree_->Branch("T_Event_HLT_IsoMu24_2p1",&T_Event_HLT_IsoMu24_2p1,"T_Event_HLT_IsoMu24_2p1/I");

    
    mytree_->Branch("T_Muon_Eta", "std::vector<float>", &T_Muon_Eta);
    mytree_->Branch("T_Muon_Phi", "std::vector<float>", &T_Muon_Phi);
    mytree_->Branch("T_Muon_Energy", "std::vector<float>", &T_Muon_Energy);
    mytree_->Branch("T_Muon_Et", "std::vector<float>", &T_Muon_Et);
    mytree_->Branch("T_Muon_Pt", "std::vector<float>", &T_Muon_Pt);
    mytree_->Branch("T_Muon_Px", "std::vector<float>", &T_Muon_Px);
    mytree_->Branch("T_Muon_Py", "std::vector<float>", &T_Muon_Py);
    mytree_->Branch("T_Muon_Pz", "std::vector<float>", &T_Muon_Pz);
    mytree_->Branch("T_Muon_Mass", "std::vector<float>", &T_Muon_Mass);
    mytree_->Branch("T_Muon_IsGlobalMuon", "std::vector<bool>", &T_Muon_IsGlobalMuon);
    mytree_->Branch("T_Muon_IsTrackerMuon", "std::vector<bool>", &T_Muon_IsTrackerMuon);
    mytree_->Branch("T_Muon_IsPFMuon", "std::vector<bool>", &T_Muon_IsPFMuon);
    mytree_->Branch("T_Muon_IsCaloMuon", "std::vector<bool>", &T_Muon_IsCaloMuon);
    mytree_->Branch("T_Muon_IsStandAloneMuon", "std::vector<bool>", &T_Muon_IsStandAloneMuon);
    mytree_->Branch("T_Muon_IsMuon", "std::vector<bool>", &T_Muon_IsMuon);
    mytree_->Branch("T_Muon_IsGlobalMuon_PromptTight", "std::vector<bool>", &T_Muon_IsGlobalMuon_PromptTight);
    mytree_->Branch("T_Muon_IsTrackerMuonArbitrated", "std::vector<bool>", &T_Muon_IsTrackerMuonArbitrated);
    mytree_->Branch("T_Muon_numberOfChambers", "std::vector<int>", &T_Muon_numberOfChambers);
    mytree_->Branch("T_Muon_numberOfChambersRPC", "std::vector<int>", &T_Muon_numberOfChambersRPC);
    mytree_->Branch("T_Muon_numberOfMatches", "std::vector<int>", &T_Muon_numberOfMatches);
    mytree_->Branch("T_Muon_numberOfMatchedStations", "std::vector<int>", &T_Muon_numberOfMatchedStations);
    mytree_->Branch("T_Muon_charge", "std::vector<int>", &T_Muon_charge);
    mytree_->Branch("T_Muon_TMLastStationTight", "std::vector<bool>", &T_Muon_TMLastStationTight);
    mytree_->Branch("T_Muon_globalTrackChi2", "std::vector<float>", &T_Muon_globalTrackChi2);
    mytree_->Branch("T_Muon_validMuonHits", "std::vector<int>", &T_Muon_validMuonHits);
    mytree_->Branch("T_Muon_trkKink", "std::vector<float>", &T_Muon_trkKink);
    mytree_->Branch("T_Muon_trkNbOfTrackerLayers", "std::vector<int>", &T_Muon_trkNbOfTrackerLayers);
    mytree_->Branch("T_Muon_trkNbOfValidTrackeHits", "std::vector<int>", &T_Muon_trkNbOfValidTrackeHits);
    mytree_->Branch("T_Muon_trkValidPixelHits", "std::vector<int>", &T_Muon_trkValidPixelHits);
    mytree_->Branch("T_Muon_trkError", "std::vector<float>", &T_Muon_trkError);
    mytree_->Branch("T_Muon_dB", "std::vector<float>", &T_Muon_dB);
    mytree_->Branch("T_Muon_dzPV", "std::vector<float>", &T_Muon_dzPV);
    mytree_->Branch("T_Muon_isoR03_emEt", "std::vector<float>", &T_Muon_isoR03_emEt);
    mytree_->Branch("T_Muon_isoR03_hadEt", "std::vector<float>", &T_Muon_isoR03_hadEt);
    mytree_->Branch("T_Muon_isoR03_hoEt", "std::vector<float>", &T_Muon_isoR03_hoEt);
    mytree_->Branch("T_Muon_isoR03_sumPt", "std::vector<float>", &T_Muon_isoR03_sumPt);
    mytree_->Branch("T_Muon_isoR03_nTracks", "std::vector<int>", &T_Muon_isoR03_nTracks);
    mytree_->Branch("T_Muon_isoR03_nJets", "std::vector<int>", &T_Muon_isoR03_nJets);
    mytree_->Branch("T_Muon_chargedHadronIsoR04", "std::vector<float>", &T_Muon_chargedHadronIsoR04);
    mytree_->Branch("T_Muon_neutralHadronIsoR04", "std::vector<float>", &T_Muon_neutralHadronIsoR04);
    mytree_->Branch("T_Muon_photonIsoR04", "std::vector<float>", &T_Muon_photonIsoR04);
    mytree_->Branch("T_Muon_chargedHadronIsoPUR04", "std::vector<float>", &T_Muon_chargedHadronIsoPUR04);
    mytree_->Branch("T_Muon_chargedHadronIsoR03", "std::vector<float>", &T_Muon_chargedHadronIsoR03);
    mytree_->Branch("T_Muon_neutralHadronIsoR03", "std::vector<float>", &T_Muon_neutralHadronIsoR03);
    mytree_->Branch("T_Muon_photonIsoR03", "std::vector<float>", &T_Muon_photonIsoR03);
    mytree_->Branch("T_Muon_chargedHadronIsoPUR03", "std::vector<float>", &T_Muon_chargedHadronIsoPUR03);
    mytree_->Branch("T_Muon_isoRingsMVA", "std::vector<float>", &T_Muon_isoRingsMVA);
    mytree_->Branch("T_Muon_HLT_Mu17_TkMu8_Mu17Leg", "std::vector<int>", &T_Muon_HLT_Mu17_TkMu8_Mu17Leg);
    mytree_->Branch("T_Muon_HLT_Mu17_TkMu8_Mu8Leg", "std::vector<int>", &T_Muon_HLT_Mu17_TkMu8_Mu8Leg);
    mytree_->Branch("T_Muon_HLT_Mu17_Mu8_Mu17Leg", "std::vector<int>", &T_Muon_HLT_Mu17_Mu8_Mu17Leg);
    mytree_->Branch("T_Muon_HLT_Mu17_Mu8_Mu8Leg", "std::vector<int>", &T_Muon_HLT_Mu17_Mu8_Mu8Leg);
    mytree_->Branch("T_Muon_HLT_Mu17_obj", "std::vector<int>", &T_Muon_HLT_Mu17_obj);
    mytree_->Branch("T_Muon_HLT_Mu8_obj", "std::vector<int>", &T_Muon_HLT_Mu8_obj);
    mytree_->Branch("T_Muon_HLT_Mu8_Ele17_Mu8Leg", "std::vector<int>", &T_Muon_HLT_Mu8_Ele17_Mu8Leg);
    mytree_->Branch("T_Muon_HLT_Ele8_Mu17_Mu17Leg", "std::vector<int>", &T_Muon_HLT_Ele8_Mu17_Mu17Leg);
    mytree_->Branch("T_Muon_HLT_IsoMu24", "std::vector<int>", &T_Muon_HLT_IsoMu24);
    mytree_->Branch("T_Muon_HLT_IsoMu24_2p1", "std::vector<int>", &T_Muon_HLT_IsoMu24_2p1);
}

// ------------ method called once each job just after ending the event loop  ------------
void
TriggerMuon::endJob()
{
    rootFile_->Write();
    rootFile_->Close();
}

// ------------ method called when starting to processes a run  ------------
void
TriggerMuon::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TriggerMuon::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TriggerMuon::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TriggerMuon::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

void
TriggerMuon::beginEvent()
{
    T_Muon_Eta = new std::vector<float>;
	T_Muon_Phi = new std::vector<float>;
	T_Muon_Energy = new std::vector<float>;
	T_Muon_Et = new std::vector<float>;
	T_Muon_Pt = new std::vector<float>;
	T_Muon_Px = new std::vector<float>;
	T_Muon_Py = new std::vector<float>;
	T_Muon_Pz = new std::vector<float>;
	T_Muon_Mass = new std::vector<float>;
	T_Muon_IsGlobalMuon = new std::vector<bool>;
	T_Muon_IsTrackerMuon = new std::vector<bool>;
	T_Muon_IsPFMuon = new std::vector<bool>;
	T_Muon_IsCaloMuon = new std::vector<bool>;
	T_Muon_IsStandAloneMuon = new std::vector<bool>;
	T_Muon_IsMuon = new std::vector<bool>;
	T_Muon_IsGlobalMuon_PromptTight = new std::vector<bool>;
	T_Muon_IsTrackerMuonArbitrated = new std::vector<bool>;
	T_Muon_numberOfChambers = new std::vector<int>;
	T_Muon_numberOfChambersRPC = new std::vector<int>;
	T_Muon_numberOfMatches = new std::vector<int>;
	T_Muon_numberOfMatchedStations = new std::vector<int>;
	T_Muon_charge = new std::vector<int>;
	T_Muon_TMLastStationTight = new std::vector<bool>;
	T_Muon_globalTrackChi2 = new std::vector<float>;
	T_Muon_validMuonHits = new std::vector<int>;
	T_Muon_trkKink = new std::vector<float>;
	T_Muon_trkNbOfValidTrackeHits = new std::vector<int>;
	T_Muon_trkNbOfTrackerLayers = new std::vector<int>;
	T_Muon_trkValidPixelHits = new std::vector<int>;
	T_Muon_trkError = new std::vector<float>;
	T_Muon_dB = new std::vector<float>;
	T_Muon_dzPV = new std::vector<float>;
	T_Muon_isoR03_emEt = new std::vector<float>;
	T_Muon_isoR03_hadEt = new std::vector<float>;
	T_Muon_isoR03_hoEt = new std::vector<float>;
	T_Muon_isoR03_sumPt = new std::vector<float>;
	T_Muon_isoR03_nTracks = new std::vector<int>;
	T_Muon_isoR03_nJets = new std::vector<int>;
	T_Muon_isoRingsMVA = new std::vector<float>;
    T_Muon_chargedHadronIsoR04 = new std::vector<float>;
    T_Muon_neutralHadronIsoR04 = new std::vector<float>;
    T_Muon_photonIsoR04 = new std::vector<float>;
    T_Muon_chargedHadronIsoPUR04 = new std::vector<float>;
    T_Muon_chargedHadronIsoR03 = new std::vector<float>;
    T_Muon_neutralHadronIsoR03 = new std::vector<float>;
    T_Muon_photonIsoR03 = new std::vector<float>;
    T_Muon_chargedHadronIsoPUR03 = new std::vector<float>;
	T_Muon_HLT_Mu17_TkMu8_Mu17Leg = new std::vector<int>;
	T_Muon_HLT_Mu17_TkMu8_Mu8Leg = new std::vector<int>;
	T_Muon_HLT_Mu17_Mu8_Mu17Leg = new std::vector<int>;
	T_Muon_HLT_Mu17_Mu8_Mu8Leg = new std::vector<int>;
	T_Muon_HLT_Mu17_obj = new std::vector<int>;
	T_Muon_HLT_Mu8_obj = new std::vector<int>;
	T_Muon_HLT_Mu8_Ele17_Mu8Leg = new std::vector<int>;
	T_Muon_HLT_Ele8_Mu17_Mu17Leg = new std::vector<int>;
	T_Muon_HLT_IsoMu24 = new std::vector<int>;
	T_Muon_HLT_IsoMu24_2p1 = new std::vector<int>;
}

void
TriggerMuon::endEvent()
{
    //Muons
	delete T_Muon_Eta;
	delete T_Muon_Phi;
	delete T_Muon_Energy;
	delete T_Muon_Et;
	delete T_Muon_Pt;
	delete T_Muon_Px;
	delete T_Muon_Py;
	delete T_Muon_Pz;
	delete T_Muon_Mass;
	delete T_Muon_IsGlobalMuon;
	delete T_Muon_IsTrackerMuon;
	delete T_Muon_IsPFMuon;
	delete T_Muon_IsCaloMuon;
	delete T_Muon_IsStandAloneMuon;
	delete T_Muon_IsMuon;
	delete T_Muon_IsTrackerMuonArbitrated;
	delete T_Muon_IsGlobalMuon_PromptTight;
	delete T_Muon_numberOfChambers;
	delete T_Muon_numberOfChambersRPC;
	delete T_Muon_numberOfMatches;
	delete T_Muon_numberOfMatchedStations;
	delete T_Muon_charge;
	delete T_Muon_TMLastStationTight;
	delete T_Muon_globalTrackChi2;
	delete T_Muon_validMuonHits;
	delete T_Muon_trkKink;
	delete T_Muon_trkNbOfTrackerLayers;
    delete T_Muon_trkNbOfValidTrackeHits;
	delete T_Muon_trkValidPixelHits;
	delete T_Muon_trkError;
	delete T_Muon_dB;
	delete T_Muon_dzPV;
	delete T_Muon_isoR03_emEt;
	delete T_Muon_isoR03_hadEt;
	delete T_Muon_isoR03_hoEt;
	delete T_Muon_isoR03_sumPt;
	delete T_Muon_isoR03_nTracks;
	delete T_Muon_isoR03_nJets;
	delete T_Muon_isoRingsMVA;
	delete T_Muon_chargedHadronIsoR04;
	delete T_Muon_neutralHadronIsoR04;
	delete T_Muon_photonIsoR04;
	delete T_Muon_chargedHadronIsoPUR04;
    delete T_Muon_chargedHadronIsoR03;
	delete T_Muon_neutralHadronIsoR03;
	delete T_Muon_photonIsoR03;
	delete T_Muon_chargedHadronIsoPUR03;
	delete T_Muon_HLT_Mu17_TkMu8_Mu17Leg;
	delete T_Muon_HLT_Mu17_TkMu8_Mu8Leg;
	delete T_Muon_HLT_Mu17_Mu8_Mu17Leg;
	delete T_Muon_HLT_Mu17_Mu8_Mu8Leg;
	delete T_Muon_HLT_Mu17_obj;
	delete T_Muon_HLT_Mu8_obj;
	delete T_Muon_HLT_Mu8_Ele17_Mu8Leg;
	delete T_Muon_HLT_Ele8_Mu17_Mu17Leg;
	delete T_Muon_HLT_IsoMu24;
	delete T_Muon_HLT_IsoMu24_2p1;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerMuon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

float TriggerMuon::deltaR(float phi1, float phi2, float eta1, float eta2)
{
    float dphi=deltaPhi(phi1,phi2);
    float deta=fabs(eta1-eta2);
    float dr = sqrt(dphi*dphi+ deta*deta);
    return dr;
}

float TriggerMuon::deltaPhi(float phi1, float phi2)
{
    float dphi;
    if(phi1<0) phi1+=2*TMath::Pi();
    if(phi2<0) phi2+=2*TMath::Pi();
    dphi=fabs(phi1-phi2);
    if(dphi>2*TMath::Pi()) dphi-=2*TMath::Pi();
    if(dphi>TMath::Pi()) dphi=2*TMath::Pi()-dphi;
    return dphi;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerMuon);
