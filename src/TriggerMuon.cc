#include "../interface/TriggerMuon.h"



TriggerMuon::TriggerMuon(const edm::ParameterSet& iConfig)

{
    
    isMC_                   = iConfig.getParameter<bool>("isMC");
    doPFPATmatching_        = iConfig.getParameter<bool>("doPFPATmatching");
    
    triggerResultsLabel_    = iConfig.getParameter<edm::InputTag>("TriggerResults");
    triggerSummaryLabel_    = iConfig.getParameter<edm::InputTag>("HLTTriggerSummaryAOD");
    muonProducers_			= iConfig.getParameter<vtag>("muonProducer");
    
    outputFile_   = iConfig.getParameter<std::string>("outputFile");
    rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");
    
    
    HLT_name.push_back("HLT_Mu17_Mu8_v");//0
    HLT_name.push_back("HLT_Mu17_TkMu8_v");//1
    
    
    HLT_triggerObjects.push_back("hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17");// 0 -> DoubleMu17Mu8_Mu17
    HLT_triggerObjects.push_back("hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8");// 1-> DoubleMu17Mu8_Mu8
    HLT_triggerObjects.push_back("hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17");// 2-> DoubleMu17TkMu8_Mu17leg
    HLT_triggerObjects.push_back("hltDiMuonGlbFiltered17TrkFiltered8");// 3-> DoubleMu17TkMu8_TkMu8leg

    

    

    

    
 
    

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
    
    edm::Handle <reco::GenParticleCollection> genParticles;
    
    // read the PAT muons
    edm::Handle<pat::MuonCollection > patMuons;
    if (doPFPATmatching_) iEvent.getByLabel( "selectedPatMuonsPFlow", patMuons );
    

    
    if (isMC_){// get the gen infos
        edm::Handle<GenEventInfoProduct> genEvent;
        iEvent.getByLabel("generator", genEvent);
        iEvent.getByLabel( "genParticles", genParticles );
    }
    
    
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
    
    //cout << "on va filler le HLT" << endl;
    
    T_Event_HLT_Mu17_Mu8 =         triggerResults->accept(theBitCorr[0]);
    T_Event_HLT_Mu17_TkMu8 =         triggerResults->accept(theBitCorr[1]);
    

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
    //cout << "event="  << iEvent.id().run() << endl;
    T_Event_RunNumber = iEvent.id().run();
    T_Event_EventNumber = iEvent.id().event();
    T_Event_LuminosityBlock = iEvent.id().luminosityBlock();

    
    int nbMuons = patMuons->size();
    for (int k = 0 ; k < nbMuons ; k++){// loop on the muons in the event
        const pat::Muon* muon = &((*patMuons)[k]);
      //  cout << "muon PT =" << muon->pt() << endl; 
        T_Muon_Eta->push_back(muon->eta());
        T_Muon_Phi->push_back(muon->phi());
        T_Muon_IsGlobalMuon->push_back(muon->isGlobalMuon());
        T_Muon_IsPFMuon->push_back(muon->isPFMuon());
        T_Muon_IsTrackerMuon->push_back(muon->isTrackerMuon());
        T_Muon_Energy->push_back(muon->energy());
        T_Muon_Et->push_back(muon->et());
        T_Muon_Pt->push_back(muon->pt());
        T_Muon_Px->push_back(muon->px());
        T_Muon_Py->push_back(muon->py());
        T_Muon_Pz->push_back(muon->pz());
        T_Muon_Mass->push_back(muon->mass());
        
    
        
        int pass_HLT_Mu17_Mu8_Mu17Leg = 0;
        int pass_HLT_Mu17_Mu8_Mu8Leg = 0;
        int pass_HLT_Mu17_TkMu8_Mu17Leg = 0;
        int pass_HLT_Mu17_TkMu8_Mu8Leg = 0;


        
        for (size_t t = 0 ; t < selectedObjects.size() ; t++){
            float HLTdeltaR = deltaR(muon->phi(), selectedObjects[t].phi(), muon->eta(), selectedObjects[t].eta());
            float relatDeltaPt = fabs(selectedObjects[t].pt()-muon->pt())/muon->pt();
            if (HLTdeltaR < 0.3){
                if (theHLTcorr[t] == 0)  pass_HLT_Mu17_Mu8_Mu17Leg = 1;
                if (theHLTcorr[t] == 1)  pass_HLT_Mu17_Mu8_Mu8Leg = 1;
                if (theHLTcorr[t] == 2)  pass_HLT_Mu17_TkMu8_Mu17Leg = 1;
                if (theHLTcorr[t] == 3)  pass_HLT_Mu17_TkMu8_Mu8Leg = 1;
            }
        }
        T_Muon_HLT_Mu17_Mu8_Mu17Leg->push_back(pass_HLT_Mu17_Mu8_Mu17Leg);
        T_Muon_HLT_Mu17_Mu8_Mu8Leg->push_back(pass_HLT_Mu17_Mu8_Mu8Leg);
        T_Muon_HLT_Mu17_TkMu8_Mu17Leg->push_back(pass_HLT_Mu17_TkMu8_Mu17Leg);
        T_Muon_HLT_Mu17_TkMu8_Mu8Leg->push_back(pass_HLT_Mu17_TkMu8_Mu8Leg);
        
        
        /// try the matching with the PAT PF muons
        bool foundAElectronPfMatch = false;
        if (doPFPATmatching_){
            for( size_t iMuon = 0; iMuon < patMuons->size(); ++iMuon ) {
                float drPF = deltaR(muon->phi(), patMuons->at( iMuon ).phi(), muon->eta(),patMuons->at( iMuon ).eta());
                if (drPF>0.1) continue;
                T_Muon_isMatchWithPAT->push_back(1);
                T_Muon_PATpt->push_back(patMuons->at( iMuon ).pt());
                T_Muon_PATeta->push_back(patMuons->at( iMuon ).eta());
                T_Muon_PATphi->push_back(patMuons->at( iMuon ).phi());
                T_Muon_PATenergy->push_back(patMuons->at( iMuon ).energy());
                T_Muon_PATpx->push_back(patMuons->at( iMuon ).px());
                T_Muon_PATpy->push_back(patMuons->at( iMuon ).py());
                T_Muon_PATpz->push_back(patMuons->at( iMuon ).pz());
            }
        }
        if (!(foundAElectronPfMatch)) T_Muon_isMatchWithPAT->push_back(0);

    }
    // now save the gen particles 
    if (isMC_){
        // now do the matching with the gen particles :
        int nbOfGen = genParticles->size();
        float minDiff= 100;
        int iteDiff = -1000;
        int motherID = 0;
        for (int m = 0 ; m < nbOfGen ; m++){//loop on the gen particles
            const reco::GenParticle & p = (*genParticles)[m];
            if (!(p.status()==1)) continue;
            if (fabs(p.pdgId())!=13) continue; // found a stable gen muon
            const reco::Candidate * theLocalCandidate = &p;
            bool hasMother = (theLocalCandidate->numberOfMothers()>0);
            const reco::Candidate * theMother;
            while (hasMother) {//check if the muon has a Z in its mother particles
                theMother = theLocalCandidate->mother();
                theLocalCandidate = theMother;
                hasMother = (theLocalCandidate->numberOfMothers()>0);
                motherID = theMother->pdgId();
                if ((theMother->pdgId()==23)||(theMother->pdgId()==22)) break;
            }
            //save the gen particle
            T_Gen_Muon_Px->push_back(p.px());
            T_Gen_Muon_Py->push_back(p.py());
            T_Gen_Muon_Pz->push_back(p.pz());
            T_Gen_Muon_Energy->push_back(p.energy());
            T_Gen_Muon_MCpart->push_back(1);
            T_Gen_Muon_PDGid->push_back(p.pdgId());
            T_Gen_Muon_status->push_back(p.status());
            if (p.numberOfMothers()>0) T_Gen_Muon_MotherID->push_back(motherID);
            else T_Gen_Muon_MotherID->push_back(-1);

        }
    }
    //save now HLT objects !
    for (size_t t = 0 ; t < selectedObjects.size() ; t++){
        T_TrigObj_Pt->push_back(selectedObjects[t].pt());
        T_TrigObj_Eta->push_back(selectedObjects[t].eta());
        T_TrigObj_Phi->push_back(selectedObjects[t].phi());
        T_TrigObj_FilterIndex->push_back(theHLTcorr[t]);
        
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
   /* mytree_->Branch("T_Event_HLT_IsoMu24",&T_Event_HLT_IsoMu24,"T_Event_HLT_IsoMu24/I");
    mytree_->Branch("T_Event_HLT_IsoMu24_2p1",&T_Event_HLT_IsoMu24_2p1,"T_Event_HLT_IsoMu24_2p1/I");
    mytree_->Branch("T_Event_HLT_Mu17_TkMu8_NoDZ",&T_Event_HLT_Mu17_TkMu8_NoDZ,"T_Event_HLT_Mu17_TkMu8_NoDZ/I");
    mytree_->Branch("T_Event_HLT_Mu13_Mu8_NoDZ",&T_Event_HLT_Mu13_Mu8_NoDZ,"T_Event_HLT_Mu13_Mu8_NoDZ/I");
    mytree_->Branch("T_Event_HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned",&T_Event_HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned,"T_Event_HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned/I");
    mytree_->Branch("T_Event_HLT_DiCentralPFJet30_PFMET80_BTagCSV07",&T_Event_HLT_DiCentralPFJet30_PFMET80_BTagCSV07,"T_Event_HLT_DiCentralPFJet30_PFMET80_BTagCSV07/I");
    mytree_->Branch("T_Event_HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80",&T_Event_HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80,"T_Event_HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80/I");
    mytree_->Branch("T_Event_HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets",&T_Event_HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets,"T_Event_HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets/I");
    mytree_->Branch("T_Event_HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets",&T_Event_HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets,"T_Event_HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets/I");
    mytree_->Branch("T_Event_HLT_L1ETM30",&T_Event_HLT_L1ETM30,"T_Event_HLT_L1ETM30/I");
    mytree_->Branch("T_Event_HLT_MET120_HBHENoiseCleaned",&T_Event_HLT_MET120_HBHENoiseCleaned,"T_Event_HLT_MET120_HBHENoiseCleaned/I");
    mytree_->Branch("T_Event_HLT_MET200_HBHENoiseCleaned",&T_Event_HLT_MET200_HBHENoiseCleaned,"T_Event_HLT_MET200_HBHENoiseCleaned/I");
    mytree_->Branch("T_Event_HLT_MET200",&T_Event_HLT_MET200,"T_Event_HLT_MET200/I");
    mytree_->Branch("T_Event_HLT_MET300_HBHENoiseCleaned",&T_Event_HLT_MET300_HBHENoiseCleaned,"T_Event_HLT_MET300_HBHENoiseCleaned/I");
    mytree_->Branch("T_Event_HLT_MET300",&T_Event_HLT_MET300,"T_Event_HLT_MET300/I");
    mytree_->Branch("T_Event_HLT_MET400_HBHENoiseCleaned",&T_Event_HLT_MET400_HBHENoiseCleaned,"T_Event_HLT_MET400_HBHENoiseCleaned/I");
    mytree_->Branch("T_Event_HLT_MET400",&T_Event_HLT_MET400,"T_Event_HLT_MET400/I");
    mytree_->Branch("T_Event_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95",&T_Event_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95,"T_Event_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95/I");
    mytree_->Branch("T_Event_HLT_PFMET150",&T_Event_HLT_PFMET150,"T_Event_HLT_PFMET150/I");
    mytree_->Branch("T_Event_HLT_PFMET180",&T_Event_HLT_PFMET180,"T_Event_HLT_PFMET180/I");*/

    
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
    mytree_->Branch("T_Muon_chargedHadronIsoR04", "std::vector<float>", &T_Muon_chargedHadronIsoR04);
    mytree_->Branch("T_Muon_neutralHadronIsoR04", "std::vector<float>", &T_Muon_neutralHadronIsoR04);
    mytree_->Branch("T_Muon_photonIsoR04", "std::vector<float>", &T_Muon_photonIsoR04);
    mytree_->Branch("T_Muon_chargedHadronIsoPUR04", "std::vector<float>", &T_Muon_chargedHadronIsoPUR04);
    mytree_->Branch("T_Muon_chargedHadronIsoR03", "std::vector<float>", &T_Muon_chargedHadronIsoR03);
    mytree_->Branch("T_Muon_neutralHadronIsoR03", "std::vector<float>", &T_Muon_neutralHadronIsoR03);
    mytree_->Branch("T_Muon_photonIsoR03", "std::vector<float>", &T_Muon_photonIsoR03);
    mytree_->Branch("T_Muon_chargedHadronIsoPUR03", "std::vector<float>", &T_Muon_chargedHadronIsoPUR03);
    
    if (doPFPATmatching_){
        mytree_->Branch("T_Muon_isMatchWithPAT", "std::vector<int>", &T_Muon_isMatchWithPAT);
        mytree_->Branch("T_Muon_PATpt", "std::vector<float>", &T_Muon_PATpt);
        mytree_->Branch("T_Muon_PATeta", "std::vector<float>", &T_Muon_PATeta);
        mytree_->Branch("T_Muon_PATphi", "std::vector<float>", &T_Muon_PATphi);
        mytree_->Branch("T_Muon_PATenergy", "std::vector<float>", &T_Muon_PATenergy);
        mytree_->Branch("T_Muon_PATpx", "std::vector<float>", &T_Muon_PATpx);
        mytree_->Branch("T_Muon_PATpy", "std::vector<float>", &T_Muon_PATpy);
        mytree_->Branch("T_Muon_PATpz", "std::vector<float>", &T_Muon_PATpz);
    }

    
    mytree_->Branch("T_Muon_HLT_Mu17_Mu8_Mu17Leg", "std::vector<bool>", &T_Muon_HLT_Mu17_Mu8_Mu17Leg);
    mytree_->Branch("T_Muon_HLT_Mu17_Mu8_Mu8Leg", "std::vector<bool>", &T_Muon_HLT_Mu17_Mu8_Mu8Leg);
    mytree_->Branch("T_Muon_HLT_Mu17_TkMu8_Mu17Leg", "std::vector<bool>", &T_Muon_HLT_Mu17_TkMu8_Mu17Leg);
    mytree_->Branch("T_Muon_HLT_Mu17_TkMu8_Mu8Leg", "std::vector<bool>", &T_Muon_HLT_Mu17_TkMu8_Mu8Leg);

    

    mytree_->Branch("T_Gen_Muon_Px", "std::vector<float>", &T_Gen_Muon_Px);
    mytree_->Branch("T_Gen_Muon_Py", "std::vector<float>", &T_Gen_Muon_Py);
    mytree_->Branch("T_Gen_Muon_Pz", "std::vector<float>", &T_Gen_Muon_Pz);
    mytree_->Branch("T_Gen_Muon_Energy", "std::vector<float>", &T_Gen_Muon_Energy);
    mytree_->Branch("T_Gen_Muon_MCpart", "std::vector<int>", &T_Gen_Muon_MCpart);
    mytree_->Branch("T_Gen_Muon_PDGid", "std::vector<int>", &T_Gen_Muon_PDGid);
    mytree_->Branch("T_Gen_Muon_status", "std::vector<int>", &T_Gen_Muon_status);
    mytree_->Branch("T_Gen_Muon_MotherID", "std::vector<int>", &T_Gen_Muon_MotherID);
    mytree_->Branch("T_Gen_Muon_deltaR", "std::vector<float>", &T_Gen_Muon_deltaR);

    
    
    mytree_->Branch("T_TrigObj_Pt", "std::vector<float>", &T_TrigObj_Pt);
    mytree_->Branch("T_TrigObj_Eta", "std::vector<float>", &T_TrigObj_Eta);
    mytree_->Branch("T_TrigObj_Phi", "std::vector<float>", &T_TrigObj_Phi);
    mytree_->Branch("T_TrigObj_FilterIndex", "std::vector<int>", &T_TrigObj_FilterIndex);

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

    T_Muon_chargedHadronIsoR04 = new std::vector<float>;
    T_Muon_neutralHadronIsoR04 = new std::vector<float>;
    T_Muon_photonIsoR04 = new std::vector<float>;
    T_Muon_chargedHadronIsoPUR04 = new std::vector<float>;
    T_Muon_chargedHadronIsoR03 = new std::vector<float>;
    T_Muon_neutralHadronIsoR03 = new std::vector<float>;
    T_Muon_photonIsoR03 = new std::vector<float>;
    T_Muon_chargedHadronIsoPUR03 = new std::vector<float>;
    
    T_Muon_isMatchWithPAT = new std::vector<int>;
    T_Muon_PATpt = new std::vector<float>;
    T_Muon_PATeta = new std::vector<float>;
    T_Muon_PATphi = new std::vector<float>;
    T_Muon_PATenergy = new std::vector<float>;
    T_Muon_PATpx = new std::vector<float>;
    T_Muon_PATpy = new std::vector<float>;
    T_Muon_PATpz = new std::vector<float>;

    
    
    T_Muon_HLT_Mu17_Mu8_Mu17Leg = new std::vector<bool>;
    T_Muon_HLT_Mu17_Mu8_Mu8Leg = new std::vector<bool>;
    T_Muon_HLT_Mu17_TkMu8_Mu17Leg = new std::vector<bool>;
    T_Muon_HLT_Mu17_TkMu8_Mu8Leg = new std::vector<bool>;
    
    T_Gen_Muon_Px = new std::vector<float>;
    T_Gen_Muon_Py = new std::vector<float>;
    T_Gen_Muon_Pz = new std::vector<float>;
    T_Gen_Muon_Energy = new std::vector<float>;
    T_Gen_Muon_MCpart = new std::vector<int>;
    T_Gen_Muon_PDGid = new std::vector<int>;
    T_Gen_Muon_status = new std::vector<int>;
    T_Gen_Muon_MotherID = new std::vector<int>;
    T_Gen_Muon_deltaR = new std::vector<float>;


    T_TrigObj_Pt = new std::vector<float>;
    T_TrigObj_Eta = new std::vector<float>;
    T_TrigObj_Phi = new std::vector<float>;
    T_TrigObj_FilterIndex = new std::vector<int>;

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
	delete T_Muon_chargedHadronIsoR04;
	delete T_Muon_neutralHadronIsoR04;
	delete T_Muon_photonIsoR04;
	delete T_Muon_chargedHadronIsoPUR04;
    delete T_Muon_chargedHadronIsoR03;
	delete T_Muon_neutralHadronIsoR03;
	delete T_Muon_photonIsoR03;
	delete T_Muon_chargedHadronIsoPUR03;

    delete T_Muon_isMatchWithPAT;
    delete T_Muon_PATpt;
    delete T_Muon_PATeta;
    delete T_Muon_PATphi;
    delete T_Muon_PATenergy;
    delete T_Muon_PATpx;
    delete T_Muon_PATpy;
    delete T_Muon_PATpz;

    
    
    delete T_Muon_HLT_Mu17_Mu8_Mu17Leg;
    delete T_Muon_HLT_Mu17_Mu8_Mu8Leg;
    delete T_Muon_HLT_Mu17_TkMu8_Mu17Leg;
    delete T_Muon_HLT_Mu17_TkMu8_Mu8Leg;
    
	delete T_Gen_Muon_Px;
	delete T_Gen_Muon_Py;
	delete T_Gen_Muon_Pz;
	delete T_Gen_Muon_Energy;
	delete T_Gen_Muon_MCpart;
	delete T_Gen_Muon_PDGid;
	delete T_Gen_Muon_status;
	delete T_Gen_Muon_MotherID;
	delete T_Gen_Muon_deltaR;
    
    
	delete T_TrigObj_Pt;
	delete T_TrigObj_Eta;
	delete T_TrigObj_Phi;
	delete T_TrigObj_FilterIndex;
    
    

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
