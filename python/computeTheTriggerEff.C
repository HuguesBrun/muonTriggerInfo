
TChain *chain = new TChain("eventsTree");


//variables in the tree
std::vector<float> *T_Muon_Pt;
std::vector<float> *T_Muon_Eta;
std::vector<float> *T_Gen_Muon_Energy;

std::vector<bool> *T_Muon_IsGlobalMuon;
std::vector<bool> *T_Muon_IsTrackerMuon;
std::vector<bool> *T_Muon_IsPFMuon;
std::vector<int> *T_Muon_isMatchWithPAT;

std::vector<int>*T_Gen_Muon_PDGid;
std::vector<int>*T_Gen_Muon_MotherID;

std::vector<bool> *T_Muon_HLT_Mu17_Mu8_Mu17Leg;
std::vector<bool> *T_Muon_HLT_Mu17_Mu8_Mu8Leg;
std::vector<bool> *T_Muon_HLT_Mu17_TkMu8_Mu17Leg;
std::vector<bool> *T_Muon_HLT_Mu17_TkMu8_Mu8Leg;


int T_Event_HLT_Mu17_Mu8;
int T_Event_HLT_Mu17_TkMu8;
// output file :
TFile *outputFile = new TFile("outputFile.root","RECREATE");


computeTheTriggerEff(){
    chain->Add("muonTriggerTree.root");
    chain->SetBranchAddress("T_Muon_Pt",&T_Muon_Pt);
    chain->SetBranchAddress("T_Muon_Eta",&T_Muon_Eta);
    
    chain->SetBranchAddress("T_Muon_IsGlobalMuon",&T_Muon_IsGlobalMuon);
    chain->SetBranchAddress("T_Muon_IsTrackerMuon",&T_Muon_IsTrackerMuon);
    chain->SetBranchAddress("T_Muon_IsPFMuon",&T_Muon_IsPFMuon);
    chain->SetBranchAddress("T_Muon_isMatchWithPAT",&T_Muon_isMatchWithPAT);
    
    chain->SetBranchAddress("T_Gen_Muon_Energy",&T_Gen_Muon_Energy);
    chain->SetBranchAddress("T_Gen_Muon_PDGid",&T_Gen_Muon_PDGid);
    chain->SetBranchAddress("T_Gen_Muon_MotherID",&T_Gen_Muon_MotherID);

    chain->SetBranchAddress("T_Muon_HLT_Mu17_Mu8_Mu17Leg",&T_Muon_HLT_Mu17_Mu8_Mu17Leg);
    chain->SetBranchAddress("T_Muon_HLT_Mu17_Mu8_Mu8Leg",&T_Muon_HLT_Mu17_Mu8_Mu8Leg);
    chain->SetBranchAddress("T_Muon_HLT_Mu17_TkMu8_Mu17Leg",&T_Muon_HLT_Mu17_TkMu8_Mu17Leg);
    chain->SetBranchAddress("T_Muon_HLT_Mu17_TkMu8_Mu8Leg",&T_Muon_HLT_Mu17_TkMu8_Mu8Leg);

    
    chain->SetBranchAddress("T_Event_HLT_Mu17_Mu8",&T_Event_HLT_Mu17_Mu8);
    chain->SetBranchAddress("T_Event_HLT_Mu17_TkMu8",&T_Event_HLT_Mu17_TkMu8);

    
    //histos :
    TH1F *h_Mu17_Mu8 = new TH1F("h_Mu17_Mu8","",2,0,2);
    TH1F *h_Mu17_TkMu8 = new TH1F("h_Mu17_TkMu8","",2,0,2);
    
    TH1F *h_Mu17_Mu8_withMatch = new TH1F("h_Mu17_Mu8_withMatch","",2,0,2);
    TH1F *h_Mu17_TkMu8_withMatch = new TH1F("h_Mu17_TkMu8_withMatch","",2,0,2);
    
    int nbEntries = chain->GetEntries();
    
    cout << "nbEntries=" << nbEntries << endl;
    int nbEventPassingGenCut=0;
    int nbEventWith2SelectedMuon=0;
    int nbEventPassingMu17Mu8=0;
    int nbEventPassingMu17TkMu8=0;
    for (int i = 0 ; i<nbEntries ; i++){
        if (i%1000==0) cout << "event number " << i << endl;
  
        chain->GetEntry(i);
   
    
        
        // check if the event is a Z->mumu event
        int nbGenMuons = T_Gen_Muon_Energy->size();
        int nbPairMuonsFromZ = 0;
        for (int j = 0 ; j < nbGenMuons ; j++){
            if (!(fabs(T_Gen_Muon_PDGid->at(j))==13 && T_Gen_Muon_MotherID->at(j)==23)) continue; //gen muons from a Z
            for (int k = (j+1) ; k < nbGenMuons ; k++){
                if (!(fabs(T_Gen_Muon_PDGid->at(k))==13 && T_Gen_Muon_MotherID->at(k)==23)) continue; //gen muons from a Z
                if ((T_Gen_Muon_PDGid->at(k)*T_Gen_Muon_PDGid->at(j))>0) continue; // two muon of opposite sign :)
                nbPairMuonsFromZ++;

            }
        }
       // cout << "nb nbZ genmuons=" << nbPairMuonsFromZ << endl;
        //if (!(nbPairMuonsFromZ>=1)) continue; // if not Z->mu mu gen event then go to next event
        nbEventPassingGenCut++;
        
        
        if (T_Muon_Pt->size()<2) continue; //need at least 2 muons in the event...
        int nbMuons = T_Muon_Pt->size();
        int nbLooseMuons = 0;
         std::vector<int> refLooseMuons = 0;
        for (int j = 0 ; j < nbMuons ; j++){
            if (!((T_Muon_Pt->at(j)>20)&&(fabs(T_Muon_Eta->at(j))<2.4))) continue;
            if (!((T_Muon_IsGlobalMuon->at(j)||T_Muon_IsTrackerMuon->at(j)))) continue;
            refLooseMuons.push_back(j);
            nbLooseMuons++;
        }
        //cout << "nb Loose muons=" << nbLooseMuons << endl;
        if (!(nbLooseMuons>=2)) continue;
        nbEventWith2SelectedMuon++;
        // fill the trigger paths histos ! 
        h_Mu17_Mu8->Fill(T_Event_HLT_Mu17_Mu8);
        h_Mu17_TkMu8->Fill(T_Event_HLT_Mu17_TkMu8);
        nbEventPassingMu17Mu8+=T_Event_HLT_Mu17_Mu8;
        nbEventPassingMu17TkMu8+=T_Event_HLT_Mu17_TkMu8;
        
        
        //fill the trigger paths after match:
        for (int j = 0 ; j < refLooseMuons.size(); j++){
            for (int k = j+1 ; k< refLooseMuons.size(); k++){
                    int passMu17Mu8andMatched = (T_Event_HLT_Mu17_Mu8&&(
                                                                        (T_Muon_HLT_Mu17_Mu8_Mu17Leg->at(refLooseMuons.at(j))&&T_Muon_HLT_Mu17_Mu8_Mu8Leg->at(refLooseMuons.at(k)))
                                                                        ||
                                                                        (T_Muon_HLT_Mu17_Mu8_Mu17Leg->at(refLooseMuons.at(k))&&T_Muon_HLT_Mu17_Mu8_Mu8Leg->at(refLooseMuons.at(j)))
                                                                        )
                                                 );
                    h_Mu17_Mu8_withMatch->Fill(passMu17Mu8andMatched);
        
                    int passMu17TkMu8andMatched = (T_Event_HLT_Mu17_TkMu8&&(
                                                                            (T_Muon_HLT_Mu17_TkMu8_Mu17Leg->at(refLooseMuons.at(j))&&T_Muon_HLT_Mu17_TkMu8_Mu8Leg->at(refLooseMuons.at(k)))
                                                                            ||
                                                                            (T_Muon_HLT_Mu17_TkMu8_Mu17Leg->at(refLooseMuons.at(k))&&T_Muon_HLT_Mu17_TkMu8_Mu8Leg->at(refLooseMuons.at(j)))
                                                                            )
                                                   );
                    h_Mu17_TkMu8_withMatch->Fill(passMu17TkMu8andMatched);
            }
        }
    }
    
    
    cout << "Mu17_Mu8=" << h_Mu17_Mu8->GetMean() << " +- " << h_Mu17_Mu8->GetMeanError() << endl;
    cout << "Mu17_TkMu8=" << h_Mu17_TkMu8->GetMean() << " +- " << h_Mu17_TkMu8->GetMeanError() << endl;
    
    cout << "after matching:" << endl;
    cout << "Mu17_Mu8=" << h_Mu17_Mu8_withMatch->GetMean() << " +- " << h_Mu17_Mu8_withMatch->GetMeanError() << endl;
    cout << "Mu17_TkMu8=" << h_Mu17_TkMu8_withMatch->GetMean() << " +- " << h_Mu17_TkMu8_withMatch->GetMeanError() << endl;
    
    
    cout << "nbEventAfterGEN=" << nbEventPassingGenCut << endl;
    cout << "nbEvent2selectedMuons=" << nbEventWith2SelectedMuon << endl;
    cout << "nbEventPassing Mu17Mu8=" << nbEventPassingMu17Mu8 << endl;
    cout << "nbEventPassing Mu17TkMu8=" << nbEventPassingMu17TkMu8 << endl;
    
    outputFile->Write();
    outputFile->Close();

}




