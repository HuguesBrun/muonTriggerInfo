import FWCore.ParameterSet.Config as cms

savePatInTree=True;

process = cms.Process("EX")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

# try to add the PAT PF sequences in the analyser...
## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
jetAlgo="AK5"
if savePatInTree: usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=True, postfix=postfix)






process.GlobalTag.globaltag = 'START53_V7A::All'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.MessageLogger.cerr.FwkReport.reportEvery = 1
#
# the MC global Tag : START53_V7A
# the RERECO of 2012 data (Jan22 reRECO) for run ABC FT_53_V21_AN3
# input
#

process.source = cms.Source(
                            "PoolSource",
                            fileNames = cms.untracked.vstring(
                                                              #      'file:skimedFile/pickevents.root'
                                                              'file:/sps/cms/hbrun/CMSSW_5_3_10_forNewSims/src/files/runDepMC/MCDY_runDep_1.root'
                                                             ),
                            secondaryFileNames = cms.untracked.vstring(),
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )



process.triggerMuon = cms.EDAnalyzer('TriggerMuon',
                                     isMC               = cms.bool(True),
                                     doPFPATmatching    = cms.bool(False),
                                     muonProducer 		= cms.VInputTag(cms.InputTag("muons")),
                                     TriggerResults          = cms.InputTag("TriggerResults", "", "HLT"),
                                     HLTTriggerSummaryAOD    = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                                     primaryVertexInputTag = cms.InputTag("goodOfflinePrimaryVertices"),
                                     outputFile		        = cms.string("muonTriggerTree.root"),
)

process.goodVertexFilter = cms.EDFilter("VertexSelector",
                                        src = cms.InputTag("offlinePrimaryVertices"),
                                        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"),
                                        filter = cms.bool(True),
                                        )
process.noScraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )

process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
                                                  filter = cms.bool(False),
                                                  src = cms.InputTag("offlinePrimaryVertices"),
                                                  filterParams = cms.PSet(
                                                                          maxZ = cms.double(24.0),
                                                                          minNdof = cms.double(4.0),
                                                                          maxRho = cms.double(2.0)
                                                                          )
                                                  )


process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")

process.triggerResultsFilter.triggerConditions = cms.vstring( 'HLT_Mu17_Mu8_v*','HLT_Mu17_TkMu8_v*','HLT_Mu17_v*','HLT_Mu24_eta2p1_v*','HLT_Mu17_TkMu8_NoDZ_v*','HLT_Mu13_Mu8_NoDZ_v*')
#process.triggerResultsFilter.triggerConditions = cms.vstring('HLT_Mu13_Mu8_NoDZ_v*','HLT_Mu17_TkMu8_NoDZ_v*')
process.triggerResultsFilter.l1tResults = ''
process.triggerResultsFilter.throw = False
process.triggerResultsFilter.hltResults = cms.InputTag( "TriggerResults", "", "HLT" )

if savePatInTree:
    # top projections in PF2PAT:
    getattr(process,"pfNoPileUp"+postfix).enable = True
    getattr(process,"pfNoMuon"+postfix).enable = True
    getattr(process,"pfNoElectron"+postfix).enable = True
    
    # tau considered as jet
    getattr(process,"pfNoTau"+postfix).enable = False
    getattr(process,"pfNoJet"+postfix).enable = False
    
    # verbose flags for the PF2PAT modules
    getattr(process,"pfNoMuon"+postfix).verbose = False
    
    #ask the analyzer to
    getattr(process,"triggerMuon").doPFPATmatching = True



#all pfMuons considered as isolated
process.pfIsolatedMuonsPFlow.isolationCut = cms.double(99999.)
process.pfSelectedMuonsPFlow.cut = cms.string('pt>5')
process.pfMuonsFromVertexPFlow.d0Cut = cms.double(99)
process.pfMuonsFromVertexPFlow.dzCut = cms.double(99)
process.pfMuonsFromVertexPFlow.d0SigCut = cms.double(9999999.)
process.pfMuonsFromVertexPFlow.dzSigCut = cms.double(9999999.)
#all pfElectrons considered as isolated
process.pfIsolatedElectronsPFlow.combinedIsolationCut = cms.double(99999.)

adaptPFIsoMuons(process,process.pfIsolatedMuonsPFlow,"PFlow", "03")


if savePatInTree:
    #sequence with PF
    process.p = cms.Path(process.goodVertexFilter * process.noScraping * process.goodOfflinePrimaryVertices * getattr(process,"patPF2PATSequence"+postfix)*process.triggerMuon)
else:
    #sequence without PF
    process.p = cms.Path(process.goodVertexFilter * process.noScraping * process.triggerMuon)

from CommonTools.ParticleFlow.PF2PAT_EventContent_cff import PF2PATStudiesEventContent
process.out.outputCommands =  PF2PATStudiesEventContent.outputCommands


