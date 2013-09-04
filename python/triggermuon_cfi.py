import FWCore.ParameterSet.Config as cms

process = cms.Process("EX")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START53_V7A::All'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.MessageLogger.cerr.FwkReport.reportEvery = 1
#
# the MC global Tag : START53_V7A
# the RERECO of 2012 data (Jan22 reRECO) for run ABC FT_53_V21_AN3
# input
#

process.source = cms.Source(
                            "PoolSource",
                            fileNames = cms.untracked.vstring(
								'file:/sps/cms/hbrun/CMSSW_5_3_7_myCode/src/files/MC_DY_1.root',
								'file:/sps/cms/hbrun/CMSSW_5_3_7_myCode/src/files/MC_DY_2.root',
								'file:/sps/cms/hbrun/CMSSW_5_3_7_myCode/src/files/MC_DY_3.root',
								'file:/sps/cms/hbrun/CMSSW_5_3_7_myCode/src/files/MC_DY_4.root',
								'file:/sps/cms/hbrun/CMSSW_5_3_7_myCode/src/files/MC_DY_5.root',
#								'file:/sps/cms/hbrun/CMSSW_5_3_10_forNewSims/src/fileRunD/reco_file.root'
                                                              #								'file:/sps/cms/hbrun/CMSSW_5_3_7_myCode/src/dataFile_runA/theFile.root'
#                                                              'file:/sps/cms/hbrun/CMSSW_5_3_10_forNewSims/src/files/runDepMC/MCDY_runDep_1.root',
 #                                                             'file:/sps/cms/hbrun/CMSSW_5_3_10_forNewSims/src/files/runDepMC/MCDY_runDep_2.root',
  #                                                            'file:/sps/cms/hbrun/CMSSW_5_3_10_forNewSims/src/files/runDepMC/MCDY_runDep_3.root',
   #                                                          'file:/sps/cms/hbrun/CMSSW_5_3_10_forNewSims/src/files/runDepMC/MCDY_runDep_4.root',
    #                                                          'file:/sps/cms/hbrun/CMSSW_5_3_10_forNewSims/src/files/runDepMC/MCDY_runDep_5.root',
                                                              ),
                            secondaryFileNames = cms.untracked.vstring(),
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )



process.triggerMuon = cms.EDAnalyzer('TriggerMuon',
                                     isMC               = cms.bool(True),
                                     muonProducer 		= cms.VInputTag(cms.InputTag("muons")),
                                     TriggerResults          = cms.InputTag("TriggerResults", "", "HLT"),
                                     HLTTriggerSummaryAOD    = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
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


process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")

process.triggerResultsFilter.triggerConditions = cms.vstring( 'HLT_Mu17_Mu8_v*','HLT_Mu17_TkMu8_v*','HLT_Mu17_v*','HLT_Mu24_eta2p1_v*','HLT_Mu17_TkMu8_NoDZ_v*','HLT_Mu13_Mu8_NoDZ_v*')
#process.triggerResultsFilter.triggerConditions = cms.vstring('HLT_Mu13_Mu8_NoDZ_v*','HLT_Mu17_TkMu8_NoDZ_v*')
process.triggerResultsFilter.l1tResults = ''
process.triggerResultsFilter.throw = False
process.triggerResultsFilter.hltResults = cms.InputTag( "TriggerResults", "", "HLT" )

#process.p = cms.Path(process.triggerResultsFilter+process.goodVertexFilter + process.noScraping+process.triggerMuon)
process.p = cms.Path(process.goodVertexFilter + process.noScraping+process.triggerMuon)





