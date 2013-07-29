import FWCore.ParameterSet.Config as cms

process = cms.Process("EX")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START53_V7A::All'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

process.MessageLogger.cerr.FwkReport.reportEvery = 1
#
# the MC global Tag : START53_V7A
# the RERECO of 2012 data (Jan22 reRECO) for run ABC FT_53_V21_AN3
# input
#

process.source = cms.Source(
                            "PoolSource",
                            fileNames = cms.untracked.vstring(
                                                              'file:/sps/cms/hbrun/CMSSW_5_3_10_forNewSims/src/files/runDepMC/MCDY_runDep_1.root',
                                                              'file:/sps/cms/hbrun/CMSSW_5_3_10_forNewSims/src/files/runDepMC/MCDY_runDep_2.root',
                                                              'file:/sps/cms/hbrun/CMSSW_5_3_10_forNewSims/src/files/runDepMC/MCDY_runDep_3.root',
                                                              'file:/sps/cms/hbrun/CMSSW_5_3_10_forNewSims/src/files/runDepMC/MCDY_runDep_4.root',
                                                              'file:/sps/cms/hbrun/CMSSW_5_3_10_forNewSims/src/files/runDepMC/MCDY_runDep_5.root',
                                                              ),
                            secondaryFileNames = cms.untracked.vstring(),
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )



process.triggerMuon = cms.EDAnalyzer('TriggerMuon',
                                     muonProducer 		= cms.VInputTag(cms.InputTag("muons")),
                                     TriggerResults          = cms.InputTag("TriggerResults", "", "HLT"),
                                     HLTTriggerSummaryAOD    = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                                     outputFile		        = cms.string("muonTriggerTree.root"),
)

process.p = cms.Path(process.triggerMuon)