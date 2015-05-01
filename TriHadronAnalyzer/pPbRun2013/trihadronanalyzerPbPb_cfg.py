import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Geometry.CaloEventSetup.CaloTopology_cfi')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('CondCore.DBCommon.CondDBCommon_cfi');
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')

process.options = cms.untracked.PSet(
                                     makeTriggerResults = cms.untracked.bool(True),
                                     wantSummary = cms.untracked.bool(True)
                                     )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TriHadronAnalysis = cms.EDAnalyzer('TriHadronAnalyzer',
                                          verbose = cms.untracked.bool(True),
                                          vertexSrc = cms.InputTag("hiSelectedVertex"),
                                          trackSrc = cms.InputTag("generalTracks"),
                                          etaMin = cms.double(-3.0),
                                          etaMax = cms.double(3.0),
                                          ptMin = cms.double(0.4),
                                          vertexZMax = cms.double(15.0),
                                          qualityString = cms.string("highPurity"),
                                          etaBins = cms.vdouble(-2.4, 0.0, 2.4),
                                          vzBins = cms.vdouble(-15.0, -13.0, -11.0, -9.0, -7.0, -5.0, -3.0, -1.0, 1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0),
                                          NptBins = cms.vdouble(1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0),
                                          cutMultMin = cms.double(0.0),
                                          cutMultMax = cms.double(1000.0),
                                          cutDzErrMax = cms.untracked.double(3.0),
                                          cutDxyErrMax = cms.untracked.double(3.0),
                                          cutPtErrMax = cms.untracked.double(0.1),
                                          etaMinTrg = cms.double(-2.4),
                                          etaMaxTrg = cms.double(2.4),
                                          etaMinAsso1 = cms.double(-2.4),
                                          etaMaxAsso1 = cms.double(2.4),
                                          etaMinAsso2 = cms.double(-2.4),
                                          etaMaxAsso2 = cms.double(2.4),
                                          ptMinTrg = cms.double(3.0),
                                          ptMaxTrg = cms.double(4.0),
                                          ptMinAsso1 = cms.double(1.0),
                                          ptMaxAsso1 = cms.double(2.0),
                                          ptMinAsso2 = cms.double(1.0),
                                          ptMaxAsso2 = cms.double(2.0),
                                          bkgFactor = cms.untracked.int32(20)
                                          )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("TriHadronFirstAttempt.root")
                                   )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.GlobalTag.globaltag = 'GR_R_53_LV6::All'

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *

process.HeavyIonGlobalParameters = cms.PSet(
                                            centralityVariable = cms.string("HFtowers"),
                                            nonDefaultGlauberModel = cms.string(""),
                                            centralitySrc = cms.InputTag("hiCentrality")
                                            )

#process.option = cms.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltSingleTrigger = process.hltHighLevel.clone()
process.hltSingleTrigger.HLTPaths = ["HLT_HIMinBiasHfOrBSC_v*"]
process.hltSingleTrigger.andOr = cms.bool(True)
process.hltSingleTrigger.throw = cms.bool(False)

process.TriHadronAnalysisMult100 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(100),
                                                                 cutMultMax = cms.double(120)
                                                                 )

process.TriHadronAnalysisMult130 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(1200),
                                                                 cutMultMax = cms.double(1400)
                                                                 )

process.TriHadronAnalysisMult160 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(1400),
                                                                 cutMultMax = cms.double(1600)
                                                                 )

process.TriHadronAnalysisMult190 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(1600),
                                                                 cutMultMax = cms.double(1800)
                                                                 )

process.Mult100 = cms.Path(process.hltSingleTrigger *
                           process.PAcollisionEventSelection *
                           process.TriHadronAnalysisMult100
                           )

process.Mult130 = cms.Path(process.hltSingleTrigger *
                           process.PAcollisionEventSelection *
                           process.TriHadronAnalysisMult130
                           )

process.Mult160 = cms.Path(process.hltSingleTrigger *
                           process.PAcollisionEventSelection *
                           process.TriHadronAnalysisMult160
                           )

process.Mult190 = cms.Path(process.hltSingleTrigger *
                           process.PAcollisionEventSelection *
                           process.TriHadronAnalysisMult190
                           )

process.schedule = cms.Schedule(process.Mult100,process.Mult130,process.Mult160,process.Mult190)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                                              '/store/hidata/HIRun2011/HIMinBiasUPC/RECO/14Mar2014-v2/00000/0008E152-62AE-E311-BB69-FA163E565820.root'
                                                              )
                            )

















