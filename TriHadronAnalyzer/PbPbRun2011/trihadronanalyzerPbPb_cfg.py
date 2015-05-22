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
                                          trackSrc = cms.InputTag("hiGeneralTracks"),
                                          etaMin = cms.double(-3.0),
                                          etaMax = cms.double(3.0),
                                          ptMin = cms.double(0.4),
                                          vertexZMax = cms.double(15.0),
                                          qualityString = cms.string("highPurity"),
                                          etaBins = cms.vdouble(-2.4, 0.0, 2.4),
                                          vzBins = cms.vdouble(-15.0, -13.0, -11.0, -9.0, -7.0, -5.0, -3.0, -1.0, 1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0),
                                          NptBins = cms.vdouble(1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0),
                                          cutMultMin = cms.double(0),
                                          cutMultMax = cms.double(39),
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
                                          ptMaxTrg = cms.double(10.0),
                                          ptMinAsso1 = cms.double(2.0),
                                          ptMaxAsso1 = cms.double(3.0),
                                          ptMinAsso2 = cms.double(2.0),
                                          ptMaxAsso2 = cms.double(3.0),
                                          bkgFactor = cms.untracked.int32(10)
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

process.TriHadronAnalysisMult0_2 = process.TriHadronAnalysis.clone(
                                                                cutMultMin = cms.double(0),
                                                                cutMultMax = cms.double(1)
                                                                )

process.TriHadronAnalysisMult0_10 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(0),
                                                                 cutMultMax = cms.double(4)
                                                                 )

process.TriHadronAnalysisMult10_20 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(4),
                                                                 cutMultMax = cms.double(8)
                                                                 )

process.TriHadronAnalysisMult20_30 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(8),
                                                                 cutMultMax = cms.double(12)
                                                                 )

process.TriHadronAnalysisMult30_40 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(12),
                                                                 cutMultMax = cms.double(16)
                                                                 )

process.TriHadronAnalysisMult40_50 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(16),
                                                                 cutMultMax = cms.double(20)
                                                                 )

process.TriHadronAnalysisMult50_60 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(20),
                                                                 cutMultMax = cms.double(24)
                                                                 )

process.TriHadronAnalysisMult60_70 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(24),
                                                                 cutMultMax = cms.double(28)
                                                                 )

#process.TriHadronAnalysisMult70_80 = process.TriHadronAnalysis.clone(
#cutMultMin = cms.double(28),
#                                                                cutMultMax = cms.double(32)
#                                                                )

process.TriHadronAnalysisMultMinBias = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(0),
                                                                 cutMultMax = cms.double(39)
                                                                 )


process.Mult010 = cms.Path(process.hltSingleTrigger *
                           process.collisionEventSelection *
                           process.TriHadronAnalysisMult0_10
                           )

process.Mult1020 = cms.Path(process.hltSingleTrigger *
                           process.collisionEventSelection *
                           process.TriHadronAnalysisMult10_20
                           )

process.Mult2030 = cms.Path(process.hltSingleTrigger *
                           process.collisionEventSelection *
                           process.TriHadronAnalysisMult20_30
                           )

process.Mult3040 = cms.Path(process.hltSingleTrigger *
                           process.collisionEventSelection *
                           process.TriHadronAnalysisMult30_40
                           )

process.Mult4050 = cms.Path(process.hltSingleTrigger *
                           process.collisionEventSelection *
                           process.TriHadronAnalysisMult40_50
                           )

process.Mult5060 = cms.Path(process.hltSingleTrigger *
                           process.collisionEventSelection *
                           process.TriHadronAnalysisMult50_60
                           )

process.Mult6070 = cms.Path(process.hltSingleTrigger *
                            process.collisionEventSelection *
                            process.TriHadronAnalysisMult60_70
                            )

process.Mult002 = cms.Path(process.hltSingleTrigger *
                            process.collisionEventSelection *
                            process.TriHadronAnalysisMult0_2
                            )

process.MultMinBias = cms.Path(process.hltSingleTrigger *
                            process.collisionEventSelection *
                            process.TriHadronAnalysisMultMinBias
                            )

process.schedule = cms.Schedule(process.Mult010,process.Mult1020,process.Mult2030,process.Mult3040,
                                process.Mult4050,process.Mult5060,process.Mult6070,process.Mult002,
                                process.MultMinBias)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                                              '/store/hidata/HIRun2011/HIMinBiasUPC/RECO/14Mar2014-v2/00000/0018A8E7-F9AF-E311-ADAB-FA163E565820.root'
                                                              )
                            )

















