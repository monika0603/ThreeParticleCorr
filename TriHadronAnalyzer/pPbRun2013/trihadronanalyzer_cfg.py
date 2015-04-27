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

process.options = cms.untracked.PSet(
                                     makeTriggerResults = cms.untracked.bool(True),
                                     wantSummary = cms.untracked.bool(True)
                                     )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TriHadronAnalysis = cms.EDAnalyzer('TriHadronAnalyzer',
                                          verbose = cms.untracked.bool(True),
                                          vertexSrc = cms.InputTag("offlinePrimaryVerticesWithBS"),
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
                                          etaMinAsso = cms.double(0),
                                          etaMaxAsso = cms.double(1.2),
                                          ptMinTrg = cms.double(3.0),
                                          ptMaxTrg = cms.double(10.0),
                                          ptMinAsso = cms.double(3.0),
                                          ptMaxAsso = cms.double(10.0)
                                          )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("pPb_Pbp_CombinedDiHadronCorrFunc_NoCutetaTrg_EtaAsso_0_1_2_Pt3_10.root")
                                   )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.GlobalTag.globaltag = 'GR_P_V42_AN3::All'

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
                                            centralityVariable = cms.string("HFtowersPlusTrunc"),
                                            nonDefaultGlauberModel = cms.string(""),
                                            centralitySrc = cms.InputTag("pACentrality"),
                                            pPbRunFlip = cms.untracked.uint32(99999999)
                                            )

process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('AwaySideAnalysis.AwayAnalyzer.PAPileUpVertexFilter_cff')
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltSingleTrigger = process.hltHighLevel.clone()
process.hltSingleTrigger.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]

process.hltMult100 = process.hltHighLevel.clone()
process.hltMult100.HLTPaths = ["HLT_PAPixelTracks_Multiplicity100_v1",
                               "HLT_PAPixelTracks_Multiplicity100_v2"]

process.hltMult130 = process.hltHighLevel.clone()
process.hltMult130.HLTPaths = ["HLT_PAPixelTracks_Multiplicity130_v1",
                               "HLT_PAPixelTracks_Multiplicity130_v2"]

process.hltMult160 = process.hltHighLevel.clone()
process.hltMult160.HLTPaths = ["HLT_PAPixelTracks_Multiplicity160_v1",
                               "HLT_PAPixelTracks_Multiplicity160_v2"]

process.hltMult190 = process.hltHighLevel.clone()
process.hltMult190.HLTPaths = ["HLT_PAPixelTracks_Multiplicity190_v1",
                               "HLT_PAPixelTracks_Multiplicity190_v2"]

process.hltMult100.andOr = cms.bool(True)
process.hltMult100.throw = cms.bool(False)

process.hltMult130.andOr = cms.bool(True)
process.hltMult130.throw = cms.bool(False)

process.hltMult160.andOr = cms.bool(True)
process.hltMult160.throw = cms.bool(False)

process.hltMult190.andOr = cms.bool(True)
process.hltMult190.throw = cms.bool(False)

process.TriHadronAnalysisMult100 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(120),
                                                                 cutMultMax = cms.double(150)
                                                                 )

process.TriHadronAnalysisMult130 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(150),
                                                                 cutMultMax = cms.double(185)
                                                                 )

process.TriHadronAnalysisMult160 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(185),
                                                                 cutMultMax = cms.double(220)
                                                                 )

process.TriHadronAnalysisMult190 = process.TriHadronAnalysis.clone(
                                                                 cutMultMin = cms.double(220),
                                                                 cutMultMax = cms.double(260)
                                                                 )

process.Mult100 = cms.Path(process.hltMult100 *
                           process.PAcollisionEventSelection *
                           process.siPixelRecHits *
                           #     process.pileupVertexFilterCutGplus *
                           process.pACentrality *
                           process.TriHadronAnalysisMult100
                           )

process.Mult130 = cms.Path(process.hltMult130 *
                           process.PAcollisionEventSelection *
                           process.siPixelRecHits *
                           #     process.pileupVertexFilterCutGplus *
                           process.pACentrality *
                           process.TriHadronAnalysisMult130
                           )

process.Mult160 = cms.Path(process.hltMult160 *
                           process.PAcollisionEventSelection *
                           process.siPixelRecHits *
                           #     process.pileupVertexFilterCutGplus *
                           process.pACentrality *
                           process.TriHadronAnalysisMult160
                           )

process.Mult190 = cms.Path(process.hltMult190 *
                           process.PAcollisionEventSelection *
                           process.siPixelRecHits *
                           #      process.pileupVertexFilterCutGplus *
                           process.pACentrality *
                           process.TriHadronAnalysisMult190
                           )

process.schedule = cms.Schedule(process.Mult100,process.Mult130,process.Mult160,process.Mult190)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                                              '/store/hidata/HIRun2013/PAHighPt/RECO/PromptReco-v1/000/211/631/00000/FEDE0B60-3F75-E211-8FE3-003048D2BC5C.root'
                                                              )
                            )

















